## Simple example of a Bayesian Inference approach to polynomial regression on 
## height data ----

## libraries ----
library(ggplot2)
library(stringr)
library(data.table)
library(stats)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(HDInterval)

## load data and prepare ----

# file can also be found at 
repo_data_path <- "~/Desktop/statistical_rethinking/data/Howell1.csv"
d <- fread(repo_data_path)

# plot raw data
ggplot(d) +
  geom_point(aes(weight,height,shape='a')) +
  scale_shape_discrete(solid = F)

# sandardize weight data
d[,"s_weight" := (weight - mean(weight))/sd(weight)]

# define polynomial posterior
poly_post <- function(pars,
                      dat,
                      sd_prior = function(sd){dunif(sd,0,50)},
                      alpha_prior = function(alpha){dnorm(alpha,178,100)},
                      b1_prior = function(b1){dnorm(b1,0,30)},
                      b2_prior = function(b2){dnorm(b2,0,30)}) {
  
  # extract/calculate relevant parameters and data
  sd <- pars[1]
  alpha <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  height <- dat[,height]
  weight <- dat[,s_weight]
  mu_vec <- alpha + b1*weight + b2*weight^2

  # calculate and return posterior
  posterior <- sum(log(dnorm(height, mean = mu_vec, sd = sd))) + 
                         log(alpha_prior(alpha)) +
                         log(b1_prior(b1)) +
                         log(b2_prior(b2)) +
                         log(sd_prior(sd))
  return(posterior)
}

## function to fit Gaussian to posterior
fit_map <- function(posterior_fun = NA,
                    sample_data = NA,
                    start = NA) {
  
  ## fit multivariate gaussian to posterior
  # find mode
  map_out <- optim(par = start,
                   fn = posterior_fun,
                   "dat" = sample_data,
                   control = list("fnscale" = -1))
  
  u <- map_out$par

  # find covariance matrix at mode
  v <- ginv(-hessian(func = posterior_fun,
                     x = map_out$par,
                     "dat" = sample_data))
  
  # return the mean and covariance
  return(list("u" = u, "sigma" = v))
  
}

# fit map
map_out <- fit_map(posterior_fun = poly_post,
                   sample_data = d,
                   start = c(5,150,20,-8))
u <- map_out$u
v <- map_out$sigma

# function to get height sample confidence at each sample weight
calc_confidence <- function(dat,
                            map_params,
                            n_samples = 1e4,
                            credMass = 0.89) {
  
  # create table of n_samples from map, and table of all unique weights
  map_samples <- data.table(mvrnorm(n = n_samples, mu = map_params$u, Sigma = map_params$v))
  setnames(map_samples,c("sigma","a","b1","b2"))
  conf_dat <- data.table("w" = dat[,unique(s_weight)])
  conf_dat[, c("low_c","high_c") := .(0,0)]
  
  # loop through each row in conf_dat, and calculate HDI
  for (ind in 1:nrow(conf_dat)) {
    h_dist <- map_samples[,rnorm(a,mean = a + b1*conf_dat[ind,w] + b2*conf_dat[ind,w]^2, sd = sigma)]
    h_hdi <- hdi(h_dist, credMass = credMass)
    conf_dat[ind, c("low_c", "high_c") := .(h_hdi[1], h_hdi[2])]
  }
  return(conf_dat)
}

# calculate confidences
h_sample_conf <- calc_confidence(dat = d,
                                 map_params = list("u" = u,
                                                   "v" = v))

# visualize best line with confidence intervals
best_line <- data.table("s_weight" = d[order(weight),unique(s_weight)])
best_line[, "mean_h" := u[2] + u[3]*s_weight + u[4]*s_weight^2]

ggplot(d) +
  geom_point(aes(s_weight,height,shape='a'),color = "blue") +
  scale_shape_discrete(solid = F) +
  geom_line(data = best_line,aes(s_weight,mean_h)) +
  geom_ribbon(data = h_sample_conf,
              aes(x = w, ymin = low_c, ymax = high_c),alpha=0.2)
