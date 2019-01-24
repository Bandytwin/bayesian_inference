## Simple example of a Bayesian Inference approach to gaussian linear regression
## using height and weight data

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

# subset out children
d2 <- d[age >= 18]

# visualize heights
ggplot(d2) + geom_density(aes(height))

## function for calculating Gaussian posterior ----
gauss_posterior <- function(pars,
                            dat,
                            sd_prior = function(sd){return(1)},
                            mu_prior = function(mu){return(1)}) {
  mu <- pars[1]
  sd <- pars[2]
  posterior <- sum(log(dnorm(dat,mean = mu,sd = sd) * sd_prior(sd) * mu_prior(mu)))
  return(posterior)
}

## EXAMPLE 1 - Model height data ----

# define mean_prior
mu_prior <- function(mu) {
  return(dnorm(mu, mean=154, sd=20))
}

# define sd_prior
sd_prior <- function(sd) {
  return(dunif(sd, min=0, max=50))
}

# visualize posterior
post_viz_data <- data.table("mu" = rnorm(100,0,1) + d2[,mean(height)],
                            "sd" = runif(100,4,10),
                            "post" = rep(0,100))
for (ind in 1:nrow(post_viz_data)) {
  pars <- c(post_viz_data[ind,mu], post_viz_data[ind,sd])
  post_val <- gauss_posterior(pars,
                              dat = d2[,height])
  post_viz_data[ind,post := post_val]
}
ggplot(post_viz_data) + 
  geom_point(aes(mu,sd,colour=post)) +
  geom_point(data = post_viz_data[which.max(post)],aes(mu,sd),color="red")

## fit multivariate gaussian to posterior

# find mode
map_out <- optim(par = c(154,d2[,sd(height)]),
                 fn = gauss_posterior,
                 "dat" = d2[,height],
                 "sd_prior" = sd_prior,
                 "mu_prior" = mu_prior,
                 control = list("fnscale" = -1))

# calculate variance matrix at mode
v <- ginv(-hessian(func = gauss_posterior,
                   x = map_out$par,
                   "dat" = d2[,height],
                   "sd_prior" = sd_prior,
                   "mu_prior" = mu_prior))

# visualize posterior
post_samps <- data.table(mvrnorm(n = 1e4, mu = map_out$par,
                      Sigma = v))

ggplot(post_samps) + geom_point(aes(V1,V2),alpha=0.2)

# # sample around the mode
# multi_fit_samps <- data.table(expand.grid("mu_samps" = map_out$par[1] + seq(-0.2,0.2,0.02),
#                                           "sd_samps" = map_out$par[2] + seq(-0.2,0.2,0.02)))
# multi_fit_samps[, "post_prob" := 0]
# for (ind in 1:nrow(multi_fit_samps)){
#   pars <- multi_fit_samps[ind,c(mu_samps,sd_samps)]
#   multi_fit_samps[ind,
#                   "post_prob" := gauss_posterior(pars = pars,
#                                                  dat = d2[,height],
#                                                  sd_prior = sd_prior,
#                                                  mu_prior = mu_prior)]
# }
# 
# # function to fit multivariate
# mvtnorm_likelihood <- function(pars,
#                                dat) {
#   # pars is vector with the following elements
#   #   pars[1] mean of first column in dat
#   #   pars[2] mean of second column in dat
#   #   pars[3] var of first column in dat
#   #   pars[4] is covariances of columns
#   #   pars[5] var of second column in dat
#   
#   mean_vec <- c(pars[1],pars[2])
#   cov_mat <- t(matrix(c(pars[3],pars[4],pars[4],pars[5]),ncol=2))
#   # check if sigma is valid, otherwise return -1
#   if (any(abs(pars[4]) > abs(diag(cov_mat)))) {
#     return(-1)
#   } else {
#     approx_post <- c()
#     for (ind in 1:nrow(dat)) {
#       approx_post <- c(approx_post, log(dmvnorm(x = as.numeric(dat[ind,1:2]),
#                                                 mean = mean_vec,
#                                                 sigma = cov_mat)))
#     }
#     return(cor(approx_post,dat[,post_prob]))
#   }
# }
# 
# # fit multivariate to posterior
# mvt_optim_out <- optim(par = c(map_out$par[1],map_out$par[2],c(1,0.5,1)),
#                        fn = mvtnorm_likelihood,
#                        "dat" = multi_fit_samps,
#                        control = list("fnscale" = -1))
# 
# # extract parameters describing distribution
# mvt_params <- list("mean" = mvt_optim_out$par[1:2],
#                    "sigma" = t(matrix(c(mvt_optim_out$par[3:4],
#                                         mvt_optim_out$par[4:5]),
#                                       ncol = 2)))
# 
# # sample and visualize posterior
# post_samp <- rmvnorm(1e4, mean = mvt_params$mean, sigma = mvt_params$sigma)
# ggplot(data.table("mu_samp" = post_samp[,1], "sd_samp" = post_samp[,2])) + 
#   geom_point(aes(mu_samp,sd_samp),alpha = 0.2)

## EXAMPLE 2 - Linear model of height based on weight ----

# look into relationship between height and weight
ggplot(d2) + geom_point(aes(weight,height))

# define alpha_prior
alpha_prior <- function(alpha) {
  return(dnorm(alpha, mean=156, sd=100))
}

# define beta_prior
beta_prior <- function(beta) {
  return(dnorm(beta, mean=0, sd=10))
}

# define sd_prior
sd_prior <- function(sd) {
  return(dunif(sd, min=0, max=50))
}

# define gaussian linear model posterior
linear_gauss_post <- function(pars,
                              dat,
                              sd_prior = function(sd){return(1)},
                              alpha_prior = function(alpha){return(1)},
                              beta_prior = function(beta){return(1)}) {
  
  # extract/calculate relevant parameters and data
  alpha <- pars[1]
  beta <- pars[2]
  sd <- pars[3]
  height <- dat[,height]
  weight <- dat[,weight]
  mu_vec <- alpha + beta*weight
  
  # calculate and return posterior
  posterior <- sum(log(dnorm(height, mean = mu_vec, sd = sd) * 
                         alpha_prior(alpha) * 
                         beta_prior(beta) *
                         sd_prior(sd)))
  return(posterior)
}

# visualize posterior with flat priors
post_viz_data <- data.table("alpha" = rnorm(100,0,1),
                            "beta" = rnorm(100,0,1) + 4,
                            "sd" = runif(100,0,10),
                            "post" = rep(0,100))
for (ind in 1:nrow(post_viz_data)) {
  pars <- c(post_viz_data[ind,alpha], 
            post_viz_data[ind,beta],
            post_viz_data[ind,sd])
  post_val <- linear_gauss_post(pars,
                                dat = d2)
  post_viz_data[ind,post := post_val]
}
# alpha vs beta
ggplot(post_viz_data) + 
  geom_point(aes(alpha,beta,colour=post)) +
  geom_point(data = post_viz_data[which.max(post)],aes(alpha,beta),color="red")
ggplot(post_viz_data) + 
  geom_point(aes(beta,sd,colour=post)) +
  geom_point(data = post_viz_data[which.max(post)],aes(beta,sd),color="red")
  
## fit multivariate gaussian to posterior

# find mode
map_out <- optim(par = c(0,4,d2[,sd(height)]),
                 fn = linear_gauss_post,
                 "dat" = d2,
                 "sd_prior" = sd_prior,
                 "alpha_prior" = alpha_prior,
                 "beta_prior" = beta_prior,
                 control = list("fnscale" = -1))

u2 <- map_out$par

# find covariance matrix at mode
v2 <- ginv(-hessian(func = linear_gauss_post,
                   x = map_out$par,
                   "dat" = d2,
                   "sd_prior" = sd_prior,
                   "alpha_prior" = alpha_prior,
                   "beta_prior" = beta_prior))

# convert to correlation
c2 <- cov2cor(v2)

# # sample around the mode
# multi_fit_samps <- data.table(expand.grid("alpha_samps" = map_out$par[1] + seq(-20,20,2),
#                               "beta_samps" = map_out$par[2] + seq(-0.5,0.5,0.1),
#                               "sd_samps" = map_out$par[3] + seq(-0.5,0.5,0.1)))
# 
# # function to fit multivariate
# lin_mvtnorm_likelihood <- function(pars,
#                                dat) {
#   # pars is vector with the following elements
#   #   pars[1] mean of first column in dat
#   #   pars[2] mean of second column in dat
#   #   pars[3] var of first column in dat
#   #   pars[4] is covariances of columns
#   #   pars[5] var of second column in dat
#   
#   mean_vec <- c(pars[1],pars[2],pars[3])
#   cov_mat <- t(matrix(c(pars[4],pars[5],pars[6],
#                         pars[5],pars[7],pars[8],
#                         pars[6],pars[8],pars[9]),ncol=3))
#   
#   # make sure cov_mat is positive semi-definite
#   
#   
#   likelihood <- 0
#   for (ind in 1:nrow(dat)) {
#     likelihood <- likelihood + log(dmvnorm(x = as.numeric(dat[ind]),
#                                            mean = mean_vec,
#                                            sigma = cov_mat))
#   }
#   return(likelihood)
# }
# 
# # fit multivariate to posterior
# mvt_optim_out <- optim(par = c(map_out$par[1],map_out$par[2],map_out$par[3],
#                                c(1,0.5,0.5,
#                                  1,0.5,1)),
#                        fn = lin_mvtnorm_likelihood,
#                        "dat" = multi_fit_samps,
#                        control = list("fnscale" = -1))
# 
# # extract parameters describing distribution
# pars <- mvt_optim_out$par
# mvt_params <- list("mean" = mvt_optim_out$par[1:3],
#                    "sigma" = t(matrix(c(pars[4],pars[5],pars[6],
#                                         pars[5],pars[7],pars[8],
#                                         pars[6],pars[8],pars[9]),ncol=3)))
# 
# # convert covariance to correlation
# cor_mat <- cov2cor(mvt_params$sigma)

## EXAMPLE 3 - Posterior confidence intervals for average height ----
##  This section along with the next one  will be used to prototype a function 
##  for calculating confidence intervals.

# sample the posterior from Example 2
post_samps <- data.table(mvrnorm(n = 1e4, mu = u2, Sigma = v2))
setnames(post_samps,c("a","b","sigma"))

# show variance of height for person with weight 50
post_samps[, "h_at_50" := a + b*50]
ggplot(post_samps) +
  geom_density(aes(h_at_50)) +
  xlab("mu|weight=50")

# calculate 89% highest density interval when weight = 50
hdi_50 <- hdi(post_samps[,"h_at_50"], credMass = 0.89)

# function to get confidence at each sample weight
calc_confidence <- function(dat,
                            map_params,
                            n_samples = 1e4,
                            credMass = 0.89) {
  
  # create table of n_samples from map, and table of all unique weights
  map_samples <- data.table(mvrnorm(n = n_samples, mu = u2, Sigma = v2))
  setnames(map_samples,c("a","b","sigma"))
  conf_dat <- data.table("w" = dat[,unique(weight)])
  conf_dat[, c("low_c","high_c") := .(0,0)]
  
  # loop through each row in conf_dat, and calculate HDI
  for (ind in 1:nrow(conf_dat)) {
    h_dist <- map_samples[,a + b*conf_dat[ind,w]]
    h_hdi <- hdi(h_dist, credMass = credMass)
    conf_dat[ind, c("low_c", "high_c") := .(h_hdi[1], h_hdi[2])]
  }
  return(conf_dat)
}

# calculate confidences
conf <- calc_confidence(dat = d2,
                        map_params = list("u" = u2,
                                          "v" = v2))

# plot best line plus confidence intervals
best_line <- data.table("x" = c(d2[,min(weight)],d2[,max(weight)]),
                        "y" = u2[1] + u2[2]*c(d2[,min(weight)],d2[,max(weight)]))
ggplot(data = d2) + 
  geom_point(aes(weight,height)) +
  geom_ribbon(data = conf,
              aes(x = w, ymin = low_c, ymax = high_c),alpha=0.2) +
  geom_line(data=best_line,aes(x,y),linetype = "dashed")

## EXAMPLE 4 - Posterior confidence intervals for actual heights