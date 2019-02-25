## Example of multivariate modeling using Mariage Data
## 

## Libraries ----
library(ggplot2)
library(stringr)
library(data.table)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(HDInterval)
library()

## Define MAP subfunction ----
fit_map <- function(posterior_fun = NA,
                    sample_data = NA,
                    start = NA) {
  # poster_fun := output of build_posterior() <- also defined in this file
  # sample data := data to use for finding the MAP estimator
  # start := defines the parameters to be searched over
  
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

## Example 1 - Predicting Divorce Rate ----
## Load and Pre-Process Data ----

# set path and load
repo_data_path <- "~/Desktop/statistical_rethinking/data"
d <- fread(str_c(repo_data_path,"/WaffleDivorce.csv"))

# standardize columns
d[, c("Marriage.s","MedianAgeMarriage.s") :=
    .((Marriage - mean(Marriage))/sd(Marriage),
      (MedianAgeMarriage - mean(MedianAgeMarriage))/sd(MedianAgeMarriage))]
           
## Define Posterior Log-Likelihood ----
multi_post <- function(pars,
                      dat) {
  
  # get parameters
  a <- pars[1]
  bR <- pars[2]
  bA <- pars[3]
  sigma <- pars[4]
  
  # get variables from data
  divorce <- dat[, Divorce]
  m <- dat[, Marriage.s]
  ma <- dat[, MedianAgeMarriage.s]
  
  # calculate secondary parameters
  mu <- a + bR*m + bA*ma
  
  # calculate priors
  a_prior <- dnorm(a, 10, 10)
  bR_prior <- dnorm(bR, 0, 1)
  bA_prior <- dnorm(bA, 0, 1)
  sigma_prior <- dunif(sigma, 0, 10)

  # calculate and return log likelihood
  log_like <- sum(log(dnorm(divorce, mean = mu, sd = sigma))) +
    log(a_prior) +
    log(bR_prior) +
    log(bA_prior) + 
    log(sigma_prior)
  
  return(log_like)
}

## Get MAP Approximation ----

# fit map
map_out <- fit_map(posterior_fun = multi_post,
                   sample_data = d,
                   start = c(10,0,0,1))

# get metrics
u <- data.table("Param" = c("a","bR","bA","sigma"),
                "Mean" = round(map_out$u,digits = 2),
                "StdDev" = round(sqrt(diag(map_out$sigma)),digits = 2))

Sig <- map_out$sigma

## Visualizations associated with multivariate models ----

# plot of parameter estimates and CIs
u[, "5.5%" := round(qnorm(0.055, mean = Mean,sd = StdDev),digits = 2)]
u[, "94.5%" := round(qnorm(0.945, mean = Mean, sd = StdDev), digits = 2)]
ggplot(u) + 
  geom_point(aes(Mean,Param),shape = 1,size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(y = Param, yend = Param, x = `5.5%`, xend = `94.5%`)) +
  xlab("Parameter Value") +
  ylab("Parameter Name") +
  ggtitle("89% Confidence Intervals on Parameter Values")

# add columns for each parameter's contribution
model_out <- d[, c("Divorce","Marriage.s","MedianAgeMarriage.s")]
model_out[,"a_contrib" := u[Param == 'a', Mean]]
model_out[,"sigma_contrib" := u[Param == 'sigma', Mean]]
model_out[,"bR_contrib" := u[Param == 'bR', Mean] * Marriage.s]
model_out[,"bA_contrib" := u[Param == 'bA', Mean] * MedianAgeMarriage.s]

# create counterfactual plots for bA and bR
p_bA <- ggplot(model_out, 
               aes(bA_contrib, Divorce - (a_contrib + bR_contrib))) +
  geom_point() +
  xlab("Prediction based on Age of Marriage") + 
  ylab("Divorce | intercept, Marriage Rate")
p_bR <- ggplot(model_out, 
               aes(bR_contrib, Divorce - (a_contrib + bA_contrib))) +
  geom_point() +
  xlab("Prediction based on Marriage Rate") + 
  ylab("Divorce | intercept, Age of Marriage")

## Example 2 - Simple example to identify spurious relationships ----
n <- 100
x_real <- rnorm(n)
x_spur <- rnorm(n, mean = x_real)
y <- rnorm(n, x_real)
dat <- data.table("n" = n,
                  "y" = y,
                  "x_real" = x_real,
                  "x_spur" = x_spur)

# show correlations between both regressors and y
p1 <- ggplot(dat,aes(x_real,y)) +
  geom_point() +
  geom_smooth() +
  xlab("x Real") +
  ylab("y")
p2 <- ggplot(dat,aes(x_spur,y)) +
  geom_point() +
  geom_smooth() +
  xlab("x Spur") +
  ylab("y")
grid.arrange(p1,p2,nrow=1)

# define model
model <- list("base_par_names" = c("a","b_spur","b_real","sig"),
              "base_priors" = alist(expression(dnorm(x = a,10,10)),
                                   expression(dnorm(b_spur,0,1)),
                                   expression(dnorm(b_real,0,1)),
                                   expression(dunif(sig,0,10))),
              "vars" = c("y","x_spur","x_real"),
              "secondary_pars" = alist(y_mu <- a + b_spur*x_spur + b_real*x_real),
              "likelihood" = expression(dnorm(y, y_mu, sig))
)

# define posterior
gen_posterior <- function(base_par_vals,
                          base_par_names,
                          base_priors,
                          vars,
                          secondary_pars,
                          likelihood = alist(""),
                          dat) {
  
  # Description of input arguments:
  #   base_par_vals := (vector) value of all parameters being treated as stochastic 
  #                variables. Will be named using base_par_names
  #   base_par_names := (vector)
  #   base_priors := (alist) list of priors associated with the base pars
  #                  i.e. "dnorm(par1, 0.1, 0.5)"
  #   vars := (vector) names of variables to extract from dat
  #   secondary_pars := (alist) parameters that are calculated directly from
  #                     either 'vars' or 'base_pars'
  #                     i.e. "a <- par1 + par2"
  #   likelihood := (string) defining the likelihood
  #   dat := (data.table) data for calculating posterior
  
  # store parameter values
  for (ind in 1:length(base_par_names)) {
    assign(base_par_names[ind],base_par_vals[ind])
  }
  
  # calculate prior contribution to likelihood
  prior_contrib <- sum(log(sapply(base_priors, eval)))
  
  # get variables from data
  for (v_name in vars) {
    assign(paste(v_name),dat[,get(v_name)])
  }  
  
  # calculate secondary parameters
  for (par in secondary_pars) {
    eval(par)
  }
  browser()
  # calculate and return log likelihood
  log_like <- sum(log(eval(likelihood))) + prior_contrib
  
  return(log_like)
}

# fit MAP
fit_map <- function(posterior_fun = gen_posterior,
                    model = NA,
                    sample_data = NA,
                    start = NA) {
  # poster_fun := output of build_posterior() <- also defined in this file
  # sample data := data to use for finding the MAP estimator
  # start := defines the parameters to be searched over
  
  # create list of args to optim
  optim_args <- list("par" = start,
                     "fn" = gen_posteiror,
                     "dat" = sample_data,
                     "control" = list("fnscale" = -1),
                     model)
  
  ## fit multivariate gaussian to posterior
  # find mode
  map_out <- do.call(optim,optim_args)
  
  u <- map_out$par
  
  # find covariance matrix at mode
  v <- ginv(-hessian(func = linear_gauss_post,
                     x = map_out$par,
                     "dat" = sample_data))
  
  # return the mean and covariance
  return(list("u" = u, "sigma" = v))
  
}

## Example 3 - Simple example to identify masked relationships ---- 