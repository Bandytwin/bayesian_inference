## This file defines functions to perform MAP on a posterior with an arbitrary 
##  number of parameters and common probability distributions. 
##  
##  A lot of this work was prototyped in statistical_inference/linear_model_height.R
##    and then functionalized here
##  
##  Author: Andrew Cowley (this is not directly from the textbook)
##  
##  Function Args:
##    parameter_dists := (list) describes the model by defining the liklihood
##                        along with all the prior distributions (see examples
##                        below for list formatting)
##          samp_data := (vector) samples to use for inference--when combined
##                        with parameters this should be enough to call the 
##                        relvant density funciton (i.e. dbinom() for binomial
##                        requires both a value x and the sample size)
##              start := (named list) start point of gradient ascent, if NA then
##                        the mean of the priors is used.
##                        
##  Function Output (list):
##    

## libraries ----
library(data.table)
library(stringr)
library(lubridate)
library(numDeriv)

## test with binomial data ----
# actual_p <- 0.2
# sample_size <- 10
# binom_data <- list("size" = sample_size,
#                    "x" = rbinom(n = 1, size = sample_size, prob = actual_p))
# 
# # create lists for each relevant distribution
# parameter_dists <- list("likelihood" = list("dist" = "binomial",
#                                            "args" = list("prob" = "binom_prob",
#                                                          "size" = 1)),
#                         "prior_prob" = list("dist" = "uniform",
#                                             "args" = list("max" = 1, "min" = 0)))
# 
# calc_prob(parameter_dists,x = list("prob" = 0.4),data=binom_data)
# 
# ## test with normal data ----
# actual_mean = 1
# actual_sd = 0.5
# sample_size = 10
# normal_data = rnorm(sample_size, mean = actual_mean, sd = actual_sd)
# 
# # define gaussian linear model posterior
# linear_gauss_post <- function(pars,
#                               dat,
#                               sd_prior = function(sd){dunif(sd, min=0, max=50)},
#                               mean_prior = function(mean){dnorm(mean, mean=0, sd=10)}) {
#   
#   # extract/calculate relevant parameters and data
#   mu <- pars[1]
#   sd <- pars[2]
#   
#   # calculate and return posterior
#   posterior <- sum(log(dnorm(dat, mean = mu, sd = sd) * 
#                          mean_prior(mu) *
#                          sd_prior(sd)))
#   return(posterior)
# }
# 
# # fit map
# map <- fit_map(posterior_fun = linear_gauss_post,
#                sample_data = normal_data,
#                start = c(0,1))

## build posterior function given specifications ----
posterior_func <- funciton(pars,
                           dat,
                           formulas) {
  
  
  
  # extract/calculate relevant parameters and data
  sd <- pars[1]
  alpha <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  height <- dat[,height]
  weight <- dat[,s_weight]
  mu_vec <- alpha + b1*weight + b2*weight^2
  
  # define posterior function in global environment
  new_fun <- function(pars,
                      dat) {
    
    
    
  }
  
    
  # assign to global environment
  assign(function_name, new_fun, envir = .GlobalEnv)
  
  # calculate and return posterior
  posterior <- sum(log(dnorm(height, mean = mu_vec, sd = sd))) + 
    log(alpha_prior(alpha)) +
    log(b1_prior(b1)) +
    log(b2_prior(b2)) +
    log(sd_prior(sd))
  return(posterior)
}

## fit MAP   ----
fit_map <- function(posterior_fun = NA,
                    sample_data = NA,
                    start = NA) {
  # poster_fun := output of build_posterior() <- also defined in this file
  # sample data := data to use for finding the MAP estimator
  # start := defines the parameters to be searched over
  
  ## fit multivariate gaussian to posterior
  # find mode
  map_out <- optim(par = start,
                   fn = linear_gauss_post,
                   "dat" = sample_data,
                   control = list("fnscale" = -1))
  
  u <- map_out$par
  
  # find covariance matrix at mode
  v <- ginv(-hessian(func = linear_gauss_post,
                     x = map_out$par,
                     "dat" = sample_data))
  
  # return the mean and covariance
  return(list("u" = u, "sigma" = v))
  
}

## Function to recursively generate distribution samples ----
calc_sample <- function(parameter,
                        dist_list) {
  
  # loop through the arguments for the distribution function and build up values
  arg_actuals <- list("n" = 1)
  for (arg_name in names(dist_list[[parameter]]$args)) {
    arg <- dist_list[[parameter]]$args[[arg_name]]
    if (dist_list[[parameter]]$args[[arg_name]] %in% names(dist_list)) {
      input_val <- calc_sample(parameter = arg,
                               dist_list = dist_list)
    } else {
      input_val <- arg
    }
    arg_actuals[[arg_name]] <- input_val
  }

  # sample using the arg_actuals
  if (parameter_dists[[parameter]]$dist == "uniform") {
    return(do.call(runif,arg_actuals))
  } else if (parameter_dists[[parameter]]$dist == "binomial") {
    return(do.call(rbinom,arg_actuals))
  } else if (parameter_dists[[parameter]]$dist == "normal") {
    return(do.call(rnorm,arg_actuals))
  }
  
}

## function to calculate un-normalized point probability given a set of
##  parameters, data, and prior + likelihood distributions ----
calc_prob <- function(dist_list,
                      x,
                      data) {
  
  # check that prior distributions exists for each parameter
  prior_names <- list()
  for (i_param in names(x)){
    if (str_c("prior_",i_param) %in% names(dist_list)){
      prior_names[[i_param]] <- str_c("prior_",i_param)
    } else {
      stop("Provide a prior for each parameter")
    }
  }
  
  # check that a value for each parameter is provided
  if (length(prior_names) != length(x)) {
    stop("Provide a value for each parameter to use for calculating prior 
         probability")
  }

  # loop through each parameter and calculate prior probability
  prob_list <- list()
  for (i_param in names(prior_names)) {
    i_dist <- dist_list[[prior_names[[i_param]]]]
    prob_func_args <- list("x" = x[[i_param]])
    
    # loop through each parameter and build up arg list for probability funcitons
    for (i_arg in names(i_dist$args)) {
      i_arg_val <- i_dist$args[[i_arg]]
      
      # if the parameter has a distribution then use x, 
      #   otherwise use provided value
      if (is.character(i_arg_val)) {
        prob_func_args[[i_arg]] <- x[[i_param]]
      } else {
        prob_func_args[[i_arg]] <- i_arg_val
      }
    }
    
    # calculate probability from prior
    if (i_dist$dist == "uniform") {
      i_prob <- do.call(dunif,prob_func_args)
    } else if (i_dist$dist == "binomial") {
      i_prob <- do.call(dbinom,prob_func_args)
    } else if (i_dist$dist == "normal") {
      i_prob <- do.call(dnorm,prob_func_args)
    }
      
    prob_list[[i_param]] <- i_prob
  }
  
  # calculate liklihood using data
  if (dist_list$likelihood$dist == "uniform") {
    prob_list["likelihood"] <- do.call(dunif,c(data,x))
  } else if (dist_list$likelihood$dist == "binomial") {
    prob_list["likelihood"] <- do.call(dbinom,c(data,x))
  } else if (dist_list$likelihood$dist == "normal") {
    prob_list["likelihood"] <- do.call(dnorm,c(data,x))
  }
  
  return(prob_list)
  
}