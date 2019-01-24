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
actual_p <- 0.2
sample_size <- 10
binom_data <- list("size" = sample_size,
                   "x" = rbinom(n = 1, size = sample_size, prob = actual_p))

# create lists for each relevant distribution
parameter_dists <- list("likelihood" = list("dist" = "binomial",
                                           "args" = list("prob" = "binom_prob",
                                                         "size" = 1)),
                        "prior_prob" = list("dist" = "uniform",
                                            "args" = list("max" = 1, "min" = 0)))

calc_prob(parameter_dists,x = list("prob" = 0.4),data=binom_data)

## test with normal data ----
actual_mean = 1
actual_sd = 0.5
sample_size = 10
normal_data = rnorm(sample_size, mean = actual_mean, sd = actual_sd)

## (Main Function) fit_map function ----
fit_map <- function(parameter_dists = NA,
                    samp_data = NA,
                    start = NA) {
  
  # check that parameter dists is specified
  if (!is.list(parameter_dists)) {
    stop("Must provide parameter distribtuions")
  }
  
  # check if start is specified, and if it isn't then initialize
  if (!is.list(start)) {
    
  }
  
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

## Function to calculate confidence intervals given data and map fit output ----