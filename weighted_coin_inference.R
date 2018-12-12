## Description ----
## Perform Bayesian Inference on a weighted coin
##  - either provie a weight (as a binomial probability value) or a sample
##    measurement as a vector of 1's and 0's
##  - can choose either a uniform prior or a cosine prior (not sure if that's
##    technically the correct name, but you get the picture...)
##  
##  - returns posterior distribution approximations using MAP and
##    grid approximations
##    
## Generalized Function Definition for MAP approximation

## Libraries ----
library(stats)
library(numDeriv)
library(ggplot2)

## Function for MAP approximation ----
MAP_approximation <- function(map_dat,
                              posterior_func = NA,
                              show_plt = T) {
  ## INPUT
  ## post_dat : data_table with two columns:
  ##    param_val - value of parameter
  ##    posterior_p - posterior probability of the associate paramater value 
  ##    
  ## posterior_func : function for the posterior distribution, if provided then  
  ##    used to find the value and location of the posterior mode. Otherwise the
  ##    mode is estimated using mat_dat.
  ##    
  ## OUTPUT
  ## var_list : a list of variables that define the fitted parabola using vertex
  ##            notation (i.e. y = a*(x - v_x)^2 + v_y)
  ##    a - the fitted second derivative that minimizes error
  ##    v_x - the x location of the vertex
  ##    v_y - the y location of the vertex
  ## The posterior MAP approximation is then calculated as follows:
  ##    MAP approx = exp(y) = exp(a*(x - v_x)^2 + v_y)) 
  
  # add column for log(posterior)
  map_dat[, "log_posterior" := log(posterior_p)]
  
  # find mode of log-posterior
  if (!is.na(posterior_func)) {
    p_mode <- optimize(posterior_func, interval=c(0,1), maximum=TRUE)$maximum
  } else{
    p_mode <- map_dat[which.max(posterior_p), param_val]
  }
  
  # define parameters for apporximating parabola vertex
  v_x <- p_mode
  if (!is.na(posterior_func)) {
    v_y <- log(posterior_func(p_mode))
  } else {
    v_y <- map_dat[param_val == p_mode, log_posterior]
  }
  
  # fit quadratic
  map_dat[, c("v_x","v_y") := .(v_x,v_y)]
  nls_mod <- nls(log_posterior ~ d_two*(param_val - v_x)^2 + v_y, 
                 data = map_dat[is.finite(log_posterior)],
                 start = list("d_two" = -40))
  d_two <- as.numeric(nls_mod$m$getPars())

  # show how good of an approximation map is
  if (show_plt) {
    
    # add predictions to map_dat
    map_dat[is.finite(log_posterior), "map_approx" := predict(nls_mod)]

    # show raw fit on the logarithm of the posterior
    p2 <- ggplot(data = map_dat) +
      geom_line(aes(param_val,log_posterior,colour="Actual"),
                size=1.2) +
      geom_line(aes(param_val,map_approx,colour="Approx"),
                size=1.2,
                linetype='dashed')
    print(p2)
    
    # show the fit on the posterior
    p1 <- ggplot(data = map_dat) +
      geom_line(aes(param_val,exp(log_posterior),colour="Actual"),
                size = 1.2) +
      geom_line(aes(param_val,exp(map_approx),colour="Approx"),
                size=1.2,
                linetype="dashed")
    print(p1)
  }
  
  output_list <- list("a" = d_two,
                      "v_x" = v_x,
                      "v_y" = v_y)
  return(output_list)
}

## Function Definition for Binary Inference ----
weighted_coin_inference <- function(coin_p = 0.5,
                                    n_samps = 10,
                                    sample_flips = NA,
                                    prior = "uniform",
                                    show_plt = T) {
  
  # make sure that either coin_p or sample_flips is NA
  if (!is.na(coin_p) & !all(!is.na(sample_flips))) {
    stop("Either 'coin_p' or 'sample_flips' must be NA")
  }
  
  # set n_samps if sample_flips is used
  if (all(!is.na(sample_flips))) {
    n_samps <- length(sample_flips)
  }
  
  # create sample vector if coin_p is used
  if (!is.na(coin_p)) {
    if (is.na(n_samps)) {
      stop("If 'coin_p' is used then n_samps must be defined")
    }
    sample_flips <- rbinom(n = n_samps, size = 1, prob = coin_p)
  }
  
  # define prior
  if (prior == "uniform"){
    prior_func <- function(p_hat) {return(1)}
  } else if (prior == "cosine") {
    prior_func <- function(p_hat) {return(sin(p_hat * pi))}
  } else {
    stop("Only supported prior values are 'uniform' and 'cosine'")
  }
  
  # define likelihood
  likelihood_func <- function(p_hat) {
    n_ones <- sum(sample_flips)
    n_trials <- length(sample_flips)
    likelihood <- dbinom(x = n_ones, size = n_trials, prob = p_hat)
    return(likelihood)
  }
  
  # define function for posterior
  posterior_func <- function(p_hat) {
    posterior <- likelihood_func(p_hat) * prior_func(p_hat)
    return(posterior)
  }
  
  # approximate posterior using grid approximation
  p_grid <- seq(0,1,0.01)
  posterior_grid <- sapply(p_grid,posterior_func)
  posterior_grid <- posterior_grid / sum(posterior_grid)
  
  # create data_table for input to map approximation
  map_dat <- data.table("param_val" = p_grid,
                        "posterior_p" = posterior_grid)
  
  # call function that performs map approximation
  map_params <- MAP_approximation(map_dat = map_dat,
                                  posterior_func = NA)

  # perform quadratic approximation
  a <- map_params$a
  v_x <- map_params$v_x
  v_y <- map_params$v_y
  posterior_map <- exp(a*(p_grid - v_x)^2 + v_y)
  
  # build data table to return
  plot_data <- data.table("Coin Weight (P)" = p_grid,
                          "Grid Approx. Density" = posterior_grid,
                          "MAP Approx. Density" = posterior_map)
  
  # return
  return(list("posterior_data" = plot_data,
              "map_param" = map_params))
  
}

## Example from Textbook to validate ----
coin_p <- NA
n_samps <- NA
sample_flips <- rep(c(1,0,1,1,1,0,1,0,1),1)
prior <- "uniform"
show_plt <- T

inference_out <- weighted_coin_inference(coin_p = coin_p,
                                         n_samps = n_samps,
                                         sample_flips = sample_flips,
                                         prior = prior,
                        show_plt = show_plt)
