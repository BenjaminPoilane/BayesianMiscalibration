### DESCRIPTION ###############################################################
#   This script provides :
#     - functions which compute the posteriors of a location-scale model when 
#       the likelihood is a Gaussian distribution, a Laplace distribution or
#       a Student distribution with fixed degrees of freedom.
#     - functions to compute the Hessian matrices of the log-posteriors of 
#       these models rescaled by - 1 / n.
#     - a function to generate observations according to a Laplace distribution
###############################################################################

### normal() ##################################################################
# normal() returns unstandardized log-posterior of a normal model.
# likelihood is :
#         dnorm(x, mu, sigma)
# and parameters are theta = c(mu, log(sigma)), prior should be on these param-
# eters.
####
# x     must be a vector
# prior must be a function of R ^ 2 * k to R ^ k with an option log = TRUE to 
#       return log-prior.
# log   must be a logical. If log is TRUE, unstandardized log-posterior is ret-
#       urned, else if log is FALSE, unstandardized posterior is returned
####
# returns a function taking as argument a matrix 2 * k and returning a vector 
# of k values corresponding to (log-)posterior of the k-columns of argument.
#########
normal <- function(x, prior, 
                   log = TRUE) {
  n <- length(x)
  logpost <- function(theta) {
    theta  <- as.matrix(theta)
    
    mu        <- theta[1, ] # location parameter
    logsigma  <- theta[2, ] # log-sd of normal dist. 
    
    k <- dim(theta)[2]   # number of theta values where to evaluate logpost()
    
    sigma <- exp(logsigma)
    logprior <- prior(theta, log = TRUE)
    
    x_mat  <- t(x) %x% rep(1, k)
    
    mu_mat    <- mu    %x% t(rep(1, n))
    sigma_mat <- sigma %x% t(rep(1, n))
    
    # std residuals
    res <- (x_mat - mu_mat) / sigma_mat
    # log-likelihood
    loglikelihood <- -0.5 * rowSums(res ^ 2) - n * log(sigma)
    
    return(loglikelihood + logprior)
  }
  
  if (log) {
    return(logpost)
  } else {
    post <- function(theta) { 
      return(exp(logpost))
    }
    return(post)
  }
}

### laplace() #################################################################
# laplace() returns unstandardized log-posterior of a normal model.
# likelihood is :
#         dnorm(x, mu, sigma)
# and parameters are theta = c(mu, log(sigma)), prior should be on these param-
# eters.
####
# x     must be a vector
# prior must be a function of R ^ 2 * k to R ^ k with an option log = TRUE to 
#       return log-prior.
# log   must be a logical. If log is TRUE, unstandardized log-posterior is ret-
#       urned, else if log is FALSE, unstandardized posterior is returned
####
# returns a function taking as argument a matrix 2 * k and returning a vector 
# of k values corresponding to (log-)posterior of the k-columns of argument.
#########
laplace <- function(x, prior,
                    log = TRUE) {
  n <- length(x)
  logpost <- function(theta) {
    theta  <- as.matrix(theta)
    
    mu        <- theta[1, ] # location parameter
    logsigma  <- theta[2, ] # log-sd of normal dist.
    
    k <- dim(theta)[2]   # number of theta values where to evaluate logpost()
    
    sigma <- exp(logsigma)
    
    # logarithm of prior evaluated at all thetas
    logprior <- prior(theta, log = TRUE)
    
    x_mat  <- t(x) %x% rep(1, k)
    
    mu_mat    <- mu    %x% t(rep(1, n))
    sigma_mat <- sigma %x% t(rep(1, n))
    
    # std residuals
    res <- (x_mat - mu_mat) / sigma_mat
    
    loglikelihood <- rowSums(dexp(abs(res), rate = sqrt(2), log = TRUE) - log(sigma_mat))
    
    return(loglikelihood + logprior)
  }
  if (log) {
    return(logpost)
  } else {
    post <- function(theta) {
      return(exp(logpost))
    }
    return(post)
  }
}

### student_fixed_df() ########################################################
# student_fixed_df() returns unstandardized log-posterior of a student model : 
#         X = mu + sigma * t  , with t ~ student(df) where df is fixed.
# likelihood is :
#         dt((x - mu) / sigma, df) / sigma
# and parameters are theta = c(mu, log(sigma)), prior should be on 
# these parameters.
####
# x     must be a vector
# prior must be a function of R ^ 2 * k to R ^ k with an option log = TRUE to 
#       return log-prior.
# dgf   must be a positive real number corresponding to the degrees of freedom.
# log   must be a logical. If log is TRUE, unstandardized log-posterior is ret-
#       urned, else if log is FALSE, unstandardized posterior is returned.
####
# returns a function taking as argument a matrix 2 * k and returning a vector 
# of k values corresponding to (log-)posterior of the k-columns of argument.
#########
student_fixed_df <- function(x, prior, dgf, 
                             log = TRUE) {
  n <- length(x)
  logpost <- function(theta) {
    theta  <- as.matrix(theta)
    
    mu        <- theta[1, ] # location parameter
    logsigma  <- theta[2, ] # log-sd of normal dist. 
    
    k <- dim(theta)[2]   # number of theta values where to evaluate logpost()
    
    dgf <- rep(dgf, k)
    sigma <- exp(logsigma)
    
    # logarithm of prior evaluated at all thetas
    logprior <- prior(theta, log = TRUE)
    
    x_mat  <- t(x) %x% rep(1, k)
    
    mu_mat    <- mu    %x% t(rep(1, n))
    sigma_mat <- sigma %x% t(rep(1, n))
    dgf_mat   <- dgf   %x% t(rep(1, n))
    
    # std residuals
    res <- (x_mat - mu_mat) / sigma_mat
    
    loglik1 <- -(dgf + 1) / 2 * rowSums(log(1 + res ^ 2 / dgf_mat))
    loglik2 <- n * (lgamma(dgf / 2 + 0.5) - lgamma(dgf / 2) - 0.5 * log(pi * dgf))
    loglikelihood <- loglik1 + loglik2 - n * logsigma
    
    return(loglikelihood + logprior)
  }
  if (log) {
    return(logpost)
  } else {
    post <- function(theta) { 
      return(exp(logpost(theta)))
    }
    return(post)
  }
}

### hessian_normal() ##########################################################
# hessian_normal() returns -1 / n times Hessian of log-posterior of normal mod-
# el with prior having hessian equal to hessian_prior, with respect to paramet-
# ers theta = (mu, log(sigma))
####
# x             must be a vector of size n
# theta         must be a vector of size 2
# hessian prior must be a 2 * 2 matrix
####
# returns a 2 * 2 matrix
#########
hessian_normal <- function(x, theta, hessian_prior) {
  mu    <- theta[1]
  sigma <- exp(theta[2])
  
  # sufficient statistics
  mx <- mean(x)
  sx <- mean((x - mx) ^ 2)
  
  # Hessian of log-likelihood of normal post
  H <- matrix(nrow = 2, ncol = 2, 
              dimnames = list(c("mu", "logsigma"), c("mu", "logsigma")))
  H["mu", "mu"]    <- -1
  H["mu", "logsigma"] <- H["logsigma", "mu"] <- -2 * (mx - mu)
  H["logsigma", "logsigma"] <- -2 * (sx + (mx - mu) ^ 2)
  
  
  return(-H / sigma ^ 2 - hessian_prior)
}

### hessian_student_fixed_df() ################################################
# hessian_student_fixed_df() returns -1 / n times Hessian of log-posterior of 
# student model with prior having hessian equal to hessian_prior, with respect 
# to parameters :
#         theta = (mu, log(sigma), log(dgf))
####
# x             must be a vector of size n
# theta         must be a vector of size 3
# dgf           must be a numeric
# hessian prior must be a 3 * 3 matrix
####
# returns a 2 * 2 matrix
#########
hessian_student_fixed_df <- function(x, theta, hessian_prior, dgf) {
  n <- length(x)
  mu    <- theta[1]
  sigma <- exp(theta[2])
  
  # Hessian of log-likelihood of student model
  H <- matrix(nrow = 2, ncol = 2, 
              dimnames = list(c("mu", "logsigma"), 
                              c("mu", "logsigma")))
  
  res <- (x - mu) / sigma
  
  denom <- (dgf + res ^ 2)
  
  H["mu", "mu"] <- -(dgf + 1) / sigma ^ 2 * mean((dgf - res ^ 2) / denom ^ 2)
  H["logsigma", "logsigma"] <- -2 * (dgf + 1) * dgf * mean((res / denom) ^ 2)
  
  H["mu", "logsigma"] <- -(dgf + 1) * dgf * 2 * mean(res / denom ^ 2) / sigma
  H["logsigma", "mu"] <- H["mu", "logsigma"]
  
  return(-H - hessian_prior / n)
}

### hessian_laplace() #########################################################
# hessian_laplace() returns -1 / n times Hessian of log-posterior of Laplace 
# model with prior having hessian equal to hessian_prior, with respect to para-
# meters theta = (mu, log(sigma))
####
# x             must be a vector of size n
# theta         must be a vector of size 2
# hessian prior must be a 2 * 2 matrix
####
# returns a 2 * 2 matrix
#########
hessian_laplace <- function(x, theta, hessian_prior) {
  mu    <- theta[1]
  sigma <- exp(theta[2])
  
  res = (x - mu)/ sigma
  
  # Hessian of log-likelihood of normal post
  H <- matrix(nrow = 2, ncol = 2,
              dimnames = list(c("mu", "logsigma"), c("mu", "logsigma")))
  H["mu", "mu"]    <- -2 * 1 # as.numeric(sum(res == 0))
  H["mu", "logsigma"] <- H["logsigma", "mu"] <- -sum(sign(res))
  H["logsigma", "logsigma"] <- 1 - 2 * sum(abs(res))
  
  return(-H / sigma ^ 2 - hessian_prior)
}


### rl() ######################################################################
# rl() returns a vector of n real values distributed according to a Laplace 
# distribution
####
# n must be a positive integer
####
# returns a vector of size n
#########
rl <- function(n) {
  return((2 * rbinom(n, size = 1, prob = 0.5) - 1)* rexp(n, rate = sqrt(2)))
}
