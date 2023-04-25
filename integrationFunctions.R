require(statmod)

### DESCRIPTION ###############################################################
#   This script provides functions for marginalization of posteriors using 
#   Gauss-Hermite quadrature as well as a function to compute quantiles of a 
#   distribution given its density.
###############################################################################

### gh_regular_grid() #########################################################
# gh_regular_grid() returns the weights and nodes to perform Gauss-Hermite qua-
# drature on R ^ dimension, i.e. to approximate
#              integral over R ^ d of f(x) * phi(x) dx 
# where phi is a standard multivariate normal density. 
# There will be n_nodes_1d ^ dimension
####
# requires library(statmod)
####
# dimension   must be a positive integer.
# n_nodes_1d  must be a positive integer.
####
# returns a list with elements nodes : a n_nodes_1d ^ dimension * dimesion mat-
# rix, and weights a vector of size n_nodes_1d ^ dimension.
#########
gh_regular_grid <- function(dimension, n_nodes_1d) {
  
  # nodes and weighs for 1 dimension
  gq <- gauss.quad(n_nodes_1d, kind = "hermite")
  weights_1D <- gq$weights
  nodes_1D   <- sqrt(2) * gq$nodes
  
  #  cartesian product of nodes, and product of weights.
  ln <- lw <- list()
  for (i in 1:dimension) {
    ln[[i]] <- nodes_1D
    lw[[i]] <- weights_1D
  }
  nodes <- expand.grid(ln)
  weights <- exp(rowSums(log(expand.grid(lw))))
  
  return(list(weights = weights, nodes = nodes))
}

######### marginal_mu() #######################################################
# marginal_mu() comupte maginal distribution of mu from logposterior over para-
# meters (mu, lambda) where lambda is vector of R ^ (d - 1) grouping all param-
# eters different from mu. Integration is approximated via Gauss-Hermite quadr-
# ature with n_gh ^ (d - 1) nodes, centered around map[-1] where map is max a 
# posteriori. hessian should be 1 / n * hessian of the negative logposterior at
# map. 
####
# logpost   should be a function taking a d * n matrix and returning a vector 
#           of size n corresponding to log-posterior evaluated at all columns.
# map       should be a vector of size d maximizing logpost over R ^ d.
# hessian   should be d * d matrix : Hessian(logpost) at map
# n_gh      should be a positive integer not taken into account if rect_grid is
#           specified
# rect_grid should be a list with elements nodes : matrix with (d - 1) columns
#           and weights a vector
# mu        should be a vector of values where to evaluate marginal posterior
# log       should be a logical. If log = TRUE, log-marginal posterior is retu-
#           rned, marginal posterior is returned otherwise
####
# returns a vector of same size than mu.
#########
marginal_mu <- function(logpost, map, hessian, mu, n,
                        n_gh = 20, rect_grid = NULL, log = TRUE) {
  # Number of parameter on wich log-posterior f is defined
  d <- length(map)[1]
  
  # Gauss-Hermite quadrature grid on regular d - 1 space
  if (is.null(rect_grid)) {
    rect_grid <- gh_regular_grid(d - 1, n_gh)  
  }
  nodes   <- rect_grid$nodes
  weights <- rect_grid$weights
  
  # transformation of nodes
  H <- hessian * n
  H_lambda  <- as.matrix(H[-1, -1])
  nodes_0   <- solve(chol(H_lambda), t(nodes))
  
  # center of grid
  mu_map <- map[1]
  lambda_map <- map[-1]
  
  H_mu <- H[1, 1]
  H_mu_lambda <- H[1, -1]
  H_lambda_mu <- H[-1, 1]
  
  # lambda_0 <- lambda_map %x% t(rep(1, n_gh ^ (d - 1)))
  
  M <- solve(H_lambda) %*% H_lambda_mu
  
  lpm <- logmarginal <- vector("numeric", length(mu))
  
  for (i in 1:length(mu)) {
    # center of grid
    lambda_c <- lambda_map - (M %*% (mu[i] - mu_map))
    
    # grid over lambda space
    lambda   <- lambda_c %x% t(rep(1, n_gh ^ (d - 1))) + nodes_0
    
    # grid over lambda with fixed mu value added : 
    theta <- rbind(rep(mu[i] - mu_map, n_gh ^ (d - 1)), lambda - lambda_map)
    
    # evaluation of order two remainder off logpost at grid
    f_1 <- logpost(rbind(rep(mu[i], n_gh ^ (d - 1)), lambda)) 
    f_2 <- 0.5 * colSums(theta * (H %*% theta))
    f_remain <- f_1 + f_2
    
    # reste[i, ] <- f_remain
    f_remain_max <- max(f_remain)
    
    # Gauss-Hermite quadrature on lambda dependent terms
    logmarginal[i] <- log(sum(weights * exp(f_remain - f_remain_max)))
    
    # gaussian term of mu only :
    # inverse variance
    H1_mu <- (H_mu - H_mu_lambda %*% solve(H_lambda) %*% H_lambda_mu)
    # quadratic term
    log_phi_mu <- -0.5 * (mu[i] - mu_map) ^ 2 * H1_mu
    lpm[i] <- -0.5 * (mu[i] - mu_map) ^ 2 * H1_mu
    # total log-marginal with mu-only terms added.
    logmarginal[i] <- logmarginal[i] + log_phi_mu + f_remain_max
  }
  if (log) {
    return(logmarginal)
  } else {
    return(exp(logmarginal))
  }
}

### marginal_mu_fct() #########################################################
# marginal_mu_fct() returns a function of mu which compute the marginal of
# the logpost with respect to parameters different than mu evaluated at mu.
# It uses the marginal_mu() function.
####
# logpost   should be a function taking a multidimensional argument, call d the
#           dimension, whose first dimension corresponds to the parameter whose 
#           marginal distribution the function will return.
# map       must be a vector of size d, corresponding to the (d-dimensional) 
#           maximum of logpost.
# hessian   must be a d * d matrix corresponding to the hessian matrix of 
#           logpost evaluated at map.
# n         must be a positive integer (the sample size of data that yield 
#           logpost)
# n_gh      must be a positive integer corresponding to the number of nodes used
#           for the Gauss-Hermite quadrature.
# rect_grid should be a list with elements nodes : matrix with (d - 1) columns
#           and weights a vector.
# log       should be a logical. If log = TRUE, log-marginal posterior is retu-
#           rned, marginal posterior is returned otherwise
####
# Returns a function which takes a vector mu and returns a vector of the same 
# corresponding to the marginal (log-)posterior evaluated at mu.
#########
marginal_mu_fct <- function(logpost, map, hessian, n,
                        n_gh = 20, rect_grid = NULL, log = TRUE) {
  f <- function(mu) {
    return(marginal_mu(logpost, map, hessian, mu, 
                       n, n_gh, rect_grid, log))
  }
  return(f)
}

### empirical_quantiles() #####################################################
# empirical_quantiles() returns the quantiles associated with density function
# fmu. 
####
# fmu       should be a vector of size l corresponding to a density evaluated 
#           at mu_tmp.
# mu_tmp    should be a vector of size l or regularly spaced values.
# quantiles should be a vector of probabilities in which to evaluate the 
#           quantile function (its variable name is poorly chosen...)
# M         should be either NULL or a lower-triangular matrix with 1 on its
#           lower half. It can be provided to speed up a little the execution 
#           of the function.
####
# returns a vector of same size than mu_tmp and fmu.
#########
empirical_quantiles <- function(fmu, mu_tmp, quantiles, M = NULL) {
  if (sum(is.na(fmu)) > 0) {
    fmu[is.na(fmu)] <- 0    
  }
  l <- length(mu_tmp)
  if (is.null(M)) {
    M <- diag(x = 1 , nrow = l)
    for (i in 1:l) {
      M[i, 1:i] <- 1
    }
  }
  # Cumulative distribution function based on fmu
  F_mu <- M %*% fmu  / sum(fmu, na.rm = T)
  
  z <- vector("numeric", length(quantiles))
  i <- 0
  for (q in quantiles) {
    i <- i + 1
    if (q < min(F_mu)) {
      z[i] <- min(mu_tmp)
    } else {
      z[i] <- max(mu_tmp[which(F_mu < q)])}
  }
  return(z)
}