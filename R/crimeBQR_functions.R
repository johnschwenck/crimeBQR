# Functions to perform Ada-Is Bayesian Spatial Quantile Regression


# Score Function according to paper
get_score <- function(y, X, tau, beta){  ## outputs px1 vector
  u <- y - X%*%beta
  psi <- ifelse(u < 0, tau - 1, tau)
  return(t(X)%*%psi)
}


# Spatial Covariance Matrix
get_spatial_covar_mat <- function(locs){  ## outputs nxn spatial covariance matrix
  dist <- rdist(locs)
  return(Matern(dist, range = .2, nu = 1.5))
}


# Likelihood Function
get_likelihood <- function(y, X, tau, beta, locs, C){  ## Outputs 1x1 scalar that is the likelihood
  score <- get_score(y, X, tau, beta)
  X_centered <- X - colMeans(X)
  coef <- -1/(2*length(y))
  kernel <- exp(coef*t(score)%*%t(X_centered)%*%solve(get_spatial_covar_mat(locs))%*%X_centered%*%score)
  return(C*kernel[1,1])
}


# Updated Parameters
get_updated_params <- function(y, X, tau, beta, Sigma, locs, draw){  ## outputs list with mu (which is a p-length vector of estimated means) and S (which is the estimated pxp covariance matrix of the features) based on the given draw of values
  n <- length(y)
  p <- ncol(X)
  
  w <- get_likelihood(y, X, tau, beta, locs, 1)*(1/(2*n)^p)/dmvnorm(draw, beta, Sigma)
  
  mu_hat <- apply(w*draw, 2, sum)/sum(w)
  
  S <- matrix(nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      S[i,j] <- sum(w*(draw[,i] - mu_hat[i])*(draw[,j] - mu_hat[j]))/(sum(w)*(n - 1))
    }
  }
  temp <- list(mu_hat, S)
  names(temp) <- c("mu", "S")
  return(temp)
}


# Ada IS Algorithm for a single quantile
adaIS_singleQuantile <- function(y, X, tau, C, locs, M, num_reps){  ## ## outputs list with mu (which is a p-length vector of estimated means) and S (which is the estimated pxp covariance matrix of the features) for the posterior distribution of our betas based on the IS algorithm
  p <- ncol(X)
  n <- length(y)
  beta <- coef(lm(y ~ X))[-1]
  S <- cov(X)
  
  initial_draw <- rmvnorm(M, beta, S)
  
  params <- get_updated_params(y, X, tau, beta, S, locs, initial_draw)
  
  for(i in 1:num_reps){
    draws <- rmvnorm(M, params$mu, params$S)
    params <- get_updated_params(y, X, tau, params$mu, params$S, locs, draws)
  }
  return(params)
}


# 
