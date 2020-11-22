# Functions to perform Ada-Is Bayesian Spatial Quantile Regression

library(fields)
library(RandomFields)
library(mvtnorm)
library(Matrix)

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
get_updated_params <- function(y, X, tau, beta, Sigma, locs, draw){  
  ## outputs list with mu (which is a p-length vector of estimated means) and S (which is the estimated pxp covariance matrix of the features) based on the given draw of values
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


##############################################################################


# Ada IS Algorithm for multiple quantiles
# tau: an n x 1 vector of tau values
adaIS_multipleQuantile <- function(y, X, tau, C, locs, M, num_reps, gam){
  
  if(!is.vector(tau)){
    stop('tau must be an n x 1 length vector')
  }
  
  tau <- as.vector(unlist(tau)) # in case input is a list
  
  #tm <- length(tau)
  p <- ncol(X)
  mu_list <- list()
  S_list <- list()
  
  # Step 1:
  for(i in 1:m){
    #print(adaIS_singleQuantile(y, X, tau[i], C, locs, M, num_reps)) # for diagnosis purposes only
    step1 <- adaIS_singleQuantile(y, X, tau[i], C, locs, M, num_reps)
    mu_list[i] <- step1$mu
    S_list[i] <- step1$S
  }
  
  # Step 2:
  # Matrix of mu vectors combined into one matrix
  mu_mat <- matrix( unlist(mu_list) , ncol = p)
  
  # Block diagonal matrix with S_list elements along the diagonal
  S_mat <- bdiag(S_list)
  
  for(i in 1:length(tau)){
    for(j in 1:length(tau)){
      
      if (i < j) { # Only compute upper triangle for speed
        out <- gam * (min(tau[i], tau[j]) - tau[i] * tau[j]) * 
          solve( (tau[i] * (1 - tau[i]) * solve(tmp[[i]]) + tau[j] * (1 - tau[j]) * solve(tmp[[j]]) ) / 2)
        
        S_mat[ ( (i-1) * p + 1 ) : (i * p), ( (j-1) * p + 1 ) : (j * p) ] <- out
        S_mat[ ( (j-1) * p + 1 ) : (j * p), ( (i-1) * p + 1 ) : (i * p) ] <- out
        
      }
      
    }
  }
  
  out <- list(mu_mat, as.matrix(S_mat))
  names(out) <- c("Mu_mat", "S_mat")
  
  # Step 3: 
  
  
  
  # Step 4:
  
  
  
  return(out)
}
