## Import packages
library(mvtnorm)
library(fields)
library(GpGp)

## Define basic methods 

augmented_feature_mat <- function(X, locs, knots, phi, sigma_sq){
  C_knots <- sigma_sq*exp(-rdist(knots, knots)/phi)
  half_inv_C_knots <- chol(solve(C_knots))
  basis_X <- sigma_sq*exp(-rdist(locs, knots)/phi)%*%half_inv_C_knots
  return(cbind(X, basis_X)[,-1])
}

get_score<- function(y, X, tau, beta){  ## See item (2) below
  n <- nrow(X)
  temp <- t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)[1,1]
  for(i in 2:n){
    temp <- temp + t(X[i,,drop=FALSE])*ifelse(y[i] - X[i,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)[1,1]
  }
  return(temp)
}

weight <- function(X, tau){  ## See item (4) below
  n <- nrow(X)
  temp <- X[1,]%*%t(X[1,])
  for(i in 2:n){
    temp <- temp + X[i,]%*%t(X[i,])
  }
  coef <- n/(tau*(1 - tau))
  return(coef*solve(temp))
}

get_likelihood <- function(y, X, tau, beta, C){ ## See item (5) below
  score <- get_score(y, X, tau, beta)
  coef <- -1/(2*length(y))
  weight_X <- weight(X, tau) 
  kernel <- exp(coef*(t(score)%*%weight_X%*%score)[1,1])
  return(C*kernel)
}

## Simulate some data 
n <- 1500
locs <- as.matrix(cbind(runif(n, 0, 10), runif(n, 0, 10)))
sigma_sq <- 2
phi <- 5
nu <- .5
tau_sq <- .5

w <- fast_Gp_sim(c(sigma_sq, phi, nu, 0), "matern_isotropic", locs, 30)

X <- as.matrix(cbind(rep(1, n), locs[,1]/10, locs[,2]/10))
beta <- as.matrix(c(1, 2, 3))
p <- length(beta)
y <- rnorm(n, X%*%beta + w, sqrt(tau_sq))

X_augmented <- augmented_feature_mat(X, locs, expand.grid(seq(0.1,9.9,length=5),seq(0.1,9.9,length=5)), phi, sigma_sq)

b <- coef(lm(y ~ X_augmented))[-1]

## B is an Mxp matrix consisting of M simulated p-length vecotrs of model coefficients
M <- 100
B <- rmvnorm(M, b, cov(X_augmented))

## Approach 1: Using approach defined in paper with augmented feature matrix
l <- get_likelihood(y, X_augmented, .5, b, 1)
l

importance_weights <- numeric(nrow(B))
for(i in 1:ncol(B)){
  importance_weights[i] <- get_likelihood(y, X_augmented, .5, B[i,], 1)*(1/(2*nrow(X))^2)/dmvnorm(B[i,], b, cov(X_augmented))
}
importance_weights

## Approach 2: Exponentiate quotient of exponents 
importanct_weights <- numeric(nrow(B))
for(i in 1:nrow(B)){
  s <- get_score(y, X_augmented, .5, B[i,])
  W <- weight(X_augmented, .5)
  numer <- (-1/(2*nrow(B)))*t(s)%*%W%*%s
  denom <- -.5*t(B[i,] - colMeans(B))%*%solve(round(cov(X_augmented), 7))%*%(B[i,] - colMeans(B))
}
importance_weights

## Approach 3: Scaling basis functions by their L^2 norms 
augmented_feature_mat <- function(X, locs, knots, phi, sigma_sq){
  C_knots <- sigma_sq*exp(-rdist(knots, knots)/phi)
  half_inv_C_knots <- chol(solve(C_knots))
  basis_X <- sigma_sq*exp(-rdist(locs, knots)/phi)%*%half_inv_C_knots
  for(i in 1:ncol(basis_X)){
    basis_X[,i] <- basis_X[,i]*sqrt(sum(basis_X[,i]^2)) 
  }
  return(cbind(X, basis_X)[,-1])
}
l <- get_likelihood(y, X_augmented, .5, b, 1)
l

## 3A: Approach 1 using feature matrix from Approach 3
importance_weights <- numeric(nrow(B))
for(i in 1:ncol(B)){
  importance_weights[i] <- get_likelihood(y, X_augmented, .5, B[i,], 1)*(1/(2*nrow(X))^2)/dmvnorm(B[i,], b, cov(X_augmented))
}
importance_weights

## 3B: Approach 2 using feature matrix from Approach 3
importanct_weights <- numeric(nrow(B))
for(i in 1:nrow(B)){
  s <- get_score(y, X_augmented, .5, B[i,])
  W <- weight(X_augmented, .5)
  numer <- (-1/(2*nrow(B)))*t(s)%*%W%*%s
  denom <- -.5*t(B[i,] - colMeans(B))%*%solve(round(cov(X_augmented), 7))%*%(B[i,] - colMeans(B))
}
importance_weights
