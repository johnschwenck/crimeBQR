## Import packages
library(fields)
library(RandomFields)
library(mvtnorm)
library(ggplot2)
library(matrixcalc)
library(gstat)
library(GpGp)

## METHODS
## Define basic methods 
augmented_feature_mat <- function(X, locs, knots, phi, sigma_sq){
  C_knots <- sigma_sq*exp(-rdist(knots, knots)/phi)
  half_inv_C_knots <- chol(solve(C_knots))
  basis_X <- sigma_sq*exp(-rdist(locs, knots)/phi)%*%half_inv_C_knots
  return(cbind(X, basis_X)[,-1])
}

get_score<- function(y, X, tau, beta){ 
  n <- nrow(X)
  temp <- t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)[1,1]
  for(i in 2:n){
    temp <- temp + t(X[i,,drop=FALSE])*ifelse(y[i] - X[i,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)[1,1]
  }
  return(temp)
}

weight <- function(X, tau){  
  n <- nrow(X)
  temp <- X[1,]%*%t(X[1,])
  for(i in 2:n){
    temp <- temp + X[i,]%*%t(X[i,])
  }
  coef <- n/(tau*(1 - tau))
  return(coef*solve(temp))
}

get_likelihood <- function(y, X, tau, beta, C){ 
  score <- get_score(y, X, tau, beta)
  coef <- -1/(2*length(y))
  weight_X <- weight(X, tau) 
  kernel <- exp(coef*(t(score)%*%weight_X%*%score)[1,1])
  return(C*kernel)
}

## Define methods used for Importance Sampling
get_updated_params <- function(y, X, tau, beta, Sigma, locs, draws){  
  n <- length(y)
  p <- ncol(X)
  
  is_weights <- numeric(nrow(draws))
  for(i in 1:nrow(draws)){
    is_weights[i] <- get_likelihood(y, X, tau, draws[i,], locs, 1)*(1/(2*n)^p)/dmvnorm(draws[i,], beta, cov(X))
    print(i)
    print(is_weights[i])
  }
  
  mu_hat <- apply(is_weights*draws, 2, sum)/sum(is_weights)
  
  S <- matrix(nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      S[i,j] <- sum(is_weights*draws[,i]*draws[,j])/sum(is_weights) - mu_hat[i]*mu_hat[j]
    }
  }
  temp <- list(mu_hat, S)
  names(temp) <- c("mu", "S")
  return(temp)
}

adaIS_singleQuantile <- function(y, X, tau, C, locs, M, num_reps){ 
  p <- ncol(X)
  n <- length(y)
  beta <- coef(lm(y ~ X))[-1]
  S <- cov(X) 
  
  a_1 <- coef(lm(y ~ X))[-1]
  a_0 <- quantile(y - X%*%a_1, tau)
  a <- c(a_0, a_1)
  X <- cbind(rep(1, n), X)
  temp <- X[1,]%*%t(X[1,])
  for(i in 2:n){
    temp <- temp + X[i,]%*%t(X[i,])
  }
  S_0 <- 1*tau*(1 - tau)*solve((1/n)*temp)
  
  params <- list(a, S_0)
  names(params) <- c("mu", "S")
  draw <- rmvnorm(M, params$mu, params$S)
  
  for(i in 2:num_reps){ 
    draw <- rmvnorm(M, params$mu, params$S)  
    params <- get_updated_params(y, X, tau, params$mu, params$S, locs, draw)
  }
  return(params)
}

## SIMULATION
n <- 10000
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

training_indices <- sample(1:n, round(.75*n), replace = FALSE)

y_train <- y[training_indices]
X_train <- X[training_indices,]
locs_train <- locs[training_indices,]

y_test <- y[-training_indices]
X_test <- X[-training_indices,]
locs_test <- locs[-training_indices,]

training_data <- as.data.frame(cbind(y_train, X_train[,-1], locs_train))
names(training_data) <- c("y_train", "X1_train", "X2_train", "long_train", "lat_train")

ggplot() + geom_point(data = training_data, aes(x = long_train, y = lat_train, col = y_train)) +
  scale_colour_gradient(low = "blue", high = "red") +
  labs(title = "Simulated Response Values",
       x = "Longitude",
       y = "Latitude",
       col = "Value")

## TEST METHODS
## Output is commented out to enhance readability
tau <- .5

X_augmented <- augmented_feature_mat(X_train, locs_train, expand.grid(seq(0.1,9.9,length=5), seq(0.1,9.9,length=5)), 5, 2)
## X_augmented[1:10,]

beta <- coef(lm(y_train ~ X_augmented))[-1]
score <- get_score(y_train, X_augmented, tau, beta)
## score

weight_mat <- weight(X_augmented, tau)
## weight_mat

lik <- get_likelihood(y_train, X_augmented, tau, beta, 1)
## lik

## B is an Mxp matrix consisting of M simulated p-length vectors of model coefficients
M <- 100
Sigma <- cov(X_augmented)
B <- rmvnorm(M, b, Sigma)

## pars <- get_updated_params(y_train, X_augmented, tau, beta, Sigma, locs_train, B)
## pars

## post_dist_params <- adaIS_singleQuantile(y_train, X_augmented, tau, 1, locs_train, 5000, 10)
## get_updated_params() isn't working, so we can't expect adaIS_singleQuantile() to work either

## ISSUES AND ATTEMPTED SOLUTIONS
## Approach 1: Using approach defined in paper with augmented feature matrix
importance_weights_1 <- numeric(nrow(B))
for(i in 1:ncol(B)){
  importance_weights_1[i] <- get_likelihood(y_train, X_augmented, tau, B[i,], 1)*(1/(2*nrow(X_augmented))^2)/dmvnorm(B[i,], b, Sigma)
}
importance_weights_1

## Approach 2: Exponentiate quotient of exponents 
importance_weights_2 <- numeric(nrow(B))
for(i in 1:nrow(B)){
  s <- get_score(y_train, X_augmented, tau, B[i,])
  W <- weight(X_augmented, .5)
  numer <- (-1/(2*nrow(X_augmented)))*t(s)%*%W%*%s
  denom <- -.5*t(B[i,] - colMeans(B))%*%solve(round(Sigma, 7))%*%(B[i,] - colMeans(B))
  importance_weights_2[i] <- exp(numer/denom)
}
importance_weights_2

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
l <- get_likelihood(y_train, X_augmented, tau, beta, 1)
l

## 3A: Approach 1 using feature matrix from Approach 3
importance_weights_3a <- numeric(nrow(B))
for(i in 1:ncol(B)){
  importance_weights_3a[i] <- get_likelihood(y_train, X_augmented, tau, B[i,], 1)*(1/(2*nrow(X_augmented))^2)/dmvnorm(B[i,], beta, cov(X_augmented))
}
importance_weights_3a

## 3B: Approach 2 using feature matrix from Approach 3
importance_weights_3b <- numeric(nrow(B))
for(i in 1:nrow(B)){
  s <- get_score(y_train, X_augmented, .5, B[i,])
  W <- weight(X_augmented, .5)
  numer <- (-1/(2*nrow(X_augmented)))*t(s)%*%W%*%s
  denom <- -.5*t(B[i,] - colMeans(B))%*%solve(round(Sigma, 7))%*%(B[i,] - colMeans(B))
  importance_weights_3b[i] <- exp(numer/denom)
}
importance_weights_3b
