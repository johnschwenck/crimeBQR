# Tests for R

library(fields)
library(RandomFields)
library(mvtnorm)

source('R/crimeBQR_functions.R')

# Test Basic functionality
n <- 100
locs <- cbind(runif(n, 0, 10), runif(n, 0, 10))
m1 <- RMexp(var = 2, scale = 1.5) + 
  RMnugget(var = .1) + 
  RMtrend(mean = 1)
sim_vals <- RFsimulate(m1, x = locs[,1], y = locs[,2])
test_X <- cbind(rnorm(n), rnorm(n), rnorm(n), rnorm(n))
test_beta <- c(1, 2, 3, 4)
test_y <- test_beta[2]*locs[,1] + test_beta[3]*locs[,2] + sim_vals$variable1
test_tau <- .5
test_tau <- c(.025, .05, .5, .95, .975)

get_score(test_y, test_X, test_tau, test_beta)
get_likelihood(test_y, test_X, test_tau, test_beta, locs, 1)



# Test IS for single quantile
m <- 10000
S_0 <- cov(test_X)
initial_draw <- rmvnorm(m, test_beta, S_0)

get_updated_params(test_y, test_X, test_tau, test_beta, cov(test_X), locs, initial_draw)

adaIS_singleQuantile(test_y, test_X, test_tau, 1, locs, m, 22)




# Test IS for multiple quantiles

# Toy example of building joint covariance block diagonal matrix
# The number of S correspond to how many quantile levels (i.e. t_1 = .05, t_2 = .5, t_3 = .9, etc --> S1, S2, S3)
# each S is composed of a p x p matrix of predictor variables (i.e. a 4 x 4 matrix implies 4 predictors)
a <- 3
s1 <- matrix(rnorm(a^2), ncol = a)
s2 <- matrix(rnorm(a^2), ncol = a)
s3 <- matrix(rnorm(a^2), ncol = a)
tmp = list(s1, s2, s3)

tau <- c(.05, .25, .5)
p = ncol(s1)
n = length(tmp)
gamma_tmp = 1 

# mat_tmp = matrix(ncol = n*p, nrow = n*p)
# for(i in seq(1:n)){
#   mat_tmp[ ((i*p) - (p-1)):(p*i), ((i*p) - (p-1)):(p*i)] <- tmp[[i]]
# }
# mat_tmp

# use this
bmat <- bdiag(tmp)

for(i in 1:length(tau)){
  for(j in 1:length(tau)){
    
    if (i < j) { # Only compute upper triangle for speed
    out <- gamma_tmp * (min(tau[i], tau[j]) - tau[i] * tau[j]) * 
        solve( (tau[i] * (1 - tau[i]) * solve(tmp[[i]]) + tau[j] * (1 - tau[j]) * solve(tmp[[j]]) ) / 2)
    
    bmat[ ( (i-1) * p + 1 ) : (i * p), ( (j-1) * p + 1 ) : (j * p) ] <- out
    bmat[ ( (j-1) * p + 1 ) : (j * p), ( (i-1) * p + 1 ) : (i * p) ] <- out
    
    }
    
  }
}

bmat



