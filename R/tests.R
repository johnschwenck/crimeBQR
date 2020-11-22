# Tests for R
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

get_score(test_y, test_X, test_tau, test_beta)
get_likelihood(test_y, test_X, test_tau, test_beta, locs, 1)



# Test IS for single quantile
m <- 10000
S_0 <- cov(test_X)
initial_draw <- rmvnorm(m, test_beta, S_0)

get_updated_params(test_y, test_X, test_tau, test_beta, cov(test_X), locs, initial_draw)

adaIS_singleQuantile(test_y, test_X, test_tau, 1, locs, m, 22)