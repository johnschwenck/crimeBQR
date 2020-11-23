paper_score <- function(y, X, tau, beta){
  n <- nrow(X)
  temp <- t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)
  for(i in 2:n){
    temp <- temp + t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)
  }
  return(temp)
}

n <- nrow(X)
temp <- t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)
for(i in 2:n){
  temp <- temp + t(X[1,,drop=FALSE])*ifelse(y[1] - X[1,,drop=FALSE]%*%as.matrix(beta) < 0, tau  - 1, tau)
}

paper_score(y_train, X_train, .5, beta = c(1, 1, 1, 1))

ifelse(y_train[1] - X_train[1,,drop=FALSE]%*%as.matrix(beta) < 0, .5  - 1, .5)
X_train[1,,drop=FALSE]%*%as.matrix(c(1, 1, 1, 1))
