#############################################
### Testing code for regularized conquer ####
#############################################

rm(list = ls())
Rcpp::sourceCpp("~/Documents/Software/conquerReg/src/conquerReg.cpp")

library(quantreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(caret)
library(rqPen)
library(tikzDevice)
library(ggplot2)

exam = function(betaMat, betaEst) {
  p = length(betaMat)
  TPR = sum(betaMat != 0 & betaEst != 0) / sum(betaMat != 0)
  FPR = sum(betaMat == 0 & betaEst != 0) / sum(betaMat == 0)
  FDR = 0
  if (sum(betaEst != 0) > 0) {
    FDR = sum(betaMat == 0 & betaEst != 0) / sum(betaEst != 0)
  }
  err = norm(betaMat - betaEst, "2")
  return (list("TPR" = TPR, "FPR" = FPR, "FDR" = FDR, "err" = err))
}


n = 500
p = 1000
s = 20
kfolds = 5
nsim = 200
tau = 0.5
beta0 = qt(tau, 2)
h = 0.5 * (log(p) / n)^(1/4)
Sigma = toeplitz(0.5^(0:(p - 1)))
M = 1
time = TPR = FPR = FDR = error = matrix(0, 4, M)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(2021 + m)
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, 2)
  ## Homo
  beta = c(runif(s, 1, 1.5), rep(0, p - s))
  betaMat = beta
  Y = X %*% beta + err
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = c(runif(s - 1, 1, 1.5), rep(0, p - s))
  #betaMat = c(beta0, beta)
  #Y = X[, 1] * err + X[, -1] %*% beta
  
  folds = createFolds(Y, kfolds, FALSE)
  U = matrix(runif(nsim * n), nsim, n)
  pivot = tau - (U <= tau)
  lambda0 = quantile(rowMaxs(abs(pivot %*% scale(X))), 0.9) / n
  lambdaSeq = seq(0.5, 1.5, length.out = 50) * lambda0
  
  ### lasso
  start = Sys.time()
  fit.lasso = cvGaussLasso(X, Y, lambdaSeq, folds, tau, kfolds, h)
  end = Sys.time()
  time[1, m] = as.numeric(difftime(end, start, units = "secs"))
  beta.lasso = as.numeric(fit.lasso$coeff)
  result = exam(betaMat, beta.lasso[-1])
  TPR[1, m] = result$TPR
  FPR[1, m] = result$FPR
  FDR[1, m] = result$FDR
  error[1, m] = result$err
  
  ### elastic net
  start = Sys.time()
  fit.elastic = cvGaussElastic(X, Y, lambdaSeq, folds, tau, alpha = 0.5, kfolds, h)
  end = Sys.time()
  time[2, m] = as.numeric(difftime(end, start, units = "secs"))
  beta.elastic = as.numeric(fit.elastic$coeff)
  result = exam(betaMat, beta.elastic[-1])
  TPR[2, m] = result$TPR
  FPR[2, m] = result$FPR
  FDR[2, m] = result$FDR
  error[2, m] = result$err
  
  ### group lasso
  beta = c(rep(1.3, 5), rep(1.5, 5), rep(1.2, 5), rep(1.4, 5), rep(0, p - s))
  betaMat = beta
  Sigma = matrix(0, p, p)
  Sigma[1:5, 1:5] = toeplitz(0.5^(0:(5 - 1)))
  Sigma[6:10, 6:10] = toeplitz(0.5^(0:(5 - 1)))
  Sigma[11:15, 11:15] = toeplitz(0.5^(0:(5 - 1)))
  Sigma[16:20, 16:20] = toeplitz(0.5^(0:(5 - 1)))
  Sigma[21:p, 21:p] = toeplitz(0.5^(0:(p - s - 1)))
  X = mvrnorm(n, rep(0, p), Sigma)
  Y = X %*% beta + err
  group = c(0, rep(0, 5), rep(1, 5), rep(2, 5), rep(3, 5), rep(4, p - s))
  start = Sys.time()
  fit.group.lasso = cvGaussGroupLasso(X, Y, lambdaSeq, folds, tau, kfolds, group = group, G = 5, h)
  end = Sys.time()
  time[3, m] = as.numeric(difftime(end, start, units = "secs"))
  beta.group.lasso = as.numeric(fit.group.lasso$coeff)
  result = exam(betaMat, beta.group.lasso[-1])
  TPR[3, m] = result$TPR
  FPR[3, m] = result$FPR
  FDR[3, m] = result$FDR
  error[3, m] = result$err
  
  ### sparse group lasso
  beta[c(1, 6, 11, 16)] = 0
  betaMat = beta
  Y = X %*% beta + err
  start = Sys.time()
  fit.sparse.group.lasso = cvGaussSparseGroupLasso(X, Y, lambdaSeq, folds, tau, kfolds, group = group, G = 5, h)
  end = Sys.time()
  time[4, m] = as.numeric(difftime(end, start, units = "secs"))
  beta.sparse.group.lasso = as.numeric(fit.sparse.group.lasso$coeff)
  result = exam(betaMat, beta.sparse.group.lasso[-1])
  TPR[4, m] = result$TPR
  FPR[4, m] = result$FPR
  FDR[4, m] = result$FDR
  error[4, m] = result$err
  
  setTxtProgressBar(pb, m / M)
}

time
TPR
FPR
FDR
error



















### fused lasso
p = 50
n = 400
kfolds = 5
nsim = 200
tau = 0.5
beta0 = qt(tau, 2)
h = 0.5 * (log(p) / n)^(1/4)
Sigma = toeplitz(0.5^(0:(p - 1)))
X = mvrnorm(n, rep(0, p), Sigma)
err = rt(n, 2)
beta = c(rep(1.5, 15), rep(2, 15), rep(-1.5, 20))
betaMat = beta
Y = X %*% beta + err

folds = createFolds(Y, kfolds, FALSE)

U = matrix(runif(nsim * n), nsim, n)
pivot = tau - (U <= tau)
lambda0 = quantile(rowMaxs(abs(pivot %*% scale(X))), 0.9) / n
lambdaSeq = seq(0.5, 1.5, length.out = 50) * lambda0

start = Sys.time()
fit.fused.lasso = cvGaussFusedLasso(X, Y, lambdaSeq, folds, tau, kfolds, h)
end = Sys.time()
as.numeric(difftime(end, start, units = "secs"))
beta.fused.lasso = as.numeric(fit.fused.lasso$coeff)
exam(betaMat, beta.fused.lasso[-1])
rbind(betaMat, beta.fused.lasso[-1])
