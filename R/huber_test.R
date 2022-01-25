rm(list = ls())
Rcpp::sourceCpp("~/Documents/Software/twoStepES/src/es.cpp")

library(esreg)
library(MASS)
library(MultiRNG)
library(matrixStats)
library(tikzDevice)
library(ggplot2)

exam = function(beta, betaHat) {
  return (sqrt(mean((betaHat - beta)^2)))
}

test_ols = function(Z, Y) {
  n = nrow(Z)
  p = ncol(Z) - 1
  gradOld = rep(0, p + 1)
  gradNew = rep(0, p + 1)
  return (as.numeric(l2Reg(Z, Y, gradOld, gradNew, 1 / n)))
}

test_Huber = function(Z, Y) {
  n = nrow(Z)
  p = ncol(Z) - 1
  gradOld = rep(0, p + 1)
  gradNew = rep(0, p + 1)
  der = rep(0, n)
  return (as.numeric(huberReg(Z, Y, der, gradOld, gradNew, n, 1 / n)))
}

n = 4000
p = 20
Sigma = toeplitz(0.5^(0:(p - 1)))
X = mvrnorm(n, rep(0, p), Sigma)
Z = cbind(1, scale(X))
beta = runif(p + 1, 0, 2)
Y = Z %*% beta + rnorm(n)

beta1 = test_ols(Z, Y)
beta2 = test_Huber(Z, Y)
rbind(beta1, beta2)
c(exam(beta, beta1), exam(beta, beta2))

