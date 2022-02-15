### test code: fixed scale

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

## ab example with fixed scale
n = 6000
alpha = 0.05
lambda = 0.5
p = n * alpha / 50
Sigma = toeplitz(0.5^(0:(p - 1)))
qr_norm = qnorm(alpha)
integrand = function(x) {qnorm(x)}
inte = integrate(integrand, lower = 0, upper = alpha)
es_norm = inte$value / alpha
df = 3
qr_t = qt(alpha, df)
integrand = function(x) {qt(x, df)}
inte = integrate(integrand, lower = 0, upper = alpha)
es_t = inte$value / alpha
M = 20
coef1 = coef2 = time = matrix(0, M, 6)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  X = mvrnorm(n, rep(0, p), Sigma)
  ## Hetero
  #err = rnorm(n)
  err = rt(n, df)
  gamma = runif(p + 1, 0, 2)
  eta = runif(1, 0, 2)
  X[, 1] = abs(X[, 1])
  Y = cbind(1, X) %*% gamma + (X[, 1] * eta) * err
  beta_qr = gamma 
  beta_qr[2] = beta_qr[2] + eta * qr_t
  beta_es = gamma 
  beta_es[2] = beta_es[2] + eta * es_t
  
  #start = Sys.time()
  #fit1 = esreg(Y ~ X, alpha = alpha)
  #end = Sys.time()
  #time[i, 1] = as.numeric(difftime(end, start, units = "secs"))
  #coef1[i, 1] = exam(beta_qr, fit1$coefficients_q)
  #coef2[i, 1] = exam(beta_es, fit1$coefficients_e)
  
  start = Sys.time()
  fit2 = oracle(X, Y, beta_qr, alpha = alpha)
  end = Sys.time()
  time[i, 2] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 2] = exam(beta_qr, fit2$beta)
  coef2[i, 2] = exam(beta_es, fit2$theta)
  
  start = Sys.time()
  fit3 = twoStep(X, Y, alpha = alpha)
  end = Sys.time()
  time[i, 3] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 3] = exam(beta_qr, fit3$beta)
  coef2[i, 3] = exam(beta_es, fit3$theta)
  
  start = Sys.time()
  fit4 = twoStepLambda(X, Y, lambda = lambda, alpha = alpha)
  end = Sys.time()
  time[i, 4] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 4] = exam(beta_qr, fit3$beta)
  coef2[i, 4] = exam(beta_es, fit3$theta)
  
  start = Sys.time()
  fit5 = twoStepRob(X, Y, alpha = alpha)
  end = Sys.time()
  time[i, 5] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 5] = exam(beta_qr, fit4$beta)
  coef2[i, 5] = exam(beta_es, fit4$theta)
  
  start = Sys.time()
  fit6 = twoStepRobLambda(X, Y, lambda = lambda, alpha = alpha)
  end = Sys.time()
  time[i, 6] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 6] = exam(beta_qr, fit4$beta)
  coef2[i, 6] = exam(beta_es, fit4$theta)
  
  setTxtProgressBar(pb, i / M)
}

rbind(colMeans(time), colMeans(coef1), colMeans(coef2))

