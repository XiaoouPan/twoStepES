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
n = 4000
alpha = 0.05
p = n * alpha / 40
Sigma = toeplitz(0.5^(0:(p - 1)))
qr0 = qnorm(alpha)
integrand = function(x) {qnorm(x)}
#df = 4
#qr0 = qt(alpha, df)
#integrand = function(x) {qt(x, df)}
inte = integrate(integrand, lower = 0, upper = alpha)
es0 = inte$value / alpha
M = 20
coef1 = coef2 = time = matrix(0, M, 4)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  X = mvrnorm(n, rep(0, p), Sigma)
  #err = rt(n, df)
  err = rnorm(n)
  ## Homo 
  beta = runif(p, 0, 2)
  Y = X %*% beta + err
  beta_qr = c(qr0, beta)
  beta_es = c(es0, beta)
  ## Hetero
  #gamma = runif(p + 1, 1, 2)
  #eta = runif(p + 1, 1, 2)
  #Y = cbind(1, X) %*% gamma + (cbind(1, X) %*% eta) * err
  #beta_qr = gamma + eta * qr0
  #beta_es = gamma + eta * es0
  
  start = Sys.time()
  fit1 = esreg(Y ~ X, alpha = alpha)
  end = Sys.time()
  time[i, 1] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 1] = exam(beta_qr, fit1$coefficients_q)
  coef2[i, 1] = exam(beta_es, fit1$coefficients_e)
  
  start = Sys.time()
  fit2 = twoStepNonstd(X, Y, alpha = alpha)
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
  fit4 = twoStepRob(X, Y, alpha = alpha)
  end = Sys.time()
  time[i, 4] = as.numeric(difftime(end, start, units = "secs"))
  coef1[i, 4] = exam(beta_qr, fit4$beta)
  coef2[i, 4] = exam(beta_es, fit4$theta)
  
  setTxtProgressBar(pb, i / M)
}

rbind(colMeans(time), colMeans(coef1), colMeans(coef2))

