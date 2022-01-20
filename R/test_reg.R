#####################################################
### Testing code for two-step expected shortfall ####
#####################################################

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

nseq = seq(1000, 10000, by = 1000)
pseq = floor(nseq / 100)
l = length(nseq)
alpha = 0.05
qr0 = qt(alpha, 2)
integrand = function(x) {qt(x, 2)}
inte = integrate(integrand, lower = 0, upper = alpha)
es0 = inte$value / alpha
M = 20
qr1 = qr2 = qr3 = matrix(0, M, l)
es1 = es2 = es3 = matrix(0, M, l)
time1 = time2 = time3 = matrix(0, M, l)

pb = txtProgressBar(style = 3)
for (j in 1:l) {
  n = nseq[j]
  p = pseq[j]
  h = ((p + log(n)) / n)^(0.4)
  Sigma = toeplitz(0.5^(0:(p - 1)))
  for (i in 1:M) {
    set.seed((j - 1) * M + i)
    X = mvrnorm(n, rep(0, p), Sigma)
    err = rt(n, 2)
    ## Homo 
    beta = runif(p, 0, 2)
    Y = X %*% beta + err
    beta_qr = c(qr0, beta)
    beta_es = c(es0, beta)
    ## Hetero
    #X[, 1] = abs(X[, 1])
    #beta = runif(p - 1, 0, 2)
    #Y = X[, 1] * err + X[, -1] %*% beta + rnorm(n)
    
    start = Sys.time()
    fit1 = esreg(Y ~ X, alpha = alpha)
    end = Sys.time()
    time1[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr1[i, j] = exam(beta_qr, fit1$coefficients_q)
    es1[i, j] = exam(beta_es, fit1$coefficients_e)
    
    start = Sys.time()
    fit2 = twoStep(X, Y, alpha = alpha)
    end = Sys.time()
    time2[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr2[i, j] = exam(beta_qr, fit2$beta)
    es2[i, j] = exam(beta_es, fit2$theta)
    
    start = Sys.time()
    fit3 = twoStepRob(X, Y, alpha = alpha)
    end = Sys.time()
    time3[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr3[i, j] = exam(beta_qr, fit3$beta)
    es3[i, j] = exam(beta_es, fit3$theta)

    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}

rbind(colMeans(time1), colMeans(time2), colMeans(time3))
rbind(colMeans(qr1), colMeans(qr2), colMeans(qr3))
rbind(colMeans(es1), colMeans(es2), colMeans(es3))
