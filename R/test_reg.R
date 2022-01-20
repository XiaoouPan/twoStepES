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

## a toy example
n = 400
p = 50
Sigma = toeplitz(0.5^(0:(p - 1)))
alpha = 0.05
M = 50
coef1 = coef2 = time = matrix(0, M, 4)

pb = txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  X = mvrnorm(n, rep(0, p), Sigma)
  err = rt(n, 2)
  ## Homo 
  beta = runif(p, 0, 2)
  Y = X %*% beta + err
  beta_qr = c(qt(alpha, 2), beta)
  integrand = function(x) {qt(x, 2)}
  inte = integrate(integrand, lower = 0, upper = alpha)
  beta_es = c(inte$value / alpha, beta)
  ## Hetero
  #X[, 1] = abs(X[, 1])
  #beta = runif(p - 1, 0, 2)
  #Y = X[, 1] * err + X[, -1] %*% beta + rnorm(n)
  
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


nseq = seq(500, 2000, by = 500)
pseq = floor(nseq / 50)
l = length(nseq)
alpha = 0.1
M = 50
coef1 = coef2 = coef3 = coef4 = matrix(0, M, l)

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
    beta_qr = c(qt(alpha, 2), beta)
    integrand = function(x) {qt(x, 2)}
    inte = integrate(integrand, lower = 0, upper = alpha)
    beta_es = c(inte$value / alpha, beta)
    ## Hetero
    #X[, 1] = abs(X[, 1])
    #beta = runif(p - 1, 0, 2)
    #Y = X[, 1] * err + X[, -1] %*% beta + rnorm(n)
    
    start = Sys.time()
    fit1 = esreg(Y ~ X, alpha = alpha)
    end = Sys.time()
    time[i, j] = as.numeric(difftime(end, start, units = "secs"))
    coef2[i, j] = mean(sqrt(colSums((list$coefficients - betaMat)^2)))
    eff2[i, j] = mean(abs(list$coefficients[indEff, ] - betaMat[indEff, ]))
    
    start = Sys.time()
    list = smqrGaussProc(X, Y, tauSeq)
    end = Sys.time()
    time3[i, j] = as.numeric(difftime(end, start, units = "secs"))
    coef3[i, j] = mean(sqrt(colSums((list$coeff - betaMat)^2)))
    eff3[i, j] = mean(abs(list$coeff[indEff, ] - betaMat[indEff, ]))
    

    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}
