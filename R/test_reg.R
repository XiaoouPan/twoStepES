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
X = mvrnorm(n, rep(0, p), Sigma)
err = rt(n, 2)
beta = runif(p, 0, 2)
Y = X %*% beta + err
alpha = 0.1
beta_qr = c(qt(alpha, 2), beta)
integrand = function(x) {qt(x, 2)}
inte = integrate(integrand, lower = 0, upper = alpha)
beta_es = c(inte$value / alpha, beta)


fit1 = esreg(Y ~ X, alpha = alpha)
fit2 = twoStepNonstd(X, Y, alpha = alpha)
fit3 = twoStep(X, Y, alpha = alpha)
fit4 = twoStepRob(X, Y, alpha = alpha)

c(exam(beta_qr, fit1$coefficients_q), exam(beta_es, fit1$coefficients_e))
c(exam(beta_qr, fit2$beta), exam(beta_es, fit2$theta))
c(exam(beta_qr, fit3$beta), exam(beta_es, fit3$theta))
c(exam(beta_qr, fit4$beta), exam(beta_es, fit4$theta))


nseq = seq(2000, 10000, by = 2000)
pseq = floor(nseq / 100)
l = length(nseq)
tauSeq = seq(0.2, 0.8, by = 0.05)
nTau = length(tauSeq)
beta0 = qt(tauSeq, 2)
M = 50
coef1 = coef2 = coef3 = coef4 = coef5 = coef6 = coef7 = matrix(0, M, l)
eff1 = eff2 = eff3 = eff4 = eff5 = eff6 = eff7 = matrix(0, M, l)
time1 = time2 = time3 = time4 = time5 = time6 = time7 = matrix(0, M, l)

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
    #beta = runif(p, 0, 2)
    #betaMat = rbind(beta0, matrix(beta, p, nTau))
    #Y = X %*% beta + err
    #indEff = 1
    ## Hetero
    X[, 1] = abs(X[, 1])
    beta = runif(p - 1, 0, 2)
    betaMat = rbind(qnorm(tauSeq), beta0, matrix(beta, p - 1, nTau))
    Y = X[, 1] * err + X[, -1] %*% beta + rnorm(n)
    indEff = 2
    
    start = Sys.time()
    list = rq(Y ~ X, tau = tauSeq, method = "fn")
    end = Sys.time()
    time2[i, j] = as.numeric(difftime(end, start, units = "secs"))
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
