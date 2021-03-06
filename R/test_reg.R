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

nseq = seq(4000, 20000, by = 2000)
alpha = 0.05
pseq = floor(nseq * alpha / 50)
l = length(nseq)
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
lambda = 0.2
qr1 = qr2 = qr3 = qr4 = qr5 = matrix(0, M, l)
es1 = es2 = es3 = es0 = es4 = es5 = matrix(0, M, l)
time1 = time2 = time3 = time4 = time5 = matrix(0, M, l)

pb = txtProgressBar(style = 3)
for (j in 1:l) {
  n = nseq[j]
  p = pseq[j]
  Sigma = toeplitz(0.5^(0:(p - 1)))
  for (i in 1:M) {
    set.seed((j - 1) * M + i)
    X = mvrnorm(n, rep(0, p), Sigma)
    ## Hetero
    #err = rnorm(n)
    err = rt(n, df)
    gamma = runif(p + 1, 0, 2)
    eta = runif(1, 0, 2)
    X = abs(X)
    Y = cbind(1, X) %*% gamma + (X[, 1] * eta) * err
    beta_qr = gamma
    beta_qr[2] = beta_qr[2] + eta * qr_t
    beta_es = gamma
    beta_es[2] = beta_es[2] + eta * es_t
    
    fit0 = oracle(X, Y, beta_qr, alpha = alpha)
    es0[i, j] = exam(beta_es, fit0$theta)
    
    #start = Sys.time()
    #fit1 = esreg(Y ~ X, alpha = alpha)
    #end = Sys.time()
    #time1[i, j] = as.numeric(difftime(end, start, units = "secs"))
    #qr1[i, j] = exam(beta_qr, fit1$coefficients_q)
    #es1[i, j] = exam(beta_es, fit1$coefficients_e)
    
    start = Sys.time()
    fit2 = twoStep(X, Y, alpha = alpha)
    end = Sys.time()
    time2[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr2[i, j] = exam(beta_qr, fit2$beta)
    es2[i, j] = exam(beta_es, fit2$theta)
    
    start = Sys.time()
    fit3 = twoStepLambda(X, Y, lambda = lambda, alpha = alpha)
    end = Sys.time()
    time3[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr3[i, j] = exam(beta_qr, fit3$beta)
    es3[i, j] = exam(beta_es, fit3$theta)
    
    start = Sys.time()
    fit4 = twoStepRob(X, Y, alpha = alpha)
    end = Sys.time()
    time4[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr4[i, j] = exam(beta_qr, fit4$beta)
    es4[i, j] = exam(beta_es, fit4$theta)
    
    start = Sys.time()
    fit5 = twoStepRobLambda(X, Y, lambda = lambda, alpha = alpha)
    end = Sys.time()
    time5[i, j] = as.numeric(difftime(end, start, units = "secs"))
    qr5[i, j] = exam(beta_qr, fit5$beta)
    es5[i, j] = exam(beta_es, fit5$theta)

    setTxtProgressBar(pb, ((j - 1) * M + i) / (l * M))
  }
}

#write.csv(rbind(time1, time2, time3), "~/Dropbox/SQR/ES/Simulation/time_t4.csv")
#write.csv(rbind(qr1, qr2, qr3), "~/Dropbox/SQR/ES/Simulation/qr_t4.csv")
#write.csv(rbind(es1, es2, es3, es0), "~/Dropbox/SQR/ES/Simulation/es_t4.csv")

### Estimation error: QR 
mean1 = colMeans(qr1, na.rm = TRUE)
mean2 = colMeans(qr2, na.rm = TRUE)
mean3 = colMeans(qr3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), 
             rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed robust method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed robust method}"))

ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error of QR") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.65, 0.82), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))



data = as.matrix(read.csv("~/Dropbox/SQR/ES/Simulation/es_t3.csv"))[, -1]
es1 = data[1:50, ]
es2 = data[51:100, ]
es3 = data[101:150, ]
es0 = data[151:200, ]

### Estimation error: ES 
mean0 = colMeans(es0, na.rm = TRUE)
#mean1 = colMeans(es1, na.rm = TRUE)
mean2 = colMeans(es2, na.rm = TRUE)
mean3 = colMeans(es3, na.rm = TRUE)
mean4 = colMeans(es4, na.rm = TRUE)
mean5 = colMeans(es5, na.rm = TRUE)
#dat = rbind(cbind(nseq, mean0), cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3), cbind(nseq, mean4), cbind(nseq, mean5))
dat = rbind(cbind(nseq, mean0), cbind(nseq, mean2), cbind(nseq, mean3), cbind(nseq, mean4), cbind(nseq, mean5))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
#dat$type = c(rep("\\texttt{Oracle}", l), rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), 
#             rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed penalized method}", l), 
#             rep("\\texttt{Proposed robust method}", l), rep("\\texttt{Proposed penalized robust method}", l))
#dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed penalized method}", 
#                                       "\\texttt{Proposed robust method}", "\\texttt{Proposed penalized robust method}", "\\texttt{Oracle}"))
dat$type = c(rep("\\texttt{Oracle}", l), rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed penalized method}", l), 
             rep("\\texttt{Proposed robust method}", l), rep("\\texttt{Proposed penalized robust method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Proposed method}", "\\texttt{Proposed penalized method}", 
                                       "\\texttt{Proposed robust method}", "\\texttt{Proposed penalized robust method}", "\\texttt{Oracle}"))


setwd("~/Dropbox/SQR/ES/Simulation")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error of ES") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.63, 0.80), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

### Standard deviation: ES 
mean0 = colSds(es0, na.rm = TRUE)
mean1 = colSds(es1, na.rm = TRUE)
mean2 = colSds(es2, na.rm = TRUE)
mean3 = colSds(es3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean0), cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Oracle}", l), rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), 
             rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed robust method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed robust method}", "\\texttt{Oracle}"))

setwd("~/Dropbox/SQR/ES/Simulation")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Standard deviation of ES") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.63, 0.82), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


## Time plot

data = as.matrix(read.csv("~/Dropbox/SQR/ES/Simulation/time_t3.csv"))[, -1]
time1 = data[1:50, ]
time2 = data[51:100, ]
time3 = data[101:150, ]

mean1 = colMeans(time1, na.rm = TRUE)
mean2 = colMeans(time2, na.rm = TRUE)
mean3 = colMeans(time3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed robust method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed robust method}"))

setwd("~/Dropbox/ES")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #scale_linetype_manual(values = c("dashed", "twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Elapsed time (in seconds)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.3, 0.8), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
  #    legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #    axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

