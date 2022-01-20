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
#qr0 = qt(alpha, 2)
#integrand = function(x) {qt(x, 2)}
#inte = integrate(integrand, lower = 0, upper = alpha)
#es0 = inte$value / alpha
qr0 = qnorm(alpha)
integrand = function(x) {qnorm(x)}
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
    #err = rt(n, 2)
    err = rnorm(n)
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

write.csv(rbind(time1, time2, time3), "~/Dropbox/ES/Simulation/time_normal.csv")
write.csv(rbind(qr1, qr2, qr3), "~/Dropbox/ES/Simulation/qr_normal.csv")
write.csv(rbind(es1, es2, es3), "~/Dropbox/ES/Simulation/es_normal.csv")

### Estimation error: quantile 
mean1 = colMeans(qr1, na.rm = TRUE)
mean2 = colMeans(qr2, na.rm = TRUE)
mean3 = colMeans(qr3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed Huberized method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed Huberized method}"))

setwd("~/Dropbox/ES")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #scale_linetype_manual(values = c("dashed", "twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error of quantile") +
  #theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  theme(legend.position = c(0.6, 0.75), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


### Estimation error: ES 
mean1 = colMeans(es1, na.rm = TRUE)
mean2 = colMeans(es2, na.rm = TRUE)
mean3 = colMeans(es3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed Huberized method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed Huberized method}"))

setwd("~/Dropbox/ES")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #scale_linetype_manual(values = c("dashed", "twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Estimation error of ES") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  #theme(legend.position = c(0.6, 0.8), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
  #      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
  #      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


## Time plot
mean1 = colMeans(time1, na.rm = TRUE)
mean2 = colMeans(time2, na.rm = TRUE)
mean3 = colMeans(time3, na.rm = TRUE)
dat = rbind(cbind(nseq, mean1), cbind(nseq, mean2), cbind(nseq, mean3))
dat = as.data.frame(dat)
colnames(dat) = c("size", "coef")
dat$type = c(rep("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", l), rep("\\texttt{Proposed method}", l), rep("\\texttt{Proposed Huberized method}", l))
dat$type = factor(dat$type, levels = c("\\texttt{Dimitriadis} \\& \\texttt{Bayer}", "\\texttt{Proposed method}", "\\texttt{Proposed Huberized method}"))

setwd("~/Dropbox/ES")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = size, y = coef)) +
  geom_line(aes(y = coef, color = type, linetype = type), size = 3) + 
  #scale_linetype_manual(values = c("dashed", "twodash", "solid")) +
  #geom_ribbon(aes(y = coef, ymin = low, ymax = upp, fill = type), alpha = 0.3)
  theme_bw() + xlab("Sample size") + ylab("Elapsed time (in seconds)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
#theme(legend.position = c(0.6, 0.8), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"),
#      legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
#      axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

