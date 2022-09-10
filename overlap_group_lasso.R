## ----setup, include=FALSE--------------------------------------------------------------------
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(eval = T)
knitr::opts_chunk$set(warning = F)
knitr::opts_chunk$set(error = F)
knitr::opts_chunk$set(message = F)
knitr::opts_chunk$set(tidy = 'styler')
knitr::opts_chunk$set(comment = ">")
library(tidyverse)


## --------------------------------------------------------------------------------------------
library(glmnet)
library(gglasso)

# for reproducibility
set.seed(1234)

# support of g is the union of groups 4 and 5 --> (beta_25,...,beta_42)
g <- sapply(0:9, function(i) 1:10+8*i)
g
beta <- rnorm(18 + 1 )
int <- beta[1]
b <-  c(rep(0,24), beta[2:19], rep(0, 40))
suppb <- rep(0, length(b)); suppb[which(b!=0)] <- 1
c(int, b)


## --------------------------------------------------------------------------------------------
Nsim <- 50
p <- 82
n <- 50
# covariates and sigma = |E[Xb+int]| calculated via monte carlo
x <- matrix(rnorm(n*p), nrow = n)
sig <- mean(replicate(1e4, abs(mean(matrix(rnorm(n*p), nrow = n)%*%b + rep(1,n)*int))))
sig


## --------------------------------------------------------------------------------------------
# model matrix and groups for overlapping group lasso
g <- sapply(0:9, function(i) 1:10+8*i)
gseq <- c(g)
xtil <- x[,gseq]
gg <- rep(1:10, each = 10)


## --------------------------------------------------------------------------------------------
# lambda sequence for every model:
# the grid is equally spaced on the log_2 scale -> if values are visualized
# as equally spaced, then the visualization is on the log_2 scale
lgrid <- 2^(seq(-4, 0.2, length = 100))

# simulation lasso n = 50
freq <- matrix(0, nrow = p, ncol = 100)
perccorrect <- numeric(100)
for(i in 1:Nsim){
  sel <- matrix(1, p, 100)
  y <- int + x%*%b + rnorm(n, sd = sig)
  m <-  glmnet(x, y, lambda = lgrid)
  which0 <- apply(m$beta, 2, function(x) which(x==0))
  for(j in 1:100) sel[which0[[j]],j] <- 0
  freq <- freq + sel
  correct <- apply(sel, 2, function(x) all(x==suppb))
  perccorrect <- perccorrect + correct
}
lasso50 <- freq/Nsim
lasso50perc <- perccorrect/Nsim
lasso50perc

# simulation overlap group lasso n = 50
freq <- matrix(0, nrow = p, ncol = 100)
perccorrect <- numeric(100)
for(i in 1:Nsim){
  sel <- matrix(0, p, 100)
  y <- int + x%*%b + rnorm(n, sd = sig)
  m <- gglasso(xtil, y, group = gg, loss = "ls", lambda = lgrid)
  which0 <- apply(m$beta, 2, function(x) which(x==0))
  which_non0 <- lapply(which0, function(x) c(g)[-x])
  for(j in 1:100) sel[which_non0[[j]],j] <- 1
  freq <- freq + sel
  correct <- apply(sel, 2, function(x) all(x==suppb))
  perccorrect <- perccorrect + correct
}
oglasso50 <- freq/Nsim
oglasso50perc <- perccorrect/Nsim
oglasso50perc
length(oglasso50perc[oglasso50perc>.9])/length(oglasso50perc)


n <- 100
x <- matrix(rnorm(n*p), nrow = n)
sig <- abs(mean(replicate(1e4, mean(matrix(rnorm(n*p), nrow = n)%*%b + rep(1,n)*int))))
xtil <- x[,gseq]

# simulation lasso n = 100
freq <- matrix(0, nrow = p, ncol = 100)
perccorrect <- numeric(100)
for(i in 1:Nsim){
  sel <- matrix(1, p, 100)
  y <- int + x%*%b + rnorm(n, sd = sig)
  m <-  glmnet(x, y, lambda = lgrid)
  which0 <- apply(m$beta, 2, function(x) which(x==0))
  for(j in 1:100) sel[which0[[j]],j] <- 0
  freq <- freq + sel
  correct <- apply(sel, 2, function(x) all(x==suppb))
  perccorrect <- perccorrect + correct
}
lasso100 <- freq/Nsim
lasso100perc <- perccorrect/Nsim
lasso100perc

# simulation overlap group lasso n = 100
freq <- matrix(0, nrow = p, ncol = 100)
perccorrect <- numeric(100)
for(i in 1:Nsim){
  sel <- matrix(0, p, 100)
  y <- int + x%*%b + rnorm(n, sd = sig)
  m <- gglasso(xtil, y, group = gg, loss = "ls", lambda = lgrid)
  which0 <- apply(m$beta, 2, function(x) which(x==0))
  which_non0 <- lapply(which0, function(x) c(g)[-x])
  for(j in 1:100) sel[which_non0[[j]],j] <- 1
  freq <- freq + sel
  correct <- apply(sel, 2, function(x) all(x==suppb))
  perccorrect <- perccorrect + correct
}
oglasso100 <- freq/Nsim
oglasso100perc <- perccorrect/Nsim
oglasso100perc
length(oglasso100perc[oglasso100perc>.9])/length(oglasso100perc)
length(oglasso100perc[oglasso100perc>.95])/length(oglasso100perc)


## ---- fig.height=7, fig.width=8, echo=F------------------------------------------------------
par(mfrow = c(2,2))
par(mar = c(5,4.5,4,2)+0.1)

image(t(lasso50[nrow(lasso50):1,ncol(lasso50):1]), col=grey(seq(0, 1, length = 256)),
      xlab = expression(log[2](lambda)), ylab = expression(hat(beta)), axes = F)
axis(2, at = seq(0.01, 1, length = 5), labels = seq(nrow(lasso50)-2, 0, by = -20))
box()

image(t(oglasso50[nrow(oglasso50):1,ncol(oglasso50):1]), col=grey(seq(0, 1, length = 256)),
      xlab = expression(log[2](lambda)), ylab = expression(hat(beta)), axes = F)
axis(2, at = seq(0.01, 1, length = 5), labels = seq(nrow(oglasso50)-2, 0, by = -20))
box()


image(t(lasso100[nrow(lasso100):1,ncol(lasso100):1]), col=grey(seq(0, 1, length = 256)),
      xlab = expression(log[2](lambda)), ylab = expression(hat(beta)), axes = F)
axis(2, at = seq(0.01, 1, length = 5), labels = seq(nrow(lasso100)-2, 0, by = -20))
box()

image(t(oglasso100[nrow(oglasso100):1,ncol(oglasso100):1]), col=grey(seq(0, 1, length = 256)),
      xlab = expression(log[2](lambda)), ylab = expression(hat(beta)), axes = F)
axis(2, at = seq(0.01, 1, length = 5), labels = seq(nrow(oglasso100)-2, 0, by = -20))
box()

par(mfrow = c(1,1))