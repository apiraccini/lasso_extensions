## --------------------------------------------------------------------------------------------
library(glinternet)

# define simulation parameters
Nsim <- 100
n <- 800
p <- 500

# for reproducibility
set.seed(24)


## --------------------------------------------------------------------------------------------
# define model matrix
snr <- 1
xmat <- data.frame(replicate(p, factor(sample(0:2, n, replace = T))), stringsAsFactors = T)
names(xmat) <- paste0("v",1:p,"_")
xxmat <- apply(xmat, 2, as.numeric)

# define true effects indices and coefficients
truemat <- xmat[,1:10]
df <- model.matrix(~ -1 + v1_ + v2_ + v3_ + v4_ + v5_ + v6_ + v7_ + v8_ + v9_ + v10_ +
                     v1_:v2_ + v3_:v4_ + v5_:v6_ + v7_:v8_ + v9_:v10_ +
                     v1_:v3_ + v2_:v4_ + v5_:v7_ + v6_:v8_ + v7_:v10_, data = truemat,
                   contrasts.arg = lapply(truemat, contrasts, contrasts=FALSE))
mainmat <- matrix(1:10, nrow = 10)
intermat <- matrix(c(1,2,3,4,5,6,7,8,9,10,1,3,2,4,5,7,6,8,7,10),
                   nrow = 10, byrow = T)
mainmat; intermat
betas <- rnorm(ncol(df))


## --------------------------------------------------------------------------------------------
get_fdr <- function(){
  #cat("Sto simulando...\n")
  y <- df%*%betas + rnorm(n, 0, sd = sqrt(var(df%*%betas)/snr))
  m <- glinternet(xxmat, y, numLevels = rep(3,p), numToFind = 10)
  ninter <- sapply(m$activeSet, function(x) NROW(x$catcat))
  ind <- max(which(sapply(ninter, function(i) i <=10)==T))
  found.main <- m$activeSet[[ind]]$cat
  found.int <- m$activeSet[[ind]]$catcat
  match.main <- duplicated(rbind(mainmat, found.main))[-(1:nrow(mainmat))]
  match.int <- duplicated(rbind(intermat, found.int))[-(1:nrow(intermat))]
  match <- c(match.main, match.int)
  fdr <- 1-mean(match)
  t(c(NROW(m$activeSet[[ind]]$catcat), fdr))
}

# obtain results
ris <- replicate(Nsim, get_fdr())


## --------------------------------------------------------------------------------------------
# arrange results into a dataframe
ris <- t(ris[1,,])
df <- data.frame(fdr = ris[,2], niter = factor(ris[,1]))
df <- df[order(df$niter),]
#df

# obtain mean fdr and sd for every number of interactions found
m <- tapply(df$fdr, df$niter, mean)
s <- tapply(df$fdr, df$niter, sd)
n <- tapply(df$fdr, df$niter, length)
lo <- m - 2*s/sqrt(n)
up <- m + 2*s/sqrt(n)

# plot results
library(tidyverse)
toplot <- data.frame(nint = as.numeric(levels(df$niter)),
                     fdr = m, lo = pmax(lo,0), up = up)
myplot <- toplot %>% ggplot(aes(x = nint, y = fdr)) +
  geom_line(size = 2) +
  geom_point(shape = 15, size = 3) +
  geom_ribbon(aes(ymin = lo, ymax = up), alpha = 0.5) +
  scale_x_continuous(breaks=1:10) +
  labs(x = "Number of interactions found", y = "False discovery rate")
myplot


## --------------------------------------------------------------------------------------------
library(tidyverse)
library(pROC)

spam <- read.table("./spambase/spambase.data", header = F, sep = ",") # Spambase data taken from the UCI Machine Learning Repository
spam <- spam %>% rename(y = V58)
spam[,-58] <- apply(spam[,-58], 2, function(x) log(1+x))

# train-test splitting
tt <- rep("train", NROW(spam))
set.seed(123)
tt[sample(1:NROW(spam), 1536)] = "test"
table(tt)

xmat <- model.matrix(y~.-1, data = spam[tt=="train",])
y <- spam$y[tt=="train"]
xmat.t <- model.matrix(y~.-1, data = spam[tt=="test",])
y.t <- spam$y[tt=="test"]


## --------------------------------------------------------------------------------------------
# define error metrics

# classification threshold
prop.table(table(spam[,"y"]))
#S <- prop.table(table(spam[,"y"]))[2]
S <- 0.5

get_err <- function(pred, real = y.t){
  require(pROC)
  real <- as.numeric(real)
  missc <- mean(real != as.numeric(pred>S))
  auc <- as.numeric(roc(real, pred)$auc)
  crosse <- -mean(real*log(pred)+(1-real)*log(1-pred))
  out <- data.frame(missc = missc, auc = auc, crosse = crosse)
  t(out)
}


## --------------------------------------------------------------------------------------------
## models

# gbm
library(gbm)
m.gbm <- gbm(y~., distribution = "bernoulli", data = spam[tt=="train",],
             n.trees = 20000, shrinkage = 0.005, interaction.depth = 2,
             bag.fraction = 1, cv.folds = 10)
nt <- gbm.perf(m.gbm, method = "cv")

p.gbm <- predict(m.gbm, newdata = spam[tt=="test",], 
                 type = "response", n.trees = nt)
p.gbm.train <- predict(m.gbm, newdata = spam[tt=="train",], 
                       type = "response", n.trees = nt)
err <- data.frame(gbm = get_err(p.gbm))
err.train <- data.frame(gbm = get_err(p.gbm.train, y))
#t(err)


## --------------------------------------------------------------------------------------------
# glinternet
library(glinternet)
m.glinternet <- glinternet.cv(xmat, y,
                              numLevels = rep(1,ncol(spam)-1),
                              lambdaMinRatio =  0.005,
                              nLambda = 40,
                              nFolds = 10,
                              family = "binomial",
                              verbose = F,
                              numCores = 5)
plot(m.glinternet)
#print(m.glinternet)

p.glinternet <- predict(m.glinternet, xmat.t, type = "response", lambdaType = "lambdaHat")
p.glinternet.train <- predict(m.glinternet, xmat, type = "response", lambdaType = "lambdaHat")
err$glinternet <- get_err(p.glinternet)
err.train$glinternet <- get_err(p.glinternet.train, y)
#t(err)
nrow(m.glinternet$activeSet[[1]]$contcont) + nrow(m.glinternet$activeSet[[1]]$cont)


## --------------------------------------------------------------------------------------------
# hiernet
library(hierNet)
m.hiernet.path <- hierNet.logistic.path(xmat, y, diagonal = F, strong = F, nlam = 40, trace = 0)
hiercv <- hierNet.cv(m.hiernet.path, xmat, y, trace = 0)


## --------------------------------------------------------------------------------------------
plot(hiercv)
hiercv$nonzero[which(hiercv$lamlist==hiercv$lamhat)]
#print(hiercv)

p.hiernet <- predict(m.hiernet.path, newx = xmat.t)$prob[,which(hiercv$lamlist == hiercv$lamhat)]
p.hiernet.train <- predict(m.hiernet.path, newx = xmat)$prob[,which(hiercv$lamlist == hiercv$lamhat)]
err$hiernet <- get_err(p.hiernet)
err.train$hiernet <- get_err(p.hiernet.train, y)
#t(err)


## --------------------------------------------------------------------------------------------
# ridge
library(glmnet)
lgrid <- exp(seq(-10, 2, length = 100))
m.ridge <- cv.glmnet(xmat, y, lambda = lgrid, family = "binomial",
                     alpha = 0, nfolds = 10)

print(m.ridge)
par(mfrow = c(1,2))
plot(m.ridge)
plot(m.ridge$glmnet.fit, xvar = "lambda")
abline(v = log(c(m.ridge$lambda.min, m.ridge$lambda.1se)), lty = 2)
par(mfrow = c(1,1))

p.ridge <- predict(m.ridge, newx = xmat.t, s = "lambda.min", type = "response") %>% as.vector()
p.ridge.train <- predict(m.ridge, newx = xmat, s = "lambda.min", type = "response") %>% as.vector()
err$ridge = get_err(p.ridge)
err.train$ridge = get_err(p.ridge.train, y)
#t(err)


# lasso
m.lasso <- cv.glmnet(xmat, y, lambda.min.ratio = 0.001, family = "binomial",
                     alpha = 1, nfolds = 10)

print(m.lasso)
par(mfrow = c(1,2))
plot(m.lasso)
plot(m.lasso$glmnet.fit, xvar = "lambda")
abline(v = log(c(m.lasso$lambda.min, m.lasso$lambda.1se)), lty = 2)
par(mfrow = c(1,1))

p.lasso <- predict(m.lasso, newx = xmat.t, s = "lambda.min", type = "response") %>% as.vector()
p.lasso.train <- predict(m.lasso, newx = xmat, s = "lambda.min", type = "response") %>% as.vector()
err$lasso <- get_err(p.lasso)
err.train$lasso <- get_err(p.lasso.train, y)
#t(err)


## --------------------------------------------------------------------------------------------
# error
knitr::kable(t(err), digits = 5, col.names = c("Misclassification error", "AUC", "Cross-Entropy"))
# plot
missc <- data.frame(missc = as.numeric(c(err.train[1,],err[1,])),
                    modello = rep(c("gbm","glinternet","hierNet","ridge","lasso"), 2),
                    set = rep(c("train","test"), each = 5))
auc <- data.frame(auc = as.numeric(c(err.train[2,],err[2,])),
                  modello = rep(c("gbm","glinternet","hierNet","ridge","lasso"), 2),
                  set = rep(c("train","test"), each = 5))
crosse <- data.frame(crosse = as.numeric(c(err.train[3,],err[3,])),
                     modello = rep(c("gbm","glinternet","hierNet","ridge","lasso"), 2),
                     set = rep(c("train","test"), each = 5)) 
misscplot <- missc %>% ggplot(aes(x = modello, y = missc, fill = set)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "0-1 loss", x = "")
aucplot <- auc %>% ggplot(aes(x = modello, y = auc, fill = set)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "auc", x = "", )
crosseplot <- crosse %>% ggplot(aes(x = modello, y = crosse, fill = set)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "cross-entropy", x = "")
ggpubr::ggarrange(misscplot, aucplot, crosseplot, nrow = 3, common.legend = T, legend = "right")
