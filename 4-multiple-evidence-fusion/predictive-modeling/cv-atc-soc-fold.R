# Multiple Evidence Fusion
#
# Cross Validation by Folds Created by ATC/SOC Levels
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: November, 2013

library('randomForest')
library('foreach')
library('doMC')
registerDoMC(64)

stdize = function (x, center = TRUE, scale = TRUE) {

  nc = ncol(x)

  if (is.logical(center)) {

    if (center) {
      center = colMeans(x, na.rm = TRUE)
      x = sweep(x, 2, center)
    }

  } else if (is.numeric(center) && length(center) == nc) {

    x = sweep(x, 2, center)

  } else {

    stop("invalid 'center'")

  }

  if (is.logical(scale)) {

    if (scale) {

      scale = sqrt(colSums(sweep(x, 2, colMeans(x))^2)/(nrow(x) - 1))
      x = sweep(x, 2, scale, "/")

    }

  } else if (is.numeric(scale) && length(scale) == nc) {

    x = sweep(x, 2, scale, "/")

  } else {

    stop("invalid 'scale'")

  }

  return(x)

}

stdize2 = function (x, x1, center = TRUE, scale = TRUE) {
  # x: train
  # x1: test
  nc = ncol(x)
  if (is.logical(center)) {
    if (center) {
      center = colMeans(x, na.rm = TRUE)
      x1 = sweep(x1, 2, center)
    }
  } else if (is.numeric(center) && length(center) == nc) {
    x1 <- sweep(x1, 2, center)
  } else {
    stop("invalid 'center'")
  }

  if (is.logical(scale)) {
    if (scale) {
      scale <- sqrt(colSums(sweep(x, 2, colMeans(x))^2)/(nrow(x) - 1))
      x1 <- sweep(x1, 2, scale, "/")
    }
  } else if (is.numeric(scale) && length(scale) == nc) {
    x1 <- sweep(x1, 2, scale, "/")
  } else {
    stop("invalid 'scale'")
  }

  return(x1)

}

tinycv = function ( y, x, fold, model = c('logistic', 'rf'), ntree = 500, mtry = 4,
                    center = FALSE, scale = FALSE ) {

  foldnum = nlevels(fold)
  foldlvl = levels(fold)
  reslist = vector('list', foldnum)

  if (model == 'rf') {

    reslist <- foreach (i = 1:foldnum, .errorhandling = 'pass') %dopar% {

      idx <- which(fold == foldlvl[i])

      # splitting training/test set
      xtr = x[-idx, ]
      xte = x[idx, ]
      ytr = y[-idx]
      yte = y[idx]

      fit <- randomForest(xtr, ytr, ntree = ntree, mtry = mtry)

      pred.resp <- as.vector(predict(fit, xte, type = 'response'))
      pred.prob <- as.vector(predict(fit, xte, type = 'prob')[, 1])
      real.resp <- as.vector(yte)

      xxx <- data.frame(idx, real.resp, pred.resp, pred.prob)

    }

    names(reslist) = foldlvl

  }

  if (model == 'logistic') {

    reslist <- foreach (i = 1:foldnum, .errorhandling = 'pass') %dopar% {

      idx <- which(fold == foldlvl[i])

      if ( center == TRUE & scale == TRUE ) {
        xtr = x[-idx, ]
        xte = x[idx, ]
        ytr = y[-idx]
        yte = y[idx]
        xtr = stdize(xtr, center = TRUE, scale = TRUE)
        xte = stdize2(xtr, xte, center = TRUE, scale = TRUE)
      }
      if ( center == TRUE & scale == FALSE ) {
        xtr = x[-idx, ]
        xte = x[idx, ]
        ytr = y[-idx]
        yte = y[idx]
        xtr = stdize(xtr, center = TRUE, scale = FALSE)
        xte = stdize2(xtr, xte, center = TRUE, scale = FALSE)
      }
      if ( center == FALSE & scale == TRUE ) {
        stop('must center then scale')
      }
      if ( center == FALSE & scale == FALSE ) {
        xtr = x[-idx, ]
        xte = x[idx, ]
        ytr = y[-idx]
        yte = y[idx]
      }

      z = cbind(ytr, xtr)
      names(z)[1] = 'y'

      fit <- glm(as.formula('y ~ .'), data = z, family = binomial())

      pred.resp <- as.vector(ifelse(predict(fit, xte, type = 'response') > 0.50, '1', '0'))
      pred.prob <- as.vector(predict(fit, xte, type = 'response'))
      real.resp <- as.vector(yte)

      xxx <- data.frame(idx, real.resp, pred.resp, pred.prob)

    }

    names(reslist) = foldlvl

  }

  return(reslist)

}

# X1 = read.table('../20131108data/feature1.txt', header = FALSE, sep = '\t')
# X2 = read.table('../20131108data/feature2.txt', header = FALSE, sep = '\t')
# X3 = read.table('../20131108data/feature3.txt', header = FALSE, sep = '\t')
# X4 = read.table('../20131108data/feature4.txt', header = FALSE, sep = '\t')
# X5 = read.table('../20131108data/feature5.txt', header = FALSE, sep = '\t')
# X6 = read.table('../20131108data/feature6.txt', header = FALSE, sep = '\t')
# X7 = read.table('../20131108data/feature7.txt', header = FALSE, sep = '\t')
# X8 = read.table('../20131108data/feature8.txt', header = FALSE, sep = '\t')
# X9 = read.table('../20131108data/feature9.txt', header = FALSE, sep = '\t')
# X10 = read.table('../20131108data/feature10.txt', header = FALSE, sep = '\t')
#
# label = as.factor(rbind(as.matrix(rep(1, 24803)), as.matrix(rep(0, 24803))))

# load y x data
load('new.RData')

# load and make fold data
loofolddrug1 = as.factor(read.table('tassociation1', sep = '\t', header = FALSE)[, 1])
loofoldadr1  = as.factor(read.table('tassociation1', sep = '\t', header = FALSE)[, 2])

loofolddrug2 = as.factor(read.table('tassociation2', sep = '\t', header = FALSE)[, 1])
loofoldadr2  = as.factor(read.table('tassociation2', sep = '\t', header = FALSE)[, 2])

loofolddrug3 = as.factor(read.table('tassociation3', sep = '\t', header = FALSE)[, 1])
loofoldadr3  = as.factor(read.table('tassociation3', sep = '\t', header = FALSE)[, 2])

loofolddrug4 = as.factor(read.table('tassociation4', sep = '\t', header = FALSE)[, 1])
loofoldadr4  = as.factor(read.table('tassociation4', sep = '\t', header = FALSE)[, 2])

loofolddrug5 = as.factor(read.table('tassociation5', sep = '\t', header = FALSE)[, 1])
loofoldadr5  = as.factor(read.table('tassociation5', sep = '\t', header = FALSE)[, 2])

loofolddrug6 = as.factor(read.table('tassociation6', sep = '\t', header = FALSE)[, 1])
loofoldadr6  = as.factor(read.table('tassociation6', sep = '\t', header = FALSE)[, 2])

loofolddrug7 = as.factor(read.table('tassociation7', sep = '\t', header = FALSE)[, 1])
loofoldadr7  = as.factor(read.table('tassociation7', sep = '\t', header = FALSE)[, 2])

loofolddrug8 = as.factor(read.table('tassociation8', sep = '\t', header = FALSE)[, 1])
loofoldadr8  = as.factor(read.table('tassociation8', sep = '\t', header = FALSE)[, 2])

loofolddrug9 = as.factor(read.table('tassociation9', sep = '\t', header = FALSE)[, 1])
loofoldadr9  = as.factor(read.table('tassociation9', sep = '\t', header = FALSE)[, 2])

loofolddrug10 = as.factor(read.table('tassociation10', sep = '\t', header = FALSE)[, 1])
loofoldadr10  = as.factor(read.table('tassociation10', sep = '\t', header = FALSE)[, 2])

# get ATC and SOC fold

# load ATC dict
x = scan('xiao-drugbank-atc.csv', what = 'character', sep = ',')
cid = grep('CID', x)
atclist = vector('list', length(cid))
names(atclist) = x[cid]
for (i in 1:(length(cid) - 1)) atclist[[i]] = x[(cid[i] + 1):(cid[i+1] - 1)]
atclist[[length(cid)]] = x[(cid[length(cid)] + 1):length(x)]

# load SOC dict
soclist = strsplit(readLines('xiao-soc.csv'), ",")

# make fold
atcfold1 = rep(NA, 49606)
socfold1 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold1[i] = substr(sample(atclist[[as.integer(loofolddrug1[i])]], 1), 1, 1)
for (i in 1:49606) socfold1[i] = sample(soclist[[as.integer(loofoldadr1[i])]], 1)

atcfold1 = as.factor(atcfold1)
atcfolddict1 = levels(atcfold1)
levels(atcfold1) = 1:length(levels(atcfold1))

socfold1 = as.factor(socfold1)
socfolddict1 = levels(socfold1)
levels(socfold1) = 1:length(levels(socfold1))

atcfold2 = rep(NA, 49606)
socfold2 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold2[i] = substr(sample(atclist[[as.integer(loofolddrug2[i])]], 1), 1, 1)
for (i in 1:49606) socfold2[i] = sample(soclist[[as.integer(loofoldadr2[i])]], 1)

atcfold2 = as.factor(atcfold2)
atcfolddict2 = levels(atcfold2)
levels(atcfold2) = 1:length(levels(atcfold2))

socfold2 = as.factor(socfold2)
socfolddict2 = levels(socfold2)
levels(socfold2) = 1:length(levels(socfold2))

atcfold3 = rep(NA, 49606)
socfold3 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold3[i] = substr(sample(atclist[[as.integer(loofolddrug3[i])]], 1), 1, 1)
for (i in 1:49606) socfold3[i] = sample(soclist[[as.integer(loofoldadr3[i])]], 1)

atcfold3 = as.factor(atcfold3)
atcfolddict3 = levels(atcfold3)
levels(atcfold3) = 1:length(levels(atcfold3))

socfold3 = as.factor(socfold3)
socfolddict3 = levels(socfold3)
levels(socfold3) = 1:length(levels(socfold3))

atcfold4 = rep(NA, 49606)
socfold4 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold4[i] = substr(sample(atclist[[as.integer(loofolddrug4[i])]], 1), 1, 1)
for (i in 1:49606) socfold4[i] = sample(soclist[[as.integer(loofoldadr4[i])]], 1)

atcfold4 = as.factor(atcfold4)
atcfolddict4 = levels(atcfold4)
levels(atcfold4) = 1:length(levels(atcfold4))

socfold4 = as.factor(socfold4)
socfolddict4 = levels(socfold4)
levels(socfold4) = 1:length(levels(socfold4))

atcfold5 = rep(NA, 49606)
socfold5 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold5[i] = substr(sample(atclist[[as.integer(loofolddrug5[i])]], 1), 1, 1)
for (i in 1:49606) socfold5[i] = sample(soclist[[as.integer(loofoldadr5[i])]], 1)

atcfold5 = as.factor(atcfold5)
atcfolddict5 = levels(atcfold5)
levels(atcfold5) = 1:length(levels(atcfold5))

socfold5 = as.factor(socfold5)
socfolddict5 = levels(socfold5)
levels(socfold5) = 1:length(levels(socfold5))

atcfold6 = rep(NA, 49606)
socfold6 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold6[i] = substr(sample(atclist[[as.integer(loofolddrug6[i])]], 1), 1, 1)
for (i in 1:49606) socfold6[i] = sample(soclist[[as.integer(loofoldadr6[i])]], 1)

atcfold6 = as.factor(atcfold6)
atcfolddict6 = levels(atcfold6)
levels(atcfold6) = 1:length(levels(atcfold6))

socfold6 = as.factor(socfold6)
socfolddict6 = levels(socfold6)
levels(socfold6) = 1:length(levels(socfold6))

atcfold7 = rep(NA, 49606)
socfold7 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold7[i] = substr(sample(atclist[[as.integer(loofolddrug7[i])]], 1), 1, 1)
for (i in 1:49606) socfold7[i] = sample(soclist[[as.integer(loofoldadr7[i])]], 1)

atcfold7 = as.factor(atcfold7)
atcfolddict7 = levels(atcfold7)
levels(atcfold7) = 1:length(levels(atcfold7))

socfold7 = as.factor(socfold7)
socfolddict7 = levels(socfold7)
levels(socfold7) = 1:length(levels(socfold7))

atcfold8 = rep(NA, 49606)
socfold8 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold8[i] = substr(sample(atclist[[as.integer(loofolddrug8[i])]], 1), 1, 1)
for (i in 1:49606) socfold8[i] = sample(soclist[[as.integer(loofoldadr8[i])]], 1)

atcfold8 = as.factor(atcfold8)
atcfolddict8 = levels(atcfold8)
levels(atcfold8) = 1:length(levels(atcfold8))

socfold8 = as.factor(socfold8)
socfolddict8 = levels(socfold8)
levels(socfold8) = 1:length(levels(socfold8))

atcfold9 = rep(NA, 49606)
socfold9 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold9[i] = substr(sample(atclist[[as.integer(loofolddrug9[i])]], 1), 1, 1)
for (i in 1:49606) socfold9[i] = sample(soclist[[as.integer(loofoldadr9[i])]], 1)

atcfold9 = as.factor(atcfold9)
atcfolddict9 = levels(atcfold9)
levels(atcfold9) = 1:length(levels(atcfold9))

socfold9 = as.factor(socfold9)
socfolddict9 = levels(socfold9)
levels(socfold9) = 1:length(levels(socfold9))

atcfold10 = rep(NA, 49606)
socfold10 = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold10[i] = substr(sample(atclist[[as.integer(loofolddrug10[i])]], 1), 1, 1)
for (i in 1:49606) socfold10[i] = sample(soclist[[as.integer(loofoldadr10[i])]], 1)

atcfold10 = as.factor(atcfold10)
atcfolddict10 = levels(atcfold10)
levels(atcfold10) = 1:length(levels(atcfold10))

socfold10 = as.factor(socfold10)
socfolddict10 = levels(socfold10)
levels(socfold10) = 1:length(levels(socfold10))



# All Feature RF ATC-fold

# 1
all_atcfold_rf_1 = tinycv(y = label, x = X1, fold = atcfold1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_1, file = 'all_atcfold_rf_1.Rdata')

# 2
all_atcfold_rf_2 = tinycv(y = label, x = X2, fold = atcfold2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_2, file = 'all_atcfold_rf_2.Rdata')

# 3
all_atcfold_rf_3 = tinycv(y = label, x = X3, fold = atcfold3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_3, file = 'all_atcfold_rf_3.Rdata')

# 4
all_atcfold_rf_4 = tinycv(y = label, x = X4, fold = atcfold4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_4, file = 'all_atcfold_rf_4.Rdata')

# 5
all_atcfold_rf_5 = tinycv(y = label, x = X5, fold = atcfold5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_5, file = 'all_atcfold_rf_5.Rdata')

# 6
all_atcfold_rf_6 = tinycv(y = label, x = X6, fold = atcfold6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_6, file = 'all_atcfold_rf_6.Rdata')

# 7
all_atcfold_rf_7 = tinycv(y = label, x = X7, fold = atcfold7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_7, file = 'all_atcfold_rf_7.Rdata')

# 8
all_atcfold_rf_8 = tinycv(y = label, x = X8, fold = atcfold8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_8, file = 'all_atcfold_rf_8.Rdata')

# 9
all_atcfold_rf_9 = tinycv(y = label, x = X9, fold = atcfold9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_9, file = 'all_atcfold_rf_9.Rdata')

# 10
all_atcfold_rf_10 = tinycv(y = label, x = X10, fold = atcfold10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_atcfold_rf_10, file = 'all_atcfold_rf_10.Rdata')



# All Feature RF SOC-fold

# 1
all_socfold_rf_1 = tinycv(y = label, x = X1, fold = socfold1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_1, file = 'all_socfold_rf_1.Rdata')

# 2
all_socfold_rf_2 = tinycv(y = label, x = X2, fold = socfold2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_2, file = 'all_socfold_rf_2.Rdata')

# 3
all_socfold_rf_3 = tinycv(y = label, x = X3, fold = socfold3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_3, file = 'all_socfold_rf_3.Rdata')

# 4
all_socfold_rf_4 = tinycv(y = label, x = X4, fold = socfold4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_4, file = 'all_socfold_rf_4.Rdata')

# 5
all_socfold_rf_5 = tinycv(y = label, x = X5, fold = socfold5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_5, file = 'all_socfold_rf_5.Rdata')

# 6
all_socfold_rf_6 = tinycv(y = label, x = X6, fold = socfold6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_6, file = 'all_socfold_rf_6.Rdata')

# 7
all_socfold_rf_7 = tinycv(y = label, x = X7, fold = socfold7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_7, file = 'all_socfold_rf_7.Rdata')

# 8
all_socfold_rf_8 = tinycv(y = label, x = X8, fold = socfold8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_8, file = 'all_socfold_rf_8.Rdata')

# 9
all_socfold_rf_9 = tinycv(y = label, x = X9, fold = socfold9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_9, file = 'all_socfold_rf_9.Rdata')

# 10
all_socfold_rf_10 = tinycv(y = label, x = X10, fold = socfold10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_socfold_rf_10, file = 'all_socfold_rf_10.Rdata')
