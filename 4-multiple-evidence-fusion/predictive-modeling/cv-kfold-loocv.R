# Multiple Evidence Fusion
#
# k-Fold Cross Validation and Leave-One-Out Cross Validation
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

      xtr = x[-idx, ]
      xte = x[idx, ]
      ytr = y[-idx]
      yte = y[idx]

      fit <- randomForest(xtr, ytr, ntree = ntree, mtry = mtry)

      pred.resp <- as.vector(predict(fit, xte, type = 'response'))
      pred.prob <- as.vector(predict(fit, xte, type = 'prob')[, 1])  # 只取第一列 预测时注意对应的label是哪个!
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



# shuffle the folds randomly
set.seed(1234)
drugmapidx = sample(rep(1:10, 75)[1:746])
set.seed(5678)
adrmapidx = sample(rep(1:10, 82)[1:817])

tenfolddrug1 = rep(NA, 746)
tenfoldadr1  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug1[which(loofolddrug1 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr1[which(loofoldadr1 == i)] = adrmapidx[i]
tenfolddrug1 = as.factor(tenfolddrug1)
tenfoldadr1  = as.factor(tenfoldadr1)

tenfolddrug2 = rep(NA, 746)
tenfoldadr2  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug2[which(loofolddrug2 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr2[which(loofoldadr2 == i)] = adrmapidx[i]
tenfolddrug2 = as.factor(tenfolddrug2)
tenfoldadr2  = as.factor(tenfoldadr2)

tenfolddrug3 = rep(NA, 746)
tenfoldadr3  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug3[which(loofolddrug3 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr3[which(loofoldadr3 == i)] = adrmapidx[i]
tenfolddrug3 = as.factor(tenfolddrug3)
tenfoldadr3  = as.factor(tenfoldadr3)

tenfolddrug4 = rep(NA, 746)
tenfoldadr4  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug4[which(loofolddrug4 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr4[which(loofoldadr4 == i)] = adrmapidx[i]
tenfolddrug4 = as.factor(tenfolddrug4)
tenfoldadr4  = as.factor(tenfoldadr4)

tenfolddrug5 = rep(NA, 746)
tenfoldadr5  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug5[which(loofolddrug5 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr5[which(loofoldadr5 == i)] = adrmapidx[i]
tenfolddrug5 = as.factor(tenfolddrug5)
tenfoldadr5  = as.factor(tenfoldadr5)

tenfolddrug6 = rep(NA, 746)
tenfoldadr6  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug6[which(loofolddrug6 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr6[which(loofoldadr6 == i)] = adrmapidx[i]
tenfolddrug6 = as.factor(tenfolddrug6)
tenfoldadr6  = as.factor(tenfoldadr6)

tenfolddrug7 = rep(NA, 746)
tenfoldadr7  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug7[which(loofolddrug7 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr7[which(loofoldadr7 == i)] = adrmapidx[i]
tenfolddrug7 = as.factor(tenfolddrug7)
tenfoldadr7  = as.factor(tenfoldadr7)

tenfolddrug8 = rep(NA, 746)
tenfoldadr8  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug8[which(loofolddrug8 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr8[which(loofoldadr8 == i)] = adrmapidx[i]
tenfolddrug8 = as.factor(tenfolddrug8)
tenfoldadr8  = as.factor(tenfoldadr8)

tenfolddrug9 = rep(NA, 746)
tenfoldadr9  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug9[which(loofolddrug9 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr9[which(loofoldadr9 == i)] = adrmapidx[i]
tenfolddrug9 = as.factor(tenfolddrug9)
tenfoldadr9  = as.factor(tenfoldadr9)

tenfolddrug10 = rep(NA, 746)
tenfoldadr10  = rep(NA, 817)
for ( i in 1:746 ) tenfolddrug10[which(loofolddrug10 == i)] = drugmapidx[i]
for ( i in 1:817 ) tenfoldadr10[which(loofoldadr10 == i)] = adrmapidx[i]
tenfolddrug10 = as.factor(tenfolddrug10)
tenfoldadr10  = as.factor(tenfoldadr10)

rm(i)



# # Test Random Fold
#
# set.seed(1234)
# rnd5fold = as.factor(sample(rep(1:5, 9922)[1:49606]))
#
# # Drug Random 5-fold CV RF
#
# # 1
# drug_rndcv_rf_1 = tinycv(y = label, x = X1[, 1:11], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_rndcv_rf_1, file = 'drug_rndcv_rf_1.Rdata')
#
# # 2
# drug_rndcv_rf_2 = tinycv(y = label, x = X2[, 1:11], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_rndcv_rf_2, file = 'drug_rndcv_rf_2.Rdata')
#
# # 3
# drug_rndcv_rf_3 = tinycv(y = label, x = X3[, 1:11], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_rndcv_rf_3, file = 'drug_rndcv_rf_3.Rdata')
#
# # ADR Random 5-fold CV RF
#
# # 1
# adr_rndcv_rf_1 = tinycv(y = label, x = X1[, 12:18], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_rndcv_rf_1, file = 'adr_rndcv_rf_1.Rdata')
#
# # 2
# adr_rndcv_rf_2 = tinycv(y = label, x = X2[, 12:18], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_rndcv_rf_2, file = 'adr_rndcv_rf_2.Rdata')
#
# # 3
# adr_rndcv_rf_3 = tinycv(y = label, x = X3[, 12:18], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_rndcv_rf_3, file = 'adr_rndcv_rf_3.Rdata')
#
# # Node Random 5-fold CV RF
#
# # 1
# node_rndcv_rf_1 = tinycv(y = label, x = X1[, c(1:8, 12:15)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(node_rndcv_rf_1, file = 'node_rndcv_rf_1.Rdata')
#
# # 2
# node_rndcv_rf_2 = tinycv(y = label, x = X2[, c(1:8, 12:15)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(node_rndcv_rf_2, file = 'node_rndcv_rf_2.Rdata')
#
# # 3
# node_rndcv_rf_3 = tinycv(y = label, x = X3[, c(1:8, 12:15)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(node_rndcv_rf_3, file = 'node_rndcv_rf_3.Rdata')
#
# # Network Random 5-fold CV RF
#
# # 1
# network_rndcv_rf_1 = tinycv(y = label, x = X1[, c(9:11, 16:19)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(network_rndcv_rf_1, file = 'network_rndcv_rf_1.Rdata')
#
# # 2
# network_rndcv_rf_2 = tinycv(y = label, x = X2[, c(9:11, 16:19)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(network_rndcv_rf_2, file = 'network_rndcv_rf_2.Rdata')
#
# # 3
# network_rndcv_rf_3 = tinycv(y = label, x = X3[, c(9:11, 16:19)], fold = rnd5fold, model = 'rf', center = FALSE, scale = FALSE)
# save(network_rndcv_rf_3, file = 'network_rndcv_rf_3.Rdata')
#
# # All Random 5-fold CV RF
#
# # 1
# all_rndcv_rf_1 = tinycv(y = label, x = X1, fold = rnd5fold, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
# save(all_rndcv_rf_1, file = 'all_rndcv_rf_1.Rdata')
#
# # 2
# all_rndcv_rf_2 = tinycv(y = label, x = X2, fold = rnd5fold, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
# save(all_rndcv_rf_2, file = 'all_rndcv_rf_2.Rdata')
#
# # 3
# all_rndcv_rf_3 = tinycv(y = label, x = X3, fold = rnd5fold, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
# save(all_rndcv_rf_3, file = 'all_rndcv_rf_3.Rdata')
#
# # Test Random Fold



# # Test feauture selection and didn't change much
# test_tenfoldcv_rf_1 = tinycv(y = label, x = X1[, c(14, 2)], fold = tenfolddrug, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_1, file = 'test_tenfoldcv_rf_1.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_1)
# tmp = tmp[order(tmp$idx), ]
# auc1 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_2 = tinycv(y = label, x = X1[, c(14, 2, 1)], fold = tenfolddrug, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_2, file = 'test_tenfoldcv_rf_2.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_2)
# tmp = tmp[order(tmp$idx), ]
# auc2= auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_3 = tinycv(y = label, x = X1[, c(14, 2, 1, 12)], fold = tenfolddrug, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_3, file = 'test_tenfoldcv_rf_3.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_3)
# tmp = tmp[order(tmp$idx), ]
# auc3 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_4 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_4, file = 'test_tenfoldcv_rf_4.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_4)
# tmp = tmp[order(tmp$idx), ]
# auc4 = auc(label, tmp$pred.prob)
#
# # best one
# test_tenfoldcv_rf_5 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_5, file = 'test_tenfoldcv_rf_5.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_5)
# tmp = tmp[order(tmp$idx), ]
# auc5 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_6 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_6, file = 'test_tenfoldcv_rf_6.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_6)
# tmp = tmp[order(tmp$idx), ]
# auc6 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_7 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_7, file = 'test_tenfoldcv_rf_7.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_7)
# tmp = tmp[order(tmp$idx), ]
# auc7 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_8 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_8, file = 'test_tenfoldcv_rf_8.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_8)
# tmp = tmp[order(tmp$idx), ]
# auc8 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_9 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_9, file = 'test_tenfoldcv_rf_9.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_9)
# tmp = tmp[order(tmp$idx), ]
# auc9 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_10 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5, 4)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_10, file = 'test_tenfoldcv_rf_10.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_10)
# tmp = tmp[order(tmp$idx), ]
# auc10 = auc(label, tmp$pred.prob)
#
# test_tenfoldcv_rf_11 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5, 4, 8)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_11, file = 'test_tenfoldcv_rf_11.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_11)
# tmp = tmp[order(tmp$idx), ]
# auc11 = auc(label, tmp$pred.prob)
#
# aucall = c(auc1, auc2, auc3, auc4, auc5, auc6, auc7, auc8, auc9, auc10, auc11)
# plot(aucall)
#
#
# test_tenfoldcv_rf_adr_1 = tinycv(y = label, x = X1[, c(14, 2)], fold = tenfoldadr, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_1, file = 'test_tenfoldcv_rf_adr_1.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_1)
# tmp = tmp[order(tmp$idx), ]
# auc1 = auc(label, tmp$pred.prob)
# print(auc1)
#
# test_tenfoldcv_rf_adr_2 = tinycv(y = label, x = X1[, c(14, 2, 1)], fold = tenfoldadr, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_2, file = 'test_tenfoldcv_rf_adr_2.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_2)
# tmp = tmp[order(tmp$idx), ]
# auc2= auc(label, tmp$pred.prob)
# print(auc2)
#
# test_tenfoldcv_rf_adr_3 = tinycv(y = label, x = X1[, c(14, 2, 1, 12)], fold = tenfoldadr, model = 'rf', mtry = 2, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_3, file = 'test_tenfoldcv_rf_adr_3.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_3)
# tmp = tmp[order(tmp$idx), ]
# auc3 = auc(label, tmp$pred.prob)
# print(auc3)
#
# test_tenfoldcv_rf_adr_4 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_4, file = 'test_tenfoldcv_rf_adr_4.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_4)
# tmp = tmp[order(tmp$idx), ]
# auc4 = auc(label, tmp$pred.prob)
# print(auc4)
#
# test_tenfoldcv_rf_adr_5 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_5, file = 'test_tenfoldcv_rf_adr_5.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_5)
# tmp = tmp[order(tmp$idx), ]
# auc5 = auc(label, tmp$pred.prob)
# print(auc5)
#
# test_tenfoldcv_rf_adr_6 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_6, file = 'test_tenfoldcv_rf_adr_6.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_6)
# tmp = tmp[order(tmp$idx), ]
# auc6 = auc(label, tmp$pred.prob)
# print(auc6)
#
# test_tenfoldcv_rf_adr_7 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_7, file = 'test_tenfoldcv_rf_adr_7.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_7)
# tmp = tmp[order(tmp$idx), ]
# auc7 = auc(label, tmp$pred.prob)
# print(auc7)
#
# test_tenfoldcv_rf_adr_8 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_8, file = 'test_tenfoldcv_rf_adr_8.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_8)
# tmp = tmp[order(tmp$idx), ]
# auc8 = auc(label, tmp$pred.prob)
# print(auc8)
#
# test_tenfoldcv_rf_adr_9 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_9, file = 'test_tenfoldcv_rf_adr_9.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_9)
# tmp = tmp[order(tmp$idx), ]
# auc9 = auc(label, tmp$pred.prob)
# print(auc9)
#
# test_tenfoldcv_rf_adr_10 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5, 4)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_10, file = 'test_tenfoldcv_rf_adr_10.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_10)
# tmp = tmp[order(tmp$idx), ]
# auc10 = auc(label, tmp$pred.prob)
# print(auc10)
#
# test_tenfoldcv_rf_adr_11 = tinycv(y = label, x = X1[, c(14, 2, 1, 12, 13, 15, 6, 7, 3, 5, 4, 8)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# save(test_tenfoldcv_rf_adr_11, file = 'test_tenfoldcv_rf_adr_11.Rdata')
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_11)
# tmp = tmp[order(tmp$idx), ]
# auc11 = auc(label, tmp$pred.prob)
# print(auc11)
#
# aucall = c(auc1, auc2, auc3, auc4, auc5, auc6, auc7, auc8, auc9, auc10, auc11)
# plot(aucall)
#
# test_tenfoldcv_rf_adr_12 = tinycv(y = label, x = X1[, c(2, 1, 12, 13, 15, 6, 7, 3, 5, 4, 8)], fold = tenfoldadr, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_12)
# tmp = tmp[order(tmp$idx), ]
# auc12 = auc(label, tmp$pred.prob)
# print(auc12)
#
# test_tenfoldcv_rf_adr_12 = tinycv(y = label, x = X1[, c(2, 1, 12, 13, 15)], fold = tenfolddrug, model = 'rf', mtry = 3, center = FALSE, scale = FALSE)
# tmp = do.call(rbind, test_tenfoldcv_rf_adr_12)
# tmp = tmp[order(tmp$idx), ]
# auc12 = auc(label, tmp$pred.prob)
# print(auc12)



# # Didn't use 8 drug and 4 ADR here because these will involve many problems
# # Drug 10-fold CV RF
#
# # 1
# drug_tenfoldcv_rf_1 = tinycv(y = label, x = X1[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_1, file = 'drug_tenfoldcv_rf_1.Rdata')
#
# # 2
# drug_tenfoldcv_rf_2 = tinycv(y = label, x = X2[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_2, file = 'drug_tenfoldcv_rf_2.Rdata')
#
# # 3
# drug_tenfoldcv_rf_3 = tinycv(y = label, x = X3[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_3, file = 'drug_tenfoldcv_rf_3.Rdata')
#
# # 4
# drug_tenfoldcv_rf_4 = tinycv(y = label, x = X4[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_4, file = 'drug_tenfoldcv_rf_4.Rdata')
#
# # 5
# drug_tenfoldcv_rf_5 = tinycv(y = label, x = X5[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_5, file = 'drug_tenfoldcv_rf_5.Rdata')
#
# # 6
# drug_tenfoldcv_rf_6 = tinycv(y = label, x = X6[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_6, file = 'drug_tenfoldcv_rf_6.Rdata')
#
# # 7
# drug_tenfoldcv_rf_7 = tinycv(y = label, x = X7[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_7, file = 'drug_tenfoldcv_rf_7.Rdata')
#
# # 8
# drug_tenfoldcv_rf_8 = tinycv(y = label, x = X8[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_8, file = 'drug_tenfoldcv_rf_8.Rdata')
#
# # 9
# drug_tenfoldcv_rf_9 = tinycv(y = label, x = X9[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_9, file = 'drug_tenfoldcv_rf_9.Rdata')
#
# # 10
# drug_tenfoldcv_rf_10 = tinycv(y = label, x = X10[, 1:8], fold = tenfolddrug, model = 'rf', center = FALSE, scale = FALSE)
# save(drug_tenfoldcv_rf_10, file = 'drug_tenfoldcv_rf_10.Rdata')
#
#
#
# # ADR_3 10-fold CV RF
#
# # 1
# adr_3_tenfoldcv_rf_1 = tinycv(y = label, x = X1[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_1, file = 'adr_3_tenfoldcv_rf_1.Rdata')
#
# # 2
# adr_3_tenfoldcv_rf_2 = tinycv(y = label, x = X2[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_2, file = 'adr_3_tenfoldcv_rf_2.Rdata')
#
# # 3
# adr_3_tenfoldcv_rf_3 = tinycv(y = label, x = X3[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_3, file = 'adr_3_tenfoldcv_rf_3.Rdata')
#
# # 4
# adr_3_tenfoldcv_rf_4 = tinycv(y = label, x = X4[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_4, file = 'adr_3_tenfoldcv_rf_4.Rdata')
#
# # 5
# adr_3_tenfoldcv_rf_5 = tinycv(y = label, x = X5[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_5, file = 'adr_3_tenfoldcv_rf_5.Rdata')
#
# # 6
# adr_3_tenfoldcv_rf_6 = tinycv(y = label, x = X6[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_6, file = 'adr_3_tenfoldcv_rf_6.Rdata')
#
# # 7
# adr_3_tenfoldcv_rf_7 = tinycv(y = label, x = X7[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_7, file = 'adr_3_tenfoldcv_rf_7.Rdata')
#
# # 8
# adr_3_tenfoldcv_rf_8 = tinycv(y = label, x = X8[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_8, file = 'adr_3_tenfoldcv_rf_8.Rdata')
#
# # 9
# adr_3_tenfoldcv_rf_9 = tinycv(y = label, x = X9[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_9, file = 'adr_3_tenfoldcv_rf_9.Rdata')
#
# # 10
# adr_3_tenfoldcv_rf_10 = tinycv(y = label, x = X10[, c(12, 13, 15)], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_3_tenfoldcv_rf_10, file = 'adr_3_tenfoldcv_rf_10.Rdata')
#
#
#
# # ADR_4 10-fold CV RF
#
# # 1
# adr_4_tenfoldcv_rf_1 = tinycv(y = label, x = X1[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_1, file = 'adr_4_tenfoldcv_rf_1.Rdata')
#
# # 2
# adr_4_tenfoldcv_rf_2 = tinycv(y = label, x = X2[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_2, file = 'adr_4_tenfoldcv_rf_2.Rdata')
#
# # 3
# adr_4_tenfoldcv_rf_3 = tinycv(y = label, x = X3[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_3, file = 'adr_4_tenfoldcv_rf_3.Rdata')
#
# # 4
# adr_4_tenfoldcv_rf_4 = tinycv(y = label, x = X4[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_4, file = 'adr_4_tenfoldcv_rf_4.Rdata')
#
# # 5
# adr_4_tenfoldcv_rf_5 = tinycv(y = label, x = X5[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_5, file = 'adr_4_tenfoldcv_rf_5.Rdata')
#
# # 6
# adr_4_tenfoldcv_rf_6 = tinycv(y = label, x = X6[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_6, file = 'adr_4_tenfoldcv_rf_6.Rdata')
#
# # 7
# adr_4_tenfoldcv_rf_7 = tinycv(y = label, x = X7[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_7, file = 'adr_4_tenfoldcv_rf_7.Rdata')
#
# # 8
# adr_4_tenfoldcv_rf_8 = tinycv(y = label, x = X8[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_8, file = 'adr_4_tenfoldcv_rf_8.Rdata')
#
# # 9
# adr_4_tenfoldcv_rf_9 = tinycv(y = label, x = X9[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_9, file = 'adr_4_tenfoldcv_rf_9.Rdata')
#
# # 10
# adr_4_tenfoldcv_rf_10 = tinycv(y = label, x = X10[, 12:15], fold = tenfoldadr, model = 'rf', center = FALSE, scale = FALSE)
# save(adr_4_tenfoldcv_rf_10, file = 'adr_4_tenfoldcv_rf_10.Rdata')



# Node 10-fold CV RF DrugFold

# 1
node_tenfoldcv_drugfold_rf_1 = tinycv(y = label, x = X1[, c(1:8, 12:15)], fold = tenfolddrug1, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_1, file = 'node_tenfoldcv_drugfold_rf_1.Rdata')

# 2
node_tenfoldcv_drugfold_rf_2 = tinycv(y = label, x = X2[, c(1:8, 12:15)], fold = tenfolddrug2, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_2, file = 'node_tenfoldcv_drugfold_rf_2.Rdata')

# 3
node_tenfoldcv_drugfold_rf_3 = tinycv(y = label, x = X3[, c(1:8, 12:15)], fold = tenfolddrug3, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_3, file = 'node_tenfoldcv_drugfold_rf_3.Rdata')

# 4
node_tenfoldcv_drugfold_rf_4 = tinycv(y = label, x = X4[, c(1:8, 12:15)], fold = tenfolddrug4, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_4, file = 'node_tenfoldcv_drugfold_rf_4.Rdata')

# 5
node_tenfoldcv_drugfold_rf_5 = tinycv(y = label, x = X5[, c(1:8, 12:15)], fold = tenfolddrug5, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_5, file = 'node_tenfoldcv_drugfold_rf_5.Rdata')

# 6
node_tenfoldcv_drugfold_rf_6 = tinycv(y = label, x = X6[, c(1:8, 12:15)], fold = tenfolddrug6, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_6, file = 'node_tenfoldcv_drugfold_rf_6.Rdata')

# 7
node_tenfoldcv_drugfold_rf_7 = tinycv(y = label, x = X7[, c(1:8, 12:15)], fold = tenfolddrug7, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_7, file = 'node_tenfoldcv_drugfold_rf_7.Rdata')

# 8
node_tenfoldcv_drugfold_rf_8 = tinycv(y = label, x = X8[, c(1:8, 12:15)], fold = tenfolddrug8, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_8, file = 'node_tenfoldcv_drugfold_rf_8.Rdata')

# 9
node_tenfoldcv_drugfold_rf_9 = tinycv(y = label, x = X9[, c(1:8, 12:15)], fold = tenfolddrug9, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_9, file = 'node_tenfoldcv_drugfold_rf_9.Rdata')

# 10
node_tenfoldcv_drugfold_rf_10 = tinycv(y = label, x = X10[, c(1:8, 12:15)], fold = tenfolddrug10, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_drugfold_rf_10, file = 'node_tenfoldcv_drugfold_rf_10.Rdata')



# Node 10-fold CV RF DrugADR

# 1
node_tenfoldcv_adrfold_rf_1 = tinycv(y = label, x = X1[, c(1:8, 12:15)], fold = tenfoldadr1, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_1, file = 'node_tenfoldcv_adrfold_rf_1.Rdata')

# 2
node_tenfoldcv_adrfold_rf_2 = tinycv(y = label, x = X2[, c(1:8, 12:15)], fold = tenfoldadr2, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_2, file = 'node_tenfoldcv_adrfold_rf_2.Rdata')

# 3
node_tenfoldcv_adrfold_rf_3 = tinycv(y = label, x = X3[, c(1:8, 12:15)], fold = tenfoldadr3, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_3, file = 'node_tenfoldcv_adrfold_rf_3.Rdata')

# 4
node_tenfoldcv_adrfold_rf_4 = tinycv(y = label, x = X4[, c(1:8, 12:15)], fold = tenfoldadr4, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_4, file = 'node_tenfoldcv_adrfold_rf_4.Rdata')

# 5
node_tenfoldcv_adrfold_rf_5 = tinycv(y = label, x = X5[, c(1:8, 12:15)], fold = tenfoldadr5, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_5, file = 'node_tenfoldcv_adrfold_rf_5.Rdata')

# 6
node_tenfoldcv_adrfold_rf_6 = tinycv(y = label, x = X6[, c(1:8, 12:15)], fold = tenfoldadr6, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_6, file = 'node_tenfoldcv_adrfold_rf_6.Rdata')

# 7
node_tenfoldcv_adrfold_rf_7 = tinycv(y = label, x = X7[, c(1:8, 12:15)], fold = tenfoldadr7, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_7, file = 'node_tenfoldcv_adrfold_rf_7.Rdata')

# 8
node_tenfoldcv_adrfold_rf_8 = tinycv(y = label, x = X8[, c(1:8, 12:15)], fold = tenfoldadr8, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_8, file = 'node_tenfoldcv_adrfold_rf_8.Rdata')

# 9
node_tenfoldcv_adrfold_rf_9 = tinycv(y = label, x = X9[, c(1:8, 12:15)], fold = tenfoldadr9, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_9, file = 'node_tenfoldcv_adrfold_rf_9.Rdata')

# 10
node_tenfoldcv_adrfold_rf_10 = tinycv(y = label, x = X10[, c(1:8, 12:15)], fold = tenfoldadr10, model = 'rf', center = FALSE, scale = FALSE)
save(node_tenfoldcv_adrfold_rf_10, file = 'node_tenfoldcv_adrfold_rf_10.Rdata')



# Network 10-fold CV RF DrugFold

# 1
network_tenfoldcv_drugfold_rf_1 = tinycv(y = label, x = X1[, c(9:11, 16:19)], fold = tenfolddrug1, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_1, file = 'network_tenfoldcv_drugfold_rf_1.Rdata')

# 2
network_tenfoldcv_drugfold_rf_2 = tinycv(y = label, x = X2[, c(9:11, 16:19)], fold = tenfolddrug2, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_2, file = 'network_tenfoldcv_drugfold_rf_2.Rdata')

# 3
network_tenfoldcv_drugfold_rf_3 = tinycv(y = label, x = X3[, c(9:11, 16:19)], fold = tenfolddrug3, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_3, file = 'network_tenfoldcv_drugfold_rf_3.Rdata')

# 4
network_tenfoldcv_drugfold_rf_4 = tinycv(y = label, x = X4[, c(9:11, 16:19)], fold = tenfolddrug4, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_4, file = 'network_tenfoldcv_drugfold_rf_4.Rdata')

# 5
network_tenfoldcv_drugfold_rf_5 = tinycv(y = label, x = X5[, c(9:11, 16:19)], fold = tenfolddrug5, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_5, file = 'network_tenfoldcv_drugfold_rf_5.Rdata')

# 6
network_tenfoldcv_drugfold_rf_6 = tinycv(y = label, x = X6[, c(9:11, 16:19)], fold = tenfolddrug6, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_6, file = 'network_tenfoldcv_drugfold_rf_6.Rdata')

# 7
network_tenfoldcv_drugfold_rf_7 = tinycv(y = label, x = X7[, c(9:11, 16:19)], fold = tenfolddrug7, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_7, file = 'network_tenfoldcv_drugfold_rf_7.Rdata')

# 8
network_tenfoldcv_drugfold_rf_8 = tinycv(y = label, x = X8[, c(9:11, 16:19)], fold = tenfolddrug8, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_8, file = 'network_tenfoldcv_drugfold_rf_8.Rdata')

# 9
network_tenfoldcv_drugfold_rf_9 = tinycv(y = label, x = X9[, c(9:11, 16:19)], fold = tenfolddrug9, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_9, file = 'network_tenfoldcv_drugfold_rf_9.Rdata')

# 10
network_tenfoldcv_drugfold_rf_10 = tinycv(y = label, x = X10[, c(9:11, 16:19)], fold = tenfolddrug10, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_drugfold_rf_10, file = 'network_tenfoldcv_drugfold_rf_10.Rdata')



# Network 10-fold CV RF FoldADR

# 1
network_tenfoldcv_adrfold_rf_1 = tinycv(y = label, x = X1[, c(9:11, 16:19)], fold = tenfoldadr1, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_1, file = 'network_tenfoldcv_adrfold_rf_1.Rdata')

# 2
network_tenfoldcv_adrfold_rf_2 = tinycv(y = label, x = X2[, c(9:11, 16:19)], fold = tenfoldadr2, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_2, file = 'network_tenfoldcv_adrfold_rf_2.Rdata')

# 3
network_tenfoldcv_adrfold_rf_3 = tinycv(y = label, x = X3[, c(9:11, 16:19)], fold = tenfoldadr3, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_3, file = 'network_tenfoldcv_adrfold_rf_3.Rdata')

# 4
network_tenfoldcv_adrfold_rf_4 = tinycv(y = label, x = X4[, c(9:11, 16:19)], fold = tenfoldadr4, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_4, file = 'network_tenfoldcv_adrfold_rf_4.Rdata')

# 5
network_tenfoldcv_adrfold_rf_5 = tinycv(y = label, x = X5[, c(9:11, 16:19)], fold = tenfoldadr5, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_5, file = 'network_tenfoldcv_adrfold_rf_5.Rdata')

# 6
network_tenfoldcv_adrfold_rf_6 = tinycv(y = label, x = X6[, c(9:11, 16:19)], fold = tenfoldadr6, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_6, file = 'network_tenfoldcv_adrfold_rf_6.Rdata')

# 7
network_tenfoldcv_adrfold_rf_7 = tinycv(y = label, x = X7[, c(9:11, 16:19)], fold = tenfoldadr7, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_7, file = 'network_tenfoldcv_adrfold_rf_7.Rdata')

# 8
network_tenfoldcv_adrfold_rf_8 = tinycv(y = label, x = X8[, c(9:11, 16:19)], fold = tenfoldadr8, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_8, file = 'network_tenfoldcv_adrfold_rf_8.Rdata')

# 9
network_tenfoldcv_adrfold_rf_9 = tinycv(y = label, x = X9[, c(9:11, 16:19)], fold = tenfoldadr9, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_9, file = 'network_tenfoldcv_adrfold_rf_9.Rdata')

# 10
network_tenfoldcv_adrfold_rf_10 = tinycv(y = label, x = X10[, c(9:11, 16:19)], fold = tenfoldadr10, model = 'rf', center = FALSE, scale = FALSE)
save(network_tenfoldcv_adrfold_rf_10, file = 'network_tenfoldcv_adrfold_rf_10.Rdata')



# # All Features 10-fold CV LogReg FoldDrug
#
# # 1
# all_folddrug_tenfoldcv_logreg_1 = tinycv(y = label, x = X1, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_1, file = 'all_folddrug_tenfoldcv_logreg_1.Rdata')
#
# # 2
# all_folddrug_tenfoldcv_logreg_2 = tinycv(y = label, x = X2, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_2, file = 'all_folddrug_tenfoldcv_logreg_2.Rdata')
#
# # 3
# all_folddrug_tenfoldcv_logreg_3 = tinycv(y = label, x = X3, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_3, file = 'all_folddrug_tenfoldcv_logreg_3.Rdata')
#
# # 4
# all_folddrug_tenfoldcv_logreg_4 = tinycv(y = label, x = X4, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_4, file = 'all_folddrug_tenfoldcv_logreg_4.Rdata')
#
# # 5
# all_folddrug_tenfoldcv_logreg_5 = tinycv(y = label, x = X5, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_5, file = 'all_folddrug_tenfoldcv_logreg_5.Rdata')
#
# # 6
# all_folddrug_tenfoldcv_logreg_6 = tinycv(y = label, x = X6, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_6, file = 'all_folddrug_tenfoldcv_logreg_6.Rdata')
#
# # 7
# all_folddrug_tenfoldcv_logreg_7 = tinycv(y = label, x = X7, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_7, file = 'all_folddrug_tenfoldcv_logreg_7.Rdata')
#
# # 8
# all_folddrug_tenfoldcv_logreg_8 = tinycv(y = label, x = X8, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_8, file = 'all_folddrug_tenfoldcv_logreg_8.Rdata')
#
# # 9
# all_folddrug_tenfoldcv_logreg_9 = tinycv(y = label, x = X9, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_9, file = 'all_folddrug_tenfoldcv_logreg_9.Rdata')
#
# # 10
# all_folddrug_tenfoldcv_logreg_10 = tinycv(y = label, x = X10, fold = tenfolddrug, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_folddrug_tenfoldcv_logreg_10, file = 'all_folddrug_tenfoldcv_logreg_10.Rdata')



# All Features 10-fold CV RF FoldDrug

# 1
all_folddrug_tenfoldcv_rf_1 = tinycv(y = label, x = X1, fold = tenfolddrug1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_1, file = 'all_folddrug_tenfoldcv_rf_1.Rdata')

# 2
all_folddrug_tenfoldcv_rf_2 = tinycv(y = label, x = X2, fold = tenfolddrug2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_2, file = 'all_folddrug_tenfoldcv_rf_2.Rdata')

# 3
all_folddrug_tenfoldcv_rf_3 = tinycv(y = label, x = X3, fold = tenfolddrug3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_3, file = 'all_folddrug_tenfoldcv_rf_3.Rdata')

# 4
all_folddrug_tenfoldcv_rf_4 = tinycv(y = label, x = X4, fold = tenfolddrug4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_4, file = 'all_folddrug_tenfoldcv_rf_4.Rdata')

# 5
all_folddrug_tenfoldcv_rf_5 = tinycv(y = label, x = X5, fold = tenfolddrug5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_5, file = 'all_folddrug_tenfoldcv_rf_5.Rdata')

# 6
all_folddrug_tenfoldcv_rf_6 = tinycv(y = label, x = X6, fold = tenfolddrug6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_6, file = 'all_folddrug_tenfoldcv_rf_6.Rdata')

# 7
all_folddrug_tenfoldcv_rf_7 = tinycv(y = label, x = X7, fold = tenfolddrug7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_7, file = 'all_folddrug_tenfoldcv_rf_7.Rdata')

# 8
all_folddrug_tenfoldcv_rf_8 = tinycv(y = label, x = X8, fold = tenfolddrug8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_8, file = 'all_folddrug_tenfoldcv_rf_8.Rdata')

# 9
all_folddrug_tenfoldcv_rf_9 = tinycv(y = label, x = X9, fold = tenfolddrug9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_9, file = 'all_folddrug_tenfoldcv_rf_9.Rdata')

# 10
all_folddrug_tenfoldcv_rf_10 = tinycv(y = label, x = X10, fold = tenfolddrug10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_tenfoldcv_rf_10, file = 'all_folddrug_tenfoldcv_rf_10.Rdata')



# # All Features 10-fold CV LogReg FoldADR
#
# # 1
# all_foldadr_tenfoldcv_logreg_1 = tinycv(y = label, x = X1, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_1, file = 'all_foldadr_tenfoldcv_logreg_1.Rdata')
#
# # 2
# all_foldadr_tenfoldcv_logreg_2 = tinycv(y = label, x = X2, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_2, file = 'all_foldadr_tenfoldcv_logreg_2.Rdata')
#
# # 3
# all_foldadr_tenfoldcv_logreg_3 = tinycv(y = label, x = X3, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_3, file = 'all_foldadr_tenfoldcv_logreg_3.Rdata')
#
# # 4
# all_foldadr_tenfoldcv_logreg_4 = tinycv(y = label, x = X4, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_4, file = 'all_foldadr_tenfoldcv_logreg_4.Rdata')
#
# # 5
# all_foldadr_tenfoldcv_logreg_5 = tinycv(y = label, x = X5, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_5, file = 'all_foldadr_tenfoldcv_logreg_5.Rdata')
#
# # 6
# all_foldadr_tenfoldcv_logreg_6 = tinycv(y = label, x = X6, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_6, file = 'all_foldadr_tenfoldcv_logreg_6.Rdata')
#
# # 7
# all_foldadr_tenfoldcv_logreg_7 = tinycv(y = label, x = X7, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_7, file = 'all_foldadr_tenfoldcv_logreg_7.Rdata')
#
# # 8
# all_foldadr_tenfoldcv_logreg_8 = tinycv(y = label, x = X8, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_8, file = 'all_foldadr_tenfoldcv_logreg_8.Rdata')
#
# # 9
# all_foldadr_tenfoldcv_logreg_9 = tinycv(y = label, x = X9, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_9, file = 'all_foldadr_tenfoldcv_logreg_9.Rdata')
#
# # 10
# all_foldadr_tenfoldcv_logreg_10 = tinycv(y = label, x = X10, fold = tenfoldadr, model = 'logistic', center = TRUE, scale = TRUE)
# save(all_foldadr_tenfoldcv_logreg_10, file = 'all_foldadr_tenfoldcv_logreg_10.Rdata')



# All Features 10-fold CV RF FoldADR

# 1
all_foldadr_tenfoldcv_rf_1 = tinycv(y = label, x = X1, fold = tenfoldadr1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_1, file = 'all_foldadr_tenfoldcv_rf_1.Rdata')

# 2
all_foldadr_tenfoldcv_rf_2 = tinycv(y = label, x = X2, fold = tenfoldadr2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_2, file = 'all_foldadr_tenfoldcv_rf_2.Rdata')

# 3
all_foldadr_tenfoldcv_rf_3 = tinycv(y = label, x = X3, fold = tenfoldadr3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_3, file = 'all_foldadr_tenfoldcv_rf_3.Rdata')

# 4
all_foldadr_tenfoldcv_rf_4 = tinycv(y = label, x = X4, fold = tenfoldadr4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_4, file = 'all_foldadr_tenfoldcv_rf_4.Rdata')

# 5
all_foldadr_tenfoldcv_rf_5 = tinycv(y = label, x = X5, fold = tenfoldadr5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_5, file = 'all_foldadr_tenfoldcv_rf_5.Rdata')

# 6
all_foldadr_tenfoldcv_rf_6 = tinycv(y = label, x = X6, fold = tenfoldadr6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_6, file = 'all_foldadr_tenfoldcv_rf_6.Rdata')

# 7
all_foldadr_tenfoldcv_rf_7 = tinycv(y = label, x = X7, fold = tenfoldadr7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_7, file = 'all_foldadr_tenfoldcv_rf_7.Rdata')

# 8
all_foldadr_tenfoldcv_rf_8 = tinycv(y = label, x = X8, fold = tenfoldadr8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_8, file = 'all_foldadr_tenfoldcv_rf_8.Rdata')

# 9
all_foldadr_tenfoldcv_rf_9 = tinycv(y = label, x = X9, fold = tenfoldadr9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_9, file = 'all_foldadr_tenfoldcv_rf_9.Rdata')

# 10
all_foldadr_tenfoldcv_rf_10 = tinycv(y = label, x = X10, fold = tenfoldadr10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_tenfoldcv_rf_10, file = 'all_foldadr_tenfoldcv_rf_10.Rdata')



# All Features LOOCV RF FoldDrug

# 1
all_folddrug_loocv_rf_1 = tinycv(y = label, x = X1, fold = loofolddrug1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_1, file = 'all_folddrug_loocv_rf_1.Rdata')

# 2
all_folddrug_loocv_rf_2 = tinycv(y = label, x = X2, fold = loofolddrug2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_2, file = 'all_folddrug_loocv_rf_2.Rdata')

# 3
all_folddrug_loocv_rf_3 = tinycv(y = label, x = X3, fold = loofolddrug3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_3, file = 'all_folddrug_loocv_rf_3.Rdata')

# 4
all_folddrug_loocv_rf_4 = tinycv(y = label, x = X4, fold = loofolddrug4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_4, file = 'all_folddrug_loocv_rf_4.Rdata')

# 5
all_folddrug_loocv_rf_5 = tinycv(y = label, x = X5, fold = loofolddrug5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_5, file = 'all_folddrug_loocv_rf_5.Rdata')

# 6
all_folddrug_loocv_rf_6 = tinycv(y = label, x = X6, fold = loofolddrug6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_6, file = 'all_folddrug_loocv_rf_6.Rdata')

# 7
all_folddrug_loocv_rf_7 = tinycv(y = label, x = X7, fold = loofolddrug7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_7, file = 'all_folddrug_loocv_rf_7.Rdata')

# 8
all_folddrug_loocv_rf_8 = tinycv(y = label, x = X8, fold = loofolddrug8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_8, file = 'all_folddrug_loocv_rf_8.Rdata')

# 9
all_folddrug_loocv_rf_9 = tinycv(y = label, x = X9, fold = loofolddrug9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_9, file = 'all_folddrug_loocv_rf_9.Rdata')

# 10
all_folddrug_loocv_rf_10 = tinycv(y = label, x = X10, fold = loofolddrug10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_folddrug_loocv_rf_10, file = 'all_folddrug_loocv_rf_10.Rdata')



# All Features LOOCV RF FoldADR

# 1
all_foldadr_loocv_rf_1 = tinycv(y = label, x = X1, fold = loofoldadr1, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_1, file = 'all_foldadr_loocv_rf_1.Rdata')

# 2
all_foldadr_loocv_rf_2 = tinycv(y = label, x = X2, fold = loofoldadr2, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_2, file = 'all_foldadr_loocv_rf_2.Rdata')

# 3
all_foldadr_loocv_rf_3 = tinycv(y = label, x = X3, fold = loofoldadr3, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_3, file = 'all_foldadr_loocv_rf_3.Rdata')

# 4
all_foldadr_loocv_rf_4 = tinycv(y = label, x = X4, fold = loofoldadr4, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_4, file = 'all_foldadr_loocv_rf_4.Rdata')

# 5
all_foldadr_loocv_rf_5 = tinycv(y = label, x = X5, fold = loofoldadr5, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_5, file = 'all_foldadr_loocv_rf_5.Rdata')

# 6
all_foldadr_loocv_rf_6 = tinycv(y = label, x = X6, fold = loofoldadr6, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_6, file = 'all_foldadr_loocv_rf_6.Rdata')

# 7
all_foldadr_loocv_rf_7 = tinycv(y = label, x = X7, fold = loofoldadr7, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_7, file = 'all_foldadr_loocv_rf_7.Rdata')

# 8
all_foldadr_loocv_rf_8 = tinycv(y = label, x = X8, fold = loofoldadr8, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_8, file = 'all_foldadr_loocv_rf_8.Rdata')

# 9
all_foldadr_loocv_rf_9 = tinycv(y = label, x = X9, fold = loofoldadr9, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_9, file = 'all_foldadr_loocv_rf_9.Rdata')

# 10
all_foldadr_loocv_rf_10 = tinycv(y = label, x = X10, fold = loofoldadr10, model = 'rf', mtry = 7, center = FALSE, scale = FALSE)
save(all_foldadr_loocv_rf_10, file = 'all_foldadr_loocv_rf_10.Rdata')
