# Multiple Evidence Fusion
#
# AUC/Precision Plots
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: November, 2013

library('ROCR')

setwd('~/MEF/result/')

load('all_atcfold_rf_1.Rdata')
load('all_atcfold_rf_2.Rdata')
load('all_atcfold_rf_3.Rdata')
load('all_atcfold_rf_4.Rdata')
load('all_atcfold_rf_5.Rdata')
load('all_atcfold_rf_6.Rdata')
load('all_atcfold_rf_7.Rdata')
load('all_atcfold_rf_8.Rdata')
load('all_atcfold_rf_9.Rdata')
load('all_atcfold_rf_10.Rdata')

load('all_socfold_rf_1.Rdata')
load('all_socfold_rf_2.Rdata')
load('all_socfold_rf_3.Rdata')
load('all_socfold_rf_4.Rdata')
load('all_socfold_rf_5.Rdata')
load('all_socfold_rf_6.Rdata')
load('all_socfold_rf_7.Rdata')
load('all_socfold_rf_8.Rdata')
load('all_socfold_rf_9.Rdata')
load('all_socfold_rf_10.Rdata')

aucmat_atcfold = matrix(0, nrow = 10, ncol = 14)
aucmat_socfold = matrix(0, nrow = 10, ncol = 24)
precmat_atcfold = matrix(0, nrow = 10, ncol = 14)
precmat_socfold = matrix(0, nrow = 10, ncol = 24)

# compute auc/prec here

for (i in 1:10) {
  for (j in 1:14) {
    eval(parse(text = paste0('pred = prediction( 1 - all_atcfold_rf_', i, '[[', j, ']][, 4], all_atcfold_rf_', i, '[[', j, ']][, 2])')))
    perf1 = performance(pred, 'auc')
    aucmat_atcfold[i, j] = perf1@y.values[[1]]
    perf2 = performance(pred, 'prec')
    k = which.min(abs(perf2@x.values[[1]] - 0.5))
    precmat_atcfold[i, j] = perf2@y.values[[1]][k]
  }
}

for (i in 1:10) {
  for (j in 1:24) {
    eval(parse(text = paste0('pred = prediction( 1 - all_socfold_rf_', i, '[[', j, ']][, 4], all_socfold_rf_', i, '[[', j, ']][, 2])')))
    perf1 = performance(pred, 'auc')
    aucmat_socfold[i, j] = perf1@y.values[[1]]
    perf2 = performance(pred, 'prec')
    k = which.min(abs(perf2@x.values[[1]] - 0.5))
    precmat_socfold[i, j] = perf2@y.values[[1]][k]
  }
}

# auc/prec barplot for 14 ATC levels + 24 SOC levels with error bars

library('ggplot2')

# get ATC and SOC fold names

# to get the fold level names here, have to re-make the fold ...
# we assume that the 10 fold level orders are the same, so use one is enough
loofolddrug = as.factor(read.table('tassociation1', sep = '\t', header = FALSE)[, 1])
loofoldadr = as.factor(read.table('tassociation1', sep = '\t', header = FALSE)[, 2])

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
atcfold = rep(NA, 49606)
socfold = rep(NA, 49606)
# look up each drug in atc dict, randomly sample one atc code, get first char
for (i in 1:49606) atcfold[i] = substr(sample(atclist[[as.integer(loofolddrug[i])]], 1), 1, 1)
for (i in 1:49606) socfold[i] = sample(soclist[[as.integer(loofoldadr[i])]], 1)

atcfold = as.factor(atcfold)
atcfolddict = levels(atcfold)

socfold = as.factor(socfold)
socfolddict = levels(socfold)

atccode = read.table('atc-code.txt', sep = '\t', header = TRUE)
soccode = read.table('soc-code.txt', sep = '\t', header = TRUE)

atcnames = rep(NA, 14)
socnames = rep(NA, 24)
for ( i in 1:14 ) atcnames[i] = as.character(atcfolddict[i])
for ( i in 1:24 ) socnames[i] = as.character(soccode[which(soccode[, 1] == socfolddict[i]), 3])

colnames(precmat_atcfold) = atcnames
colnames(precmat_socfold) = socnames



df1 = data.frame(
  'ATC.Group' = factor(atcnames, levels = atcnames),
  'AUC Score' = colMeans(aucmat_atcfold),
  se = apply(aucmat_atcfold, FUN = sd, MARGIN = 2)
)

df2 = data.frame(
  'SOC.Group' = factor(socnames, levels = socnames),
  'AUC Score' = colMeans(aucmat_socfold),
  se = apply(aucmat_socfold, FUN = sd, MARGIN = 2)
)

df3 = data.frame(
  'ATC.Group' = factor(atcnames, levels = atcnames),
  'Precision' = colMeans(precmat_atcfold),
  se = apply(precmat_atcfold, FUN = sd, MARGIN = 2)
)

df4 = data.frame(
  'SOC.Group' = factor(socnames, levels = socnames),
  'Precision' = colMeans(precmat_socfold),
  se = apply(precmat_socfold, FUN = sd, MARGIN = 2)
)

require(RColorBrewer)
pal24 = c(brewer.pal(9, 'Set1'), setdiff(brewer.pal(8, 'Set2'), "#FFD92F"), setdiff(brewer.pal(11, 'Set3'), c("#FFFFB3", "#FDB462", '#D9D9D9')))

# define the top and bottom of the errorbars
limits = aes(ymax = AUC.Score + se, ymin = AUC.Score - se)

p1 = ggplot(df1, aes(fill = ATC.Group, y = AUC.Score, x = ATC.Group)) + coord_cartesian(ylim = c(0.9, 1))
dodge = position_dodge(width = 0.9)
p1 + geom_bar(position = dodge) + scale_fill_manual(values = pal24) + geom_errorbar(limits, position = dodge, width = 0.25)

ggsave(file = 'atc_auc_barplot.pdf', scale = 2.5)
ggsave(file = 'atc_auc_barplot.png', scale = 2.5)

p2 = ggplot(df2, aes(fill = SOC.Group, y = AUC.Score, x = SOC.Group)) + coord_cartesian(ylim = c(0.9, 1))
dodge = position_dodge(width = 0.9)
p2 + geom_bar(position = dodge) + scale_fill_manual(values = pal24) + geom_errorbar(limits, position = dodge, width = 0.25)

ggsave(file = 'soc_auc_barplot.pdf', scale = 2.5)
ggsave(file = 'soc_auc_barplot.png', scale = 2.5)

# Define the top and bottom of the errorbars
limits = aes(ymax = Precision + se, ymin = Precision - se)

p3 = ggplot(df3, aes(fill = ATC.Group, y = Precision, x = ATC.Group)) + coord_cartesian(ylim = c(0.7, 1))
dodge = position_dodge(width = 0.9)
p3 + geom_bar(position = dodge) + scale_fill_manual(values = pal24) + geom_errorbar(limits, position = dodge, width = 0.25)

ggsave(file = 'atc_prec_barplot.pdf', scale = 2.5)
ggsave(file = 'atc_prec_barplot.png', scale = 2.5)

p4 = ggplot(df4, aes(fill = SOC.Group, y = Precision, x = SOC.Group)) + coord_cartesian(ylim = c(0.7, 1))
dodge = position_dodge(width = 0.9)
p4 + geom_bar(position = dodge) + scale_fill_manual(values = pal24) + geom_errorbar(limits, position = dodge, width = 0.25)

ggsave(file = 'soc_prec_barplot.pdf', scale = 2.5)
ggsave(file = 'soc_prec_barplot.png', scale = 2.5)



# write csv

AUC.Score = colMeans(aucmat_atcfold)
se = apply(aucmat_atcfold, FUN = sd, MARGIN = 2)
AUC.Score = format(round(AUC.Score, 4), nsmall = 4)
se = format(round(se, 4), nsmall = 4)

x = data.frame(
  ATC.Group = atcnames,
  AUC.Score = paste(AUC.Score, '+/-', se)
)

write.csv(x, 'atc-auc.csv', row.names = FALSE, quote = FALSE)

AUC.Score = colMeans(aucmat_socfold)
se = apply(aucmat_socfold, FUN = sd, MARGIN = 2)
AUC.Score = format(round(AUC.Score, 4), nsmall = 4)
se = format(round(se, 4), nsmall = 4)

x = data.frame(
  SOC.Group = socnames,
  AUC.Score = paste(AUC.Score, '+/-', se)
)

write.csv(x, 'soc-auc.csv', row.names = FALSE, quote = FALSE)

Precision = colMeans(precmat_atcfold)
se = apply(precmat_atcfold, FUN = sd, MARGIN = 2)
Precision = format(round(Precision, 4), nsmall = 4)
se = format(round(se, 4), nsmall = 4)

x = data.frame(
  ATC.Group = atcnames,
  Precision = paste(Precision, '+/-', se)
)

write.csv(x, 'atc-prec.csv', row.names = FALSE, quote = FALSE)

Precision = colMeans(precmat_socfold)
se = apply(precmat_socfold, FUN = sd, MARGIN = 2)
Precision = format(round(Precision, 4), nsmall = 4)
se = format(round(se, 4), nsmall = 4)

x = data.frame(
  SOC.Group = socnames,
  Precision = paste(Precision, '+/-', se)
)

write.csv(x, 'soc-prec.csv', row.names = FALSE, quote = FALSE)
