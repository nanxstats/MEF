# Multiple Evidence Fusion
#
# Compute Connectivity Map (cMap) based Drug-Drug Similarity
#
# The following functions require R packages "cMap2data" and "proxy",
# which could be acquired from Bioconductor and CRAN, respectively.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 24, 2013

library('cMap2data')
library('proxy')

x = read.table('drugkegg.txt', sep = '\t')
data(druglabels)

drugnames = tolower(as.character(x$V4))

intersect( drugnames, tolower(colnames(drugRL) ))
# total matched drugs: 343

data(drugRL)
colnames(drugRL) = tolower(colnames(drugRL))

featuremat = matrix(0, nrow = nrow(drugRL), ncol = length(drugnames))

for ( i in 1:length(drugnames)) {
  if ( drugnames[i] %in% colnames(drugRL) ) {
    featuremat[, i] = drugRL[, which(drugnames[i] == colnames(drugRL))]
  }
}

featuremat = featuremat[c(1:250, (nrow(featuremat) - 250 + 1):(nrow(featuremat))), ]
featuremat = t(featuremat)

# cMap profile for 746 drugs
write.table(featuremat, 'cmap-profile.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# calculate similarity matrix

geneexpsimmat = as.matrix(proxy::simil(featuremat, method = 'cosine'))
diag(geneexpsimmat) = 1

# replace 0/0=1 with 0
geneexpsimmat[which(geneexpsimmat == 1)] = 0
diag(geneexpsimmat) = 1

mat4digit = format(round(geneexpsimmat, 4), nsmall = 4)

write.table(mat4digit, 'cmapsimmat.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
