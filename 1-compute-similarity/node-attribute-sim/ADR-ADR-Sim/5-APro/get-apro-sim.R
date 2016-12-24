# Multiple Evidence Fusion
#
# Compute APro-based ADR-ADR Similarity
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Oct 1, 2013

adr = read.csv('ADR.txt', header = FALSE, stringsAsFactors = FALSE)

# ADR UMLS code -- protein ID data that has been curated by removing
# status == 'wrong', 'false positive' and 'opposite phenotype' cases.
# Original data from supplementary material 2 of Kuhn, et al.
# 'Systematic identification of proteins that elicit drug side effects'.
predrel = read.csv('se-protein.csv', header = TRUE, stringsAsFactors = FALSE)

high = as.character(unique(predrel$UMLS.code.of.SE))
real = as.character(unique(adr$V1))
length(intersect(high, real))  # 377 in common

predrelunique = aggregate(predrel[-1], by = list(predrel$UMLS.code.of.SE), c)

# remove @inh @act

predrellist = vector('list', nrow(predrelunique))
for ( i in 1:length(predrellist) ) {
  predrellist[[i]] = strsplit(paste(unlist(predrelunique[i, 2]),
                                    collapse = ' '), ' ')[[1]]
}
names(predrellist) = as.character(predrelunique[, 1])

predrellist1 = vector('list', nrow(predrelunique))
for ( i in 1:length(predrellist) ) {
  predrellist1[[i]] = unique(gsub('@act', '', gsub('@inh', '', predrellist[[i]])))
}
names(predrellist1) = as.character(predrelunique[, 1])

actinhdiff = data.frame(V1 = sapply(predrellist, length),
                        V2 = sapply(predrellist1, length),
                        V3 = sapply(predrellist, length) -
                          sapply(predrellist1, length))

mat = matrix(0.0, ncol = 817L, nrow = 817L)

for ( i in 1:817 ) {
  tmp1 = which(adr[i, ] == names(predrellist1))
  for ( j in i:817 ) {
    tmp2 = which(adr[j, ] == names(predrellist1))
    if ( length(tmp1) == 0L | length(tmp2) == 0L ) {
      mat[i, j] = 0
    } else {
      mat[i, j] = length(intersect(predrellist1[[i]], predrellist1[[j]]))/
        length(union(predrellist1[[i]], predrellist1[[j]]))
    }
  }
}

mat[lower.tri(mat)] = t(mat)[lower.tri(t(mat))]
diag(mat) = 1

mat4digit = format(round(mat, 4), nsmall = 4)
write.table(mat4digit, 'aprosim.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
