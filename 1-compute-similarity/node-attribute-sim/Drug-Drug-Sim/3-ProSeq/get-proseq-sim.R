# Multiple Evidence Fusion
#
# Compute Protein Sequence Based Drug-Drug Similarity
#
# The following functions require R packages "foreach", "doMC" and "Biostrings".
# The R data file proseq.RData contains the protein sequences for the
# drug targets (retrieved from UniProt).
#
# This task was completed using the high-performance computing server of
# CBDD Group, Central South University. It took about 7 days to complete
# with 64 parallel tasks.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 05, 2013

load('proseq.RData')

library('foreach')
library('doMC')
library('Biostrings')
registerDoMC(64)
data(BLOSUM62)

#' Function to compute sequence alignment based similarity values
#' between two protein lists in parallel
#'
#' @param twoid length-2 integer vector specifying the indexes in a list,
#' whose components stored the protein sequences
#'
#' @return The similarity matrix between to protein sequence lists.

getPairSim = function (twoid) {

  id1 = twoid[1]
  id2 = twoid[2]

  n2 = length(protmap[[id2]])

  if (all(protmap[[id1]] == '') | all(protmap[[id2]] == '')) {

    mat = matrix(0)

  } else {

    id1good = which(protmap[[id1]] != '')
    id2good = which(protmap[[id2]] != '')

    mat = matrix(0L, nrow = length(id1good), ncol = length(id2good))

    for (i in 1:length(id1good)) {
      for (j in 1:length(id2good)) {

        s1 = try(AAString(protmap[[id1]][id1good][i]), silent = TRUE)
        s2 = try(AAString(protmap[[id2]][id2good][j]), silent = TRUE)
        s12 = try(pairwiseAlignment(s1, s2, type = "local",
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE),
                  silent = TRUE)
        s11 = try(pairwiseAlignment(s1, s1, type = "local",
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE),
                  silent = TRUE)
        s22 = try(pairwiseAlignment(s2, s2, type = "local",
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE),
                  silent = TRUE)

        if ( is.numeric(s12) == FALSE | is.numeric(s11) == FALSE | is.numeric(s22) == FALSE ) {
          mat[i, j] = 0.0
        } else if ( abs(s11) < .Machine$double.eps | abs(s22) < .Machine$double.eps ) {
          mat[i, j] = 0.0
        } else {
          mat[i, j] = s12/sqrt(s11 * s22)
        }

      }
    }

  }

  return(mat)

}

# generate lower matrix index
idx = combn(1:length(protmap), 2)

# then use foreach parallelization
# the input is all pair combination
seqsimlist = vector('list', ncol(idx))

seqsimlist <- foreach (i = 1:length(seqsimlist), .errorhandling = 'pass') %dopar% {
  xxx <- getPairSim(rev(idx[, i]))
}

seqsimMaxtotal = sapply(seqsimlist, max)
seqsimAvgtotal = sapply(seqsimlist, mean)
seqsimAvgmax = sapply(seqsimlist, function(x) mean(c(apply(x, 1, max), apply(x, 2, max))) )

# convert list to matrix
seqsimmatMaxtotal = matrix(0, length(protmap), length(protmap))
seqsimmatAvgtotal = matrix(0, length(protmap), length(protmap))
seqsimmatAvgmax = matrix(0, length(protmap), length(protmap))

for (i in 1:length(seqsimlist)) {
  seqsimmatMaxtotal[idx[2, i], idx[1, i]] = seqsimMaxtotal[i]
}

for (i in 1:length(seqsimlist)) {
  seqsimmatAvgtotal[idx[2, i], idx[1, i]] = seqsimAvgtotal[i]
}

for (i in 1:length(seqsimlist)) {
  seqsimmatAvgmax[idx[2, i], idx[1, i]] = seqsimAvgmax[i]
}

seqsimmatMaxtotal[upper.tri(seqsimmatMaxtotal)] = t(seqsimmatMaxtotal)[upper.tri(t(seqsimmatMaxtotal))]
seqsimmatAvgtotal[upper.tri(seqsimmatAvgtotal)] = t(seqsimmatAvgtotal)[upper.tri(t(seqsimmatAvgtotal))]
seqsimmatAvgmax[upper.tri(seqsimmatAvgmax)] = t(seqsimmatAvgmax)[upper.tri(t(seqsimmatAvgmax))]

diag(seqsimmatMaxtotal) = rep(1, length(seqsimlist))
diag(seqsimmatAvgtotal) = rep(1, length(seqsimlist))
diag(seqsimmatAvgmax) = rep(1, length(seqsimlist))

seqsimmatMaxtotal4digit = format(round(seqsimmatMaxtotal, 4), nsmall = 4)
seqsimmatAvgtotal4digit = format(round(seqsimmatAvgtotal, 4), nsmall = 4)
seqsimmatAvgmax4digit = format(round(seqsimmatAvgmax, 4), nsmall = 4)

write.table(seqsimmatMaxtotal4digit, 'seqsimmat-maxtotal.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(seqsimmatAvgtotal4digit, 'seqsimmat-avgtotal.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(seqsimmatAvgmax4digit, 'seqsimmat-avgmax.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
