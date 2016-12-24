# Multiple Evidence Fusion
#
# Compute PPI-based Drug-Drug Similarity
#
# The following functions require R packages "igraph", "foreach" and "doMC".
# The R data file pro-geneid.Rdata contains the Entrez Gene IDs for the
# drug targets (retrieved from UniProt).
#
# This task was completed using the high-performance computing server of
# CBDD Group, Central South University. It took about 2 days to complete
# with 64 parallel tasks.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 14, 2013

library('igraph')
library('foreach')
library('doMC')
registerDoMC(64)

load('pro-geneid.Rdata')  # uniprot id to gene id mapping

biogrid = read.table(gzfile('BIOGRID-ORGANISM-Homo_sapiens-3.2.104.tab2.tar.gz'),
                     sep = '\t', header = FALSE, fill = TRUE,
                     stringsAsFactors = FALSE, quote = '')

biogrid = biogrid[, 2:3]
biogrid[, 1] = as.character(biogrid[, 1])
biogrid[, 2] = as.character(biogrid[, 2])

g = graph.data.frame(biogrid, directed = FALSE)
nodes = unique(as.vector(as.matrix(biogrid)))

# A * e ^ -(|P1 - P2|)
# |P1 - P2| = shortest path

A = 0.9 * exp(1)

#' Function to compute PPI-based similarity values between two protein
#' lists in parallel
#'
#' @param twoid length-2 integer vector specifying the indexes in a list,
#' whose components stored the Gene IDs of the proteins
#'
#' @return The similarity matrix between to protein sequence lists.

ppiSim = function (twoid) {

  id1 = twoid[1]
  id2 = twoid[2]

  if (all(geneid[[id1]] == '') | all(geneid[[id2]] == '')) {
    mat = matrix(0)
  } else {
    mat = matrix(0L, nrow = length(geneid[[id1]]), ncol = length(geneid[[id2]]))
    for ( i in 1:length(geneid[[id1]]) ) {
      for ( j in 1:length(geneid[[id2]]) ) {
        gid1 = as.character(geneid[[id1]][i])
        gid2 = as.character(geneid[[id2]][j])
        if (gid1 == gid2) {
          mat[i, j] = 1
        } else if ( (gid1 %in% nodes) & (gid2 %in% nodes) ) {
          spath = length(get.shortest.paths(g, from = gid1, to = gid2,
                                            output = 'epath')[[1]])
          mat[i, j] = A * ( exp(1)^(-spath) )
        } else {
          mat[i, j] = 0
        }
      }
    }
  }

  return(mat)

}

# generate lower matrix index
idx = combn(1:length(geneid), 2)

# then use foreach parallelization
# input is all pair combination

ppisimlist = vector('list', ncol(idx))

ppisimlist <- foreach (i = 1:ncol(idx), .errorhandling = 'pass') %dopar% {
  xxx <- ppiSim(rev(idx[, i]))
}

ppisimAvgtotal = sapply(ppisimlist, mean)

# convert list to matrix
ppisimmatAvgtotal = matrix(0, length(geneid), length(geneid))

for (i in 1:length(ppisimlist)) ppisimmatAvgtotal[idx[2, i], idx[1, i]] = ppisimAvgtotal[i]
ppisimmatAvgtotal[upper.tri(ppisimmatAvgtotal)] = t(ppisimmatAvgtotal)[upper.tri(t(ppisimmatAvgtotal))]
diag(ppisimmatAvgtotal) = 1
ppisimmatAvgtotal6digit = format(round(ppisimmatAvgtotal, 6), nsmall = 6)

write.table(ppisimmatAvgtotal6digit, 'ppisimmat.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
