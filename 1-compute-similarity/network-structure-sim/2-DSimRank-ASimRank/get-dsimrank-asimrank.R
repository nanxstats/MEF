# Multiple Evidence Fusion
#
# Compute Drug-Based and ADR-Based Simrank Similarity
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Aug 23, 2013

adradjmat = as.matrix(read.table('AdrAdjMatrix.txt', header = FALSE, sep = '\t'))
drugadjmat = as.matrix(read.table('DrugAdjMatrix.txt', header = FALSE, sep = '\t'))
adjmat = as.matrix(read.table('AdjMatrix.txt', header = FALSE, sep = '\t'))

library('compiler')  # enable JIT for getsim()
library('foreach')
library('doMC')
registerDoMC(64)

getsim = function (adjmat, simmatold, node, C) {

  Ia = sum(adjmat[, node[1L] ])
  Ib = sum(adjmat[, node[2L] ])

  if ( Ia != 0L & Ib != 0L ) {

    idxi = which( adjmat[, node[1]] == 1L )
    idxj = which( adjmat[, node[2]] == 1L )
    sumij = sum(simmatold[as.matrix(expand.grid(idxi, idxj))])

    rab = (C * sumij) / (Ia * Ib)

  } else {

    rab = 0

  }

  return(rab)

}

getsim = cmpfun(getsim)

parasim = function (adjmat, C = 0.8, maxiter = 100, eps = 1e-6, silent = FALSE) {

  simmatold = diag(nrow(adjmat))
  simmat = diag(9999L, nrow(adjmat))  # to ensure there is a simmat when i = 1

  idx = combn(1:nrow(adjmat), 2)

  for ( i in 1:maxiter ) {

    error = max(abs(simmat - simmatold))
    if (!silent) print(paste('Iter:', i, 'Error:', format(round(error, 8), nsmall = 8)))

    if ( all(abs(simmat - simmatold) <= eps) == TRUE ) {

      break

    } else {

      if ( i == 1L ) simmat = diag(nrow(adjmat))

      simmatold = simmat

      simmat = diag(nrow(adjmat))

      simlist <- foreach(k = 1:ncol(idx), .errorhandling = 'pass') %dopar% {
        xxx <- getsim(adjmat, simmatold, rev(idx[, k]), C)
      }

      for (j in 1:length(simlist)) {
        simmat[idx[2, j], idx[1, j]] = simlist[[j]]
      }

      simmat[upper.tri(simmat)] = t(simmat)[upper.tri(t(simmat))]

    }

  }

  return(simmat)

}

adrres10 = parasim(adradjmat, C = 0.8, maxiter = 10, eps = 1e-4, silent = FALSE)
adrres4digit10 = format(round(adrres10, 4), nsmall = 4)
write.table(adrres4digit10, 'adr-simranksimmat-10.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

drugres10 = parasim(drugadjmat, C = 0.8, maxiter = 10, eps = 1e-4, silent = FALSE)
drugres4digit10 = format(round(drugres10, 4), nsmall = 4)
write.table(drugres4digit10, 'drug-simranksimmat-10.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
