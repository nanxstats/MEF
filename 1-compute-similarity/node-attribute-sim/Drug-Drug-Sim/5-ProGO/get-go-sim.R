# Multiple Evidence Fusion
#
# Compute GO-based Drug-Drug Similarity
#
# The following functions require R packages "GOSemSim" and "parallel".
# The R data file pro-geneid.Rdata contains the Entrez Gene IDs for the
# drug targets (retrieved from UniProt).
#
# This task was completed using the high-performance computing server of
# CBDD Group, Central South University. It took about 2 days to complete
# with 64 parallel tasks.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 01, 2013

library('GOSemSim')
library('parallel')

load('pro-geneid.Rdata')  # uniprot id to gene id mapping

#' Function to compute the pairwise GO semantic similarity values
#' between two drug targets in parallel
#'
#' @param twoid length-2 integer vector specifying the indexes in a list,
#' whose components stored the Entrez Gene IDs for each drug target
#' @param goont subontologies: "MF", "BP", and "CC"
#' @param gomeasure similarity measure: "Lin"
#'
#' @return This function will write the similarity values to a txt file,
#' which would be loaded and transformed into a similarity matrix later.

goPairSim = function (twoid, goont, gomeasure) {

  id1 = twoid[1]
  id2 = twoid[2]

  if (all(geneid[[id1]] == '') | all(geneid[[id2]] == '')) {

    mat = matrix(0)

  } else {

    id1good = 1:length(geneid[[id1]])
    id2good = 1:length(geneid[[id2]])

    mat = matrix(0L, nrow = length(id1good), ncol = length(id2good))

    for (i in 1:length(id1good)) {
      for (j in 1:length(id2good)) {

        gid1 = as.character(geneid[[id1]][id1good][i])
        gid2 = as.character(geneid[[id2]][id2good][j])

        res = try(suppressWarnings(geneSim(gid1, gid2, organism = 'human',
                                           combine = 'BMA', ont = goont,
                                           measure = gomeasure)), silent = TRUE)

        if ( is.list(res) ) {
          if ( is.numeric( res$geneSim ) ) {  # eliminate the NA's
            mat[i, j] = res$geneSim
          } else {
            mat[i, j] = 0
          }
        } else {
          mat[i, j] = 0
        }

      }
    }

    mat[is.na(mat)] = 0  # to excessively eliminate NA's

  }

  write(paste(c(paste(rev(twoid), collapse = ','),
                as.character(mean(mat))), collapse = ','),
        paste0(goont, '-', gomeasure, '-go-simmat.txt'), append = TRUE)

}

# calculate GO-based similarity (MF)

theont = 'MF'
themeasure = 'Lin'

idx = combn(1:length(geneid), 2)  # generate lower matrix index

cl = makeCluster(64)

clusterEvalQ(cl, {
  require('GOSemSim')
  load('pro-geneid.Rdata')
  bigdrug = which(sapply(geneid, length) > 10.5)

  for (i in bigdrug) {
    set.seed(1001)
    geneid[[i]] = sample(geneid[[i]], 10L)
  }

  idx = combn(1:length(geneid), 2)
})

parallel::parCapply(cl = cl, x = idx, FUN = goPairSim, theont, themeasure)
stopCluster(cl)

# calculate GO-based similarity (BP)

theont = 'BP'
themeasure = 'Lin'

idx = combn(1:length(geneid), 2)  # generate lower matrix index

cl = makeCluster(64)

clusterEvalQ(cl, {
  require('GOSemSim')
  load('pro-geneid.Rdata')
  bigdrug = which(sapply(geneid, length) > 3.5)

  for (i in bigdrug) {
    set.seed(1001)
    geneid[[i]] = sample(geneid[[i]], 3L)
  }

  idx = combn(1:length(geneid), 2)
})

parallel::parCapply(cl = cl, x = idx, FUN = goPairSim, theont, themeasure)
stopCluster(cl)

# calculate GO-based similarity (CC)

theont = 'CC'
themeasure = 'Lin'

idx = combn(1:length(geneid), 2)  # generate lower matrix index

cl = makeCluster(64)

clusterEvalQ(cl, {
  require('GOSemSim')
  load('pro-geneid.Rdata')
  bigdrug = which(sapply(geneid, length) > 10.5)

  for (i in bigdrug) {
    set.seed(1001)
    geneid[[i]] = sample(geneid[[i]], 10L)
  }

  idx = combn(1:length(geneid), 2)
})

parallel::parCapply(cl = cl, x = idx, FUN = goPairSim, theont, themeasure)
stopCluster(cl)

# load the output result and save as matrices

# MF

txtavgtotal = read.csv(paste0(getwd(), 'MF-Lin-go-simmat.txt'), header = FALSE)
gosimmatAvgtotal = matrix(0, 746, 746)

for (i in 1:nrow(txtavgtotal)) {
  gosimmatAvgtotal[txtavgtotal[i, 1], txtavgtotal[i, 2]] = txtavgtotal[i, 3]
}

gosimmatAvgtotal[upper.tri(gosimmatAvgtotal)] =
  t(gosimmatAvgtotal)[upper.tri(t(gosimmatAvgtotal))]
diag(gosimmatAvgtotal) = 1
gosimmatAvgtotal4digit = format(round(gosimmatAvgtotal, 4), nsmall = 4)
write.table(gosimmatAvgtotal4digit, paste0(getwd(), 'gosimmat-MF-Lin.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# BP

txtavgtotal = read.csv(paste0(getwd(), 'BP-Lin-go-simmat.txt'), header = FALSE)
gosimmatAvgtotal = matrix(0, 746, 746)

for (i in 1:nrow(txtavgtotal)) {
  gosimmatAvgtotal[txtavgtotal[i, 1], txtavgtotal[i, 2]] = txtavgtotal[i, 3]
}

gosimmatAvgtotal[upper.tri(gosimmatAvgtotal)] =
  t(gosimmatAvgtotal)[upper.tri(t(gosimmatAvgtotal))]
diag(gosimmatAvgtotal) = 1
gosimmatAvgtotal4digit = format(round(gosimmatAvgtotal, 4), nsmall = 4)
write.table(gosimmatAvgtotal4digit, paste0(getwd(), 'gosimmat-BP-Lin.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# CC

txtavgtotal = read.csv(paste0(getwd(), 'CC-Lin-go-simmat.txt'), header = FALSE)
gosimmatAvgtotal = matrix(0, 746, 746)

for (i in 1:nrow(txtavgtotal)) {
  gosimmatAvgtotal[txtavgtotal[i, 1], txtavgtotal[i, 2]] = txtavgtotal[i, 3]
}

gosimmatAvgtotal[upper.tri(gosimmatAvgtotal)] =
  t(gosimmatAvgtotal)[upper.tri(t(gosimmatAvgtotal))]
diag(gosimmatAvgtotal) = 1
gosimmatAvgtotal4digit = format(round(gosimmatAvgtotal, 4), nsmall = 4)
write.table(gosimmatAvgtotal4digit, paste0(getwd(), 'gosimmat-CC-Lin.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
