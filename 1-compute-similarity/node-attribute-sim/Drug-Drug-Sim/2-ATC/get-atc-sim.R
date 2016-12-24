# Multiple Evidence Fusion
#
# Compute ATC-Based Drug-Drug Similarity
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Aug 23, 2013

x = scan('drugbank-atc.csv', what = 'character', sep = ',')

# CID (drugs) position
cid = grep('CID', x)

atclist = vector('list', length(cid))
names(atclist) = x[cid]

# fill the ATC code(s) for each drug
for (i in 1:(length(cid) - 1)) atclist[[i]] = x[(cid[i] + 1):(cid[i+1] - 1)]
atclist[[length(cid)]] = x[(cid[length(cid)] + 1):length(x)]

# reduce all ATC codes to a single vector
atcvec = unlist(atclist)
names(atcvec) = NULL

# split the vector for further extraction
atcsplit = strsplit(atcvec, split = '')

# combine the second and fifth level number
for (i in 1:length(atcsplit)) atcsplit[[i]][2] = paste0(atcsplit[[i]][2], atcsplit[[i]][3])
for (i in 1:length(atcsplit)) atcsplit[[i]][6] = paste0(atcsplit[[i]][6], atcsplit[[i]][7])
for (i in 1:length(atcsplit)) atcsplit[[i]] = atcsplit[[i]][c(-3, -7)]

# cumulative combine level number
for (i in 1:length(atcsplit)) atcsplit[[i]][2] = paste0(atcsplit[[i]][1], atcsplit[[i]][2])
for (i in 1:length(atcsplit)) atcsplit[[i]][3] = paste0(atcsplit[[i]][2], atcsplit[[i]][3])
for (i in 1:length(atcsplit)) atcsplit[[i]][4] = paste0(atcsplit[[i]][3], atcsplit[[i]][4])
for (i in 1:length(atcsplit)) atcsplit[[i]][5] = paste0(atcsplit[[i]][4], atcsplit[[i]][5])

# construct level 1 list
lvl1code = sapply(atcsplit, '[[', 1)
lvl1codeunique = unique(lvl1code)
lvl1list = vector('list', length(lvl1codeunique))
names(lvl1list) = lvl1codeunique
for (i in 1:length(lvl1list)) lvl1list[[i]] = atcvec[which(lvl1code == lvl1codeunique[i])]

# construct level 2 list
lvl2code = sapply(atcsplit, '[[', 2)
lvl2codeunique = unique(lvl2code)
lvl2list = vector('list', length(lvl2codeunique))
names(lvl2list) = lvl2codeunique
for (i in 1:length(lvl2list)) lvl2list[[i]] = atcvec[which(lvl2code == lvl2codeunique[i])]

# construct level 3 list
lvl3code = sapply(atcsplit, '[[', 3)
lvl3codeunique = unique(lvl3code)
lvl3list = vector('list', length(lvl3codeunique))
names(lvl3list) = lvl3codeunique
for (i in 1:length(lvl3list)) lvl3list[[i]] = atcvec[which(lvl3code == lvl3codeunique[i])]

# construct level 4 list
lvl4code = sapply(atcsplit, '[[', 4)
lvl4codeunique = unique(lvl4code)
lvl4list = vector('list', length(lvl4codeunique))
names(lvl4list) = lvl4codeunique
for (i in 1:length(lvl4list)) lvl4list[[i]] = atcvec[which(lvl4code == lvl4codeunique[i])]

# compute information content for each level

atctotal = length(atcvec)

lvl1info = -log(sapply(lvl1list, length)/atctotal)
lvl2info = -log(sapply(lvl2list, length)/atctotal)
lvl3info = -log(sapply(lvl3list, length)/atctotal)
lvl4info = -log(sapply(lvl4list, length)/atctotal)

#' Function to find most informative ancestor and return the IC value
#'
#' @param atc1 ATC code 1
#' @param atc2 ATC code 2
#'
#' @return The IC value between two ATC codes.

atc2sim = function (atc1, atc2) {

  atc1split = strsplit(atc1, '')[[1]]
  atc2split = strsplit(atc2, '')[[1]]
  # combine the second and fifth level number
  atc1split[2] = paste0(atc1split[2], atc1split[3])
  atc1split[6] = paste0(atc1split[6], atc1split[7])
  atc1split = atc1split[c(-3, -7)]
  atc2split[2] = paste0(atc2split[2], atc2split[3])
  atc2split[6] = paste0(atc2split[6], atc2split[7])
  atc2split = atc2split[c(-3, -7)]
  lvlcompare = (atc1split == atc2split)[-5] == TRUE
  if (lvlcompare[1] == FALSE) {
    lvlnum = 0L
  } else if (lvlcompare[2] == FALSE) {
    lvlnum = 1L
  } else if (lvlcompare[3] == FALSE) {
    lvlnum = 2L
  } else if (lvlcompare[4] == FALSE) {
    lvlnum = 3L
  } else {
    lvlnum = 4L
  }

  # Retrieve the IC value based on lvlnum and atc
  if (lvlnum == 0L) {
    sim = 0  # if in level1 and still not match, then similarity should be 0
  } else if (lvlnum == 1L) {
    sim = lvl1info[atc1split[1]]
  } else if (lvlnum == 2L) {
    sim = lvl2info[paste0(atc1split[1], atc1split[2])]
  } else if (lvlnum == 3L) {
    sim = lvl3info[paste0(atc1split[1], atc1split[2], atc1split[3])]
  } else if (lvlnum == 4L) {
    sim = lvl4info[paste0(atc1split[1], atc1split[2], atc1split[3], atc1split[4])]
  }
  names(sim) = NULL

  return(sim)

}

#' Function to calculate the similarity value between two CIDs
#'
#' @param cid1 CID 1
#' @param cid2 CID 2
#'
#' @return The similarity value between two CIDs
#' (each CID could corresponse to more than one ATC codes).

cid2sim = function (cid1, cid2) {

  cid1atcs = atclist[[cid1]]
  cid2atcs = atclist[[cid2]]
  allcombn = expand.grid(cid1atcs, cid2atcs)

  sims = rep(NA, nrow(allcombn))

  for (i in 1:nrow(allcombn)) {

    sims[i] = atc2sim(atc1 = as.character(allcombn[i, 1]),
                      atc2 = as.character(allcombn[i, 2]))

  }

  return(max(sims))

}

#' Function to calculate the similarity matrix for several CIDs
#'
#' cids - the CID vector
#'
#' @return The similarity matrix for several CIDs
#' (each CID could corresponse to more than one ATC codes).

cidSimMatrix = function (cids) {

  d = length(cids)
  if (length(cids) < 1.5) stop('Must provide more than one cid')
  simmat = diag(d)
  # for diagnal elements, set the similarity = its IC value
  maxic = -log(1/atctotal)
  diag(simmat) = rep(maxic, d)
  # Generate lower matrix index
  idxmat = combn(1:d, 2)

  for (i in 1:ncol(idxmat)) {
    simmat[idxmat[2, i], idxmat[1, i]] = cid2sim(cids[idxmat[2, i]],
                                                 cids[idxmat[1, i]])
  }

  simmat[upper.tri(simmat)] = t(simmat)[upper.tri(t(simmat))]

  return(simmat/maxic)

}

cids = x[cid]
atcsimmat = cidSimMatrix(cids)

atcsimmat4digit = format(round(atcsimmat, 4), nsmall = 4)
write.table(atcsimmat4digit, 'atcsimmat.txt', sep = '\t',
            quote = FALSE, col.names = FALSE, row.names = FALSE)
