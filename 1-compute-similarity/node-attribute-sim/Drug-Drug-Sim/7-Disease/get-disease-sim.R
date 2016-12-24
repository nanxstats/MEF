# Multiple Evidence Fusion
#
# Compute Disease-based Drug-Drug Similarity
#
# The following functions require R packages "foreach" and "doMC".
#
# This task was completed using the high-performance computing server of
# CBDD Group, Central South University. It took about several hours to
# complete with 64 parallel tasks.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 22, 2013

library('foreach')
library('doMC')
registerDoMC(64)

dbid = as.character(read.table('drugkegg.txt', sep = '\t')[, 2])
chem = read.csv('CTD_chemicals.csv', comment.char = "#")
dise = read.csv('CTD_chemicals_diseases.csv', comment.char = "#")
omim = read.table('MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.mat',
                  sep = '\t', row.names = 1, header = FALSE)
colnames(omim) = row.names(omim)

chem = chem[, c(2, 9)]
dise = dise[, c(2, 9)]

chem = chem[which(chem[, 2] != ''), ]
chem[, 1] = gsub('MESH:', '', chem[, 1])
dise = dise[which(dise[, 2] != ''), ]

# expand duplicated DrugBank IDs
chemtmp = vector('list', nrow(chem))
for (i in 1:nrow(chem)) chemtmp[[i]] = as.character(chem[i, 2])
names(chemtmp) = as.character(chem[, 1])
chemtmp = sapply(chemtmp, strsplit, split = '\\|')

tmp = unlist(chemtmp)
chem2 = data.frame(DrugBankIDs = tmp, ChemicalID = names(tmp))
chem2[, 2] = substring(chem2[, 2], 1, 7)

# expand duplicated disease OMIM IDs
disetmp = vector('list', nrow(dise))
for (i in 1:nrow(dise)) disetmp[[i]] = as.character(dise[i, 2])
names(disetmp) = as.character(dise[, 1])
disetmp = sapply(disetmp, strsplit, split = '\\|')
tmp = unlist(disetmp)
dise2 = data.frame(OmimIDs = tmp, ChemicalID = names(tmp))
dise2[, 2] = substring(dise2[, 2], 1, 7)
dise2[, 1] = as.character(dise2[, 1])

# map from DrugBank id to mesh id
dbmesh = rep('', 746)
for (i in 1:746) try(dbmesh[i] <- chem2[which(dbid[i] == chem2[, 1]), 2], silent = TRUE)

# map from mesh id to OMIM ids
meshomim = vector('list', 746)
names(meshomim) = dbid

for (i in 1:746) {
  if ( dbmesh[i] != '' ) {
    try(meshomim[[i]] <- as.character(dise2[which(dbmesh[i] == dise2[, 2]), 1]),
        silent = TRUE)
  }
}

for (i in 1:746) {
  if ( length(meshomim[[i]]) > 0.5 ) {
    meshomim[[i]] = meshomim[[i]][ ( meshomim[[i]] %in% colnames(omim) )]
  }
}

# calculate similarity matrix
idx = combn(1:length(meshomim), 2)

getmat = function (idx1, idx2) {

  if ( length(meshomim[[idx1]]) > 0.5 & length(meshomim[[idx2]]) > 0.5 ) {

    mymat = matrix(0, nrow = length(meshomim[[idx1]]), ncol = length(meshomim[[idx2]]))

    for ( i in 1:length(meshomim[[idx1]]) ) {
      for ( j in 1:length(meshomim[[idx2]]) ) {
        mymat[i, j] = omim[ meshomim[[idx1]][i], meshomim[[idx2]][j] ]
      }
    }

  } else {

    mymat = matrix(0)

  }

  return(mymat)

}


diseasesimlist = vector('list', ncol(idx))

diseasesimlist <- foreach (i = 1:length(diseasesimlist), .errorhandling = 'pass') %dopar% {
  xxx <- getmat(idx[2, i], idx[1, i])
}

diseasesimAvgtotal = sapply(diseasesimlist, mean)
diseasesimmatAvgtotal = matrix(0, 746, 746)

for (i in 1:length(diseasesimlist)) {
  diseasesimmatAvgtotal[idx[2, i], idx[1, i]] = diseasesimAvgtotal[i]
}

diseasesimmatAvgtotal[upper.tri(diseasesimmatAvgtotal)] =
  t(diseasesimmatAvgtotal)[upper.tri(t(diseasesimmatAvgtotal))]
diag(diseasesimmatAvgtotal) = 1
diseasesimmatAvgtotal4digit = format(round(diseasesimmatAvgtotal, 4), nsmall = 4)
write.table(diseasesimmatAvgtotal4digit, 'diseasesimmat.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
