# Multiple Evidence Fusion
#
# Compute Pathway-based Drug-Drug Similarity
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Sep 10, 2013

dbid = as.character(read.table('drugkegg.txt', sep = '\t')[, 2])
chem = read.csv('CTD_chemicals.csv', comment.char = "#")
path = read.csv('CTD_chem_pathways_enriched.csv', comment.char = "#")

chem = chem[, c(2, 9)]
path = path[, c(2, 5)]

chem = chem[which(chem[, 2] != ''), ]
chem[, 1] = gsub('MESH:', '', chem[, 1])

# expand duplicated DrugBank IDs
chemtmp = vector('list', nrow(chem))
for (i in 1:nrow(chem)) chemtmp[[i]] = as.character(chem[i, 2])
names(chemtmp) = as.character(chem[, 1])
chemtmp = sapply(chemtmp, strsplit, split = '\\|')

tmp = unlist(chemtmp)
chem2 = data.frame(DrugBankIDs = tmp, ChemicalID = names(tmp))
chem2[, 2] = substring(chem2[, 2], 1, 7)

# map from DrugBank id to mesh id
dbmesh = rep('', 746)
for (i in 1:746) try(dbmesh[i] <- chem2[which(dbid[i] == chem2[, 1]), 2], silent = TRUE)

# map from mesh id to pathway ids
meshpath = vector('list', 746)
names(meshpath) = dbid

for (i in 1:746) {
  if ( dbmesh[i] != '' ) {
    try(meshpath[[i]] <- as.character(path[which(dbmesh[i] == path[, 1]), 2]),
        silent = TRUE)
  }
}

# calculate similarity matrix
mat = matrix(0.0, ncol = 746, nrow = 746)
idx = combn(1:length(meshpath), 2)

for (i in 1:ncol(idx)) {
  if ( length(meshpath[[ idx[2, i] ]]) > 0.5 & length(meshpath[[ idx[1, i] ]]) > 0.5 ) {
    C = length(intersect(meshpath[[ idx[2, i] ]], meshpath[[ idx[1, i] ]]))
    A = length(meshpath[[ idx[2, i] ]])
    B = length(meshpath[[ idx[1, i] ]])
    mat[idx[2, i], idx[1, i]] = C / (A + B - C)
  }
}

mat[upper.tri(mat)] = t(mat)[upper.tri(t(mat))]
diag(mat) = 1

mat4digit = format(round(mat, 4), nsmall = 4)

write.table(mat4digit, 'pathwaysimmat.txt', sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
