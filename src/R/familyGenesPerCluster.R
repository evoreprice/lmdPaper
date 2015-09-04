#!/usr/bin/Rscript

library(data.table)

# find tfdb
tfdbFile <- 'data/tfdb/os.Rds'
if (!file.exists(tfdbFile)) {
  cat("tfdb file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
tfdb <- readRDS(tfdbFile)

# check for mfuzz output
clusterFile <- "output/mfuzz/clusters.Rds"
if (!file.exists(clusterFile)) {
  cat("clusterFile not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
clusters <- readRDS(clusterFile)

# find expressed genes
expGenFile <- "output/expressedGenes/expressedGenesAll.Rds"
if (!file.exists(expGenFile)) {
  cat("tfdb file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
expressedGenesAll <- readRDS(expGenFile)

# think we need a long table of cluster vs gene vs family
clusters.wide <- data.table(msuId = unique(unlist(sapply(clusters, rownames))))
n <- paste("Cluster", 1:length(clusters))
clusters.wide[, (n) := lapply(clusters, function(x) msuId %in% rownames(x)), by = msuId]
clusters.wide[msuId %in% tfdb$Protein.ID, Family := tfdb[Protein.ID == msuId, Family], by = msuId]
clusters.long <- reshape2::melt(clusters.wide, id.vars = c("msuId", "Family"), variable.name = "Cluster")


reshape2::dcast(clusters.long[!is.na(Family) & value], msuId + Family ~ Cluster, fun.aggregate = sum, value.var = "value")


# cluster table
genesPerCluster <- tfdb[,.("*n*" = length(Protein.ID),
        "*n*~exp~" = sum(Protein.ID %in% expressedGenesAll)), by = Family]  

paste("Cluster", 1:length(clusters))
genesPerCluster


x <- 'ABI3VP1'
tfdb[Family == x, sum(Protein.ID %in% expressedGenesAll), by = Family]


