#!/usr/bin/Rscript

# find tfdb
tfdbFile <- 'data/tfdb/os.Rds'
if (!file.exists(tfdbFile)) {
  cat("tfdb file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for mfuzz output
mfuzzDir <- "output/mfuzz"
if (!dir.exists(mfuzzDir)) {
  cat("mfuzzDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# load up the data
tfdb <- readRDS(tfdbFile)
vstFiltered <- readRDS(paste0(mfuzzDir, "/vstFiltered.Rds"))
clusters <- readRDS(paste0(mfuzzDir, "/clusters.Rds"))

# split tfdb by family
tfdbByFamily <- split(tfdb, tfdb$Family)

# total expressed genes
numExpr <- dim(vstFiltered)[1]

# number of genes from each family that is expressed (i.e. background for hypergeom)
exprPerFamily <- sapply(tfdbByFamily, function(x)
  sum(x$Protein.ID %in% rownames(vstFiltered)))

# number of genes per family in each cluster
inClusterPerFamily <- sapply(tfdbByFamily, function(family)
  sapply(clusters, function(x)
    sum(family$Protein.ID %in% rownames(x))))
rownames(inClusterPerFamily) <- paste0("C", 1:dim(inClusterPerFamily)[1])

# number of genes in each cluster
clusterSize <- sapply(clusters, function(x) dim(x)[1])
names(clusterSize) <- paste0("C", 1:dim(inClusterPerFamily)[1])

# run hypergeometric test
hyperGeomResults <- sapply(1:length(tfdbByFamily), function(i)
  phyper(inClusterPerFamily[,i] - 1, exprPerFamily[i], numExpr - exprPerFamily[i],
         clusterSize, lower.tail = FALSE))
colnames(hyperGeomResults) <- names(tfdbByFamily)

# combine the results for adjusting p-values
resultsLong <- merge(reshape2::melt(inClusterPerFamily,
                                    varnames = c("Cluster", "Family"),
                                    value.name = "Number"),
                     reshape2::melt(hyperGeomResults,
                                    varnames = c("Cluster", "Family"),
                                    value.name = "p-value"))

# don't do z-tests if number of genes < 2
resultsLong[resultsLong$Number < 2, ]$`p-value` <- NA

# calculate adjusted p-values
resultsLong$padj <- p.adjust(resultsLong$`p-value`, method = "BH",
                             n = sum(!is.na(resultsLong$`p-value`)))

resultsLong$genesInFamily <- genesPerFamily[as.character(resultsLong$Family)]

# add the number of genes expressed (think I should skip)
resultsLong$exprGenesInFamily <- exprPerFamily[as.character(resultsLong$Family)]

# go wide
padj <- reshape2::dcast(resultsLong, Family ~ Cluster, value.var = 'padj')
rownames(padj) <- padj$Family
padj$Family <- NULL
colnames(padj) <- paste0("padj(", colnames(padj), ")")

Number <- reshape2::dcast(resultsLong, Family ~ Cluster, value.var = 'Number')
rownames(Number) <- Number$Family
Number$Family <- NULL
colnames(Number) <- paste0("n(", colnames(Number), ")")

# remove columns that don't have p-values
keep <- apply(padj, 1, function(x) !all(is.na(x)))

sigFamilies <- cbind(Number[keep, ], padj[keep,])

# add the number of genes per family
genesPerFamily <- sapply(tfdbByFamily, function(x) dim(x)[1])
sigFamilies$`Number of genes` <- genesPerFamily[rownames(sigFamilies)]

# add the number of expresed genes per family
sigFamilies$`Number of expressed genes` <- exprPerFamily[rownames(sigFamilies)]

# set the column order
sigFamilies <- sigFamilies[,c('Number of genes', 'Number of expressed genes',
                              'n(C1)', 'padj(C1)', 'n(C2)', 'padj(C2)', 'n(C3)',
                              'padj(C3)', 'n(C4)', 'padj(C4)', 'n(C5)',
                              'padj(C5)', 'n(C6)', 'padj(C6)')]

# output
saveRDS(sigFamilies, paste0(mfuzzDir, "/sigFamilies.Rds"))
