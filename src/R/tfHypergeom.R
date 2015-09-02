#!/usr/bin/Rscript

library(data.table)

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

# what about ALOGs (not TFs)
ALOG <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
          'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
          'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500', 'LOC_Os05g28040')


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
resultsLong$padj <- p.adjust(resultsLong$`p-value`, method = "BH")

# go wide
padj <- reshape2::dcast(resultsLong, Family ~ Cluster, value.var = 'padj')
rownames(padj) <- padj$Family
padj$Family <- NULL
colnames(padj) <- paste0("C~", sub(".*(\\d+).*", "\\1", colnames(padj)),"~ *p*~adj~")

Number <- reshape2::dcast(resultsLong, Family ~ Cluster, value.var = 'Number')
rownames(Number) <- Number$Family
Number$Family <- NULL
colnames(Number) <- paste0("C~", sub(".*(\\d+).*", "\\1", colnames(Number)),"~ *n*")

# remove columns that don't have p-values
keep <- apply(padj, 1, function(x) !all(is.na(x)))

sigFamilies <- cbind(Number[keep, ], padj[keep,])

# add the number of genes per family
genesPerFamily <- sapply(tfdbByFamily, function(x) dim(x)[1])
sigFamilies$`*n*` <- genesPerFamily[rownames(sigFamilies)]

# add the number of expresed genes per family
sigFamilies$`*n*~expr~` <- exprPerFamily[rownames(sigFamilies)]

# set the column order programatically
# these two are always first
generalCn <- c("*n*", "*n*~expr~")
# use data.table sorting to order...
cn <- colnames(sigFamilies)
cnDt <- data.table(cn = cn[!cn %in% generalCn])
# ... first by cluster number ...
cnDt[, cnum := as.numeric(gsub("^.*(\\d).*$", "\\1", cn))]
# ... then by n in cluster / p-val
cnDt[, pn := gsub(".*([np]).*", "\\1", cn)]
# do the sort
setkey(cnDt, "cnum", "pn")
# re-order sigFamilies
sigFamilies <- sigFamilies[, c(generalCn, cnDt[,cn])]

# output
saveRDS(sigFamilies, paste0(mfuzzDir, "/sigFamilies.Rds"))
