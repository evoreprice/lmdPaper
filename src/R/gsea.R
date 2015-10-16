#!/usr/bin/Rscript

library(gage)
library(DESeq2)
library(data.table)

# check for DESeq2 output
deseqDir <- "output/DESeq2"
if (!dir.exists(deseqDir)) {
  cat("deseqDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
ddsFile <- paste0(deseqDir, "/ddsWald.Rds")
if (!file.exists(ddsFile)) {
  cat("ddsFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
dds <- BiocGenerics::updateObject(readRDS(ddsFile))

# find tfdb
tfdbFile <- 'data/tfdb/os.Rds'
if (!file.exists(tfdbFile)) {
  cat("tfdb file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
tfdb <- readRDS(tfdbFile)
famCatFile <- "data/tfdb/famCat.Rds"
if (!file.exists(famCatFile)) {
  cat("famCat file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
famCat <- readRDS(famCatFile)

# load expressed genes
exprGenFile <- "output/expressedGenes/expressedGenesAll.Rds"
if (!file.exists(exprGenFile)) {
  cat("exprGenFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
expressedGenes <- readRDS(exprGenFile)

# one-by-one comparisons of each stage to all other stages using contrast
# vectors like this: cv <- c(0,1,1,-1,-1,0,0). Repeat for each level of stage
stage <- levels(dds$stage)

getRes <- function(x) {
  x <- make.names(x)
  # contrast vector based on results names
  cv <- resultsNames(dds)
  # resultsNames that don't start with stage aren't in the comparison
  cv[!grepl("^stage", cv)] <- 0
  # resultsNames for our stage of interest are +1.
  cv[grepl(x, cv)] <- 1
  # all other stages are for comparison. Set them to -1.
  cv[!(cv == 1 | cv == 0)] <- -1
  # make numeric for results call
  cv <- as.numeric(cv)
  return(results(dds, contrast = cv)$log2FoldChange)
}

lfcFrame <- sapply(stage, getRes)
rownames(lfcFrame) <- rownames(counts(dds))

# run gage on the expressed TFDB genes and extract test statistics
tfdb.expressed <- tfdb[Protein.ID %in% expressedGenes]

res <- gage(lfcFrame, gsets = split(tfdb.expressed$Protein.ID, tfdb.expressed$Family),
            samp = NULL, ref = NULL)
gageFrame <- res$stats
gageFrame <- gageFrame[, stage]
gageFrame <- gageFrame[complete.cases(gageFrame), ]

# get the family order for the plot. Cluster on TFs and others separately since
# they will be facetted in the plot.
tfFrame <- gageFrame[(rownames(gageFrame) %in% famCat[Category == 'TF', Family]),]
tfOrder <- rownames(tfFrame)[hclust(dist(tfFrame))$order]
otherFrame <- gageFrame[(rownames(gageFrame) %in% famCat[Category == 'Other', Family]),]
otherOrder <- rownames(otherFrame)[hclust(dist(otherFrame))$order]

famOrder <- c(tfOrder, otherOrder)

# make p-value data.frames
greaterP <- res$greater
greaterP <- greaterP[, stage]
lessP <- res$less
lessP <- lessP[,stage]

# data.table magic to bind it all together
gageTable <- data.table(gageFrame, keep.rownames = TRUE, key = 'rn')
gpTable <- data.table(greaterP, keep.rownames = TRUE, key = 'rn')
lpTable <- data.table(lessP, keep.rownames = TRUE, key = 'rn')

# melt then merge
gageTable.melt <- reshape2::melt(gageTable, id.vars = 'rn', measure.vars = stage,
               variable.name = "Stage", value.name = "Test statistic")
gpTable.melt <- reshape2::melt(gpTable, id.vars = 'rn', measure.vars = stage,
                                 variable.name = "Stage", value.name = "pg")
lpTable.melt <- reshape2::melt(lpTable, id.vars = 'rn', measure.vars = stage,
                                 variable.name = "Stage", value.name = "pl")
setkey(gageTable.melt, "rn", "Stage")
setkey(gpTable.melt, "rn", "Stage")
setkey(lpTable.melt, "rn", "Stage")

plotData <- lpTable.melt[gpTable.melt][gageTable.melt]

# pick which p-value to use
plotData[`Test statistic` < 0, pval := pl]
plotData[`Test statistic` > 0, pval := pg]

# adjust p-values, need un-data-table syntax to avoid segfault. yuck.
plotData[!is.na(pval), padj := p.adjust(pval, "BH")]

# factor for lazy ggplotting
plotData[, showPval := padj < 0.1]

# make the pvals pretty
plotData[, padj := format.pval(c(0.12, padj), digits = 2, eps = 0.001, na.form = '')[-1]]

# order the y-axis
plotData[, rn := factor(rn, levels = famOrder)]

# MAKE OUTPUT FOLDER
outDir <- "output/gsea"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(plotData, paste0(outDir, "/gseaTable.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
