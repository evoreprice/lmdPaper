#!/usr/bin/Rscript

library(data.table)

# load zhang file
zhangFile <- "data/zhangGenes.tab"
if (!file.exists(zhangFile)) {
  cat("zhangFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
zhangGenes <- data.table(read.table(zhangFile, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE))

# load expressed genes
exprGenFile <- "output/expressedGenes/expressedGenesByLibrary.Rds"
if (!file.exists(exprGenFile)) {
  cat("exprGenFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
expressedGenesByLibrary <- readRDS(exprGenFile)

# split the zhang genes expression field.
trim <- function(x) {gsub("^\\s+|\\s+$", "", x)}
zhangGenes <- zhangGenes[,.(
  geneExp = trim(unlist(strsplit(geneExp, "[,]")))
  ), by = c('msuId', 'zhangRef', 'plotLabel')]

# get the LOCs that have expression we're interested in
RM <- c("IM", "RM")
PBM <- c("BM", "PBM")
SBM <- c("BM", "SBM")
SM <- c("SM")
FM <- c("FM")
tissues <- unique(c(RM, PBM, SBM, SM, FM))
expZhangGenes <- zhangGenes[geneExp %in% tissues]

# make a truth table for the zhang genes
expZhangGenes[,RM := geneExp %in% RM]
expZhangGenes[,PBM := geneExp %in% PBM]
expZhangGenes[,SBM := geneExp %in% SBM]
expZhangGenes[,SM := geneExp %in% SM]
expZhangGenes[,FM := geneExp %in% FM]

# collapse multiple lines by msuId and zhangRef
setkey(expZhangGenes, "msuId", "zhangRef")
expZhangGenes <- expZhangGenes[, .(
  plotLabel = unique(plotLabel),
  geneExp = paste(unique(geneExp), collapse = " | "),
  RM = any(RM),
  PBM = any(PBM),
  SBM = any(SBM),
  SM = any(SM),
  FM = any(FM)
), by = c("msuId", "zhangRef")]

# add another ID column for genes with more than one reference and order it on plot labels
setkey(expZhangGenes, "msuId", "zhangRef")
expZhangGenes[, id := paste0("id", 1:dim(unique(expZhangGenes))[1])]
expZhangGenes[,id := factor(id, levels = rev(id[order(expZhangGenes[, plotLabel])]))]

# melt
expZhangGenes.melt <- reshape2::melt(
  expZhangGenes, id.vars = c("msuId", "zhangRef", "plotLabel", "id"),
  measure.vars = c("RM", "PBM", "SBM", "SM", "FM"), variable.name = "stage",
  value.name = "call.z")

# find out if we called these genes expressed (2 out of 3 libraries per stage)
egbl <- data.table(msuId = unique(unlist(expressedGenesByLibrary)), key = "msuId")
egbl <- cbind(egbl, egbl[, lapply(expressedGenesByLibrary, function(x) msuId %in% x)])

# melt and dcast to sum logicals
egbl.melt <- reshape2::melt(egbl, id.vars = "msuId", variable.name = "sample",
                            value.name = "call")
egbl.melt[, stage := sub("r\\d+", "", sample)]
egbl.melt[, stage := plyr::mapvalues(stage, from = c("n1", "n2", "n3", "n4"),
                                     to = c("RM", "PBM", "SBM", "SM"))]
geneCalls <- data.table(reshape2::dcast(egbl.melt, msuId ~ stage, value.var = "call",
                                        fun.aggregate = function(x) sum(x) > 1), key = "msuId")
# add a fake "FM" column to check if we are getting a lot of genes in the SM sample
geneCalls[, FM := SM]

# melt geneCalls and join
geneCalls.melt <- reshape2::melt(geneCalls, id.vars = "msuId",
                                 variable.name = "stage", value.name = "call.rs")
setkey(geneCalls.melt, "msuId", "stage")
setkey(expZhangGenes.melt, "msuId", "stage")

plotData.long <- geneCalls.melt[expZhangGenes.melt]

# fix NAs
plotData.long[is.na(call.rs), call.rs := FALSE]
# set levels of stage
plotData.long[, stage := factor(stage, levels = c("RM", "PBM", "SBM", "SM", "FM"))]

# add column to compare in situ vs. rnaseq. 1 = both, 2 = in situ only, 3 =
# rnaseq only
plotData.long[call.rs & call.z, compare := 1]
plotData.long[!call.rs & call.z, compare := 2]
plotData.long[call.rs & !call.z, compare := 3]


# MAKE OUTPUT FOLDER
outDir <- "output/compare"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(plotData.long, paste0(outDir, "/compare.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)