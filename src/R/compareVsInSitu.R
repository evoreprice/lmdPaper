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

# check for deseq2 output
ddsFile <- 'output/DESeq2/ddsLrt.Rds'
if (!file.exists(ddsFile)) {
  cat("ddsFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
dds <- readRDS(ddsFile)

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

# make sure we have numerical, non-NA refs for each record
if (dim(expZhangGenes[!is.numeric(zhangRef) | is.na(zhangRef)])[1] > 0) {
  sink(file = stderr())
  cat("ERROR: non-numerical or NA reference detected. Check refs.\n")
  expZhangGenes[!is.numeric(zhangRef) | is.na(zhangRef)]
  sink()
  quit(save = "no", status = 1)
}

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
expZhangGenes[, id := factor(id, levels = rev(id[order(expZhangGenes[, plotLabel])]))]

# melt
expZhangGenes.melt <- reshape2::melt(
  expZhangGenes, id.vars = c("msuId", "zhangRef", "plotLabel", "id"),
  measure.vars = c("RM", "PBM", "SBM", "SM", "FM"), variable.name = "stage",
  value.name = "call.z")

# compare this to the expressed genes truth table
expGenTT.wide <- data.table(readRDS("output/expressedGenes/expGenTT.Rds"))
expGenTT <- reshape2::melt(expGenTT.wide, id.vars = "id", variable.name = "lib", value.name = "expressed")

# add stage names from deseq2, but map "ePBM/SBM" to "SBM"
colData <- data.table(as.data.frame(GenomicRanges::colData(
  readRDS('output/DESeq2/ddsLrt.Rds'))), keep.rownames = TRUE)
expGenTT[, stage := colData[rn == as.character(lib), stage]]
expGenTT[, stage := plyr::mapvalues(stage, "ePBM/SBM", "SBM")]

# call a gene expressed if it's detected in > 1 library per stage
expGenTT[, call.rs := sum(expressed) > 1, by = c("id", "stage")]

# collapse by stage + id
setkey(expGenTT, "id", "stage")
expGenCalls <- unique(expGenTT)
expGenCalls[, c("expressed", "lib") := NULL ]

# add a fake FM column for a quick comparo
expGenCalls <-  rbind(expGenCalls, expGenCalls[stage == "SM", .(id, stage = "FM", call.rs)])

# rename and sort for easy join
expGenCalls[,msuId := id][, id := NULL]
setkey(expGenCalls, "msuId", "stage")
setkey(expZhangGenes.melt, "msuId", "stage")

# join
plotData.long <- expGenCalls[expZhangGenes.melt]

# set levels of stage
plotData.long[, stage := factor(stage, levels = c("RM", "PBM", "SBM", "SM", "FM"))]

# add column to compare in situ vs. rnaseq. × = both, + = in situ only, · =
# rnaseq only
plotData.long[call.rs & call.z, compare := "×"]
plotData.long[!call.rs & call.z, compare := "+"]
plotData.long[call.rs & !call.z, compare := "·"]

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