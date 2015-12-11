#!/usr/bin/Rscript

library(data.table)

tfdbFile <- 'data/tfdb/os.Rds'
hbGenesFile <- 'output/homeobox/hbGenes.Rds'
expressedGenesAllFile <- 'output/expressedGenes/expressedGenesAll.Rds'
vstFile <- 'output/DESeq2/vst.Rds'

# check files exist
lapply(list(tfdbFile, expressedGenesAllFile, vstFile, hbGenesFile), function(x)
  if(!file.exists(x)) {
    cat(x, " not found, exiting\n", file = stderr())
    quit(save = "no", status = 1)
  })

# vst geometric means
vst <- readRDS(vstFile)
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}
vstMeans.matrix <- sapply(levels(vst$stage), function(x)
  apply(GenomicRanges::assay(vst[, vst$stage == x]), 1, gm_mean))

# homeobox classes
hbGenes <- readRDS(hbGenesFile)

# tfdb
tfdb <- readRDS(tfdbFile)

# add tfdb Hbs missing from Jain et. al data
missingHbs <- tfdb[Family == "HB", Protein.ID[which(!(Protein.ID %in% hbGenes[, msuId]))]]
oryzr::LocToGeneName(missingHbs)
if (length(missingHbs) >0 ){
  warning(paste(
    "The following msuIds were not found in the Jain et. al data:\n",
    paste(missingHbs, collapse = ", "),
    "\nThey will be categorized as \"Not annotated\""))
  notAnn <- data.table(
    class = "Not annotated",
    msuId = missingHbs
  )
  notAnn[, symbol := oryzr::LocToGeneName(msuId)$symbols, by = msuId]
  hbGenes <- rbind(hbGenes, notAnn)
}

# expressed genes
expressedGenesAll <- readRDS(expressedGenesAllFile)

# expressed homeobox genes
expHb <- hbGenes[msuId %in% expressedGenesAll] 

# get vst read counts for expHb
hbVst <- vstMeans.matrix[expHb[, msuId],]

# scale the VST read counts
hbVstScaled <- t(scale(t(hbVst), center = TRUE))

# choose a cluster method
# dist.methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")
# hclust.methods <- c("ward.D", "ward.D2")
# for (dm in dist.methods){
#   for (hm in hclust.methods){
#     plot(hclust(dist(hbVstScaled, method = dm), method = hm), main = dm)
#   }
# }
dm <- "minkowski"
hm <- "ward.D2"

# cut the tree into 5 groups of expression
hc <- hclust(dist(hbVstScaled, method = dm), method = hm)
cuts <- cutree(hc, k = 5)

# get the order for the x-axis
geneOrder <- rownames(hbVstScaled)[hc$order]

# merge scaled vst to hb data
hbVstScaled.table <- data.table(hbVstScaled, keep.rownames = TRUE, key = 'rn')
setnames(hbVstScaled.table, 'rn', 'msuId')
setkey(expHb, 'msuId')
plotData <- expHb[hbVstScaled.table, ]

# set up plot
plotData[!is.na(symbol), symbol := paste(symbol, msuId, sep = "\n")]
plotData[is.na(symbol), symbol := msuId, by = msuId]
plotData.long <- reshape2::melt(plotData, id.vars = c('msuId', 'symbol', 'class'), variable.name = 'Stage',
                                value.name = 'Scaled reads')
plotData.long[, msuId := factor(msuId, levels = geneOrder)]
plotData.long[, symbol := factor(symbol, levels = plotData[geneOrder, symbol])]

# set up segment data
segData <- plotData[,.(class, msuId, symbol)]
idx <- 1:length(geneOrder)
names(idx) <- geneOrder
segData[, idx := idx[msuId]]
segData[, c("y", "yend") := .(
  idx - 0.5, idx + 0.5
)]
segData[, colour := plyr::mapvalues(
  class,
  from = unique(class),
  to = RColorBrewer::brewer.pal(length(unique(class)), "Set3")
  )]

# add cuts to plotdata and segdata
plotData.long[, cut := cuts[as.character(msuId)], by = msuId]
segData[, cut := cuts[as.character(msuId)], by = msuId]

# rescale y for segData relative to each facet
segData[, minY := min(y), by = cut]
segData[, yRel := y - minY + 0.5]

# MAKE OUTPUT FOLDER
outDir <- "output/homeobox"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(plotData.long, paste0(outDir, "/plotData.long.Rds"))
saveRDS(segData, paste0(outDir, "/segData.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
