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
setkey(zhangGenes, 'msuId')

# load expressed genes
exprGenFile <- "output/expressedGenes/expressedGenesByLibrary.Rds"
if (!file.exists(exprGenFile)) {
  cat("exprGenFile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
expressedGenesByLibrary <- readRDS(exprGenFile)

# split the zhang genes expression field.
trim <- function(x) {gsub("^\\s+|\\s+$", "", x)}
zhangGenes[, geneExp := toupper(trim(unlist(strsplit(geneExp, ",", fixed = TRUE)))),
           by = msuId]

# remove some columns that we don't need for this plot
zhangGenes[,c("geneName", "geneFunction", "geneEvidence"):= NULL]

# get the LOCs that have expression we're interested in
RM <- c("IM", "RM")
PBM <- c("BM", "PBM")
SBM <- c("BM", "SBM")
SM <- c("SM")
FM <- c("FM")
tissues <- unique(c(RM, PBM, SBM, SM))
expZhangGenes <- zhangGenes[geneExp %in% tissues]

# make a truth table for the zhang genes
expZhangGenes[,RMz := geneExp %in% RM]
expZhangGenes[,PBMz := geneExp %in% PBM]
expZhangGenes[,SBMz := geneExp %in% SBM]
expZhangGenes[,SMz := geneExp %in% SM]
expZhangGenes[,FMz := geneExp %in% FM]

# find out if we called these genes expressed (2 out of 3 libraries per stage)
egbl <- data.table(msuId = unique(unlist(expressedGenesByLibrary)), key = "msuId")
egbl <- cbind(egbl, egbl[, lapply(expressedGenesByLibrary, function(x) msuId %in% x)])
egbl <- egbl[, .(
  RMr = sum(n1r1,n1r3, n1r4) > 1,
  PBMr = sum(n2r1,n2r3, n2r4) > 1,
  SBMr = sum(n3r1,n3r2, n3r3) > 1,
  SMr = sum(n4r1,n4r2, n1r3) > 1
  ), by = msuId]

##################
### UP TO HERE ###
##################

# todo: fix multiple lines for zhangGene truth table (just use "OR")

# join 
isVsRs <- egbl[expZhangGenes]
isVsRs

# how is this going to be plotted?
randLab <- function(i) {paste0(sample(letters, 3, replace = FALSE), collapse = "")}
pd <- data.frame(id = 1:10,
           value1 = sample(c(TRUE, FALSE), size = 10, replace = TRUE),
           value2 = sample(c(TRUE, FALSE), size = 10, replace = TRUE),
           label = sapply(1:10, randLab))
pd.melt <- reshape2::melt(pd, id.vars = c('id', 'label'))

# hack ggplot

# mappings as normal
ggplot(pd.melt, aes(x = variable, y = id, fill = value)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  xlab(NULL) + ylab(NULL) +
  # make a custom scale specifying the labels and the breaks (breaks == columns
  # in the df). Negative expand to shrink label columns
  scale_x_discrete(labels = c("", "T1", "T2", ""),
                   limits = c("id", "value1", "value2", "label"),
                   expand = c(0,-0.4)) + 
  geom_raster() +
  # map one label to max + 0.5 and hjust to the right
  geom_text(mapping = aes(x = 3.5, y = id, label = label), hjust = -0.1) +
  # map the other label to 1.5 and hjust to the left
  geom_text(mapping = aes(x = 1.5, y = id, label = id), hjust = 1.1)


