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
  geneExp = paste(geneExp, collapse = " | "),
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
                                        fun.aggregate = function(x) sum(x) > 2), key = "msuId")
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

# test plot
colours <- RColorBrewer::brewer.pal(3, "Set1")[c(2,1,3)]
labels <- c("1" = "Detected\nin both", "2" = "Only detected\nby in situ", "3" = "Only detected by\nRNA-sequencing")
ggplot(plotData.long, aes(x = stage, y = id, colour = as.factor(compare))) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  xlab(NULL) + ylab(NULL) +
  scale_colour_manual(values = colours, na.value = NA,
                      labels = labels) +
  geom_point(size = 1) +
  scale_x_discrete(labels = c("", "RM", "PBM", "SBM", "SM", "FM", ""),
                   limits = c("id", "RM", "PBM", "SBM", "SM", "FM", "zhangRef"),
                   expand = c(0,-0.5)) +
  # map one label to 1.5 and hjust to the left
  geom_text(mapping = aes(x = 1.9, y = id, label = plotLabel),
            hjust = 1, colour = "black", size = 2, fontface = "italic") +
  # map the other to max + 0.5 and hjust to the right
  geom_text(mapping = aes(x = 6.1, y = id, label = zhangRef),
            hjust = 0, colour = "black", size = 2, fontface = "plain") 


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


