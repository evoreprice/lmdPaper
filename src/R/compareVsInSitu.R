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

# split the zhang genes expression field
trim <- function(x) {gsub("^\\s+|\\s+$", "", x)}
idVsExp <- zhangGenes[, .(
  meristem = toupper(trim(unlist(strsplit(geneExp, "[,]"))))
  ), by = msuId]
setkey(idVsExp, 'msuId')

# get the LOCs that have expression we're interested in
RM <- c("IM", "RM")
PBM <- c("BM", "PBM")
SBM <- c("BM", "SBM")
SM <- c("SM", "FM")

RMt <- idVsExp[meristem %in% RM, .(msuId = unique(msuId), RM = TRUE)]
PBMt <- idVsExp[meristem %in% PBM, .(msuId = unique(msuId), PBM = TRUE)]
SBMt <- idVsExp[meristem %in% SBM, .(msuId = unique(msuId), SBM = TRUE)]
SMt <- idVsExp[meristem %in% SM, .(msuId = unique(msuId), SM = TRUE)]

# this is nasty
merge(RMt, merge(PBMt, merge(SBMt, SMt, all = TRUE), all = TRUE), all = TRUE)



inSituResults <- idVsExp[meristem %in% tissues]

newtab <- inSituResults[, .(msuId = unique(msuId))]


inSituResults[,expressed := TRUE]



# cast to wide
reshape2::dcast(inSituResults, msuId ~ meristem, value.var = "expressed")


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


