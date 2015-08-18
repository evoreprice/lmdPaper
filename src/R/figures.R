#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
extrafont::loadfonts()

#####################
### SETUP FOR RMD ###
#####################

# https://gist.githubusercontent.com/rmflight/3858973/raw/591ee603ce0996ac407902f7dd59520b35b3702a/tableFigureFuncs.r
pasteLabel <- function(preText, inObj, objName, insLink=TRUE){
  objNum <- inObj[objName]
  
  useText <- paste(preText, objNum, sep=" ")
  if (insLink){
    useText <- paste("[", useText, "](#", objName, ")", sep="")
  }
  useText
}

# this increments the counter, and gives a name to the number so we can reference it later
incCount <- function(inObj, useName){
  nObj <- length(inObj)
  useNum <- max(inObj) + 1
  inObj <- c(inObj, useNum)
  names(inObj)[nObj+1] <- useName
  inObj
}
figCount <- c("_"=0)
tableCount <- c("_"=0)

tableCat <- function(inFrame){
  outText <- paste(names(inFrame), collapse=" | ")
  outText <- c(outText, paste(rep("---", ncol(inFrame)), collapse=" | "))
  invisible(apply(inFrame, 1, function(inRow){
    outText <<- c(outText, paste(inRow, collapse=" | "))
  }))
  return(outText)
}

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

###################
### STATS TABLE ###
###################

t_libStats <- readRDS('output/quantStats/libStats.Rds')

#############
### Mfuzz ###
#############

expressionMatrix <- readRDS('output/mfuzz/expressionMatrix.Rds')
c1 <- readRDS('output/mfuzz/c1.Rds')
memCutoff <- 0.5

# get clusters and membership
cluster <- data.table(id = names(c1$cluster), Cluster = c1$cluster, Membership = apply(c1$membership, 1, max), key = "id")

# get standardised VST counts
exprs <- data.table(Biobase::exprs(expressionMatrix), keep.rownames = TRUE, key = "rn")

# set up cluster labels
numGenes <- cluster[Membership > memCutoff, .(number = length(id)), by = Cluster]
setkey(numGenes, 'Cluster')
clustLabels <- numGenes[, paste0("Cluster ", Cluster, "\n(", number, " genes)")]

# merge clusters with expression values
plotData.wide <- cluster[exprs]

# make long
plotData <- reshape2::melt(plotData.wide[ Membership > memCutoff ], id.vars = c('id', 'Cluster', 'Membership'), variable.name = 'Stage', value.name = "Normalised, transformed read counts")

# add cluster labels
plotData[, Cluster := plyr::mapvalues(as.factor(Cluster), seq(1, length(unique(Cluster))), clustLabels)]

# fix stage label
plotData[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# set up heatscale
heatscale <- RColorBrewer::brewer.pal(n = 6, name = "YlOrRd")

# add centres
centres.wide <- data.table(c1$centers, Cluster = clustLabels)
centres <- reshape2::melt(centres.wide, id.vars = 'Cluster', variable.name = 'Stage', value.name = "Normalised, transformed read counts")
centres[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# number of clusters (for text)
c <- length(unique(c1$cluster))
maxClust <- readRDS('output/mfuzz/maxClust.Rds')

# main cluster figure
f_mfuzzClusters <- ggplot(plotData, aes(x = Stage, y = `Normalised, transformed read counts`,
                                        colour = Membership, group = id)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = grid::unit(1, "lines")) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1), colour = "black", alpha = 0.5) +
  facet_wrap(~ Cluster, nrow = 3)

# centroid dis vs. c (for SI)
f_mfuzzCentroids <- readRDS('output/mfuzz/centPlot.Rds')
# 
# # PCA (for SI)
# clustered <- unique(names(c1$cluster))
# 
# vg.d <- dist(Biobase::exprs(expressionMatrix[clustered,]))
# vg.mds <- cmdscale(vg.d, 2)
# vg.mds <- data.frame(MDS1 = vg.mds[,1], MDS2 = vg.mds[,2], cluster = c1$cluster,
#                      max.membership = apply(c1$membership,1,max))
# 
# f_mfuzzPca <- ggplot(vg.mds, aes(x = MDS1, y=MDS2, colour = factor(cluster))) +
#   theme_minimal(base_size = 8, base_family = "Helvetica") +
#   coord_fixed(ratio = 1) +
#   geom_point(aes(size=max.membership),alpha=0.5, shape = 16) +
#   scale_colour_brewer(palette = "Set1", name = "Cluster") +
#   scale_size_area(guide = FALSE)

