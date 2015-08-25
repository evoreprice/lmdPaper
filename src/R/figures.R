#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
extrafont::loadfonts()

#####################
### SETUP FOR RMD ###
#####################

# https://gist.githubusercontent.com/rmflight/3858973/raw/591ee603ce0996ac407902f7dd59520b35b3702a/tableFigureFuncs.r
pasteLabel <- function(preText, inObj, objName, insLink=FALSE){
  objNum <- inObj[objName]
  
  useText <- paste0(preText, objNum)
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
s_figCount <- c("_"=0)
s_tableCount <- c("_"=0)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

###################
### STATS TABLE ###
###################

st_libStats <- readRDS('output/quantStats/libStats.Rds')
s_tableCount <- incCount(s_tableCount, "st_libStats")

##################
### LMD FIGURE ###
##################

figCount <- incCount(figCount, "f_lmdFigure")

################
### PCA PLOT ###
################

expressedGenes <- readRDS('output/expressedGenes/expressedGenesAll.Rds')
vst <- readRDS('output/DESeq2/vstAll.Rds')

exprVst <- GenomicRanges::assay(vst)[expressedGenes,]
exprVst <- GenomicRanges::assay(vst)
pca <- prcomp(t(exprVst))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

pcaPlotData <- data.frame(
  label = toupper(rownames(pca$x)),
  PCA1 = pca$x[,1],
  PCA2 = pca$x[,2],
  Stage = GenomicRanges::colData(vst)$stage
)

# row index to help with overplotted labels
pcaPlotData$rIdx <- 1:dim(pcaPlotData)[1]
opIdx <- pcaPlotData[c('n2r3', 'n4r1'),'rIdx']

sf_pca <- ggplot(data = pcaPlotData,
                aes(x = PCA1, y = PCA2, colour = Stage, label = label)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  scale_colour_brewer(palette = "Set1", name = 'Stage') +
  guides(colour=guide_legend(title=NULL)) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  xlim(-50, 77) +
  geom_point(shape = 16, size = 5, alpha = 0.7) +
  geom_text(data = pcaPlotData[-opIdx,],
            hjust = 1.2, show.legend = FALSE) +
  geom_text(data = pcaPlotData[opIdx,],
            hjust = -0.2, show.legend = FALSE)

s_figCount <- incCount(s_figCount, "sf_pca")


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
  facet_wrap(~ Cluster, ncol = 2)

figCount <- incCount(figCount, "f_mfuzzClusters")

# centroid dis vs. c (for SI)
sf_mfuzzCentroids <- readRDS('output/mfuzz/centPlot.Rds') +
  theme_minimal(base_size = 8, base_family = "Helvetica")

s_figCount <- incCount(s_figCount, "sf_mfuzzCentroids")

# MDS for SI
vg.mds <- readRDS('output/mfuzz/vg.mds.Rds')
sf_mfuzzPca <- ggplot(vg.mds, aes(x = MDS1, y=MDS2, colour = factor(cluster))) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  coord_fixed(ratio = 1) +
  geom_point(aes(size=max.membership),alpha=0.5, shape = 16) +
  scale_colour_brewer(palette = "Set1", name = "Cluster") +
  scale_size_area(guide = FALSE)

s_figCount <- incCount(s_figCount, "sf_mfuzzPca")

#################
### HYPERGEOM ###
#################

t_hypergeom <- readRDS('output/mfuzz/sigFamilies.Rds')

# fix p-values. 
pFix <- function(x) {
  format.pval(c(0.122, x), digits = 2, eps = 0.001, na.form = "")[-1]
}
# just use a for loop, should use *apply() but too complicated
for (x in colnames(t_hypergeom)) {
  if (grepl("^\\*p", x)) {
    t_hypergeom[,x] <- pFix(t_hypergeom[,x])
  }
}

# remove NAs
t_hypergeom[is.na(t_hypergeom)] <- ""

tableCount <- incCount(tableCount, "t_hypergeom")

############
### GSEA ###
############

gsea <- readRDS('output/gsea/gseaTable.Rds')

# format some labels
setnames(gsea, old = "Test statistic", new = "Test\nstatistic")
gsea[, Stage := plyr::mapvalues(Stage, "ePBM/SBM", "ePBM/\nSBM")]

heatscale <- rev(RColorBrewer::brewer.pal(6, "RdBu"))

f_gsea <- ggplot(gsea, aes(x = Stage, y = rn, label = padj, fill = `Test\nstatistic`)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  xlab(NULL) + ylab(NULL) +
  theme(legend.key.size = grid::unit(1, "lines"),
        axis.text.x = element_text(vjust = 0.5)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(colours = heatscale) +
  geom_raster() +
  geom_text(data = gsea[showPval == TRUE], size = 2)

figCount <- incCount(figCount, "f_gsea")

#######################
### COMPARE IN SITU ###
#######################

compare <- readRDS('output/compare/compare.Rds')

# test plot
colours <- RColorBrewer::brewer.pal(3, "Set1")[c(2,1,3)]
labels <- c("1" = "Detected\nin both", "2" = "Only detected\nby in situ", "3" = "Only detected by\nRNA-sequencing")
f_reviewInSitu <- ggplot(compare, aes(x = stage, y = id, colour = as.factor(compare))) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank()) +
  xlab(NULL) + ylab(NULL) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(values = colours, na.value = NA,
                      labels = labels, name = NULL) +
  geom_point(size = 1) +
  scale_x_discrete(labels = c("", "RM", "PBM", "SBM", "SM", "FM", ""),
                   limits = c("id", "RM", "PBM", "SBM", "SM", "FM", "zhangRef"),
                   expand = c(2,0)) +
  # map one label to 1.5 and hjust to the left
  geom_text(mapping = aes(x = 1.5, y = id, label = plotLabel),
            hjust = 1, colour = "black", size = 2, fontface = "italic",
            family = "Helvetica") +
  # map the other to max + 0.5 and hjust to the right
  geom_text(mapping = aes(x = 6.5, y = id, label = zhangRef),
            hjust = 0, colour = "black", size = 2, fontface = "plain",
            family = "Helvetica") 
f_reviewInSitu
figCount <- incCount(figCount, "f_reviewInSitu")
