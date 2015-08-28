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

#######################
### COMPARE IN SITU ###
#######################

compare <- readRDS('output/compare/compare.Rds')
setkey(compare, "zhangRef")

# convert zhangRef to citekey (these papers were mined from zhang paper manually)
citekeys <- structure(c("@Ikeda:2007ja", "@Li:2011es", "@Ren:2013jv", "@Suzaki:2004ib", 
                        "@Chu:2006fw", "@Lee:2012hj", "@Yoshida:2013ff", "@Xue:2008ki", 
                        "@Ashikari:2005eg", "@Yan:2011hw", "@Li:2013iq", "@Yoshida:2012fg", 
                        "@Kurakawa:2007go", "@Lee:2012hj", "@Lee:2007cj", "@Gao:2010iz", 
                        "@IkedaKawakatsu:2012co", "@Horigome:2009gt", "@Lee:2007cj", 
                        "@Jiao:2010ft", "@Miura:2010it"),
                      .Names = c("56", "85", "118", 
                                 "129", "21", "82", "154", "148", "3", "150", "86", "153", "78", 
                                 "82", "83", "48", "59", "54", "83", "68", "100"))
compare[, Reference := citekeys[as.character(zhangRef)]]

# SI table
st_reviewInSitu <- data.table(reshape2::dcast(compare, msuId + Reference ~ stage,
                                              value.var = 'compare'))
st_reviewInSitu[, `Gene symbol` :=
                  oryzr::LocToGeneName(msuId, plotLabels = FALSE)$symbols,
                by = msuId]
setnames(st_reviewInSitu, old = "msuId", new = "MSU identifier")
setcolorder(st_reviewInSitu, neworder =
              c('Gene symbol', 'MSU identifier', 'RM', 'PBM', 'SBM', 'SM', 'FM',
                'Reference'))
setkey(st_reviewInSitu, "Gene symbol")
s_tableCount <- incCount(s_tableCount, "st_reviewInSitu")

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
cluster <- data.table(id = names(c1$cluster), Cluster = c1$cluster,
                      Membership = apply(c1$membership, 1, max), key = "id")

# get standardised VST counts
exprs <- data.table(Biobase::exprs(expressionMatrix), keep.rownames = TRUE)
exprs[,id := rn][, rn := NULL]
setkey(exprs, "id")

plotData.wide <- exprs[cluster[Membership > memCutoff]]
plotData.wide[, number := length(id), by = Cluster]
plotData.wide[, label := factor(paste0("Cluster ", Cluster, "\n(", number, " genes)"))]

# relevel clusters
centres.wide <- data.table(c1$centers)

# re-order the cluster for the plot
centres.wide[, Cluster := paste("Cluster", 1:nrow(centres.wide))]

# find the changes between RM and PBM and PBM and SM for each cluster
centres.wide[, c("n1n2", "n2n4") :=
               list(PBM - RM,
                    SM - PBM)]
# divide these changes into categories
centres.wide[, c("dn1n2", "dn2n4") :=
               list(c('dec', 'small', 'inc')[cut(n1n2, breaks = c(-Inf, -0.5, 0.5, Inf))],
                    c('dec', 'small', 'inc')[cut(n2n4, breaks = c(-Inf, -1, 1, Inf))])]               

# first, show gradual increase / decrease
centres.wide[dn1n2 == dn2n4, cOrder := c(1,2)[order(RM)]]

# next, big changes in n1n2 but small in n2n4
centres.wide[dn2n4 == 'small', cOrder := c(3,4)[order(RM)]]

# small changes in n1n2, then large change
centres.wide[dn1n2 == 'small', cOrder := c(5,6)[order(RM)]]

# complex patterns 
centres.wide[!dn1n2 == dn2n4 & !dn1n2 == "small" & !dn2n4 == "small",
             cOrder := c(7,8)[order(SM)]]

# order any leftovers on RM
if (any(is.na(centres.wide[, cOrder]))) {
  orderMax <- max(centres.wide[,cOrder], na.rm = TRUE)
  centres.wide[is.na(cOrder), cOrder := c((orderMax + 1):nrow(centres.wide))]
}

# relevel the clusters by cOrder
plotData.wide[, label :=
                factor(label, levels = levels(label)[order(centres.wide[,cOrder])])]

# add label to centres.wide
setkey(centres.wide, 'cOrder')
centres.wide[, label := plotData.wide[, levels(label)]]
centres.wide[, label := factor(label, levels = label)]

# make long
plotData <- reshape2::melt(plotData.wide,
                           id.vars = c("id", "Cluster", "Membership", "label", "number"),
                           variable.name = "Stage",
                           value.name = "Normalised, transformed read counts")

# fix stage label
plotData[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# set up heatscale
heatscale <- RColorBrewer::brewer.pal(n = 6, name = "YlOrRd")

# add centres to plot
centres <- reshape2::melt(centres.wide, id.vars = 'label',
                          measure.vars = c("RM", "PBM", "ePBM.SBM", "SM"),
                          variable.name = 'Stage',
                          value.name = "Normalised, transformed read counts")
centres[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# number of clusters (for text)
c <- length(unique(c1$cluster))
maxClust <- readRDS('output/mfuzz/maxClust.Rds')

# main cluster figure
f_mfuzzClusters <- ggplot(plotData,
                          aes(x = Stage, y = `Normalised, transformed read counts`,
                              colour = Membership, group = id)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = grid::unit(1, "lines")) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1), colour = "black", alpha = 0.5) +
  facet_wrap("label", ncol = 2)

figCount <- incCount(figCount, "f_mfuzzClusters")

# function for inserting cluster number in text
getClusterName <- function(deln1n2, deln2n4){
  as.numeric(centres.wide[dn1n2 == deln1n2 & dn2n4 == deln2n4,
                          sub(".*(\\d+).*", "\\1", Cluster)])
}

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
famCat <- readRDS("data/tfdb/famCat.Rds")

# format some labels
setnames(gsea, old = "Test statistic", new = "Test\nstatistic")
gsea[, Stage := plyr::mapvalues(Stage, "ePBM/SBM", "ePBM/\nSBM")]

# separate by TF / other proteins
setkey(famCat, 'Family')
setkey(gsea, 'rn')
gsea[, Category := famCat[gsea][,Category]]
gsea[, Category := plyr::mapvalues(Category, from = c("TF", "Other"),
                                   to = c("Transcription factors", "Other regulators"))]
gsea[, Category := factor(Category, levels = c("Transcription factors", "Other regulators"))]

heatscale <- rev(RColorBrewer::brewer.pal(6, "RdBu"))

f_gsea <- ggplot(gsea, aes(x = Stage, y = rn, label = padj, fill = `Test\nstatistic`)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  xlab(NULL) + ylab(NULL) +
  theme(legend.key.size = grid::unit(1, "lines"),
        axis.text.x = element_text(vjust = 0.5)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(colours = heatscale) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  geom_raster() +
  geom_text(data = gsea[showPval == TRUE], size = 2)

figCount <- incCount(figCount, "f_gsea")