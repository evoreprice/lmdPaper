#!/usr/bin/Rscript

library(Mfuzz)
library(ggplot2)
library(data.table)
library(xlsx)

# check for DESeq2 output
deseqDir <- "output/DESeq2"
if (!dir.exists(deseqDir)) {
  cat("deseqDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for expressed genes output
cutoffDir <- "output/expressedGenes"
if (!dir.exists(cutoffDir)) {
  cat("cutoffDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for msu annotation file
msuAnn.file <- "data/genome/os/all.locus_brief_info.7.0.tab"
if (!file.exists(msuAnn.file)) {
  cat("MSU annotation not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# make output folder
outDir <- "output/mfuzz"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# take the geometric mean of the vst expression values
vst <- readRDS(paste0(deseqDir, "/vst.Rds"))
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}
vstMeans.matrix <- sapply(levels(vst$stage), function(x)
  apply(GenomicRanges::assay(vst[, vst$stage == x]), 1, gm_mean))

# load the expressed gene list
expressedGenes <- readRDS(paste0(cutoffDir, '/expressedGenesAll.Rds'))

# only keep expressed genes
vstMeans <- data.frame(vstMeans.matrix)
vstFiltered <- vstMeans[expressedGenes,]

# get the most variable genes
vstByVar <- vstFiltered[(rev(order(apply(vstFiltered, 1, var)))),]
varGenes <- vstByVar[1:(0.1 * dim(vstByVar)[1]),]

# set up the expressionSet
pData <- data.frame(Stage = as.factor(colnames(varGenes)), row.names = colnames(varGenes))
phenoData <- new('AnnotatedDataFrame', data = pData)
vg.e <- ExpressionSet(assayData = as.matrix(varGenes), phenoData = phenoData)

# standardise
vg.s <- standardise(vg.e)

# estimate the fuzzifier
m1 <- mestimate(vg.s)

# estimate the cluster number
maxClust <- 25
centroids <- data.frame(
  x = 2:(maxClust),
  y = Dmin(vg.s, m = m1, crange = seq(2, maxClust, 1), repeats = 3, visu = FALSE)
)

# visualise dmin vs cluster number
centPlot <- ggplot(centroids, aes(x = x, y = y)) +
  theme_minimal() +
  xlab(expression(Cluster~number~"("*italic(c)*")")) +
  ylab("Minimum centroid distance") +
  geom_point() +
  stat_smooth(method = loess, se = FALSE)

# can try to find inflection points, doesn't work very well.
points <- seq(2, maxClust, length.out = 10000)
pred <- predict(loess(centroids$y ~ centroids$x), points)
infl <- c(FALSE, diff(diff(diff(pred)) > 0) != 0)
centPlotWithPoints <- centPlot +
  geom_point(data = data.frame(x = points[infl], y = pred[infl]), colour = 'red')

# we will go with 6 clusters
c <- 6
memCutoff <- 0.5

# run the clustering
set.seed(1)
c1 <- mfuzz(vg.s, c = c, m = m1)
clusters <- acore(vg.s, c1, min.acore = memCutoff)

# annotate the clusters for output

# take a copy of clusters
clusterExpr <- lapply(clusters, as.data.table)

# only output genes above cutoff
clusterExpr <- lapply(clusterExpr, function(x) x[MEM.SHIP >= memCutoff,])

# add ID to each dt
clusterExpr <- lapply(1:length(clusterExpr), function(i)
  clusterExpr[[i]][, Cluster := i])
# combine
clusterExpr <- do.call(rbind, clusterExpr)

# get gene names from LocToGeneName
setkey(clusterExpr, 'NAME')
annotatedClusters <- cbind(clusterExpr, oryzr::LocToGeneName(clusterExpr[, NAME], plotLabels = FALSE))

# get MSU annotation from all.locus_brief_info.7.0.tab
acNAME <- as.character(annotatedClusters[,NAME])
msuAnn <- data.table(read.delim(file = msuAnn.file, sep = "\t",header = TRUE,
                                fill = TRUE), key = 'locus')  
# data.table magic to get 1 annotation per LOC ID
anns <- msuAnn[acNAME, .(annotation = paste(unique(annotation), sep = ",")), by = locus]
annotatedClusters[, MSU.annotation := as.character(anns[, annotation])]

# rename and reorder columns
annotatedClusters[, gene_id := NAME][, membership := MEM.SHIP][, c('NAME', 'MEM.SHIP') := NULL]
setcolorder(annotatedClusters, c('gene_id', 'Cluster', 'membership', 'RapID',
                                 'symbols', 'names', 'MSU.annotation'))

# create a workbook
wb <- xlsx::createWorkbook()

# add sheet for each cluster
for(i in 1:max(annotatedClusters[,Cluster])){
  sheet <- xlsx::createSheet(wb, sheetName = paste("Cluster", i))
  xlsx::addDataFrame(annotatedClusters[Cluster == i], sheet = sheet,
                     showNA = FALSE, row.names = FALSE)
}

# write the wb
saveWorkbook(wb, "xlsx/annotatedClusters.xlsx")

# write the other output
saveRDS(vg.s, paste0(outDir, "/expressionMatrix.Rds"))
saveRDS(c1, paste0(outDir, "/c1.Rds"))
saveRDS(clusters, paste0(outDir, "/clusters.Rds"))
saveRDS(annotatedClusters, paste0(outDir, "/annotatedClusters.Rds"))
saveRDS(centPlot, paste0(outDir, "/centPlot.Rds"))
saveRDS(maxClust, paste0(outDir, "/maxClust.Rds"))
saveRDS(vstFiltered, paste0(outDir, "/vstFiltered.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
