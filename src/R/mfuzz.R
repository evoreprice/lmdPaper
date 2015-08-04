#!/usr/bin/Rscript

#SBATCH --job-name Rscript
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output log/shuffledTpm.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

library(Mfuzz)

# set variables
scriptName <- 'mfuzz'
outputBasename <- paste(
  scriptName,
  Sys.Date(),
  sep = "-"
)

# find results
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
tpmDir <- rev(sort(outputDirs[grep('tpm-', outputDirs)]))[1]
cutadaptDir <- rev(sort(outputDirs[grep('cutadapt', outputDirs)]))[1]
outputDirs <- list.dirs(path = cutadaptDir, full.names = TRUE, recursive = FALSE)
starDir <- rev(sort(outputDirs[grep('STAR', outputDirs)]))[1]
outputDirs <- list.dirs(path = starDir, full.names = TRUE, recursive = FALSE)
deseqDir <- rev(sort(outputDirs[grep('DESeq2', outputDirs)]))[1]
outputDirs <- list.dirs(tpmDir, full.names = TRUE, recursive = FALSE)
cutoffDir <- rev(sort(outputDirs[grep('calculateCutoffs', outputDirs)]))[1]

# take the geometric mean of the vst expression values
vst <- readRDS(paste0(deseqDir, "/vst.Rds"))
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}
vstMeans.matrix <- sapply(levels(vst$stage), function(x)
  apply(GenomicRanges::assay(vst[, vst$stage == x]), 1, gm_mean))


# select genes that were significant in the LRT from DESeq2
# sigLRT <- read.table('/home/tom/Desktop/laserdissect/finalAnalysis/agriGO/LRToutput/lrtSig.txt')$V1
# varGenes <- vstMeans.matrix[sigLRT,]
# rmean <- apply(varGenes, 1, mean)
# rmax <- apply(varGenes, 1, max)
# varGenes <- varGenes[rmean > 6 | rmax > 6,]

# filter by geometric mean of tpm at cutoff
vstMeans <- data.frame(vstMeans.matrix)

expressedGenes <- readRDS('../deseq/LRT/expressedGenes.Rds')

vstFiltered <- vstMeans[expressedGenes,]
#varGenes <- vstFiltered

#vstMeans$rowmean <- rowMeans(vstMeans)
#vstMeans$rowmax <- apply(vstMeans[,c(1:4)], 1, max)
#vstFiltered <- vstMeans[vstMeans$rowmean > 8 | vstMeans$rowmax > 8, ]
#vstFiltered <- subset(vstFiltered, select = Rachis.Meristem:Spikelet.Meristem)

# variable genes across the whole dataset
vstByVar <- vstFiltered[(rev(order(apply(vstFiltered, 1, var)))),]

varGenes <- vstByVar[1:(0.25 * dim(vstByVar)[1]),]

#write(rownames(varGenes), sep = '\n', file = '/home/tom/Desktop/laserdissect/finalAnalysis/agriGO/mfuzzOutput/background.txt')

# set up the expressionSet
colnames(varGenes) <- gsub('\\.', '\n', colnames(varGenes))
pData <- data.frame(Stage = as.factor(colnames(varGenes)), row.names = colnames(varGenes))
phenoData <- new('AnnotatedDataFrame', data = pData)
vg.e <- ExpressionSet(assayData = as.matrix(varGenes), phenoData = phenoData)

# standardise and find clustering parameters
vg.s <- standardise(vg.e)
m1 <- mestimate(vg.s)
# 3 clusters of expression
Dmin(vg.s, m = m1, crange = seq(2, 20, 1), repeats = 3, visu = TRUE)
cselection(vg.s, m = m1, crange = seq(2, 20, 2), repeats = 3, visu = TRUE)

set.seed(1)
c1 <- mfuzz(vg.s, c = 6, m = m1 * 1)
clusters <- acore(vg.s, c1)

length(unique(unlist(sapply(clusters, rownames)))) ## number of genes that got clustered

ggplotClusters(clusters = c1, expressionMatrix = exprs(vg.s), memCutoff = 0.7, pointsize = 10, ncol = 3)

g <- ggplotClusters(clusters = c1, expressionMatrix = exprs(vg.s), memCutoff = 0.7, pointsize = 10, ncol = 3) +
  theme_minimal(base_size = 12) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
g
# add cluster centres to plot
centreLines <- melt(c1$centers)
centreLines$cluster <- levels(g$data$cluster)[centreLines$Var1]
gL <- g + geom_line(data = centreLines,
                    mapping = aes(x = Var2, y = value, group = 1)) +
  labs(colour = 'Membership\n')
cairo_pdf('clustersWithCentres.pdf', width = width, height = height, family = "Verdana")
gL
dev.off()

## output for aiguablava

exprs <- (data.frame(exprs(vg.s)))
data.table::setnames(exprs,
                     old = c("Rachis.Meristem",
                             "Primary.Branch.Meristem",
                             "Secondary.Branch.Meristem",
                             "Spikelet.Meristem"),
                     new = c("RM",
                             "PBM",
                             "ePBM/\nAM",
                             "SM"))
cl <- centreLines
cl$Var2 <- plyr::revalue(cl$Var2, c(
  "Rachis\nMeristem" = "RM",
  "Primary\nBranch\nMeristem" = "PBM",
  "Secondary\nBranch\nMeristem" = "ePBM/\nAM",
  "Spikelet\nMeristem" = "SM"))

g <- ggplotClusters(clusters = c1, expressionMatrix = exprs, memCutoff = 0.7, pointsize = 24, ncol = 3) +
  theme_minimal(base_size = 24, base_family = 'Lato') +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        rect = element_rect(fill = 'transparent', colour = NA)) +
        #panel.background = element_rect(fill = 'transparent', colour = NA),
        #plot.background = element_rect(fill = 'transparent', colour = NA)) +
  geom_line(data = cl, mapping = aes(x = Var2, y = value, group = 1)) +
  labs(colour = 'Membership')

# set height for AID figures
extrafont::loadfonts()
wA <- convertUnit(unit(452.038, 'mm'), unitTo = 'in', valueOnly = TRUE)
hA <- convertUnit(unit(172.945, 'mm'), unitTo = 'in', valueOnly = TRUE)
cairo_pdf('/home/tom/Dropbox/temp/clusters.pdf', width = wA, height = hA,
          family = 'Lato', bg = 'transparent')
g + theme(
  legend.text = element_text(size = 10),
  legend.text.align = 0.5,
  legend.title = element_text(size = 14),
  axis.ticks = element_blank()
  )
dev.off()

### output for knitr

#saveRDS(gL, '/home/tom/Desktop/laserdissect/finalAnalysis/reports/fuzzyClusters.Rds')

#output for DAVID

#saveRDS(clusters, 'clusters.Rds')
length(unique(unlist(sapply(clusters, rownames)))) ## number of genes that got clustered

names(clusters) <- paste0('Cluster', 1:length(clusters))
for (cluster in names(clusters)){
  write(rownames(clusters[[cluster]]), sep = '\n', file = paste0('/home/tom/Desktop/laserdissect/finalAnalysis/agriGO/mfuzzOutput/', cluster, '.txt'))
  #write(rownames(clusters[[cluster]]), sep = '\n', file = paste0('/home/tom/Desktop/laserdissect/finalAnalysis/agriGO/LRToutput/', cluster, '.txt'))
}
write(rownames(vstMeans), sep = '\n', file = '/home/tom/Desktop/laserdissect/finalAnalysis/agriGO/mfuzzOutput/background.txt')

# MDS plot for clusters

length(names(c1$cluster))

clustered <- unique(unlist(sapply(clusters, rownames)))

vg.d <- dist(exprs(vg.s)[clustered,])
vg.mds <- cmdscale(vg.d, 2)

c1$cluster[clustered]
apply(c1$membership,1,max)[clustered]

vg.mds <- data.frame(MDS1 = vg.mds[,1], MDS2 = vg.mds[,2], cluster = c1$cluster,
                     max.membership = apply(c1$membership,1,max))
MDSgeneplot <- ggplot(vg.mds, aes(x = MDS1, y=MDS2, colour = factor(cluster))) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 10),
        legend.title = element_text(face = 'plain'),
        rect = element_blank()) +
  scale_colour_brewer(palette = "Set1", name = "Cluster")

pdf(file = 'MDS.pdf', width = 3.93701, height = 3.93701) # 1 column figure: width = 3.93701, height = 3.93701
MDSgeneplot + geom_point(aes(size=max.membership),alpha=0.1) + scale_size(range = c(3,5), guide = FALSE) + coord_fixed(ratio = 1)
dev.off()

# output figures for talk

cairo_pdf(filename = 'clusters.pdf', family = 'Verdana', width = 4.055, height = 5.482)
ggplotClusters(clusters = c1, expressionMatrix = exprs(vg.s), memCutoff = 0.7, pointsize = 12, ncol = 3) + theme_minimal() + xlab(NULL) + theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.title = element_text(face = 'plain'))
dev.off()

cairo_pdf(filename = 'clustersMDS.pdf', family = 'Verdana', width = 4.055, height = 2.505)
MDSgeneplot + geom_point(aes(size=max.membership),alpha=0.7) + scale_size(range = c(1,3), guide = FALSE) + coord_fixed(ratio = 1)
dev.off()

cairo_pdf(filename = 'clustersMDS.pdf', family = 'Verdana', width = 4.055, height = 2.505)
MDSgeneplot + geom_point(aes(size=max.membership),alpha=0.7) + scale_size_area(max_size = 1.5, guide = FALSE) + coord_fixed(ratio = 1)
dev.off()

#########################
### ANNOTATE CLUSTERS ###
#########################

load('/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/cummeRbund/annotationToAdd.Rdata')

# take a copy of clusters
clusterExpr <- lapply(clusters, as.data.table)

# only output genes above cutoff
clusterExpr <- lapply(clusterExpr, function(x) x[MEM.SHIP >= 0.7,])

# add ID to each dt
clusterExpr <- lapply(1:length(clusterExpr), function(i)
  clusterExpr[[i]][, Cluster := i])
# combine
clusterExpr <- do.call(rbind, clusterExpr)

# get gene names from LocToGeneName
setkey(clusterExpr, 'NAME')
annotatedClusters <- cbind(clusterExpr, oryzr::LocToGeneName(clusterExpr[, NAME], plotLabels = FALSE))

# get annotation from annotationToAdd
acNAME <- as.character(annotatedClusters[,NAME])
annotatedClusters[, MSU.annotation := as.character(annotationToAdd[acNAME, gene_description])]

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
saveWorkbook(wb, paste0('annotatedClusters', Sys.Date(), '.xlsx'))

##############################
### TROUBLESHOOT LABELLING ###
##############################

clusterData <- data.frame(ids = names(c1$cluster), cluster = factor(c1$cluster), Membership = apply(c1$membership, 1, max), exprs(vg.s)[as.vector(names(c1$cluster)),])


####################
### OLD ANALYSIS ###
####################

vstMeans.matrix <- sapply(levels(vst$stage), function(x) apply(assay(vst[,vst$stage == x]), 1, mean))
vstMeans <- data.frame(vstMeans.matrix)
vstMeans$rowmean <- rowMeans(vstMeans)
vstMeans$rowmax <- apply(vstMeans[,c(1:4)], 1, max)

# get highly expressed genes by filtering and take the top 10% of them by
# variance

vstFiltered <- vstMeans[vstMeans$rowmean > 8 | vstMeans$rowmax > 8, ]
vstFiltered <- subset(vstFiltered, select = Rachis.Meristem:Spikelet.Meristem)

vstByVar <- vstFiltered[(rev(order(apply(vstFiltered, 1, var)))),]
varGenes <- vstByVar[1:(0.25 * dim(vstByVar)[1]),]

#########################
### OUTPUT FOR AGRIGO ###
#########################

write(paste(rownames(varGenes), collapse = '\n'), file = 'background.txt')

write(paste(sigLRT, collapse = '\n'), file = 'sigLRT.txt')
write(paste(clusters[[2]]$NAME, collapse = '\n'), file = 'cluster2.txt')
write(paste(clusters[[4]]$NAME, collapse = '\n'), file = 'cluster4.txt')


sapply(c(1:length(clusters)), function(i) 'LOC_Os10g33780' %in% clusters[[i]]$NAME)

gL

tawawa <- melt(data.frame(exprs(vg.s['LOC_Os10g33780',])))
tawawa$variable <- gsub('\\.', '\n', tawawa$variable)

['LOC_Os10g33780',]

gL + geom_line(data = tawawa,
               aes(x = variable, y = value, group = 1), colour = 'blue')

sapply(c(1:length(clusters)), function(i) 'LOC_Os04g22870' %in% clusters[[i]]$NAME)
