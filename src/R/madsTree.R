#!/usr/bin/Rscript

library(ggdendro)
library(ggplot2)
library(data.table)

# DESeq2 results
ddsOsFile <- "output/DESeq2/ddsWald.Rds"
ddsSlFile <- "output/madsComp/ddsSl.Rds"
ddsAtFile <- "output/madsComp/ddsAt.Rds"
# clustal alignment of MADS peptides
clustalAlignFile <- 'output/madsComp/clustal/madsPeptides.aln'
# mads peptides
madsPeptidesFile <- "output/madsComp/clustal/madsPeptides.Rds"

lapply(list(ddsAtFile, ddsSlFile, ddsOsFile, clustalAlignFile, madsPeptidesFile), function(x)
  if(!file.exists(x)) {
    cat(x, " not found, exiting\n", file = stderr())
    quit(save = "no", status = 1)
  })

# get DESeq2 results
ddsOs <- readRDS(ddsOsFile)
resOs <- data.table(as.data.frame(
  DESeq2::results(ddsOs, contrast = c("stage", "SM", "RM"))),
  keep.rownames = TRUE, key = "rn")
ddsAt <- readRDS(ddsAtFile)
resAt <- data.table(as.data.frame(
  DESeq2::results(ddsAt, contrast = c("stage", "FM", "IM"))),
  keep.rownames = TRUE, key = "rn")
ddsSl <- readRDS(ddsSlFile)
resSl <- data.table(as.data.frame(
  DESeq2::results(ddsSl, contrast = c("stage", "fm", "sim"))),
  keep.rownames = TRUE, key = "rn")

# make DESeq results data.table
resAll <- rbind(resOs, resAt, resSl)

# get LFCs for aligned peptides
madsPeptides <- readRDS(madsPeptidesFile)
madsPeptides[, log2FoldChange := resAll[rn == gene_name1, log2FoldChange], by = gene_name1]
#remove LFC for outgroup
madsPeptides["LOC_Os03g14850", log2FoldChange := NA]

# make tree
clustalAlign <- seqinr::read.alignment(clustalAlignFile, format = 'clustal')
cDist <- seqinr::dist.alignment(clustalAlign, matrix = "similarity")
hc <- hclust(cDist, method = "average")
hcdata <- lapply(dendro_data(hc), data.table)

# add lfc to labels, hack species names
label(hcdata)[, log2FoldChange := madsPeptides[name == label, log2FoldChange], by = label]
label(hcdata)[, label := sub("\\(", " (", label)]

ggdendrogram(hc, rotate = TRUE)

heatscale <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
ggplot(segment(hcdata)) +
  xlab(NULL) + ylab(NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  coord_flip() +
  scale_y_reverse(limits = c(0.65, -0.12)) +
  geom_label(data = label(hcdata),
             mapping = aes(x = x, y = y, label = label, fill = log2FoldChange),
             hjust = "left", size = 2) + 
  guides(fill = guide_colourbar(title = expression(Log[2]*"-fold change"))) +
  scale_fill_gradient2(low = heatscale[1], mid = 'grey90', high = heatscale[5],
                       midpoint = 0, na.value = NA) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), lineend = "round") 

ggsave(filename = "~/test.pdf", width = 3.150 * 2, height = 3.150)

DESeq2::plotCounts(ddsAt, "AT3G02310", intgroup = "stage")
DESeq2::plotCounts(ddsSl, "Solyc05g015750.2", intgroup = "stage")
DESeq2::plotCounts(ddsOs, "LOC_Os03g11614", intgroup = "stage")
DESeq2::plotCounts(ddsOs, "LOC_Os03g54170", intgroup = "stage")
resOs["LOC_Os03g54170",]


plot(ape::as.phylo(hc))
plot(ape::as.phylo(hc), type = "unrooted")
