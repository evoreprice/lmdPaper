#!/usr/bin/Rscript

library(ggdendro)
library(ggplot2)
library(data.table)

clustalAlign <- seqinr::read.alignment('output/madsComp/clustal/madsPeptides.aln',
                                       format = 'clustal')
cDist <- seqinr::dist.alignment(clustalAlign, matrix = "similarity")
hc <- hclust(cDist, method = "average")
hcdata <- lapply(dendro_data(hc), data.table)
#segment(hcdata)[yend == 0, yend := y - 0.1]


ggdendrogram(hc, rotate = TRUE)

ggplot(segment(hcdata)) +
  theme_dendro() +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  geom_label(data = label(hcdata), mapping = aes(x = x, y = y, label = label),
             hjust = 0) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) 
  

plot(ape::as.phylo(hc))
plot(ape::as.phylo(hc), type = "unrooted")
