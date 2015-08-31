#!/usr/bin/Rscript

clustalAlign <- seqinr::read.alignment('output/madsComp/clustal/madsPeptides.aln',
                                       format = 'clustal')
ggdendro::ggdendrogram(
  
  ape::njs(cDist)


cDist <- seqinr::dist.alignment(clustalAlign)
mat <- seqinr::as.matrix.alignment(clustalAlign)
ape::boot.phylo(phylo, mat, seqinr::dist.alignment)



hclust(cDist)

hc <- hclust(seqinr::dist.alignment(clustalAlign))
