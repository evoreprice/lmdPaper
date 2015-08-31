#!/usr/bin/Rscript

clustalAlign <- seqinr::read.alignment('output/madsComp/clustal/madsPeptides.aln',
                                       format = 'clustal')
#ggdendro::ggdendrogram(

ape::as.alignment(clustalAlign)

cDist <- seqinr::dist.alignment(clustalAlign, matrix = "similarity")
mat <- as.matrix(cDist)
which(apply(mat, 2, function(x) any(is.nan(x))))


mat[, "AT1G59920"]

phy <- ape::njs(cDist)



mat <- seqinr::as.matrix.alignment(clustalAlign)
ape::boot.phylo(phylo, mat, seqinr::dist.alignment)



hclust(cDist, "average")

hc <- hclust(seqinr::dist.alignment(clustalAlign))
