#!/usr/bin/Rscript

library(data.table)

args <- commandArgs(TRUE)
outdir <- args[1]

# find the gtfLength
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
tpmDir <- rev(sort(outputDirs[grep('tpm', outputDirs)]))[1]

# import gtfLength
gtfLength <- data.table(readRDS(paste0(tpmDir, "/gtfLength.Rds")), keep.rownames = TRUE)

# import GTF
gtfFile <- "data/genome/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
gtf <- rtracklayer::import.gff(gtfFile, format = 'gtf', genome = 'Osativa_204_v7.0', asRangedData=F, feature.type="exon")

# make a dummy gff file with each 'gene' on the right chromosome but with the 
# coordinates 1 --> CDS length. This will be shuffled by bedtools shuffle so
# that the intergenic windows are the same size as acual genes.

##### NOT WORKING!!! TEST BELOW

# take ss of gtfLength

set.seed(1)
pick <- sample(1:dim(gtfLength)[1], size = 1000, replace = FALSE)

ss <- gtfLength[pick,]

ss[,.(
  seqid = unique(as.character(GenomeInfoDb::seqnames(gtf[gtf$gene_name == rn,]))),
  source = 'phytozomev10',
  type = 'CDS',
  start = 1,
  end = Length,
  score = ".",
  strand = unique(as.character(rtracklayer::strand(gtf[gtf$gene_name == rn,]))),
  phase = ".",
  attributes = paste0("ID=", rn)
), by = rn]

# Uses by = rn to collect the duplicate seqids and strands that are unhelpfully
# returned from the GRanges object, then a second fast grouping to only keep the
# first of each.
dummyGff <- gtfLength[,.(
  seqid = as.character(GenomeInfoDb::seqnames(gtf[gtf$gene_name == rn,])),
  source = 'phytozomev10',
  type = 'CDS',
  start = 1,
  end = Length,
  score = ".",
  strand = as.character(rtracklayer::strand(gtf[gtf$gene_name == rn,])),
  phase = ".",
  attributes = paste0("ID=", rn)
), by = rn][,.(seqid = seqid[1], strand = strand[1]), by = rn]


saveRDS(dummyGff, paste0(outdir, "/dummyGff.Rds"))

write.table(dummyGff, paste0(outdir, "/dummyGff.gff3"), sep = '\t', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

