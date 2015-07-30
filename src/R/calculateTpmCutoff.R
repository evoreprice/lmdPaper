#!/usr/bin/Rscript

# set variables
scriptName <- 'calculateTpmCutoff'
outputBasename <- paste(
  Sys.Date(),
  scriptName,
  sep = "-"
)

# FEATURE LENGTHS FROM GTF

# code from devon ryan:
# http://seqanswers.com/forums/showpost.php?p=129175&postcount=3

# import GTF
gtfFile <- 'data/genome/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf'
gtf <- rtracklayer::import.gff(gtfFile, format = 'gtf', genome = 'Osativa_204_v7.0', asRangedData=F, feature.type="exon")

# reduce ranges by gene_name (MSU ID), i.e. merge overlapping exons
grl <- GenomicRanges::reduce(split(gtf, elementMetadata(gtf)$gene_name))
reducedGtf <- unlist(grl, use.names = TRUE)

# add metadata
elementMetadata(reducedGtf)$gene_name <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGtf)$widths <- width(reducedGtf)

# calculate feature lengths 
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
output <- sapply(split(reducedGtf, elementMetadata(reducedGtf)$gene_name), calc_length)
gtfLength <- data.frame(Length = output, row.names = names(output))

# PARSE STAR FILES

# find the most recent cutadapt output
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
cutadaptDir <- rev(sort(outputDirs[grep('cutadapt', outputDirs)]))[1]

# find the most recent STAR output
outputDirs <- dir(path = cutadaptDir, pattern = "STAR", full.names = TRUE)
starDir <- rev(sort(outputDirs[grep('STAR', outputDirs)]))[1]

# parse the log.final.out files

# parse STAR files for the average mapped fragment length
files <- list.files(starDir, full.names = TRUE, pattern = 'final.out')

parseStarFile <- function(x) {
  lines <- readLines(x)
  return(
    as.numeric(
      # remove tab character
      gsub("\t", '', unlist(
        # split at the bar
        strsplit(
          # find the line
          lines[grep("Average mapped length", lines, fixed = TRUE)],
          split = "|", fixed = TRUE))[2], fixed = TRUE))
  )}

mu <- sapply(files, parseStarFile)
names(mu) <- gsub('.Log.final.out', '', basename(names(mu)), fixed = TRUE)

# CALCULATE TPM

# find DESeq2 output
outputDirs <- dir(path = starDir, pattern = "DESeq2", full.names = TRUE)
deseqDir <- rev(sort(outputDirs[grep('DESeq2', outputDirs)]))[1]

dds <- readRDS(paste0(deseqDir, '/ddsLrt.Rds'))
counts <- counts(dds, normalized = TRUE)

# from code at
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

# calculate efflength per sample
effLength <- sapply(colnames(counts), function(x)
  # use the gene names in the count matrix to get the effective gene length with
  # mu and gtfLength
  gtfLength[names(counts[,x]),'Length'] - mu[x] + 1
)
rownames(effLength) <- rownames(counts)

# formula for calculating the per-sample TPM using the colnames in counts
calcTpm <- function(x){
  rate <- log(counts[,x]) - log(effLength[,x])
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}

# apply formula across rownames
tpm <- sapply(colnames(counts), calcTpm)

# SAVE OUTPUT

saveRDS(tpm, )
saveRDS(gtfLength, )

# SAVE LOGS

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/", scriptName, '-', Sys.Date(), '.out')
writeLines(sInf, logLocation)

