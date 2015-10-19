#!/usr/bin/Rscript

library(rtracklayer)
library(dplyr)

# check for shuffled GTF
shuffleDir <- "output/shuffle"
if (!dir.exists(shuffleDir)) {
  cat("shuffleDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for htseq counts
htseqDir <- "output/dnaTpm/shuffledCounts"
if (!dir.exists(htseqDir)) {
  cat("htseqDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for STAR output
starDir <- "output/STAR"
if (!dir.exists(starDir)) {
  cat("starDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for DESeq2 output
deseqDir <- "output/DESeq2"
if (!dir.exists(deseqDir)) {
  cat("deseqDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for tpm output
tpmDir <- "output/tpm"
if (!dir.exists(tpmDir)) {
  cat("tpmDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# FEATURE LENGTHS FROM GTF

# code from devon ryan:
# http://seqanswers.com/forums/showpost.php?p=129175&postcount=3

# import GTF
gtfFile <- list.files(shuffleDir, pattern = 'shuffled.gff3', full.names = TRUE)
gtf <- import.gff(gtfFile, format = 'gff3', genome = 'Osativa_204_v7.0',
                  feature.type="CDS")

# reduce ranges by gene_name (MSU ID), i.e. merge overlapping exons
grl <- reduce(split(gtf, elementMetadata(gtf)$ID))
reducedGtf <- unlist(grl, use.names = FALSE)

# add metadata
elementMetadata(reducedGtf)$gene_name <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGtf)$widths <- width(reducedGtf)

# calculate feature lengths with dplyr
output <- group_by(as.data.frame(reducedGtf), gene_name) %>%
  summarize(length = sum(widths))
gtfLength <- data.frame(Length = output$length, row.names = output$gene_name)

# PARSE STAR FILES

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

countFiles <- list.files(htseqDir, pattern = 'htseq-count$', full.names = TRUE)
countMatrix <- do.call(cbind, lapply(countFiles, read.table, header = FALSE, sep = "\t", row.names = 1))
colnames(countMatrix) <- gsub("\\..*", "", basename(countFiles))

# remove n1r2 and n2r2 (degraded RNA)
countMatrix$n1r2 <- NULL
countMatrix$n2r2 <- NULL

# generate colData
sample <- colnames(countMatrix)
sampleToStage <- c(
  n1 = "RM",
  n2 = "PBM",
  n3 = "ePBM/SBM",
  n4 = "RM"
)
stage <- sampleToStage[substr(sample, 1, 2)]
colData <- data.frame(row.names = sample, stage)

# generate DESeq2 object (exclude the htseq-count diagnostic lines starting with "__")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = subset(countMatrix, !grepl("^__", rownames(countMatrix))),
                                      colData = colData,
                                      design = ~ stage)
dds <- DESeq2::estimateSizeFactors(dds)
counts <- DESeq2::counts(dds, normalized = TRUE)
  
# from code at
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

# calculate efflength per sample
effLength <- sapply(colnames(counts), function(x)
  # use the gene names in the count matrix to get the effective gene length with
  # mu and gtfLength
  gtfLength[names(counts[,x]),'Length'] - mu[x] + 1
)
rownames(counts) <- paste0("dna_", rownames(counts))
rownames(effLength) <- rownames(counts)

# combine with the counts and efflength dataframes from the real tpm
realDds <- readRDS(paste0(deseqDir, '/ddsLrt.Rds'))
realCounts <- DESeq2::counts(realDds, normalized = TRUE)

realEffLength <- readRDS(paste0(tpmDir, '/effLength.Rds'))

counts <- rbind(counts, realCounts)
effLength <- rbind(effLength, realEffLength)

# formula for calculating the per-sample TPM using the colnames in counts
calcTpm <- function(x){
  rate <- log(counts[,x]) - log(effLength[,x])
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}

# apply formula across rownames
tpm <- sapply(colnames(counts), calcTpm)

# make output folder
outDir <- "output/dnaTpm"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(tpm, paste0(outDir, "/dnaTpm.Rds"))
saveRDS(gtfLength, paste0(outDir, "/shuffledGtfLength.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
