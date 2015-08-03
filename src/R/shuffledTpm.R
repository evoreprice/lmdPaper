#!/usr/bin/Rscript

#SBATCH --job-name Rscript
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output log/shuffledTpm.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

library(rtracklayer)

# set variables
scriptName <- 'shuffledTpm'
outputBasename <- paste(
  scriptName,
  Sys.Date(),
  sep = "-"
)

# find the most recent shuffled gff3
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
shuffleDir <- rev(sort(outputDirs[grep('shuffle', outputDirs)]))[1]

# make output folder
outDir <- paste(shuffleDir, outputBasename, sep = "/")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# FEATURE LENGTHS FROM GTF

# code from devon ryan:
# http://seqanswers.com/forums/showpost.php?p=129175&postcount=3

# import GTF
gtfFile <- list.files(shuffleDir, pattern = 'shuffled.gff3', full.names = TRUE)
gtf <- import.gff(gtfFile, format = 'gff3', genome = 'Osativa_204_v7.0', asRangedData=F, feature.type="CDS")

# reduce ranges by gene_name (MSU ID), i.e. merge overlapping exons
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_name))
reducedGtf <- unlist(grl, use.names = TRUE)

# add metadata
elementMetadata(reducedGtf)$gene_name <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGtf)$widths <- width(reducedGtf)

# calculate feature lengths 
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
output <- OpenRepGrid::sapply_pb(split(reducedGtf, elementMetadata(reducedGtf)$gene_name), calc_length)
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

# read shuffled htseq-count results into DESeq to get normalised counts
outputDirs <- list.dirs(shuffleDir)
htseqDir <- rev(sort(outputDirs[grep('htseq-count', outputDirs)]))[1]

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
saveRDS(tpm, paste0(outDir, "/dnaTpm.Rds"))
saveRDS(gtfLength, paste0(outDir, "/shuffledGtfLength.Rds"))

# SAVE LOGS

sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/", scriptName, '-', Sys.Date(), '.out')
writeLines(sInf, logLocation)

