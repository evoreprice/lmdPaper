#!/usr/bin/Rscript

# set variables
scriptName <- 'runDeseq'
outputBasename <- paste(
  Sys.Date(),
  scriptName,
  sep = "-"
)

# 1. SETUP

# find the most recent cutadapt output
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
cutadaptDir <- rev(sort(outputDirs[grep('cutadapt', outputDirs)]))[1]

# find the most recent STAR output
outputDirs <- dir(path = cutadaptDir, pattern = "STAR", full.names = TRUE)
starDir <- rev(sort(outputDirs[grep('STAR', outputDirs)]))[1]

# load the quant files
starFiles <- list.files(starDir, pattern = "ReadsPerGene", full.names = TRUE)
starCounts <- do.call(cbind,
                      lapply(starFiles, read.table, header = FALSE, sep = "\t", row.names = 1,
                             colClasses = c("character", "integer", rep("NULL", 2))))
colnames(starCounts) <- gsub("\\..*", "", basename(starFiles))

# remove n1r2 and n2r2 (degraded RNA)
starCounts$n1r2 <- NULL
starCounts$n2r2 <- NULL

# generate colData
sample <- colnames(starCounts)
sampleToStage <- c(
  n1 = "RM",
  n2 = "PBM",
  n3 = "ePBM/SBM",
  n4 = "RM"
)
stage <- sampleToStage[substr(sample, 1, 2)]
batch <- as.numeric(substr(sample, 4, 4))
batch[batch <= 2] <- 'first'
batch[!batch == 'first'] <- 'second'
colData <- data.frame(row.names = sample, stage, batch)

# generate DESeq2 object (exclude the STAR diagnostic lines starting with "N")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = subset(starCounts, !grepl("^N_", rownames(starCounts))),
                                      colData = colData,
                                      design = ~ stage + batch)

# RUN DESEQ2

dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds)

# RUN LRT

ddsLrt <- DESeq2::nbinomLRT(dds, reduced = ~ stage, betaPrior = FALSE, maxit = 1000)

# RUN WALD TESTS

ddsWald <- DESeq2::nbinomWaldTest(dds)

# RUN TRANSORMATIONS

vst <- DESeq2::varianceStabilizingTransformation(ddsLrt)
rld <- DESeq2::rlogTransformation(dds)

# MAKE OUTPUT FOLDER

outDir <- paste0(starDir, "/DESeq2-", Sys.Date())
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT

saveRDS(ddsLrt, paste0(outDir, "/ddsLrt.Rds"))
saveRDS(ddsWald, paste0(outDir, "/ddsWald.Rds"))
saveRDS(vst, paste0(outDir, "/vst.Rds"))
saveRDS(rld, paste0(outDir, "/rld.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/", scriptName, '-', Sys.Date(), '.out')
writeLines(sInf, logLocation)
