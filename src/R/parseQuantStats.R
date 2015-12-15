#!/usr/bin/Rscript

library(data.table)

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

# check for rRNA/tRNA quant output
rnaStatsDir <- "output/quantStats"
if (!dir.exists(rnaStatsDir)) {
  cat("rnaStatsDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for detect_expressed_genes output
expGenDir <- "output/expressedGenes"
if (!dir.exists(expGenDir)) {
  cat("expGenDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# from lab book
inputRNA <- c(27.580,7.795,5.3,2.3, # N1
              33.700,7.900,3.3,5.4, # N2
              79.635,24.225,11.4, # N3
              28.115,17.285,56.3) # N4
RIN <- c(7.08,6.1,8,7.6, # N1
         7.05,5.4,8,8.7, # N2
         7,7.05,8.1, # N3
         7.08,7.05,7.0) # N4
area <- c(727647,511191,300886,848070, # N1
          600176,2735798,394641,439534, # N2
          1563411,1574719,1109344, # N3
          975315,758919,2347647) # N4
nb <- c(7, 7, 5, 4, # N1
        5, 4, 2, 3, # N2
        5, 5, 4, # N3
        5, 4, 4) # N4

# function to get required info from STAR output
parseStarInfo <- function(x){
  starName <- gsub('.Log.final.out', '', basename(x), fixed = TRUE)
  starLines <- readLines(x)
  ir <- as.numeric(unlist(strsplit(
    grep('Number of input reads', starLines, value = TRUE),
    split = "|\t", fixed = TRUE))[2])/1000000
  um <- as.numeric(unlist(strsplit(
    grep('Uniquely mapped reads number', starLines, value = TRUE),
    split = "|\t", fixed = TRUE))[2])/1000000
  mm <- as.numeric(unlist(strsplit(
    grep('Number of reads mapped to multiple loci', starLines, value = TRUE),
    split = "|\t", fixed = TRUE))[2])/1000000
  mp <- (um + mm) * 100 / ir
  return(c(starName, round(c(ir, um, mm, mp), 1)))
}

# get mapping information from STAR files
starFiles <- list.files(starDir, pattern = "Log.final.out", full.names = TRUE)
starInfo <- as.data.table(do.call(rbind, lapply(starFiles, parseStarInfo)))
setnames(starInfo, c('lib', 'Reads (M)', 'Uniquely mapped reads (M)',
                     'Multimapped reads (M)', 'Total mapping percentage'))

# rRNA / tRNA reads
quantFile <- paste0(rnaStatsDir, "/quantStats.csv")
if (!file.exists(quantFile)) {
  cat("quantFile not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
rnaStats <- as.data.table(read.csv(quantFile))
rnaStats[, rRNA := round(rRNA/1000000, 2)][, tRNA := round(tRNA/1000000, 2)]
setnames(rnaStats, c('lib', 'rRNA reads (M)', 'tRNA reads (M)'))

# reads in genes
ddsFile <- paste0(deseqDir, "/ddsLrt.Rds")
if (!file.exists(ddsFile)) {
  cat("ddsFile not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
dds <- readRDS(ddsFile)
readsPerLib <- as.data.table(round(colSums(DESeq2::counts(dds))/1000000, 1),
                             keep.rownames = TRUE)
setnames(readsPerLib, c('lib', 'Reads in genes (M)'))

# expressed genes
expGenFile <- paste0(expGenDir, "/expressedGenesByLibrary.Rds")
if (!file.exists(expGenFile)) {
  cat("expGenFile not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
expGen <- readRDS(expGenFile)
expGenPerLib <- as.data.table(sapply(expGen, function(x) length(unique(x))),
                              keep.rownames = TRUE)
setnames(expGenPerLib, c('lib', "Detected genes"))

# merge
setkey(starInfo, 'lib')
setkey(expGenPerLib, 'lib')
setkey(readsPerLib, 'lib')
setkey(rnaStats, 'lib')
libStats <- expGenPerLib[readsPerLib][rnaStats[starInfo]]
libStats[, "Yield (ng)" := round(inputRNA, 1)][, RIN := RIN ]
libStats[, "Dissected area (mm²)" := area/1000000] 
libStats[, "Number of panicles" := nb]

# make all columns except 'lib' numeric
libStats <- libStats[, lapply(.SD, as.numeric), by = lib]

# calculate rRNA/tRNA percentages
libStats[, "rRNA (%)" := round(`rRNA reads (M)`/`Reads (M)`*100, 1)]
libStats[, "tRNA (%)" := round(`tRNA reads (M)`/`Reads (M)`*100, 1)]

# library names in upper case!
libStats[, lib := toupper(lib)]

# stage and replicate number
libStats[, Sample := levels(GenomicRanges::colData(dds)$stage)[
  as.integer(gsub("N([[:digit:]+]).*", "\\1", lib))]]
libStats[, Replicate := as.integer(gsub(".*R(\\d+)", "\\1", lib))]

# arrange table
setcolorder(libStats, c('lib', "Sample", "Replicate", "Number of panicles",
                        "Dissected area (mm²)", 'Yield (ng)', 'RIN',
                        'Reads (M)', 'rRNA reads (M)', 'rRNA (%)',
                        'tRNA reads (M)', 'tRNA (%)',
                        'Uniquely mapped reads (M)', 'Multimapped reads (M)',
                        'Total mapping percentage', 'Reads in genes (M)',
                        'Detected genes'))
setnames(libStats, 'lib', "Library")

# save output
outDir <- rnaStatsDir
saveRDS(libStats, paste0(outDir, '/libStats.Rds'))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)

