#!/usr/bin/Rscript

#SBATCH --job-name Rscript
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output log/calculateCutoffs.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

library(ggplot2)
library(dplyr)

### SETUP ----------------------------------------------------------------------

# find results
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
shuffleDir <- rev(sort(outputDirs[grep('shuffle', outputDirs)]))[1]
tpmDir <- rev(sort(outputDirs[grep('tpm-', outputDirs)]))[1]
outputDirs <- list.dirs(shuffleDir, full.names = TRUE, recursive = FALSE)
shufTpmDir <- rev(sort(outputDirs[grep('shuffledTpm', outputDirs)]))[1]

scriptName <- 'calculateCutoffs'
outputBasename <- paste(
  scriptName,
  Sys.Date(),
  sep = "-"
)

# make output folder
outDir <- paste(tpmDir, outputBasename, sep = "/")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}


### SCRIPT ---------------------------------------------------------------------

# intergenic shuffle
dnaTpm <- as.data.frame(readRDS(paste0(shufTpmDir, '/dnaTpm.Rds')))
gtfLengthInt <- readRDS(paste0(shufTpmDir, '/shuffledGtfLength.Rds'))

# real results
gtfLength <- readRDS(paste0(tpmDir, '/gtfLength.Rds'))

# observe length distributions
gtfLengthInt$type <- "dna"
gtfLength$type <- "real"

combLength <- rbind(gtfLengthInt, gtfLength)

g <- ggplot(combLength, aes(x = log2(Length), colour = type, fill = type), alpha = 0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.5)

# add factor for separating distributions
head(dnaTpm)
dnaTpm$type <- substr(rownames(dnaTpm), 1, 4)
dnaTpm$type[!dnaTpm$type == "dna_"] <- "Real"
dnaTpm$type[dnaTpm$type == "dna_"] <- "DNA"

# add id column and melt
dnaTpm$id <- rownames(dnaTpm)
dnaTpm.plot <- reshape2::melt(dnaTpm, id.vars = c('id', 'type'), variable.name = 'lib', value.name = 'tpm')

p <- ggplot(dnaTpm.plot, aes(x = log(tpm), colour = type, fill = type), alpha = 0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.5) +
  facet_wrap(~lib)

# calculate 95th percentile per library
q95 <- as.data.frame(dplyr::group_by(dnaTpm.plot, type, lib) %>%
  dplyr::summarise(q95 = quantile(tpm, 0.95, type = 7)))
q95.plot <- q95[q95$type == 'DNA', ]
q95.plot$type <- NULL

p <- p + geom_vline(data = q95.plot, mapping = aes(xintercept = log(q95)))

# vector of cutoff values
cutoffs <- c(q95.plot$q95)
names(cutoffs) <- q95.plot$lib

# which genes have tpm higher than q95?
expressedGenesByLibrary <- lapply(levels(q95.plot$lib), function(libName)
  subset(dnaTpm.plot, type == 'Real' & lib == libName & tpm > cutoffs[libName])$id)
names(expressedGenesByLibrary) <- levels(q95.plot$lib)

expressedGenes <- length(unique(do.call(c, expressedGenesByLibrary)))

# what are the real cutoffs in the tpm calculations performed *without* the shuffled GTF?
realTpm <- readRDS(paste0(tpmDir, '/tpm.Rds'))

realTpmCutoff <- function(libName){
  min(subset(realTpm, rownames(realTpm) %in% expressedGenesByLibrary[[libName]], select = libName))
}
realTpmCutoffs <- sapply(colnames(realTpm), realTpmCutoff)

### OUTPUT ---------------------------------------------------------------------

saveRDS(g, paste0(outDir, "/gtfLengthDistributions.ggplot2.Rds"))
saveRDS(p, paste0(outDir, "/realAndShuffledTpmDistributions.ggplot2.Rds"))
saveRDS(expressedGenesByLibrary, paste0(outDir, "/expressedGenesByLibrary.Rds"))
saveRDS(expressedGenesAll, paste0(outDir, "/expressedGenesAll.Rds"))
saveRDS(realTpmCutoffs, paste0(outDir, "/tpmCutoffs.Rds"))

sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
writeLines(sInf, paste0(outDir, "/log.out"))