#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)

### SETUP ----------------------------------------------------------------------

# find results

# check for shuffled GTF
shuffleDir <- "output/shuffle"
if (!dir.exists(shuffleDir)) {
  cat("shuffleDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for tpm output
tpmDir <- "output/tpm"
if (!dir.exists(tpmDir)) {
  cat("tpmDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# check for shufTpmDir output
shufTpmDir <- "output/dnaTpm"
if (!dir.exists(shufTpmDir)) {
  cat("shufTpmDir not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# make output folder
outDir <- "output/expressedGenes"
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
dnaTpm$type <- substr(rownames(dnaTpm), 1, 4)
dnaTpm$type[!dnaTpm$type == "dna_"] <- "Genic"
dnaTpm$type[dnaTpm$type == "dna_"] <- "Intergenic"

# add id column and melt
dnaTpm$id <- rownames(dnaTpm)

## remove ChrSy and ChrUn "genes"
#rem <- c(grep(toupper('ChrSy'), toupper(dnaTpm$id), fixed = TRUE),
#grep(toupper('ChrUn'), toupper(dnaTpm$id), fixed = TRUE))
#dnaTpm <- dnaTpm[-rem,]

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
q95.plot <- q95[q95$type == 'Intergenic', ]
q95.plot$type <- NULL

p <- p + geom_vline(data = q95.plot, mapping = aes(xintercept = log(q95)))

# vector of cutoff values
cutoffs <- c(q95.plot$q95)
names(cutoffs) <- q95.plot$lib

# which genes have tpm higher than q95?
expressedGenesByLibrary <- lapply(levels(q95.plot$lib), function(libName)
  subset(dnaTpm.plot, type == 'Genic' & lib == libName & tpm > cutoffs[libName])$id)
names(expressedGenesByLibrary) <- levels(q95.plot$lib)

expressedGenes <- unique(do.call(c, expressedGenesByLibrary))

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
saveRDS(expressedGenes, paste0(outDir, "/expressedGenesAll.Rds"))
saveRDS(realTpmCutoffs, paste0(outDir, "/tpmCutoffs.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)