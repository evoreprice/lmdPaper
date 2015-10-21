#!/usr/bin/Rscript

library(data.table)

# HB annotations from Jain et al in doc format, of course
fileUrl <- "http://onlinelibrary.wiley.com/store/10.1111/j.1742-4658.2008.06424.x/asset/supinfo/FEBS_6424_sm_Table%20S1.doc?v=1&s=e522a4da4ba65322a6ef5bfd597e4e09dbfdb977"
temp <- tempfile()
download.file(fileUrl, temp)

# need a dictionary of class names
hbClasses <- c("HD-ZIP I", "HD-ZIP II", "HD-ZIP III", "HD-ZIP IV", "BLH", 
               "KNOX I", "KNOX II", "WOX", "ZF-HD", "PHD", "Unclassified")

# read the data using catdoc (http://freecode.com/projects/catdoc)
dirtyData <- system(paste("catdoc", temp), intern = TRUE)

# remove blank lines
dirtyData <- dirtyData[!dirtyData == ""]

# find lines that match each class
classHits <- sapply(hbClasses, grep, dirtyData)
matchIndex <- data.table(
  class = names(classHits),
cStart = sapply(classHits, function(x) x[1])
)
matchIndex[, cStop := c(cStart[-1] - 1, NA)]
matchIndex[class == "Unclassified", cStop := grep("References", dirtyData) - 1]

# find genes in each set of lines
retrieveGeneIds <- function(cStart, cStop) {
  genes <- grep("Os\\d+g\\d+", dirtyData[cStart:cStop], value = TRUE)
  return(unique(gsub(".*(Os\\d+g\\d+).*", "\\1", genes)))
}

# retrieve gene IDs for each class
hbGenes <- matchIndex[, .(msuId = retrieveGeneIds(cStart, cStop)), by = class]

# Jain et al. are confused about identifiers, fix
hbGenes[, msuId := paste0("LOC_", msuId)]

# add oryzabase gene names
hbGenes[, symbol := oryzr::LocToGeneName(msuId)$symbols, by = msuId]

# make output folder
outDir <- "output/homeobox"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(hbGenes, paste0(outDir, "/hbGenes.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          paste("catdoc version:", NA),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)