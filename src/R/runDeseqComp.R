#!/usr/bin/Rscript

# 1. SETUP

# check for STAR output
starDirSl <- "output/madsComp/sl/STAR"
if (!dir.exists(starDirSl)) {
  cat("starDirSl not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
starDirAt <- "output/madsComp/at/STAR"
if (!dir.exists(starDirAt)) {
  cat("starDirAt not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

# load the quant files
starFilesSl <- list.files(starDirSl, pattern = "ReadsPerGene", full.names = TRUE)
if (length(starFilesSl) == 0) {
  cat("Couldn't find STAR count files, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}
# load the quant files
starFilesAt <- list.files(starDirAt, pattern = "ReadsPerGene", full.names = TRUE)
if (length(starFilesAt) == 0) {
  cat("Couldn't find STAR count files, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}


# tomato
starCountsSl <- do.call(cbind,
                      lapply(starFilesSl, read.table, header = FALSE, sep = "\t", row.names = 1,
                             colClasses = c("character", "integer", rep("NULL", 2))))
colnames(starCountsSl) <- gsub("\\..*", "", basename(starFilesSl))
colDataSl <- data.frame(row.names = colnames(starCountsSl),
stage = sapply(colnames(starCountsSl), function(x) unlist(strsplit(x, split = "_", fixed = TRUE))[2]))

ddsSl <- DESeq2::DESeqDataSetFromMatrix(countData = subset(starCountsSl, !grepl("^N_", rownames(starCountsSl))),
                                      colData = colDataSl,
                                      design = ~ stage)
ddsSl <- DESeq2::DESeq(ddsSl)

# arabidopsis
starCountsAt <- do.call(cbind,
                        lapply(starFilesAt, read.table, header = FALSE, sep = "\t", row.names = 1,
                               colClasses = c("character", "integer", rep("NULL", 2))))
colnames(starCountsAt) <- gsub("\\..*", "", basename(starFilesAt))
colDataAt <- data.frame(row.names = colnames(starCountsAt), stage = sub("_.*", "", colnames(starCountsAt)))
ddsAt <- DESeq2::DESeqDataSetFromMatrix(countData = subset(starCountsAt, !grepl("^N_", rownames(starCountsAt))),
                                        colData = colDataAt,
                                        design = ~ stage)
ddsAt <- DESeq2::DESeq(ddsAt)

# MAKE OUTPUT FOLDER
outDir <- "output/DESeq2"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# SAVE OUTPUT
saveRDS(ddsSl, paste0(outDir, "/ddsSl.Rds"))
saveRDS(ddsAt, paste0(outDir, "/ddsAt.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)

