#!/usr/bin/Rscript

library(data.table)

# set up biomaRt for phytozome
phytozome <- biomaRt::useMart(biomart = 'phytozome_mart', dataset = "phytozome")

# query biomaRt
filters <- c("panther_id_list", "organism_id")
values <- list("PTHR11945", c("167", "204", "225"))
atts <- c("transcript_name", "gene_name1", "organism_name", "synonyms", "peptide_sequence")
martResults <- biomaRt::getBM(atts, filters, values, phytozome, uniqueRows = TRUE)

# only deal with first transcript
madsPeptides <- unique(data.table(martResults, key = c("gene_name1", "transcript_name")))
madsPeptides[, primaryTx := c(sort(transcript_name))[[1]], by = gene_name1]
madsPeptides <- madsPeptides[transcript_name == primaryTx]

# add oryzr gene symbols to synonyms
madsPeptides[organism_name == "Osativa", synonyms := oryzr::LocToGeneName(gene_name1)$symbols]

# set up names for un-named genes and add "At" or "Os" for ambiguous genes
madsPeptides[, name := synonyms ]
madsPeptides[name == '', name := NA]
madsPeptides[!is.na(name) & organism_name == "Osativa", name := paste("Os", name)]
madsPeptides[!is.na(name) & organism_name == "Athaliana", name := paste("At", name)]
madsPeptides[is.na(name), name := gene_name1]

# make a list of lines for writing
madsLines <- c(apply(madsPeptides, 1, function(x)
  c(paste0(">", x["name"]), x["peptide_sequence"], "")))

# MAKE OUTPUT FOLDER
outDir <- "output/madsComp/clustal"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

writeLines(madsLines, paste0(outDir, "/madsPeptides.fasta"))
saveRDS(madsPeptides, paste0(outDir, "/madsPeptides.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)