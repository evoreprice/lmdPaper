#!/usr/bin/Rscript

library(data.table)

# load TFDBs

osTFDBfile <- "data/tfdb/os.Rds"
if (!file.exists(osTFDBfile)) {
  cat("osTFDBfile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

atTFDBfile <- "data/tfdb/at.Rds"
if (!file.exists(atTFDBfile)) {
  cat("atTFDBfile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

slTFDBfile <- "data/tfdb/puTfdbSl.Rds"
if (!file.exists(slTFDBfile)) {
  cat("slTFDBfile file not found, exiting\n", file = stderr())
  quit(save = "no", status = 1)
}

osTfdb <- readRDS(osTFDBfile)
atTfdb <- readRDS(atTFDBfile)
slTfdb <- data.table(readRDS(slTFDBfile))
slTfdb[,c("plantTFDB_id", "data_source") := NULL]
setnames(slTfdb, names(osTfdb))

# strip transcript numbers from at and sl
atTfdb[,Protein.ID := sub("\\.\\d+$", "", Protein.ID)]
slTfdb[,Protein.ID := sub("\\.\\d+$", "", Protein.ID)]

# get MADS identifiers
osMads <- osTfdb[Family == "MADS", Protein.ID]
atMads <- atTfdb[Family == "MADS", Protein.ID]
slMads <- slTfdb[Family == "MIKC", Protein.ID]

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



writeLines(madsLines, "output/test.fasta")
