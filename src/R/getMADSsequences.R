#!/usr/bin/Rscript

library(data.table)

# interesting MADS family from arora paper
agl2os <- c('LOC_Os03g11614', 'LOC_Os03g54170', 'LOC_Os09g32948', 'LOC_Os06g06750', 'LOC_Os08g41950')

# arabidopsis genes
agl2at <- c("AT5G15800", "AT2G03710", "AT3G02310", "AT1G24260")

# set up biomaRt for phytozome
phytozome <- biomaRt::useMart(biomart = 'phytozome_mart', dataset = "phytozome")

# query biomaRt for tomato genes
filters <- c("gene_name_filter")
values <- list(c(agl2os, agl2at))
atts <- c("ortholog__dm_ortholog_organism_name", "ortholog__dm_ortholog_gene_name")
martResults <- data.table(biomaRt::getBM(atts, filters, values, phytozome),
                          key = "ortholog__dm_ortholog_organism_name")
agl2sl <- martResults["Slycopersicum", unique(ortholog__dm_ortholog_gene_name)]

# query biomart for protein sequences
filters <- c("gene_name_filter")
values <- list(c(agl2os, agl2at, agl2sl))
atts <- c("transcript_name", "gene_name1", "organism_name", "synonyms", "peptide_sequence")
peptResults <- data.table(biomaRt::getBM(atts, filters, values, phytozome, uniqueRows = TRUE),
                          key = c("gene_name1", "transcript_name"))

# only work on primary transripts
peptResults[, primaryTx := c(sort(transcript_name))[[1]], by = gene_name1]
# ignore multiple synonyms
madsPeptides <- unique(peptResults[transcript_name == primaryTx])
# add rice names
madsPeptides[organism_name == "Osativa", synonyms := oryzr::LocToGeneName(gene_name1)$symbols]

# set up names for un-named genes and add "At" or "Os" for ambiguous genes
madsPeptides[, name := synonyms ]
madsPeptides[name == '', name := NA]
madsPeptides[!is.na(name) & organism_name == "Osativa", name := paste0("Os_", name)]
madsPeptides[!is.na(name) & organism_name == "Athaliana", name := paste0("At_", name)]
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