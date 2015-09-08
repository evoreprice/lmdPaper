#!/usr/bin/Rscript

library(data.table)

# set up biomaRt for phytozome
phytozome <- biomaRt::useMart(biomart = 'phytozome_mart', dataset = "phytozome")

# query biomaRt
filters <- c("panther_id_list", "organism_id")
values <- list("PTHR11945", c("167", "204", "225"))
atts <- c("transcript_name", "gene_name1", "organism_name", "synonyms", "peptide_sequence")
martResults <- data.table(biomaRt::getBM(atts, filters, values, phytozome, uniqueRows = TRUE),
                          key = c( "gene_name1", "transcript_name", "organism_name"))
rm(phytozome)

# deal with multiple synonyms
formatSynonymList <- function(x) {
  xLen <- length(x)
  if(xLen == 1) {
    return(x)
  }
  return(paste0(x[1], " (", paste(x[2:xLen], collapse = "/"), ")"))
}
martResults[, synonyms := formatSynonymList(unlist(list(synonyms))),
            by = c("gene_name1", "transcript_name", "organism_name")]

# only work with first transcript
martResults[, primaryTx := sort(transcript_name)[1], by = gene_name1]
madsPeptides <- unique(martResults[transcript_name == primaryTx])

# add oryzr gene symbols to synonyms
madsPeptides[organism_name == "Osativa", synonyms := oryzr::LocToGeneName(gene_name1)$symbols]

# deal with multiple MADS peptides for one oryza name
gene_names <- c("LOC_Os03g03070","LOC_Os03g03100")
tagDuplicateSynonyms <- function(gene_names) {
  if(length(gene_names) == 1) {
    return(madsPeptides[gene_name1 == gene_names, synonyms])
  }
  madsPeptides[organism_name == "Osativa" & !is.na(synonyms) & gene_name1 %in% gene_names,
               paste(synonyms, toupper(letters)[1:length(synonyms)], sep = "_")]
}
madsPeptides[organism_name == "Osativa" & !is.na(synonyms),
             synonyms := tagDuplicateSynonyms(unlist(list(gene_name1))),
             by = synonyms]

# set up names for un-named genes and add "At" or "Os" for ambiguous genes
madsPeptides[, name := synonyms ]
madsPeptides[name == '', name := NA]
madsPeptides[!is.na(name) & organism_name == "Osativa", name := paste0("OS_", name)]
madsPeptides[!is.na(name) & organism_name == "Athaliana", name := paste0("AT_", name)]
madsPeptides[is.na(name), name := toupper(gene_name1)]

# make a list of lines for writing
madsLines <- c(apply(madsPeptides, 1, function(x)
  c(paste0(">", x["name"]), x["peptide_sequence"], "")))

# write the FASTA
writeLines(madsLines, "/tmp/madslines.fasta")

# run clustalo on generated FASTA
system("clustalo -i /tmp/madslines.fasta --full --force --outfmt=fa --outfile=/tmp/mads.faa")

# read clustal back in
clustalAlign <- seqinr::read.alignment("/tmp/mads.faa", format = 'fasta')

# convert to a matrix
clustal <- seqinr::as.matrix.alignment(clustalAlign)

# set sliding window size
window <- 140

# make a list of matrices of width window
getWindow <- function(clustal, i, window) {
  indEnd <- i + window - 1
  return(clustal[, i:indEnd])
}
startIdxs <- seq(1:(dim(clustal)[2] - window + 1))
clustalWindows <- lapply(startIdxs, function(i)
  getWindow(clustal, i, window))
names(clustalWindows) <- startIdxs

# make each window into an alignment
alignmentWindows <- lapply(clustalWindows, function(x)
  seqinr::as.alignment(nb = dim(x)[1], nam = rownames(x), seq = apply(x, 1, paste0, collapse = "")))

# run dist on each window
windowDists <- lapply(alignmentWindows, seqinr::dist.alignment, matrix = "similarity")

# sum distances (NaN -> 3)
sumDist <- function(x) {
  # take a copy (don't mutate original)
  y <- x
  y[is.nan(y)] <- 3
  return(sum(y))
}
distSums <- lapply(windowDists, function(x) sumDist(x))

# get the window with the lowest distance
startIdx <- as.integer(names(distSums[which.min(distSums)]))

# get the region + 10 residues each side
bestWindow <- clustal[, (startIdx - 10) : (startIdx + window + 9) ]

# make ungapped sequences
gappedDomains <- apply(bestWindow, 1, paste0, collapse = "")
names(gappedDomains) <- rownames(bestWindow)
domains <- sapply(gappedDomains, gsub, pattern = "-", replacement = "")

# chuck out "domains" with < 70 AAs
keptDomains <- domains[!sapply(domains, nchar) < 70]

# write a new fasta file
keptDomainsLines <- c()
for (x in names(keptDomains)) {
  keptDomainsLines <- c(keptDomainsLines, paste0(">", x), toupper(keptDomains[[x]]), "")
}

writeLines(keptDomainsLines, "/tmp/madsDomain.fasta")

# run clustalo on new FASTA
system("clustalo -i /tmp/madsDomain.fasta --full --force --outfmt=fa --outfile=/tmp/madsDomain.faa")

# read clustal back in
domainAlign <- seqinr::read.alignment("/tmp/madsDomain.faa", format = 'fasta')
cDist <- seqinr::dist.alignment(domainAlign, matrix = "similarity")
hc <- hclust(cDist, method = "average")
hcdata <- lapply(dendro_data(hc), data.table)
ggdendrogram(hc, rotate = TRUE)



