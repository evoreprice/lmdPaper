#!/usr/bin/Rscript

library(data.table)

# DESeq2 results
ddsOsFile <- "output/DESeq2/ddsWald.Rds"
ddsSlFile <- "output/madsComp/ddsSl.Rds"
ddsAtFile <- "output/madsComp/ddsAt.Rds"

lapply(list(ddsAtFile, ddsSlFile, ddsOsFile), function(x)
  if(!file.exists(x)) {
    cat(x, " not found, exiting\n", file = stderr())
    quit(save = "no", status = 1)
  })

# get DESeq2 results
ddsOs <- readRDS(ddsOsFile)
resOs <- data.table(as.data.frame(
  DESeq2::results(ddsOs, contrast = c("stage", "SM", "RM"))),
  keep.rownames = TRUE, key = "rn")
ddsAt <- readRDS(ddsAtFile)
resAt <- data.table(as.data.frame(
  DESeq2::results(ddsAt, contrast = c("stage", "FM", "IM"))),
  keep.rownames = TRUE, key = "rn")
ddsSl <- readRDS(ddsSlFile)
resSl <- data.table(as.data.frame(
  DESeq2::results(ddsSl, contrast = c("stage", "fm", "sim"))),
  keep.rownames = TRUE, key = "rn")

# make DESeq results data.table
resAll <- rbind(resOs, resAt, resSl)

# set up biomaRt for phytozome
phytozome <- biomaRt::useMart(biomart = 'phytozome_mart', dataset = "phytozome")

# query biomaRt
filters <- c("panther_id_list", "organism_id")
values <- list("PTHR11945", c("167", "204", "225"))
atts <- c("transcript_name", "gene_name1", "organism_name", "synonyms", "peptide_sequence")
res <- biomaRt::getBM(atts, filters, values, phytozome, uniqueRows = TRUE)
martResults <- data.table(res, key = c("gene_name1", "transcript_name", "organism_name"))
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

# add l2fc values
madsPeptides[, log2FoldChange := resAll[rn == gene_name1, log2FoldChange], by = gene_name1]

# sneak peek at number of genes
#madsPeptides[, .(length(unlist(gene_name1))), by = organism_name]

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

# what happens if i do this with the whole alignment?
#domainAlign <- clustalAlign

# clean the alignment
alignmentLength <- nchar(domainAlign$seq[[1]])
nb <- domainAlign$nb
# minimum percentages
minpcnongap <- 30
minpcident <- 5
# number of non-gaps at each position
nbNonGaps <- sapply(1:alignmentLength, function(i)
  sum(!sapply(domainAlign$seq, substr, start = i, stop = i) == "-"))
# percentage of non-gap pairs that are identical at each position
countPcNonGapIdent <- function(letterSet) {
  # remove gap letters
  letterSet <- letterSet[!letterSet == "-"]
  # return 0 if no comparisons
  if (length(letterSet) < 2) {return(0)} 
  # how many comparisons?
  nbPairs <- length(letterSet) * (length(letterSet) - 1)/2
  # do the comparisons
  identPairs <- 0
  while(length(letterSet) > 1) {
    # how many of the pairs are identical?
    identPairs <- identPairs + sum(letterSet[1] == letterSet[-1])
    # remove the current letter
    letterSet <- letterSet[-1]
  }
  return(identPairs * 100 / nbPairs)
}
pcNonGapIdent <- sapply(1:alignmentLength, function(i)
  countPcNonGapIdent(sapply(domainAlign$seq, substr, start = i, stop = i)))
# keep the letters that match our criteria
keptLetters <- which(nbNonGaps * 100 / nb > minpcnongap & pcNonGapIdent > minpcident)
# make a new alignment
cleanedMatrix <- seqinr::as.matrix.alignment(domainAlign)[,keptLetters]
nb <- dim(cleanedMatrix)[1]
nam <- rownames(cleanedMatrix)
seq <- apply(cleanedMatrix, 1, paste0, collapse = "")
names(seq) <- NULL
cleanedAlignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq)

# make a protein tree
cDist <- seqinr::dist.alignment(domainAlign, matrix = "similarity")
hc <- hclust(cDist, method = "average")
hcdata <- lapply(dendro_data(hc), data.table)

# nj tree with ape, visualise with ggtree?
# need a function to make the tree
makeNjTree <- function(domainAlign) {
  if (class(domainAlign) == "matrix") {
    # convert matrix to alignment
    nb <- dim(domainAlign)[1]
    nam <- rownames(domainAlign)
    seq <- apply(domainAlign, 1, paste0, collapse = "")
    names(seq) <- NULL
    domainAlign <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq)
  }
  cDist <- seqinr::dist.alignment(domainAlign, matrix = "similarity")
  njTree <- ape::bionj(cDist)
  return(njTree)
}

njTree <- makeNjTree(cleanedAlignment)
njBoot <- ape::boot.phylo(njTree, seqinr::as.matrix.alignment(domainAlign), makeNjTree)
ggtree(njTree, layout = "circular") + geom_text(aes(label = label, angle = angle))
t1 <- ggtree(njTree, layout = "circular") + geom_tiplab(aes(angle = angle + 90), hjust = -0.2, vjust = 0.5)
ggtree(njTree) %>% hilight(292, fill ="red", alpha = 0.5)
t1 <- ggtree(njTree, layout = "circular") %>% collapse(node = 292)
t1
annotation_clade2(t1, tip1 = 197, tip2 = 169, label = "bof")
ggsave("~/test3.pdf", width = 20, height = 20)

ggtree(njTree) + geom_tiplab() + geom_text(aes(label = node))
t <- ggtree(njTree) + geom_tiplab()
annotation_clade2(t, tip1 = 197, tip2 = 169, label = "bof")
t %>% hilight(292, fill ="red", alpha = 0.5)
ggtree(njTree) + geom_text(aes(label = node)) %>% collapse(node = 292)
ggsave("~/test2.pdf", width = 10, height = 49)

# can i make the tree with ape?
cDistAll <- seqinr::dist.alignment(clustalAlign, matrix = "similarity")
njsTree <- ape::njs(cDistAll)
njsTree <- ape::root(njsTree, interactive = TRUE)
ape::boot.phylo(njsTree)
ape::chronopl(njsTree, 0, tol = -Inf)

# similarity histogram
hist(1-(cDist^2))

hist(1-(seqinr::dist.alignment(clustalAlign)^2))

ape::is.ultrametric(njsTree)

ape::as.hclust.phylo(njsTree)


#ggdendrogram(hc, rotate = TRUE) # sneak peek
# plot(hc)
# rect.hclust(hc, k = 16)

# add L2FCs to label
label(hcdata)[, log2FoldChange := madsPeptides[name == label, log2FoldChange], by = label]

# plot (move to figures.R)
heatscale <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
ggplot(segment(hcdata)) +
  xlab(NULL) + ylab(NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  coord_flip() +
  scale_y_reverse(limits = c(0.775, -0.2)) +
  scale_x_continuous(expand = c(0,1)) +
  geom_label(data = label(hcdata),
             mapping = aes(x = x, y = y, label = label, fill = log2FoldChange),
             hjust = "left", size = 2) + 
  guides(fill = guide_colourbar(title = expression(Log[2]*"-fold change"))) +
  scale_fill_gradient2(low = heatscale[1], mid = 'grey90', high = heatscale[5],
                       midpoint = 0, na.value = NA) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), lineend = "round") 

ggsave(filename = "~/test.eps", width = 8.3, height = 11.7 * 4) 

