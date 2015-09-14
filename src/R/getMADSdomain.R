#!/usr/bin/Rscript

library(data.table)
set.seed(1)

# how many CPUs we got?
SLURM_NTASKS <- as.integer(Sys.getenv("SLURM_NTASKS"))
if(!is.na(SLURM_NTASKS)) {
  maxCpus <- SLURM_NTASKS
} else {
  maxCpus <- 1
}

# DESeq2 results
ddsOsFile <- "output/DESeq2/ddsWald.Rds"
ddsSlFile <- "output/madsComp/ddsSl.Rds"
ddsAtFile <- "output/madsComp/ddsAt.Rds"
# atTfdb
atTfdbFile <- "data/attfdb/atTfdb.Rds"

lapply(list(ddsAtFile, ddsSlFile, ddsOsFile, atTfdbFile), function(x)
  if(!file.exists(x)) {
    cat(x, " not found, exiting\n", file = stderr())
    quit(save = "no", status = 1)
  })

# MAKE OUTPUT FOLDER
outDir <- "output/madsComp/clustal"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

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

# set up names for un-named genes and add "At" or "Os" to named proteins
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
allMads.file <- paste0(outDir, "/allMads.fasta")
writeLines(madsLines, allMads.file)

# run clustalo on generated FASTA
allMads.clustal.file <- paste0(outDir, "/allMads.faa")
cmd <- paste0("clustalo --full --force --outfmt=fa", " --threads=", maxCpus,
              " -i ", allMads.file, " --outfile=", allMads.clustal.file)
system(cmd)

# read clustal back in
allMads.clustal <- seqinr::read.alignment(allMads.clustal.file, format = "fasta")

# function to clean the alignment
cleanAlignment <- function(myAlignment, minpcnongap, minpcident) {
  alignmentLength <- nchar(myAlignment$seq[[1]])
  nb <- myAlignment$nb
  # number of non-gaps at each position
  nbNonGaps <- sapply(1:alignmentLength, function(i)
    sum(!sapply(myAlignment$seq, substr, start = i, stop = i) == "-"))
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
    countPcNonGapIdent(sapply(myAlignment$seq, substr, start = i, stop = i)))
  # keep the letters that match our criteria
  keptLetters <- which(nbNonGaps * 100 / nb > minpcnongap & pcNonGapIdent > minpcident)
  # make a new alignment
  cleanedMatrix <- seqinr::as.matrix.alignment(myAlignment)[,keptLetters]
  nb <- dim(cleanedMatrix)[1]
  nam <- rownames(cleanedMatrix)
  seq <- apply(cleanedMatrix, 1, paste0, collapse = "")
  names(seq) <- NULL
  cleanedAlignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq)
  return(cleanedAlignment)
}

# minimum percentages
minpcnongap <- 30
minpcident <- 5

# clean
allMads.cleaned <- cleanAlignment(allMads.clustal, minpcnongap, minpcident)

# find the clade of the alignment with the most type ii mads genes
# lists of type I and type II names
setkey(madsPeptides, "gene_name1")
atTfdb <- readRDS(atTfdbFile)
typeI <-atTfdb[FamilyID == "MADS" & SubFamily == "TypeI", unique(toupper(LocusName))]
typeII <- atTfdb[FamilyID == "MADS" & SubFamily == "TypeII", unique(toupper(LocusName))]
t1names <- madsPeptides[typeI, unique(name[!is.na(name)])]
t2names <- madsPeptides[typeII, unique(name[!is.na(name)])]
# distance and UPGMA tree
allMads.dist <- seqinr::dist.alignment(allMads.cleaned, matrix = "similarity")
allMads.dist[is.na(allMads.dist)] <- 0
allMads.hc <- hclust(allMads.dist, method = "average")
# get the lowest cut with as many type II mads as possible
optimClades <- function(h, hc) {
  clades <- cutree(hc, h = h)
  # which clade has the most type II mads
  mikccClade <- as.integer(names(which.max(table(clades[names(clades) %in% t2names]))))
  # list of genes in this clade
  genesInClade <- unique(names(clades[clades == mikccClade]))
  # return number of type II mads minus number of type I mads in this clade
  return(sum(genesInClade %in% t2names) - sum(genesInClade %in% t1names)/10)
}
hMin <- optimize(f = optimClades, interval = c(0,1), hc = allMads.hc,
                 maximum = TRUE, tol = 10^-10)

# get the genes in this clade
clades <- cutree(allMads.hc, h = hMin$maximum)
mikcClade <- as.integer(names(which.max(table(clades[names(clades) %in% t2names]))))
genesInClade <- unique(names(clades[clades == mikcClade]))

# get the clustal matrix for these genes
allMads.matrix <- seqinr::as.matrix.alignment(allMads.clustal)
mikcMatrix <- allMads.matrix[genesInClade,]

# make ungapped sequences
gappedMikc <- apply(mikcMatrix, 1, paste0, collapse = "")
names(gappedMikc) <- rownames(mikcMatrix)
mikc <- sapply(gappedMikc, gsub, pattern = "-", replacement = "")
# chuck out proteins with < 175 AAs
minProtLength <- 200
#maxProtLength <- 400
keptMikc <- mikc[sapply(mikc, nchar) >= minProtLength]

############################################
### DISABLE SLIDING WINDOW STUFF FOR NOW ###
############################################
# # set sliding window size
# window <- 140
# 
# # make a list of matrices of width window
# getWindow <- function(clustal, i, window) {
#   indEnd <- i + window - 1
#   return(clustal[, i:indEnd])
# }
# startIdxs <- seq(1:(dim(clustal)[2] - window + 1))
# clustalWindows <- lapply(startIdxs, function(i)
#   getWindow(clustal, i, window))
# names(clustalWindows) <- startIdxs
# 
# # make each window into an alignment
# alignmentWindows <- lapply(clustalWindows, function(x)
#   seqinr::as.alignment(nb = dim(x)[1], nam = rownames(x), seq = apply(x, 1, paste0, collapse = "")))
# 
# # run dist on each window
# windowDists <- lapply(alignmentWindows, seqinr::dist.alignment, matrix = "similarity")
# 
# # sum distances (NaN -> 3)
# sumDist <- function(x) {
#   # take a copy (don't mutate original)
#   y <- x
#   y[is.nan(y)] <- 3
#   return(sum(y))
# }
# distSums <- lapply(windowDists, function(x) sumDist(x))
# 
# # get the window with the lowest distance
# startIdx <- as.integer(names(distSums[which.min(distSums)]))
# 
# # get the region + 10 residues each side
# bestWindow <- clustal[, (startIdx - 10) : (startIdx + window + 9) ]
# 
# # make ungapped sequences
# gappedDomains <- apply(bestWindow, 1, paste0, collapse = "")
# names(gappedDomains) <- rownames(bestWindow)
# domains <- sapply(gappedDomains, gsub, pattern = "-", replacement = "")
# 
# # density of domain size?
# #dLength <- data.table(domainLength = sapply(domains, nchar))
# #ggplot(dLength, aes(x = domainLength)) + stat_density() +
# #  geom_vline(xintercept = dLength[, quantile(domainLength, 0.15)], colour = "red") +
# #  geom_vline(xintercept = 70, colour = "blue")
# #domains[sapply(domains, nchar) > 70 & sapply(domains, nchar) < dLength[, quantile(domainLength, 0.15)]]
# 
# # chuck out "domains" with < 70 AAs
# keptDomains <- domains[!sapply(domains, nchar) < 70]

##############################################
### SWITCH BELOW BACK TO KEPTDOMAIN IF NEC ###
##############################################

# write a new fasta file
keptMikcLines <- c()
for (x in names(keptMikc)) {
  keptMikcLines <- c(keptMikcLines, paste0(">", x), toupper(keptMikc[[x]]), "")
}
mikcProteins.file <- paste0(outDir, "/mikcProteins.fasta")
writeLines(keptMikcLines, mikcProteins.file)

# run clustalo on new FASTA
mikc.clustal.file <- paste0(outDir, "/mikcProteins.faa")
cmd <- paste0("clustalo --iterations=10 --outfmt=fa --force", " --threads=",
              maxCpus, " -i ", mikcProteins.file, " --outfile=", mikc.clustal.file)
system(cmd)

# read clustal back in
mikc.clustal <- seqinr::read.alignment(mikc.clustal.file, format = 'fasta')

# clean
mikcCleaned <- cleanAlignment(mikc.clustal, minpcnongap, minpcident)

##############################################
### disable second MIKCc selection for now ###
##############################################

# #select MIKCc clade (plus outgroup)
# cDist <- seqinr::dist.alignment(cleanedAlignment, matrix = "similarity")
# hc <- hclust(cDist, method = "average")
# clades <- cutree(hc, h = 0.475)
# mikcc <- clades[clades == which.max(table(clades))]
# # will use AGL33 to root tree
# mikccMatrix <- cleanedMatrix[c(names(mikcc), "AT_AGL33"),]
# nb <- dim(mikccMatrix)[1]
# nam <- rownames(mikccMatrix)
# seq <- apply(mikccMatrix, 1, paste0, collapse = "")
# names(seq) <- NULL
# mikccAlignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq)

# nj tree with ape, visualise with ggtree
# need a function to make the tree
makeNjTree <- function(myAlignment, outgroup) {
  if (class(myAlignment) == "matrix") {
    # convert matrix to alignment
    nb <- dim(myAlignment)[1]
    nam <- rownames(myAlignment)
    seq <- apply(myAlignment, 1, paste0, collapse = "")
    names(seq) <- NULL
    myAlignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq)
  }
  cDist <- seqinr::dist.alignment(myAlignment, matrix = "similarity")
  njTree <- ape::bionj(cDist)
  rootedNjt <- ape::root(njTree, outgroup, resolve.root = TRUE)
  return(rootedNjt)
}
# choose an outgroup
# og <- "AT_AGL41"
og <- mikcCleaned$nam[mikcCleaned$nam %in% t1names][1]
njTree <- makeNjTree(mikcCleaned, og)

# draw a tree (move to figures)
library(ggplot2)
library(ggtree)

heatscale <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
gtree <- ggplot(njTree, aes(x = x, y = y, label = label)) +
  xlab(NULL) + ylab(NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_gradient2(low = heatscale[1], mid = 'grey90', high = heatscale[5],
                       midpoint = 0, na.value = "white") +
  geom_tree()
# add expression values as an annotation
setkey(madsPeptides, "name")
exprAnnot <- madsPeptides[unique(njTree$tip.label), .(name, log2FoldChange)]
gtree <- gtree %<+% exprAnnot
gtree <- gtree + geom_label(mapping = aes(fill = log2FoldChange), size = 2) +
  scale_y_continuous(expand = c(0,1))
# gtree <- annotation_clade(gtree, node = 151, "AGL2-like")
# gtree <- annotation_clade(gtree, node = 144, "AGL6-like")
# annotation_clade(gtree, node = 125, "SQUA-like")
gtree
ggsave(filename = paste0(outDir, "/tempTree.pdf"), width = 8.3, height = 11.7 * 2) 

# 
# # similarity histogram
# hist(1-(cDist^2))
# 
# hist(1-(seqinr::dist.alignment(clustalAlign)^2))
# 
# ape::is.ultrametric(njsTree)
# 
# ape::as.hclust.phylo(njsTree)
# 
# 
# #ggdendrogram(hc, rotate = TRUE) # sneak peek
# # plot(hc)
# # rect.hclust(hc, k = 16)
# 
# # add L2FCs to label
# label(hcdata)[, log2FoldChange := madsPeptides[name == label, log2FoldChange], by = label]
# 
# # plot (move to figures.R)
# heatscale <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
# ggplot(segment(hcdata)) +
#   xlab(NULL) + ylab(NULL) +
#   theme_minimal(base_size = 8, base_family = "Helvetica") +
#   theme(axis.text = element_blank(),
#         panel.grid = element_blank()) +
#   coord_flip() +
#   scale_y_reverse(limits = c(0.775, -0.2)) +
#   scale_x_continuous(expand = c(0,1)) +
#   geom_label(data = label(hcdata),
#              mapping = aes(x = x, y = y, label = label, fill = log2FoldChange),
#              hjust = "left", size = 2) + 
#   guides(fill = guide_colourbar(title = expression(Log[2]*"-fold change"))) +
#   scale_fill_gradient2(low = heatscale[1], mid = 'grey90', high = heatscale[5],
#                        midpoint = 0, na.value = NA) +
#   geom_segment(aes(x=x, y=y, xend=xend, yend=yend), lineend = "round") 
# 
# ggsave(filename = "~/test.eps", width = 8.3, height = 11.7 * 4) 

# OUTPUT TO SAVE
saveRDS(njTree, paste0(outDir, "/njTree.Rds"))
saveRDS(madsPeptides, paste0(outDir, "/madsPeptides.Rds"))
saveRDS(minpcnongap, paste0(outDir, "/minpcnongap.Rds"))
saveRDS(minpcident, paste0(outDir, "/minpcident.Rds"))
saveRDS(minProtLength, paste0(outDir, "/minProtLength.Rds"))
saveRDS(og, paste0(outDir, "/og.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          paste("clustalo version:", system("clustalo --version", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
