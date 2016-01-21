#!/usr/bin/Rscript

library(data.table)

# load tpm data
data.s1 <- data.table(readRDS('output/tpm/tpm.Rds'), keep.rownames = TRUE,
                      key = "rn")

# make combined annotation (this should be part of oryzr)
annotation.data <- data.table(data.s1[, oryzr::LocToGeneName(rn)],
                              keep.rownames = TRUE, key = "rn")
msu.annotation <- data.table(read.delim(
  file = "data/genome/os/all.locus_brief_info.7.0.tab", sep = "\t",
  header = TRUE, fill = TRUE), key = 'locus')  
annotation.data <- msu.annotation[annotation.data, .(
  locus,
  CGSNL.symbol = symbols,
  CGSNL.name = names,
  annotation)]
annotation.data[, TIGR.annotation := paste(unique(annotation), collapse = ","),
                by = locus]
annotation.data[, annotation := NULL]
annotation.data <- unique(annotation.data)

# fast join with tpm
data.s1 <- data.s1[annotation.data, .(
  msuId = rn,
  n1r1, n1r3, n1r4, n2r1, n2r3, n2r4, n3r1, n3r2, n3r3, n4r1, n4r2, n4r3,
  CGSNL.symbol, CGSNL.name, TIGR.annotation
)]

saveRDS(data.s1, "output/tpm/datas1.Rds")

quit(save = "no", status = 0)
