#!/usr/bin/Rscript

library(data.table)

# download AGRIS AtTFDB
temp <- tempfile()
download.file("http://arabidopsis.med.ohio-state.edu/Downloads/AtTFDB.zip", temp)
atTfdb <- data.table(read.table(unz(temp, "families_data.tbl"), header = FALSE,
                                sep = "\t")[,1:8])
setnames(atTfdb, old = names(atTfdb),
         new = c("FamilyID", "LocusName", "GeneName", "Description",
                 "RNIntegration", "SubFamily", "BindingSite", "SpecialInfo"))


# MAKE OUTPUT FOLDER
outDir <- "data/attfdb"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# output
saveRDS(atTfdb, paste0(outDir, "/atTfdb.Rds"))

