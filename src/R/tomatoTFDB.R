#!/usr/bin/Rscript

temp <- tempfile()
download.file("http://planttfdb.cbi.pku.edu.cn/download/gene_model_family/Sly",
              temp)

puTfdb <- read.table(temp, header = TRUE, stringsAsFactors = FALSE)

saveRDS(puTfdb, 'data/tfdb/puTfdbSl.Rds')

quit(save = "no", status = 0)
