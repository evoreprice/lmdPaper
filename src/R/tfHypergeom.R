#!/usr/bin/Rscript

# load tfdb

tfdbFile <- 'data/tfdb/os.Rds'

if (!file.exists(tfdbFile)) {
cat("tfdb file not found, exiting\n", file = stderr())
quit(save = "no", status = 1)
}