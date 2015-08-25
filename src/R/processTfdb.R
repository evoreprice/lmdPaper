#!/usr/bin/Rscript

library(data.table)

########################
### PROCESS THE TFDB ###
########################

# read tfdb
tfdbRaw <- data.table(read.table('data/tfdb/tfdb.tab', header = T, sep = '\t',
                                 stringsAsFactors = F))

# sort by Protein.ID and Family
setkey(tfdbRaw, 'Species', 'Family', 'Protein.ID')

# remove duplicates, store just in case
dups <- tfdbRaw[duplicated(tfdbRaw)]
tfdbRaw <- unique(tfdbRaw)

os <- unique(tfdbRaw[Species == 'Oryza sativa subsp. japonica', .(
  Protein.ID = gsub("\\..*", '', Protein.ID),
  Family
)])

at <- unique(tfdbRaw[Species == "Arabidopsis thaliana" , .(
  Protein.ID,
  Family
)])

# list of family vs. category
famCat <- unique(tfdbRaw["Oryza sativa subsp. japonica",.(Family, Category)])

# output
saveRDS(os, 'data/tfdb/os.Rds')
saveRDS(at, 'data/tfdb/at.Rds')
saveRDS(famCat, 'data/tfdb/famCat.Rds')

quit(save = "no", status = 0)
