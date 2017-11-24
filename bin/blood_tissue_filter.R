#! /usr/bin/env Rscript

# get all tables in current working directory
all_tables = list.files(".", pattern = "*.tsv")

# get the corresponding TCGA barcode in the form TCGA-**-****
barcodes = unlist(lapply(all_tables, function(table){
  splits = strsplit(table, "-")[[1]]
  paste(splits[1:3], collapse = "-")
}))

# extract tissue type from table file name
types = unlist(lapply(all_tables, function(table){
  strsplit(table, "-")[[1]][4]
}))

# get duplicated barcodes (both blood and tissue)
dups = barcodes[which(duplicated(barcodes))]

# remove the table corresponding to tissue for duplicated individuals
for(d in dups){
  tissue_table = all_tables[which(barcodes == d & grepl("11", types))]
  system(paste("rm", tissue_table, sep=" "))
}
