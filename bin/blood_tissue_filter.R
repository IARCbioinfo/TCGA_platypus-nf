#! /usr/bin/env Rscript

source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
source("blood_tissue_filter_function.R")
library(foreach)
library(doParallel)

# get the parameters
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$min_blood_QUAL)) args$min_blood_QUAL = 20 #default in platypus is 20 for PASS variant
if(is.null(args$min_tissue_QUAL)) args$min_tissue_QUAL = 20 #default in platypus is 20 for PASS variant
if(is.null(args$cpu)) args$cpu = 1

min_blood_QUAL = as.numeric(args$min_blood_QUAL)
min_tissue_QUAL = as.numeric(args$min_tissue_QUAL)

# get all tables in current working directory
all_tables = list.files(".", pattern = "*.tsv")

# get the corresponding TCGA barcode in the form TCGA-**-****
barcodes = unlist(lapply(all_tables, function(table){
  substr(table, regexpr("TCGA", table)[1], regexpr("TCGA", table)[1] + 11)
}))

# extract tissue type from table file name
tissue_type = unlist(lapply(all_tables, function(table){
  substr(table, regexpr("TCGA", table)[1] + 13 , regexpr("TCGA", table)[1] + 15)
}))

# get tables from samples for which only blood or tissue is available, and make a symlink to create a filtered file, which is required by nextflow
uniques = which(!(duplicated(barcodes) | duplicated(barcodes, fromLast = T)))
for (u in uniques) { system(paste("ln -s", all_tables[u], gsub(".tsv", "_blood_tissue_filtered.tsv", all_tables[u]) )) }

# for each of duplicated barcode, run in parallel the script to filter on blood/tissue
cl<-makeCluster(cpu)
registerDoParallel(cl)
foreach(d = which(duplicated(barcodes))) %dopar% blood_tissue_filter(barcodes[d], all_tables)
date()
stopCluster(cl)
