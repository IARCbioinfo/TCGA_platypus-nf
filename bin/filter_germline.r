#! /usr/bin/env Rscript

source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")

# get the parameters
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$min_af)) args$min_af = 0.1
if(is.null(args$min_DP)) args$min_DP = 10

min_af = as.numeric(args$min_af)
min_DP = as.numeric(args$min_DP)
VCF = args$vcf

# read the input VCF
vcf = read.table(pipe(paste("grep -v '^##' ",VCF," | sed s/^#//" )),stringsAsFactors=F,header=T,sep="\t",quote=NULL,check.names = F,fill=TRUE)
header_size = as.numeric(system(paste(paste("sed -e '/^#CHROM/,$d' ",VCF,sep="")," | wc -l",sep=""), intern = T))
con = file(VCF,"r")
header = readLines(con, n=header_size)
close(con)

# compute filtered VCF
geno_col = which(colnames(vcf) == "FORMAT") + 1
DP = get_genotype(vcf[,geno_col],vcf$FORMAT[1],"NR")
vaf = get_genotype(vcf[,geno_col],vcf$FORMAT[1],"NV") / DP
new_vcf = vcf[which(vaf >= min_af & vcf$FILTER=="PASS" & DP >= min_DP),]
colnames(new_vcf)[1]="#CHROM"

# output the new VCF
cat(header, file=paste(gsub(".vcf","",VCF),"_filter.vcf",sep=""), sep = "\n")
write.table(new_vcf, paste(gsub(".vcf","",VCF),"_filter.vcf",sep=""), sep="\t", row.names = FALSE, append = T, quote = F)
