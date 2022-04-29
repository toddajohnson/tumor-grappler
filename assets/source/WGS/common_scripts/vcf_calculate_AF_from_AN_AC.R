#!/usr/bin/env Rscript

options(width=200)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

data.file <- args[[1]]
output.file <- args[[2]]

#data.file <- "/home/tjohnson/reference/alfa/20210106/chromosome_AN_AC_text/chr22_alfa_AN_AC.txt.gz"
#output.file <- "/home/tjohnson/reference/alfa/20210106/chromosome_AF_text/chr22_alfa_AF.txt"

message(paste("Calculating AF from AC/AN fields for ", data.file, sep=""))

ac.an.dt <- fread(data.file, sep="\t")

message("Modifying column names")
old.cols <- copy(colnames(ac.an.dt))
new.cols <- gsub("# ", replacement="", x=old.cols)
new.cols <- unlist(lapply(X=new.cols, FUN=function(curr.str){
	strsplit(curr.str, split="]", fixed=TRUE)[[1]][[2]]}))
new.cols <- gsub(":", replacement=".", x=new.cols)

setnames(ac.an.dt, old=old.cols, new=new.cols)

#EUR.AC EUR.AN AFO.AC AFO.AN EAS.AC EAS.AN AFA.AC AFA.AN LAC.AC LAC.AN LEN.AC LEN.AN OAS.AC OAS.AN SAS.AC SAS.AN OTR.AC OTR.AN AFR.AC AFR.AN ASN.AC ASN.AN TOT.AC TOT.AN
#pop.ls <- c('EUR', 'AFO', 'EAS', 'AFA', 'LAC', 'LEN', 'OAS', 'SAS', 'OTR', 'AFR', 'ASN', 'TOT')
#pop.indexes.to.use <- c(length(pop.ls), 3L)

ac.an.dt[,AF_EAS:=EAS.AC/EAS.AN]
ac.an.dt[,AF_TOT:=TOT.AC/TOT.AN]
ac.an.dt[,MAF_EAS:=ifelse(AF_EAS>0.5, 1-AF_EAS, AF_EAS)]
ac.an.dt[,MAF_TOT:=ifelse(AF_TOT>0.5, 1-AF_TOT, AF_TOT)]

col.names.use <- c(TRUE, rep(FALSE, 24))
names(col.names.use) <- paste('chr', as.character(c(1:22, 'X', 'Y', 'M')), sep="")

message(paste("Writing table to ", output.file, sep=""))
fwrite(ac.an.dt[,list(CHROM, POS, ID, REF, ALT, AF_EAS, AF_TOT, MAF_EAS, MAF_TOT)],
	file=output.file, sep="\t", scipen=6, col.names=col.names.use[[ac.an.dt[1,]$CHROM]])
