#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=350)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )


library(data.table)
library(parallel)
library(ggpubr)

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "SVs", "figures")

load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]

#modified QC categories 3/11/2022
if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

if (remove.deleted.genes.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_DELETED_GENES', QCStatus, invert=TRUE),]
}

if (remove.gender.mismatch.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_GENDER_MISMATCH', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined.QCd[,low.purity:=0L]
purple.qc.dt.combined.QCd[grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),low.purity:=1L]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[which(!(low.purity==1L & tumor.depth<tumor.depth.cutoff)),]

#purple.qc.dt.combined.QCd <- rbind(
#	purple.qc.dt.combined.QCd[QCStatus%in%c('PASS'),],
#	purple.qc.dt.combined.QCd[tumor.depth>=tumor.depth.cutoff,][grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),])

purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[!tumor.sample.id%in%extra.samples.to.exclude,]

QC.pass.samples <- purple.qc.dt.combined.QCd$tumor.sample.id
length(QC.pass.samples)
# 34

#merged.SV.CN.info.dt', 'linx.SV.CN.info.dt', 'SV_CN.func.summary

clinical.data <- merge(
	x = purple.purity.dt[,list(tumor.sample.id, gender, purity, ploidy)],
	y = patient.info.dt[,list(tumor.sample.id,
		histology, histology.short, years.dialysis, age)],
	by = c('tumor.sample.id'),
	all.x = TRUE)

## clinical.data[!is.na(age),list(sample.ct=.N, mean.age=mean(age)), by=list(histology.short)]
#		histology.short sample.ct mean.age
#1:            pRCC         5 47.40000
#2:           ccRCC        18 65.61111
#3:         ACD-RCC         9 56.88889
#4:          ccpRCC         2 50.00000
#5:           chRCC         3 60.00000
#6:          uncRCC         1 63.00000

clinical.data <- clinical.data[tumor.sample.id%in%QC.pass.samples,]
clinical.data[,Histology:=factor(histology.short, levels=c('ccRCC', 'ACD-RCC', 'pRCC', 'ccpRCC', 'chRCC'))]

Histology.ls <- unique(clinical.data$Histology)
Histology.ls <- Histology.ls[order(Histology.ls)]
Histology.colors.ls <- get_palette(palette='jco', k=length(Histology.ls))
names(Histology.colors.ls) <- Histology.ls
##> Histology.colors.ls
##	ccRCC     ACD-RCC        pRCC      ccpRCC       chRCC 
##	"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" 
##	pRCC       ccRCC     ACD-RCC      ccpRCC       chRCC 
##	"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" 

## Note that passing Histology here must be as an integer index into the factor...
#clinical.data[,Histology.color:=Histology.colors.ls[as.character(Histology)]]
clinical.data[,Histology.color:=Histology.colors.ls[Histology]]

save(list=c('Histology.colors.ls', 'clinical.data'),
		file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))
