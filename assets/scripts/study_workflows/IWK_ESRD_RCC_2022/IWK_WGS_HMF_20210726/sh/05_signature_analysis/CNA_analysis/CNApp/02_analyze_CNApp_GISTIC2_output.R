#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh


#AD.threshold <- 0.1
#GISTIC2.date <- "20220202_median"

AD.threshold <- 0.3
GISTIC2.date <- "20220205_none_noChrX"

options(width=220)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

gpl_prefix <- "GRIDSS-2.12.0"

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})

library(stringr)
#library(R.matlab)

HGNC.dir <- "~/reference/HGNC"
load(file = file.path(HGNC.dir, 'hgnc_complete_set_20220131.Rdata'))

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
#HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "maftools", "figures")

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/linx"
#load(file = paste(ref.dir, "/ensembl.gene.info.Rdata", sep=""))

##hg38.genes <- readMat( file.path('~/tools/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat') )

print("Loading cytoband info")
cytoband.dir <- '~/tools/snpEff/data/GRCh38.99'
cytoband.dt <- fread(file.path(cytoband.dir, 'cytoBand.txt.gz'), col.names=c('chrom', 'start', 'end', 'cytoband', 'stain'))

chrom.map.ls <- c(1:22, 'X', 'Y', 'M')
names(chrom.map.ls) <- paste('chr', c(1:22, 'X', 'Y', 'M'), sep='')

cytoband.dt[,chr:=chrom.map.ls[chrom]]
cytoband.dt <- cytoband.dt[!is.na(chr),]
cytoband.dt[,chromband:=paste(chr, cytoband, sep='')]
cytoband.dt[,chromarm:=chrom]
cytoband.dt[grep('p', chromband),chromarm:=paste(chromarm, 'p', sep='')]
cytoband.dt[grep('q', chromband),chromarm:=paste(chromarm, 'q', sep='')]


#load( file = file.path( run.dir, "result", "maftools", "clinical_data.Rdata" ))
load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

clinical.data <- merge(
		x = clinical.data,
		y = clinical.data[,list(subtype.ct=.N), by=list(Histology)],
		by = 'Histology')

chrom.map.ls <- c(as.character(1:22), 'X', 'Y')
names(chrom.map.ls) <- as.character(1:24)

annotation.ls <- unique(as.character(clinical.data$Histology))
annotation.ls <- c('AllSamples', 'ccRCC', 'non-ccRCC', 'ACD-RCC', 'pRCC', 'ccpRCC', 'chRCC', 'other-RCC')

#save(list=c('cna_freqs.dt'),
#load( file = file.path( run.dir, "result/CNApp/merged_cna_freqs.Rdata" ))
load(file = file.path( run.dir, "result/CNApp/merged_cna_freqs.Rdata" ))
load(file = file.path(run.dir, 'result', 'GISTIC2', paste('CN.', GISTIC2.date, '.AllSamples.non-ccRCC.with.subtype.Rdata', sep="")))

status.ls <- names(resolved.CN.lesions.merged.ls)
names(status.ls) <- status.ls

#cna_arm_freqs.dt', 'cna_chrom_freqs.dt', 'cna_subcytoband_freqs.dt

resolved.CNApp.freqs.ls <- lapply(
	X = status.ls,
	FUN = function( curr.status ){
#		curr.status <- status.ls[[3]]
		print(paste("Merging ", curr.status, sep=""))
		curr.resolved <- resolved.CN.lesions.merged.ls[[curr.status]]
		curr.resolved[,chr:=as.integer(chr)]
		
		if (curr.status=='Broad chrom.'){
			curr.freqs.dt <- copy(cna_chrom_freqs.dt)
			
			curr.dt <- merge(
				x = curr.freqs.dt[,list(annotation, Type, chrom, chr, chromarm='', chromband='', CNApp.freq=freq/100)],
				y = curr.resolved,
				by = c('Type', 'chrom', 'chr'))
			
			dcast(
				data = curr.dt,
				formula = analysis.source + status + Type + chrom + chr + chromarm + chromband ~ annotation,
				value.var = "CNApp.freq")
		}else if(curr.status=='Broad arm'){
			curr.freqs.dt <- copy(cna_arm_freqs.dt)
			
			curr.dt <- merge(
				x = curr.freqs.dt[,list(annotation, Type, chrom, chr, chromarm, chromband='', CNApp.freq=freq/100)],
				y = curr.resolved,
				by = c('Type', 'chrom', 'chr', 'chromarm'))
		
			dcast(
				data = curr.dt,
				formula = analysis.source + status + Type + chrom + chr + chromarm + chromband ~ annotation,
				value.var = "CNApp.freq")
		}else if(curr.status=='Focal'){
			curr.freqs.dt <- copy(cna_subcytoband_freqs.dt)
			
			curr.dt <- merge(
				x = curr.freqs.dt[,list(annotation, Type, chrom, chr, chromarm, chromband, CNApp.freq=freq/100)],
				y = curr.resolved,
				by = c('Type', 'chrom', 'chr', 'chromarm', 'chromband'))
		
			dcast(
				data = curr.dt,
				formula = analysis.source + status + Type + chrom + chr + chromarm + chromband ~ annotation,
				value.var = "CNApp.freq",
				fill = 0,
				fun.aggregate = 'mean')
		}
	})

resolved.CNApp.freqs <- rbindlist(resolved.CNApp.freqs.ls)[order(Type, -AllSamples),]
resolved.CNApp.freqs[,region:=ifelse(status=='Broad chrom.', chrom, ifelse(status=='Broad arm', chromarm, ifelse(status=='Focal', chromband, 'UNK')))]
fwrite(resolved.CNApp.freqs[,c('Type', 'chr', 'region', 'analysis.source', 'AllSamples', 'ccRCC', 'ACD-RCC', 'pRCC', 'ccpRCC', 'chRCC'), with=FALSE],
	file.path( run.dir, "result", "CNApp", paste("CN.", GISTIC2.date, ".Amp.Del.CNApp.freqs.by.subtype.", AD.threshold, ".tsv", sep="")), sep='\t')

