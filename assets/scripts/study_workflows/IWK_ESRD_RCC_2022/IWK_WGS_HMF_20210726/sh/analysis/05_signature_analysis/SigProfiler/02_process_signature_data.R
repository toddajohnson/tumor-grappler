#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=210)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)
			library(maftools)
			library(sigminer)})

library(stringr)

## arguments that could be added to common_config.R
remove.copy.number.noise.samples <- TRUE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "SigProfiler", "figures")

sig.type.ls <- c("SBS", "DBS", "ID", "CN_S_40", "CN_S_48", "CN_W")
names(sig.type.ls) <- sig.type.ls

sig.class.ls <- c("SBS", "DBS", "ID", "CN_S", "CN_S", "CN_W")
names(sig.class.ls) <- sig.type.ls

sig.mode.ls <- c("SBS", "DBS", "ID", "copynumber", "copynumber", "copynumber")
names(sig.mode.ls) <- sig.type.ls

curr.categorical.annotation.col <- 'Histology'
curr.numerical.annotation.col <- 'years.dialysis'

load(file = file.path( run.dir, "result/sample_summaries", "clinical_data_with_colors.Rdata" ))
clinical.data$cat.annotation <- clinical.data[[curr.categorical.annotation.col]]
clinical.data$num.annotation <- clinical.data[[curr.numerical.annotation.col]]

base.cols <- c('tumor.sample.id', 'gender', 'purity', 'ploidy', 'cat.annotation', 'num.annotation')

sig.denovo.dt.long.filt.ls <- lapply(
	X = sig.type.ls[1:3],
	FUN = function(curr.type){
		curr.class <- sig.class.ls[[curr.type]]
		load(file = file.path(run.dir,
						paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/", curr.class, ".refit.activities.Rdata", sep="")))
		curr.sig.dt.cols <- copy(colnames(sig.denovo.dt))
		
		sig.classes <- curr.sig.dt.cols[which(!curr.sig.dt.cols%in%base.cols)]
		
		#paste('SBS', c('1', '2', '3', '5', '7a', '8', '13', '16', '17b', '18', '33', '93'), sep="")
		
		dt.long <- melt(
				data = sig.denovo.dt,
				id.vars = c('tumor.sample.id', 'gender', 'purity', 'ploidy', 'cat.annotation', 'num.annotation'),
				measure.vars = sig.classes,
				variable.name = 'sig.class',
				value.name = "sig.mut.ct",
				na.rm = FALSE,
				variable.factor = FALSE,
				value.factor = FALSE)
		
		dt.long[,sig.mut.total:=sum(sig.mut.ct), by=list(tumor.sample.id)]
		dt.long[,sig.mut.propn:=sig.mut.ct/sig.mut.total]
		
		dt.long[,max.sig.mut.propn:=max(sig.mut.propn), by=list(sig.class)]
		dt.long[,sig.mut.non.zero:=ifelse(sig.mut.ct>0, 1L, 0L)]
		dt.long[,sig.mut.non.zero.ct:=sum(sig.mut.non.zero), by=list(sig.class)]
		dt.long[,sig.mut.ct.median:=median(sig.mut.ct), by=list(sig.class)]

		sig.class.in.group.summary <- dt.long[,
			list(sig.mut.ct.in.group.median=median(sig.mut.ct)),
			by=list(cat.annotation, sig.class)]
		sig.class.in.group.summary[,sig.mut.ct.in.group.median.non.zero:=ifelse(sig.mut.ct.in.group.median>0, 1L, 0L)]

		dt.long <- merge(
			x = dt.long,
			y = sig.class.in.group.summary[,
				list(sig.mut.ct.in.group.median.non.zero.ct=sum(sig.mut.ct.in.group.median.non.zero)),
				by=list(sig.class)],
			by = c('sig.class'),
			all.x = TRUE)

#		dt.long[,sig.mut.ct:=ifelse(sig.mut.ct==0, 0.01, sig.mut.ct)]
#		dt.long[,sig.mut.propn:=ifelse(sig.mut.propn==0, 0.01, sig.mut.propn)]
		
#		dt.long[,log.sig.mut.ct:=log(sig.mut.ct)]
#		dt.long[,log.sig.mut.propn:=log(sig.mut.propn)]
		
		dt.long[,median.filter:=ifelse(sig.mut.ct.median>0, 1L, 0L),]
		dt.long[,group.median.filter:=ifelse(sig.mut.ct.in.group.median.non.zero.ct>=2, 1L, 0L),]

		dt.long.sig.class.filt <- dt.long[median.filter==1,]
		included.sig.classes <- sig.classes[which(sig.classes%in%unique(dt.long.sig.class.filt$sig.class))]
		dt.long.sig.class.filt[,sig.class:=factor(sig.class, levels=included.sig.classes, order=TRUE)]
		
		dt.long.group.sig.class.filt <- dt.long[group.median.filter==1,]
		included.group.sig.classes <- sig.classes[which(sig.classes%in%unique(dt.long.group.sig.class.filt$sig.class))]
		dt.long.group.sig.class.filt[,sig.class:=factor(sig.class, levels=included.group.sig.classes, order=TRUE)]
		
		
		list(
				dt.long = dt.long,
				dt.long.sig.class.filt = dt.long.sig.class.filt,
				dt.long.group.sig.class.filt = dt.long.group.sig.class.filt)
	})

sig.dt.long.filt.ls <- lapply(
	X = sig.type.ls[1:3],
	FUN = function(curr.type){
		curr.class <- sig.class.ls[[curr.type]]
		load(file = file.path(run.dir,
						paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/", curr.class, ".refit.activities.Rdata", sep="")))
		curr.sig.dt.cols <- copy(colnames(sig.dt))
		
		sig.classes <- curr.sig.dt.cols[which(!curr.sig.dt.cols%in%base.cols)]

		dt.long <- melt(
				data = sig.dt,
				id.vars = c('tumor.sample.id', 'gender', 'purity', 'ploidy', 'cat.annotation', 'num.annotation'),
				measure.vars = sig.classes,
				variable.name = 'sig.class',
				value.name = "sig.mut.ct",
				na.rm = FALSE,
				variable.factor = FALSE,
				value.factor = FALSE)
		
		dt.long[,sig.mut.total:=sum(sig.mut.ct), by=list(tumor.sample.id)]
		dt.long[,sig.mut.propn:=sig.mut.ct/sig.mut.total]
		dt.long[,max.sig.mut.propn:=max(sig.mut.propn), by=list(sig.class)]
		dt.long[,sig.mut.non.zero:=ifelse(sig.mut.ct>0, 1L, 0L)]
		dt.long[,sig.mut.non.zero.ct:=sum(sig.mut.non.zero), by=list(sig.class)]
		dt.long[,sig.mut.ct.median:=median(sig.mut.ct), by=list(sig.class)]
		
		sig.class.in.group.summary <- dt.long[,
				list(sig.mut.ct.in.group.median=median(sig.mut.ct)),
				by=list(cat.annotation, sig.class)]
		sig.class.in.group.summary[,sig.mut.ct.in.group.median.non.zero:=ifelse(sig.mut.ct.in.group.median>0, 1L, 0L)]
		
		dt.long <- merge(
				x = dt.long,
				y = sig.class.in.group.summary[,
						list(sig.mut.ct.in.group.median.non.zero.ct=sum(sig.mut.ct.in.group.median.non.zero)),
						by=list(sig.class)],
				by = c('sig.class'),
				all.x = TRUE)
		
#		dt.long[,sig.mut.ct:=ifelse(sig.mut.ct==0, 0.01, sig.mut.ct)]
#		dt.long[,sig.mut.propn:=ifelse(sig.mut.propn==0, 0.01, sig.mut.propn)]
#		
#		dt.long[,log.sig.mut.ct:=log(sig.mut.ct)]
#		dt.long[,log.sig.mut.propn:=log(sig.mut.propn)]
		
		dt.long[,median.filter:=ifelse(sig.mut.ct.median>0, 1L, 0L),]
		dt.long[,group.median.filter:=ifelse(sig.mut.ct.in.group.median.non.zero.ct>=2, 1L, 0L),]
		
		dt.long.sig.class.filt <- dt.long[median.filter==1,]
		included.sig.classes <- sig.classes[which(sig.classes%in%unique(dt.long.sig.class.filt$sig.class))]
		dt.long.sig.class.filt[,sig.class:=factor(sig.class, levels=included.sig.classes, order=TRUE)]
		
		dt.long.group.sig.class.filt <- dt.long[group.median.filter==1,]
		included.group.sig.classes <- sig.classes[which(sig.classes%in%unique(dt.long.group.sig.class.filt$sig.class))]
		dt.long.group.sig.class.filt[,sig.class:=factor(sig.class, levels=included.group.sig.classes, order=TRUE)]
		
		
		list(
				dt.long = dt.long,
				dt.long.sig.class.filt = dt.long.sig.class.filt,
				dt.long.group.sig.class.filt = dt.long.group.sig.class.filt)	})
	
#save(list=c('sig.denovo.dt.long.filt.ls', 'sig.dt.long.filt.ls', 'sig.denovo.glm.ls', 'sig.COSMIC.glm.ls'),
save(list=c('sig.denovo.dt.long.filt.ls', 'sig.dt.long.filt.ls'),
	file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/merged.refit.activities.Rdata", sep="")))
