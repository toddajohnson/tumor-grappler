#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
#module load R/4.0.2


options(width=350)

working_dir <- getwd()
#gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "/home/tjohnson/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)
	library(XML)
	library(RCurl)
	library(rlist) })

sage.version <- "2.8"

if ( length(list.dirs(paste(run.dir, "/result_summaries/", gpl_prefix, sep=""), recursive=FALSE))==0 ){
	system(paste("mkdir ", run.dir, "/result_summaries/", gpl_prefix, sep=""))
}

load( file = paste(run.dir, "/config_files/sage_run_info.Rdata", sep=""))

trim.leading <- function (x)  sub("^\\s+", "", x)
trim.trailing <- function (x) sub("\\s+$", "", x)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

germline.subject.ids.ls <- sage.germline.run.info.dt$subject.id
names(germline.subject.ids.ls) <- germline.subject.ids.ls

somatic.subject.ids.ls <- sage.somatic.run.info.dt$subject.id
names(somatic.subject.ids.ls) <- somatic.subject.ids.ls

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

saved.objects <- c()

if (import_panel_germline==TRUE){
	germline.variant.rate.details.ls <- mclapply(
		X = germline.subject.ids.ls,
		FUN = function( curr.subject.id ){
	#		curr.subject.id <- germline.subject.ids.ls[[1]]
			print( paste("Extracting variant rate detail tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[2]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant rate detail tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Chromosome, length, Variants, Variants.rate)]
		},
		mc.cores = 6)
	
	germline.variant.rate.details <- rbindlist(germline.variant.rate.details.ls)
	
	germline.variant.type.counts.ls <- mclapply(
		X = germline.subject.ids.ls,
		FUN = function( curr.subject.id ){
	#		curr.subject.id <- germline.subject.ids.ls[[1]]
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[3]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant type count tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Type, Total)]
		},
		mc.cores = 6)
	
	germline.variant.type.counts <- rbindlist(germline.variant.type.counts.ls)
	
	germline.variant.genes.cols.ls <- mclapply(
			X = germline.subject.ids.ls,
			FUN = function( curr.subject.id ){
	#		curr.subject.id <- somatic.subject.ids.ls[[1]]
				print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
				curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
				
				curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
				setnames(curr.dt, old='#GeneName', new='GeneName')
				orig.cols <- copy(colnames(curr.dt))
				curr.dt[,subject.id:=curr.subject.id]
				
				print( paste("Finished extracting variant gene table columns for ", curr.subject.id, sep=""))
				
				colnames(curr.dt[,c('subject.id', orig.cols), with=FALSE])
			},
			mc.cores = 6)
	
	germline.variant.genes.cols <- unique(unlist(germline.variant.genes.cols.ls))
	
	germline.variant.genes.ls <- mclapply(
		X = germline.subject.ids.ls,
		FUN = function( curr.subject.id ){
	#		curr.subject.id <- somatic.subject.ids.ls[[1]]
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
			
			curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
			setnames(curr.dt, old='#GeneName', new='GeneName')
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.id:=curr.subject.id]
	
			missing.cols <- germline.variant.genes.cols[which(!germline.variant.genes.cols%in%colnames(curr.dt))]
			set(curr.dt, i=NULL, j=missing.cols, value=0L)
			
			print( paste("Finished extracting variant gene table for ", curr.subject.id, sep=""))
			
			curr.dt[,c(germline.variant.genes.cols), with=FALSE]
		},
		mc.cores = 6)
	
	
	germline.variant.genes.dt <- rbindlist(germline.variant.genes.ls)
	
	saved.objects <- c(saved.objects, 'germline.variant.rate.details', 'germline.variant.type.counts', 'germline.variant.genes.dt')
	
	print("Germline gene panel variant counts")
	print( germline.variant.type.counts[Chromosome=='Total',list(subject.id, Variants, Variants.rate)] )
	
	print("Germline gene panel variant counts")
	print( germline.variant.rate.details[Chromosome=='Total',list(subject.id, Variants, Variants.rate)] )
}



if (import_genomewide_germline==TRUE){
	germline.genomewide.subject.ids.ls <- sage.germline.run.info.dt$subject.id
	names(germline.genomewide.subject.ids.ls) <- germline.genomewide.subject.ids.ls
	
	germline.genomewide.variant.rate.details.ls <- mclapply(
		X = germline.genomewide.subject.ids.ls,
		FUN = function( curr.subject.id ){
			print( paste("Extracting variant rate detail tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[2]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant rate detail tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Chromosome, length, Variants, Variants.rate)]
		},
		mc.cores = 6)
	
	germline.genomewide.variant.rate.details <- rbindlist(germline.genomewide.variant.rate.details.ls)
	
	germline.genomewide.variant.type.counts.ls <- mclapply(
		X = germline.genomewide.subject.ids.ls,
		FUN = function( curr.subject.id ){
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[3]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant type count tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Type, Total)]
		},
		mc.cores = 6)
	
	germline.genomewide.variant.type.counts <- rbindlist(germline.genomewide.variant.type.counts.ls)
	
	germline.genomewide.variant.genes.cols.ls <- mclapply(
		X = germline.genomewide.subject.ids.ls,
		FUN = function( curr.subject.id ){
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
			
			curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
			setnames(curr.dt, old='#GeneName', new='GeneName')
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant gene table columns for ", curr.subject.id, sep=""))
			
			colnames(curr.dt[,c('subject.id', orig.cols), with=FALSE])
		},
		mc.cores = 6)
	
	germline.genomewide.variant.genes.cols <- unique(unlist(germline.genomewide.variant.genes.cols.ls))
	
	germline.genomewide.variant.genes.ls <- mclapply(
		X = germline.genomewide.subject.ids.ls,
		FUN = function( curr.subject.id ){
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
			
			curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
			setnames(curr.dt, old='#GeneName', new='GeneName')
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.id:=curr.subject.id]
			
			missing.cols <- germline.genomewide.variant.genes.cols[which(!germline.genomewide.variant.genes.cols%in%colnames(curr.dt))]
			set(curr.dt, i=NULL, j=missing.cols, value=0L)
			
			print( paste("Finished extracting variant gene table for ", curr.subject.id, sep=""))
			
			curr.dt[,c(germline.genomewide.variant.genes.cols), with=FALSE]
		},
		mc.cores = 6)

	germline.genomewide.variant.genes.dt <- rbindlist(germline.genomewide.variant.genes.ls)
	
	
	germline.genomewide.AF_MAF_filt.subject.ids.ls <- sage.germline.run.info.dt$subject.id
	names(germline.genomewide.AF_MAF_filt.subject.ids.ls) <- germline.genomewide.AF_MAF_filt.subject.ids.ls
	
	germline.genomewide.AF_MAF_filt.variant.rate.details.ls <- mclapply(
			X = germline.genomewide.AF_MAF_filt.subject.ids.ls,
			FUN = function( curr.subject.id ){
				print( paste("Extracting variant rate detail tables for ", curr.subject.id, sep=""))
				curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".AF_MAF_filtered.html", sep="")
				
				tables <- readHTMLTable( doc=curr.snpEff.file)
				tables <- list.clean(tables, fun = is.null, recursive = FALSE)
				
				curr.dt <- as.data.table(tables[[2]])
				colnames(curr.dt) <- trim(colnames(curr.dt))
				colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
				curr.dt[,subject.id:=curr.subject.id]
				
				print( paste("Finished extracting variant rate detail tables for ", curr.subject.id, sep=""))
				
				curr.dt[,list(subject.id, Chromosome, length, Variants, Variants.rate)]
			},
			mc.cores = 6)
	
	germline.genomewide.AF_MAF_filt.variant.rate.details <- rbindlist(germline.genomewide.AF_MAF_filt.variant.rate.details.ls)
	
	germline.genomewide.AF_MAF_filt.variant.type.counts.ls <- mclapply(
			X = germline.genomewide.AF_MAF_filt.subject.ids.ls,
			FUN = function( curr.subject.id ){
				print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
				curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".AF_MAF_filtered.html", sep="")
				tables <- readHTMLTable( doc=curr.snpEff.file)
				tables <- list.clean(tables, fun = is.null, recursive = FALSE)
				
				curr.dt <- as.data.table(tables[[3]])
				colnames(curr.dt) <- trim(colnames(curr.dt))
				colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
				curr.dt[,subject.id:=curr.subject.id]
				
				print( paste("Finished extracting variant type count tables for ", curr.subject.id, sep=""))
				
				curr.dt[,list(subject.id, Type, Total)]
			},
			mc.cores = 6)
	
	germline.genomewide.AF_MAF_filt.variant.type.counts <- rbindlist(germline.genomewide.AF_MAF_filt.variant.type.counts.ls)
	
	germline.genomewide.AF_MAF_filt.variant.genes.cols.ls <- mclapply(
			X = germline.genomewide.AF_MAF_filt.subject.ids.ls,
			FUN = function( curr.subject.id ){
				print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
				curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".AF_MAF_filtered.genes.txt", sep="")
				
				curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
				setnames(curr.dt, old='#GeneName', new='GeneName')
				orig.cols <- copy(colnames(curr.dt))
				curr.dt[,subject.id:=curr.subject.id]
				
				print( paste("Finished extracting variant gene table columns for ", curr.subject.id, sep=""))
				
				colnames(curr.dt[,c('subject.id', orig.cols), with=FALSE])
			},
			mc.cores = 6)
	
	germline.genomewide.AF_MAF_filt.variant.genes.cols <- unique(unlist(germline.genomewide.AF_MAF_filt.variant.genes.cols.ls))
	
	germline.genomewide.AF_MAF_filt.variant.genes.ls <- mclapply(
			X = germline.genomewide.AF_MAF_filt.subject.ids.ls,
			FUN = function( curr.subject.id ){
				print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
				curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/germline_genomewide/snpEff_summary_", snpEff.db, ".AF_MAF_filtered.genes.txt", sep="")
				
				curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
				setnames(curr.dt, old='#GeneName', new='GeneName')
				orig.cols <- copy(colnames(curr.dt))
				curr.dt[,subject.id:=curr.subject.id]
				
				missing.cols <- germline.genomewide.AF_MAF_filt.variant.genes.cols[which(!germline.genomewide.AF_MAF_filt.variant.genes.cols%in%colnames(curr.dt))]
				set(curr.dt, i=NULL, j=missing.cols, value=0L)
				
				print( paste("Finished extracting variant gene table for ", curr.subject.id, sep=""))
				
				curr.dt[,c(germline.genomewide.AF_MAF_filt.variant.genes.cols), with=FALSE]
			},
			mc.cores = 6)
	
	germline.genomewide.AF_MAF_filt.variant.genes.dt <- rbindlist(germline.genomewide.AF_MAF_filt.variant.genes.ls)
	
	saved.objects <- c(saved.objects,
		'germline.genomewide.variant.rate.details', 'germline.genomewide.variant.type.counts', 'germline.genomewide.variant.genes.dt',
		'germline.genomewide.AF_MAF_filt.variant.rate.details', 'germline.genomewide.AF_MAF_filt.variant.type.counts', 'germline.genomewide.AF_MAF_filt.variant.genes.dt')
		
	print("Germline genomewide variant counts")
	print( germline.genomewide.variant.rate.details[Chromosome=='Total',list(subject.id, Variants, Variants.rate)] )
	
	print("Germline genomewide AF_MAF_filtered variant counts")
	print( germline.genomewide.AF_MAF_filt.variant.rate.details[Chromosome=='Total',list(subject.id, Variants, Variants.rate)] )
}



somatic.variant.rate.details.ls <- mclapply(
		X = somatic.subject.ids.ls,
		FUN = function( curr.subject.id ){
#		curr.subject.id <- somatic.subject.ids.ls[[1]]
			print( paste("Extracting variant rate detail tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/somatic/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[2]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant rate detail tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Chromosome, length, Variants, Variants.rate)]
		},
		mc.cores = 6)

somatic.variant.rate.details <- rbindlist(somatic.variant.rate.details.ls)

somatic.variant.type.counts.ls <- mclapply(
		X = somatic.subject.ids.ls,
		FUN = function( curr.subject.id ){
#		curr.subject.id <- somatic.subject.ids.ls[[1]]
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/somatic/snpEff_summary_", snpEff.db, ".filtered.html", sep="")
			tables <- readHTMLTable( doc=curr.snpEff.file)
			tables <- list.clean(tables, fun = is.null, recursive = FALSE)
			
			curr.dt <- as.data.table(tables[[3]])
			colnames(curr.dt) <- trim(colnames(curr.dt))
			colnames(curr.dt) <- gsub(pattern=" ", replacement=".", x=colnames(curr.dt), fixed=TRUE)
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant type count tables for ", curr.subject.id, sep=""))
			
			curr.dt[,list(subject.id, Type, Total)]
		},
		mc.cores = 6)

somatic.variant.type.counts <- rbindlist(somatic.variant.type.counts.ls)


somatic.variant.genes.cols.ls <- mclapply(
		X = somatic.subject.ids.ls,
		FUN = function( curr.subject.id ){
#		curr.subject.id <- somatic.subject.ids.ls[[1]]
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/somatic/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
			
			curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
			setnames(curr.dt, old='#GeneName', new='GeneName')
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.id:=curr.subject.id]
			
			print( paste("Finished extracting variant gene table columns for ", curr.subject.id, sep=""))
			
			colnames(curr.dt[,c('subject.id', orig.cols), with=FALSE])
		},
		mc.cores = 6)

somatic.variant.genes.cols <- unique(unlist(somatic.variant.genes.cols.ls))


somatic.variant.genes.ls <- mclapply(
		X = somatic.subject.ids.ls,
		FUN = function( curr.subject.id ){
#		curr.subject.id <- somatic.subject.ids.ls[[1]]
			print( paste("Extracting variant type count tables for ", curr.subject.id, sep=""))
			curr.snpEff.file <- paste( run.dir, "/result/mutations/", curr.subject.id, "/SAGE-", sage.version, "/somatic/snpEff_summary_", snpEff.db, ".filtered.genes.txt", sep="")
			
			curr.dt <- fread( curr.snpEff.file, skip=1, header=TRUE)
			setnames(curr.dt, old='#GeneName', new='GeneName')
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.id:=curr.subject.id]
			
			if (exists('germline.genomewide.variant.genes.cols.ls')){
				genes.cols <- germline.genomewide.variant.genes.cols
			}else if(exists('germline.variant.genes.cols')){
				genes.cols <- germline.variant.genes.cols
			}
			missing.cols <- genes.cols[which(!genes.cols%in%colnames(curr.dt))]
			set(curr.dt, i=NULL, j=missing.cols, value=0L)
			
			print( paste("Finished extracting variant gene table for ", curr.subject.id, sep=""))
			
			curr.dt[,c(genes.cols), with=FALSE]
		},
		mc.cores = 6)

somatic.variant.genes.dt <- rbindlist(somatic.variant.genes.ls)

saved.objects <- c(saved.objects, 'somatic.variant.rate.details', 'somatic.variant.type.counts', 'somatic.variant.genes.dt')


save(list=saved.objects,
	file = paste( run.dir, "/result_summaries/sage_variant_snpEff_summaries.Rdata", sep=""))

print("Somatic variant counts")
print( somatic.variant.rate.details[Chromosome=='Total',list(subject.id, Variants, Variants.rate)] )
