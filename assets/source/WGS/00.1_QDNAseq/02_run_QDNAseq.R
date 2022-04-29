#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

base.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( base.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(ACE)
	library(QDNAseq)
	library(QDNAseq.hg38)
	library(parallel)
	library(data.table)
	library(utils) })

bam.dir <- paste(base.dir, "/result/downsampled_bam", sep="")
analysis.dir <- paste(base.dir, "/result/QDNAseq_ACE", sep="")
fig.dir <- paste(analysis.dir, "/figures", sep="")
qdnaseq.analysis.dir <- paste(analysis.dir, "/qdnaseq_output", sep="")
ace.analysis.dir <- paste(analysis.dir, "/ace_output", sep="")

dir.create(analysis.dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(fig.dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(qdnaseq.analysis.dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(ace.analysis.dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

load(file = paste(base.dir, "/config_files/fastq_file_info.Rdata", sep=""))

curr.bin.size <- 100

subject.ids.ls <- unique(sample.info.dt$subject.id)
names(subject.ids.ls) <- subject.ids.ls
		
qdnaseq.rds.files <- list.files(qdnaseq.analysis.dir, pattern=glob2rx("*.segments.rds"))
qdnaseq.rds.subject.ids.ls <- substring(qdnaseq.rds.files, 1, nchar(subject.ids.ls[[1]]))
qdnaseq.rds.subject.ids.ls <- unlist(lapply(X=qdnaseq.rds.files, FUN=function(curr.file){strsplit(curr.file, split='.', fixed=TRUE)[[1]][[1]]}))

catch.output.ls <- mclapply(
	X = subject.ids.ls[which(!subject.ids.ls%in%qdnaseq.rds.subject.ids.ls)],
	FUN = function( curr.subject.id ){
		curr.tumor.id <- sample.info.dt[subject.id==curr.subject.id & sample.class.short=='T',]$sample.id
		curr.normal.id <- sample.info.dt[subject.id==curr.subject.id & sample.class.short=='R',]$sample.id
		
		print( paste("Setting bam files for ", curr.subject.id, sep=""))
		setwd(qdnaseq.analysis.dir)
		tumor_bam <- paste(bam.dir, "/", curr.tumor.id, ".downsampled.bam", sep="")
		
		if (length(curr.normal.id)==1){
			normal_bam <- paste(bam.dir, "/", curr.normal.id, ".downsampled.bam", sep="")
			bam.files <- c(tumor_bam, normal_bam)
		}else{
			bam.files <- c(tumor_bam)
		}
		
		print( paste("Running QDNAseq for ", curr.subject.id, sep=""))
		bins <- getBinAnnotations(binSize=curr.bin.size, genome="hg38")
		readCounts <- binReadCounts(bins, bamfiles=bam.files, ext="downsampled.bam")
		readCountsFiltered <- applyFilters(readCounts)
		readCountsFiltered <- estimateCorrection(readCountsFiltered)
		copyNumbers <- correctBins(readCountsFiltered)
		copyNumbersNormalized <- normalizeBins(copyNumbers)
		copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
		copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
		copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
		copyNumbersCalled <- callBins(copyNumbersSegmented)
		
		print( paste("Saving QDNAseq output for ", curr.subject.id, sep=""))
		
		saveRDS(copyNumbersSegmented,
			file = file.path(qdnaseq.analysis.dir, paste(curr.subject.id, ".segments.rds", sep="")))
		saveRDS(copyNumbersCalled,
			file = file.path(qdnaseq.analysis.dir, paste(curr.subject.id, ".calls.rds", sep="")))

		print( paste("Plotting QDNAseq output for ", curr.subject.id, sep=""))
	
		pdf(file=paste(qdnaseq.analysis.dir, "/", curr.subject.id, ".segments.pdf", sep=""), onefile=TRUE,
			width=7, height=5)
		plot(copyNumbersSegmented, cex.axis=0.6, cex.lab=1, las=2)
		dev.off()

		pdf(file=paste(qdnaseq.analysis.dir, "/", curr.subject.id, ".calls.pdf", sep=""), onefile=TRUE,
			width=7, height=5)
		plot(copyNumbersCalled, cex.axis=0.6, cex.lab=1, las=2)
		dev.off()
		
		ace.output.dir <- paste(ace.analysis.dir, "/", curr.subject.id, sep="")
		dir.create(ace.output.dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
		setwd(ace.output.dir)

		print( paste("Running ACE for ", curr.subject.id, sep=""))
		
		ploidyplotloop(copyNumbersSegmented, ace.output.dir, ploidies = c(2,3), 
			imagetype = 'pdf', method = 'RMSE', penalty = 0.5, cap = 12, 
			bottom = 0, trncname = FALSE, printsummaries = TRUE, 
			autopick = FALSE)
	}, mc.cores=18)



