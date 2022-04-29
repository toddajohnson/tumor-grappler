#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)
	library(VariantAnnotation)})

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

sample.ids.ls <- unique(sample.info.dt$sample.id)
names(sample.ids.ls) <- sample.ids.ls

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")


vbe.summary.cols <- c(
	'taxid_genus', 'name_genus', 'reads_genus_tree', 'taxid_species', 'name_species',
	'reads_species_tree', 'taxid_assigned', 'name_assigned', 'reads_assigned_tree', 'reads_assigned_direct',
	'reference', 'reference_taxid', 'reference_kmer_count', 'alternate_kmer_count', 'rname',
	'startpos', 'endpos', 'numreads', 'covbases', 'coverage',
	'meandepth', 'meanbaseq', 'meanmapq', 'integrations', 'QCStatus')

vbe.summary.colClasses <- c(
	'integer', 'character', 'integer', 'integer', 'character',
	'integer', 'integer', 'character', 'integer', 'integer',
	'character', 'integer', 'integer', 'integer', 'character',
	'integer', 'integer', 'integer', 'integer', 'double',
	'double', 'double', 'double', 'integer', 'character')

vbe.base.dir <- file.path(run.dir, 'result/VIRUSBreakend')

vbe.summaries.ls <- lapply(
	X = sample.ids.ls,
	FUN = function( curr.sample.id ){
		print(curr.sample.id)
		curr.subject.id <- unique(sample.info.dt[sample.id==curr.sample.id,]$subject.id);
		curr.subject.index <- sample.info.dt[sample.id==curr.sample.id,]$subject.index
		curr.sample.index <- sample.info.dt[sample.id==curr.sample.id,]$sample.index
		
		curr.vbe.sample.dir <- file.path(vbe.base.dir, curr.subject.id, curr.sample.id)
		
		curr.file <- file.path( curr.vbe.sample.dir, paste(curr.sample.id, '.vbe.vcf.gz.summary.tsv', sep='') )
		
		if ( file.exists(curr.file) ){
			if( file.size(curr.file)>0){
				curr.dt <- fread(
					curr.file,
					skip = 1,
					col.names = vbe.summary.cols,
					colClasses = vbe.summary.colClasses)
				
				curr.dt[,subject.id:=curr.subject.id]
				curr.dt[,sample.id:=curr.sample.id]
				curr.dt[,subject.index:=curr.subject.index]
				curr.dt[,sample.index:=curr.sample.index]
				curr.dt[,c('subject.id', 'sample.id', 'subject.index', 'sample.index', vbe.summary.cols), with=FALSE]
				curr.dt
			}else{
				NULL
			}
		}else{
			NULL
		}
	})

vbe.summaries.dt <- rbindlist(vbe.summaries.ls)

parse_ALT <- function( curr.ALT ){
	paste(as.character(curr.ALT[[1]]), collapse=",")
}


sample.ids.with.viral.insertes.ls <- unique(vbe.summaries.dt$sample.id)

vbe.sv.variant.info.ls <- mclapply(
	X = sample.ids.with.viral.insertes.ls,
	FUN = function( curr.sample.id ){
		print( curr.sample.id )
		curr.subject.id <- unique(sample.info.dt[sample.id==curr.sample.id,]$subject.id)
		curr.subject.index <- sample.info.dt[sample.id==curr.sample.id,]$subject.index
		curr.sample.index <- sample.info.dt[sample.id==curr.sample.id,]$sample.index
		curr.vbe.sample.dir <- file.path(vbe.base.dir, curr.subject.id, curr.sample.id)

		sv.destination.file <- tempfile()
		sv.vcf.file <- file.path( curr.vbe.sample.dir, paste(curr.sample.id, '.vbe.vcf', sep='') )
#		sv.tabix.file <- TabixFile(sv.vcf.file, yieldSize=100000)		
		sv.vcf <- readVcf(sv.vcf.file, "GRCh38")
			
		sv.vcf.header.dt <- as.data.table(as.data.frame(rowRanges(sv.vcf)), keep.rownames = "ID")
		sv.vcf.header.dt[,CHROM:=as.character(seqnames)]
		sv.vcf.header.dt[,chr:=chrom.map.ls[CHROM]]
		sv.vcf.header.dt[,chr:=ifelse(is.na(chr), 99L, chr)];
		sv.vcf.header.dt[,ALT:=parse_ALT(curr.ALT=ALT), by=1:nrow(sv.vcf.header.dt)]
		sv.vcf.header.dt$ALT <- unlist(sv.vcf.header.dt$ALT)
		sv.vcf.header.dt[,subject.id:=curr.subject.id]
		sv.vcf.header.dt[,sample.id:=curr.sample.id]
		sv.vcf.header.dt[,subject.index:=curr.subject.index]
		sv.vcf.header.dt[,sample.index:=curr.sample.index]
			
		sv.vcf.dt <- as.data.table(info(sv.vcf))
			
		cbind(sv.vcf.header.dt[,list(subject.id, sample.id, subject.index, sample.index, chr, CHROM, POS=start, REF, ALT, ID, QUAL, FILTER)], sv.vcf.dt)
		},
		mc.cores = 4)

vbe.sv.variant.info.dt <- rbindlist(vbe.sv.variant.info.ls)

dir.create(file.path(run.dir, "/result_summaries/VIRUSBreakend"), showWarnings = TRUE, recursive = FALSE, mode = "0777")

save( list=c("vbe.summaries.dt", "vbe.sv.variant.info.dt"),
	file = paste(run.dir, "/result_summaries/VIRUSBreakend/merged_summaries.Rdata", sep=""))

fwrite( vbe.summaries.dt,
	file = paste(run.dir, "/result_summaries/VIRUSBreakend/merged_viral_insert_summaries.csv", sep=""), sep=",")

vbe.sv.variant.info.dt[,list(subject.id, sample.id, CHROM, POS, REF, ALT, ID, MATEID, BEALN, BEID, INSRMRC)]

fwrite( vbe.sv.variant.info.dt,
		file = paste(run.dir, "/result_summaries/VIRUSBreakend/merged_host_viral_insert_SV_info.csv", sep=""), sep=",")

