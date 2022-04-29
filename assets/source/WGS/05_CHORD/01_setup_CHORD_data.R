# module use /usr/local/package/modulefiles/
# module load R/4.0.2

## Bioconductor packages required by mutSigExtractor
#library('BiocManager')
#install('BSgenome') ## Install genome parser
#install('BSgenome.Hsapiens.UCSC.hg38') ## Install the default genome
## randomForest is required by CHORD
#install('randomForest')

## Install CHORD and mutSigExtractor directly from github using devtools
#install("devtools")
#devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
#devtools::install_github('https://github.com/UMCUGenetics/CHORD/')

options(width=350)

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
	library(BSgenome)
	library(BSgenome.Hsapiens.UCSC.hg38)
	library(randomForest)
	library(mutSigExtractor)
	library(CHORD)
	library(VariantAnnotation)
	library(parallel) })

run_prefix <- gpl_prefix

base.output.dir <- paste(run.dir, "/result/CHORD", sep="")
purple.dir <- paste(run.dir, "/result/", run_prefix, sep="")

if ( file.exists(base.output.dir, recursive=FALSE)==FALSE ){
	dir.create(base.output.dir)
}

contexts.dir <- paste(base.output.dir, "/contexts", sep="")

if ( file.exists(contexts.dir, recursive=FALSE)==FALSE ){
	dir.create(contexts.dir)
}

output.dir <- paste(base.output.dir, "/output", sep="")
if ( file.exists(output.dir, recursive=FALSE)==FALSE ){
	dir.create(output.dir)
}


purple.qc.dt <- fread(paste(run.dir, "/result_summaries/", run_prefix, "/purple.qc.summary.txt", sep=""))
purple.qc.dt.passed <- purple.qc.dt[!grep("FAIL_CONTAMINATION", QCStatus, fixed=TRUE),]

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL[sample.id.T%in%purple.qc.dt.passed$tumor.sample.id,]

tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls

contexts.dt.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( tumor.sample.id ){
#		tumor.sample.id <- tumor.sample.ids.ls[[7]]
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==tumor.sample.id,]$subject.id

		out_path <- paste0(contexts.dir,'/', tumor.sample.id,'_contexts.txt')
		
		curr.vcf.snv <- paste(purple.dir, "/", curr.subject.id, "/purple/", tumor.sample.id, "/", tumor.sample.id, ".purple.somatic.vcf.gz", sep="")
		curr.vcf.sv <- paste(purple.dir, "/", curr.subject.id, "/purple/", tumor.sample.id, "/", tumor.sample.id, ".purple.sv.vcf.gz", sep="")
		
		contexts <- extractSigsChord(
			vcf.snv = curr.vcf.snv,
			vcf.indel = curr.vcf.snv,
			vcf.sv = curr.vcf.sv,
			sv.caller = 'gridss',
			output.path = out_path,
			sample.name = tumor.sample.id,
			ref.genome = BSgenome.Hsapiens.UCSC.hg38,
			verbose = T)
		fread(out_path)
	},
	mc.cores = 12)

l_contexts <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( tumor.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==tumor.sample.id,]$subject.id
		out_path <- paste0(contexts.dir,'/', tumor.sample.id,'_contexts.txt')
		read.delim(out_path, check.names=F)
	})


## Merged the contexts into a matrix
merged_contexts <- do.call(rbind, l_contexts)

## Write to output directory
write.table(merged_contexts, paste(output.dir, '/merged_contexts.txt', sep=""), sep='\t', quote=F)

chord_output <- chordPredict(merged_contexts, do.bootstrap=T, verbose=F)
write.table(chord_output, paste(output.dir, '/chord_pred.txt', sep=""), sep='\t', quote=F)

save( chord_output,
	file = paste(output.dir, '/chord_pred.Rdata', sep=""))

