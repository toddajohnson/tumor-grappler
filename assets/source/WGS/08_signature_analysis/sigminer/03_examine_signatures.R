#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

run.dir <- args[1]

options(width=200)

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

library(sigminer)
library(NMF)

purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")

sig.type.ls <- c("SBS", "DBS", "ID", "CN_S_40", "CN_S_48", "CN_W")
names(sig.type.ls) <- sig.type.ls

###
##sig.type <- sig.type.ls[[6]]
###

sig.class.ls <- c("SBS", "DBS", "ID", "CN_S", "CN_S", "CN_W")
names(sig.class.ls) <- sig.type.ls

sig.mode.ls <- c("SBS", "DBS", "ID", "copynumber", "copynumber", "copynumber")
names(sig.mode.ls) <- sig.type.ls

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL
tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

tumor.normal.pair.info[,somatic.vcf:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.somatic.vcf.gz", sep="") )]
tumor.normal.pair.info[,cnv.file:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.segment.tsv", sep="") )]

fig.dir <- file.path(run.dir, "result/sigminer/figures")

lapply(
	X = sig.type.ls,
	FUN = function(sig.type){
		sig.mode <- sig.mode.ls[[sig.type]]
		sig.class <- sig.class.ls[[sig.type]]
		
		load( file = paste(run.dir, "/result/sigminer/extracted_signatures_", sig.type, ".Rdata", sep=""))
		
		pdf( file=file.path(fig.dir, paste(sig.type, ".signature.profiles.pdf", sep="")),
			width=14, height=8, onefile=TRUE )
		
		print(paste('plotting ', sig.type, sep=''))
		
		if (sig.type%in%c('SBS', 'DBS', 'ID') ){
			load( file = paste(run.dir, "/result/sigminer/mt_tally_", sig.type, ".Rdata", sep=""))
			curr.nmf_matrix <- mt_tally$nmf_matrix
		}else if (sig.type=='CN_W'){
			load( file = paste(run.dir, "/result/sigminer/CN_W_tally.Rdata", sep=""))
			curr.nmf_matrix <- cn_tally$nmf_matrix
		}else if (sig.type=='CN_S_48'){
			load( file = paste(run.dir, "/result/sigminer/CN_S_tally.Rdata", sep=""))
			curr.nmf_matrix <- cn_tally$all_matrices$CN_48
		}else if (sig.type=='CN_S_40'){
			load( file = paste(run.dir, "/result/sigminer/CN_S_tally.Rdata", sep=""))
			curr.nmf_matrix <- cn_tally$all_matrices$CN_40
		}
		
		if (sig.type%in%c('CN_S_40', 'CN_S_48', 'CN_W')){
			load(file = file.path( run.dir, "result", "sigminer", paste("merged_", sig.class, ".Rdata", sep="") ))
			
			show_cn_group_profile( data = merged.cn, genome_build = "hg38")
			show_cn_freq_circos( data = merged.cn, genome_build = "hg38")
		}
		
		p1 <- bp_show_survey(sig_est, add_score=TRUE)
		print(p1)
		
		curr.db_type <- "human-genome"
		
		if (sig.type%in%c("SBS", "DBS")){
			curr.sig_db <- paste('latest_', sig.mode.ls[[sig.type]], '_GRCh38', sep='')
			curr.normalize <- 'row'
		}else if (sig.type=="ID"){
			curr.sig_db <- sig.mode.ls[[sig.type]]
			curr.normalize <- 'row'
		}else if (sig.type=='CN_S_48'){
			curr.sig_db <- "CNS_TCGA"
			curr.normalize <-  'row'
		}else if (sig.type=='CN_S_40'){
			curr.sig_db <- "CNS_USARC"
			curr.normalize <-  'row'
		}else if (sig.type=='CN_W'){
			curr.sig_db <- NULL
			curr.normalize <- 'feature'
		}
		
		if (sig.type%in%c("SBS", "DBS", "ID")){
			print( paste('Creating ', sig.type, ' similarity plot', sep=""))
		
			sim <- get_sig_similarity(
				Signature = obj_suggested,
				sig_db = curr.sig_db,
				db_type = curr.db_type,
				normalize = curr.normalize)
		
			# Visualize the match result
			pheatmap::pheatmap(sim$similarity)
		
			fitted.sigs <- sig_fit(
				catalogue_matrix = t(curr.nmf_matrix),
				sig = obj_suggested,
				sig_index = "ALL",
				sig_db = curr.sig_db,
				db_type = curr.db_type,
				sig.mode = sig.mode.ls[[sig.type]],
		#				return_class = "data.table",
				rel_threshold = 0.05)
		
			print(paste("Plotting profile for ", sig.type, sep=""))
			p3 <- show_sig_profile( obj_suggested,
				mode = sig.mode.ls[[sig.type]],
				paint_axis_text = FALSE,
				x_label_angle = 90)
			print(p3)
		}else if( sig.class=="CN_S" ){
			print( paste('Creating ', sig.type, ' similarity plot', sep=""))
			
			sim <- get_sig_similarity(
				Signature = obj_suggested,
				sig_db = curr.sig_db,
				db_type = curr.db_type,
				normalize = curr.normalize)
		
			# Visualize the match result
			pheatmap::pheatmap(sim$similarity)
		
			fitted.sigs <- sig_fit(
				catalogue_matrix = t(curr.nmf_matrix),
				sig = obj_suggested,
				sig_index = "ALL",
				sig_db = curr.sig_db,
				db_type = curr.db_type,
				sig.mode = sig.mode.ls[[sig.type]],
				rel_threshold = 0.05)
			
			print(paste("Plotting profile for ", sig.type, sep=""))
			p3 <- show_sig_profile( obj_suggested,
				mode = sig.mode.ls[[sig.type]],
				normalize = curr.normalize,
				method = "S",
				style = "cosmic")
			print(p3)
		}else if( sig.type=="CN_W" ){
			print( paste('Creating ', sig.type, ' similarity plot', sep=""))
			curr.normalize <- 'feature'
						
			sim <- NULL
			fitted.sigs <- NULL
			
			print(paste("Plotting profile for ", sig.type, sep=""))
			p3 <- show_sig_profile( obj_suggested,
				mode = sig.mode.ls[[sig.type]],
				normalize = curr.normalize,
				method = "W",
				style = "cosmic")
			print(p3)
		}
		p4 <- show_sig_exposure(obj_suggested, hide_samps=FALSE)
		print(p4)
		
		## if (!is.null(fitted.sigs)){
		##     print(paste("Plotting fitted signature exposure for ", sig.type, sep=""))
		##     p5 <- show_sig_fit(fitted.sigs)
		##     print(p5)
		## }
		
		dev.off()
		
		save( list=c('sig_est', 'obj_suggested', 'sim', 'fitted.sigs'),
			file = paste(run.dir, "/result/sigminer/processed_signatures.", sig.type, ".Rdata", sep=""))
	})