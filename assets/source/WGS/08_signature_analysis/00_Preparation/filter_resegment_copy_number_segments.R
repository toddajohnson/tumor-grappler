#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

min.segment.length <- as.integer(args[1])

options(width=350)
options (repos = "http://cran.ism.ac.jp")

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

needed.packages <- c("data.table", "parallel")
install.packages(setdiff(needed.packages, rownames(installed.packages()))) 

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)})

source( file.path(run.dir, "config_files/common_config.R") )

## chrX.PAR1.start <- 10001
## chrX.PAR1.end <- 2781479
## 
## chrX.PAR2.start <- 55701383
## chrX.PAR2.end <- 156030895

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

names(run.dirs.ls) <- basename(run.dirs.ls)
study.run.names <- names(run.dirs.ls)

study.stored.object.names <- c('purple.qc.dt', 'purple.purity.dt', 'purple.cnv.somatic.info')
names(study.stored.object.names) <- study.stored.object.names

study.data.ls <- lapply(
	X = run.dirs.ls,
	FUN = function( curr.run.dir ){
		source( file.path(curr.run.dir, "config_files/common_config.R"), local=TRUE )
		if(exists("version.date.for.driver.catalog.file")==TRUE){
			load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.", version.date.for.driver.catalog.file, ".Rdata", sep=""))
		}else{
			load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.Rdata", sep=""))
		}

		mget(x=study.stored.object.names)
	})
#load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))

merged.study.data.ls <- lapply(
	X = study.stored.object.names,
	FUN = function( curr.object.name ){
		rbindlist(lapply(
			X = study.run.names,
			FUN = function( curr.study.run.name ){
				curr.dt <- study.data.ls[[curr.study.run.name]][[curr.object.name]]
				curr.cols <- copy(colnames(curr.dt))
				curr.dt[,study.run.name:= curr.study.run.name]
				curr.dt[,c('study.run.name', curr.cols), with=FALSE]
		}))
})

purple.qc.dt <- merged.study.data.ls[['purple.qc.dt']]
purple.purity.dt <- merged.study.data.ls[['purple.purity.dt']]
purple.cnv.somatic.info <- merged.study.data.ls[['purple.cnv.somatic.info']]

## load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
## load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))
## load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.Rdata", sep=""))

purple.qc.dt.combined <- merge(
	x = purple.qc.dt[,list(subject.id, tumor.sample.id,
		Sex=ifelse(AmberGender=='FEMALE', 'F', ifelse(AmberGender=='MALE', 'M', 'UNK')),
		QCStatus, AmberMeanDepth=as.numeric(AmberMeanDepth))],
	y = purple.purity.dt[,list(subject.id, tumor.sample.id, purity, ploidy)],
	by = c('subject.id', 'tumor.sample.id'))

full.sample.ct <- nrow(purple.qc.dt.combined)

purple.qc.dt.combined[,tumor.depth:=purity*AmberMeanDepth]
#purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
#purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]
#
#if (remove.copy.number.noise.samples==TRUE){
#	purple.qc.dt.combined <- purple.qc.dt.combined[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
#}
#
#purple.qc.dt.combined <- rbind(
##	purple.qc.dt.combined[QCStatus%in%c('PASS', 'WARN_HIGH_COPY_NUMBER_NOISE'),],
#	purple.qc.dt.combined[QCStatus%in%c('PASS', 'WARN_GENDER_MISMATCH', 'WARN_HIGH_COPY_NUMBER_NOISE'),],
#	purple.qc.dt.combined[tumor.depth>=tumor.depth.cutoff,][grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),])
#
#purple.qc.dt.combined <- purple.qc.dt.combined[!tumor.sample.id%in%extra.samples.to.exclude,]
#QCd.sample.ct <- nrow(purple.qc.dt.combined)
#

if( exists('exclude.from.CN.analysis')==FALSE ){
	exclude.from.CN.analysis <- c()
}
samples.to.exclude <- unique(c(exclude.from.CN.analysis, extra.samples.to.exclude))

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
purple.qc.dt.combined <- purple.qc.dt.combined.QCd
purple.qc.dt.combined <- purple.qc.dt.combined[!tumor.sample.id%in%samples.to.exclude,]

QCd.sample.ct <- nrow(purple.qc.dt.combined)

print(paste(QCd.sample.ct, '/', full.sample.ct, ' samples passed QC', sep="" ))

## vep.sample.info.dt <- fread( file.path(run.dir, "result_summaries", gpl_prefix, "vep_annotation_sample_info.tsv"))
## 
## vep.tumor.sample.ids.ls <- vep.sample.info.dt[tumor.sample.id%in%purple.qc.dt.combined$tumor.sample.id,]$tumor.sample.id
## names(vep.tumor.sample.ids.ls) <- vep.tumor.sample.ids.ls

purple.cnv.segments <- purple.cnv.somatic.info[tumor.sample.id%in%purple.qc.dt.combined$tumor.sample.id,
	list(subject.id, tumor.sample.id, chromosome, start, end, seg.length=end-start, copyNumber,
		minorAlleleCopyNumber, majorAlleleCopyNumber, segmentStartSupport, segmentEndSupport, method,
		depthWindowCount, bafCount, baf)]
#		depthWindowCount, bafCount, observedBAF)]

purple.cnv.segments <- merge(
	x = purple.cnv.segments,
	y = purple.qc.dt.combined[,list(subject.id, tumor.sample.id, Sex, ploidy)],
	by = c('subject.id', 'tumor.sample.id'))

purple.cnv.segments[,baf:=ifelse(baf>1,1,baf)]
purple.cnv.segments[,NumMarkers:=ifelse(depthWindowCount>bafCount, depthWindowCount, ifelse(bafCount>0, bafCount, ifelse(seg.length>=3, 3, 1)))]
purple.cnv.segments[,copyNumber:=ifelse(copyNumber<=0, 0.005, copyNumber)]
purple.cnv.segments[,majorAlleleCopyNumber:=ifelse(majorAlleleCopyNumber<=0, 0.005, majorAlleleCopyNumber)]
purple.cnv.segments[,minorAlleleCopyNumber:=ifelse(minorAlleleCopyNumber<=0, 0.005, minorAlleleCopyNumber)]

#purple.cnv.segments[,adjusted_diploid:=ifelse(chromosome=='chrX', ifelse(Sex=='M', ifelse( (start<chrX.PAR1.end & end>chrX.PAR1.start) | (start<chrX.PAR2.end & end>chrX.PAR2.start), 2, 1), 2), ifelse(chromosome=='chrY', 1, 2))]
#purple.cnv.segments[,adjusted_ploidy:=ifelse(chromosome=='chrX', ifelse(Sex=='M', ifelse( (start<chrX.PAR1.end & end>chrX.PAR1.start) | (start<chrX.PAR2.end & end>chrX.PAR2.start), ploidy, ploidy/2), ploidy), ifelse(chromosome=='chrY', ploidy/2, ploidy))]

purple.cnv.segments[,adjusted_diploid:=ifelse(chromosome=='chrX', ifelse(Sex=='M', 1, 2), ifelse(chromosome=='chrY', 1, 2))]
purple.cnv.segments[,adjusted_ploidy:=ifelse(chromosome=='chrX', ifelse(Sex=='M', ploidy/2, ploidy), ifelse(chromosome=='chrY', ploidy/2, ploidy))]

cn.segments <- purple.cnv.segments[(depthWindowCount>0 | bafCount>0) & method!='LONG_ARM',
	list(sample=tumor.sample.id, ploidy, adjusted_diploid, adjusted_ploidy, chromosome, startpos=start, endpos=end,
		BAF=baf, seg.length,
		NumMarkers, copyNumber,
		seg.mean.diploid = ifelse(copyNumber==0.005, -5, log2(copyNumber) - log2(2)),
		seg.mean.ploidy = ifelse(copyNumber==0.005, -5, log2(copyNumber) - log2(ploidy)),
		seg.mean.adj_diploid = ifelse(copyNumber==0.005, -5, log2(copyNumber) - log2(adjusted_diploid)),
		seg.mean.adj_ploidy = ifelse(copyNumber==0.005, -5, log2(copyNumber) - log2(adjusted_ploidy)),
		majorAlleleCopyNumber, minorAlleleCopyNumber)]

cn.segments <- cn.segments[seg.length>min.segment.length,];

cn.segments.split <- split( cn.segments, by="sample", keep.by=TRUE)

n.loops <- 100L

max.dist.segm <- 100000L
percent.dist <- 2

dev.to.use <- 0.16
dev.baf <- 0.1

cn.segments.resegmented.dt.ls <- mclapply(
	X = cn.segments.split,
	FUN = function( curr.segs.dt ){
		new.file<-list()
		
		new.file[[1]] <- as.data.frame(curr.segs.dt)
		zz<-2
		new.file[[zz]] <- data.frame(sample=as.character(NA), ploidy=as.numeric(NA), adjusted_diploid=as.numeric(NA), adjusted_ploidy=as.numeric(NA),
			chromosome=as.character(NA), startpos=as.integer(NA), endpos=as.integer(NA),
			BAF=as.numeric(NA), seg.length=as.integer(NA), NumMarkers=as.numeric(NA),
			copyNumber=as.numeric(NA),
			seg.mean.diploid=as.numeric(NA), seg.mean.ploidy=as.numeric(NA),
			seg.mean.adj_diploid=as.numeric(NA), seg.mean.adj_ploidy=as.numeric(NA),
			majorAlleleCopyNumber=as.numeric(NA), minorAlleleCopyNumber=as.numeric(NA))

		for (zz in 2:n.loops) {
			
			print(zz)
			new.file[[zz]] <- data.frame(sample=as.character(NA), ploidy=as.numeric(NA), adjusted_diploid=as.numeric(NA), adjusted_ploidy=as.numeric(NA),
				chromosome=as.character(NA), startpos=as.integer(NA), endpos=as.integer(NA),
				BAF=as.numeric(NA), seg.length=as.integer(NA), NumMarkers=as.numeric(NA),
				copyNumber=as.numeric(NA),
				seg.mean.diploid=as.numeric(NA), seg.mean.ploidy=as.numeric(NA),
				seg.mean.adj_diploid=as.numeric(NA), seg.mean.adj_ploidy=as.numeric(NA),
				majorAlleleCopyNumber=as.numeric(NA), minorAlleleCopyNumber=as.numeric(NA))
			
			k<-1
			flag<-0
			
			for (i in 1:(nrow(new.file[[zz-1]])-1)) {
				if (flag==1) {flag<-0; next}
				s1<-new.file[[zz-1]][i,c('sample', 'ploidy', 'adjusted_diploid', 'adjusted_ploidy', 'chromosome', 'startpos', 'endpos',
					'BAF', 'seg.length', 'NumMarkers', 'copyNumber',
					'seg.mean.diploid', 'seg.mean.ploidy', 'seg.mean.adj_diploid', 'seg.mean.adj_ploidy',
					'majorAlleleCopyNumber', 'minorAlleleCopyNumber')]
				new.file[[zz]][k,]<-s1
				if (nrow( new.file[[zz-1]])<2) {break}
				s2<-new.file[[zz-1]][i+1,c('sample', 'ploidy', 'adjusted_diploid', 'adjusted_ploidy', 'chromosome', 'startpos', 'endpos',
					'BAF', 'seg.length', 'NumMarkers', 'copyNumber',
					'seg.mean.diploid', 'seg.mean.ploidy', 'seg.mean.adj_diploid', 'seg.mean.adj_ploidy',
					'majorAlleleCopyNumber', 'minorAlleleCopyNumber')]

				if (s1$chromosome==s2$chromosome) {
					if (s2$seg.mean.ploidy<(s1$seg.mean.ploidy+dev.to.use) & s2$seg.mean.ploidy>(s1$seg.mean.ploidy-dev.to.use) & 
						(s2$startpos-s1$endpos)<max.dist.segm  & s2$BAF<(s1$BAF+dev.baf) & s2$BAF>(s1$BAF-dev.baf) & 
						(s2$startpos-s1$endpos)<(s2$seg.length+s1$seg.length)*percent.dist) {
						new.file[[zz]][k,]<-s1
						new.file[[zz]][k,]$endpos <- s2$endpos
						new.file[[zz]][k,]$seg.length <- s2$endpos - s1$startpos + 1L
						new.file[[zz]][k,]$NumMarkers <- s1$NumMarkers + s2$NumMarkers
						new.file[[zz]][k,]$copyNumber <- (s1$copyNumber*s1$NumMarkers + s2$copyNumber*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$majorAlleleCopyNumber <- (s1$majorAlleleCopyNumber*s1$NumMarkers + s2$majorAlleleCopyNumber*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$minorAlleleCopyNumber <- (s1$minorAlleleCopyNumber*s1$NumMarkers + s2$minorAlleleCopyNumber*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$seg.mean.diploid <- (s1$seg.mean.diploid*s1$NumMarkers + s2$seg.mean.diploid*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$seg.mean.ploidy <- (s1$seg.mean.ploidy*s1$NumMarkers + s2$seg.mean.ploidy*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$seg.mean.adj_diploid <- (s1$seg.mean.adj_diploid*s1$NumMarkers + s2$seg.mean.adj_diploid*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						new.file[[zz]][k,]$seg.mean.adj_ploidy <- (s1$seg.mean.adj_ploidy*s1$NumMarkers + s2$seg.mean.adj_ploidy*s2$NumMarkers)/new.file[[zz]][k,]$NumMarkers
						
						k<-k+1
						flag<-1
						if (i==(nrow(new.file[[zz-1]])-2)) {
							new.file[[zz]][k,]<-new.file[[zz-1]][i+2,c('sample', 'ploidy', 'adjusted_diploid', 'adjusted_ploidy', 'chromosome', 'startpos', 'endpos',
								'BAF', 'seg.length', 'NumMarkers', 'copyNumber', 'seg.mean.diploid', 'seg.mean.ploidy', 'seg.mean.adj_diploid', 'seg.mean.adj_ploidy',
								'majorAlleleCopyNumber', 'minorAlleleCopyNumber')]; break
						}
						next
					} else {
						if (i==nrow(new.file[[zz-1]])-1) {
							new.file[[zz]][k+1,]<-s2
							break
						} 
					} 
				} else {
					if (i==nrow(new.file[[zz-1]])-1) {
						new.file[[zz]][k+1,]<-s2
						break
					}
				}
				k<-k+1
				
			}
			
			if (zz>2) {
				if (nrow(new.file[[zz]])==nrow(new.file[[zz-1]])) {break}
			}
			
		}
		new.file[[length(new.file)]]
	}, mc.cores = 6)

cn.segments.resegmented.dt <- rbindlist(cn.segments.resegmented.dt.ls)

cn.segments.resegmented.dt[,status:='']
cn.segments.resegmented.dt[(round(majorAlleleCopyNumber, digits=0) + round(minorAlleleCopyNumber, digits=0)) == round(copyNumber, digits=0),status:='Maj+Min == CN']
cn.segments.resegmented.dt[(round(majorAlleleCopyNumber, digits=0) + round(minorAlleleCopyNumber, digits=0)) > round(copyNumber, digits=0),status:='Maj+Min > CN']
cn.segments.resegmented.dt[(round(majorAlleleCopyNumber, digits=0) + round(minorAlleleCopyNumber, digits=0)) < round(copyNumber, digits=0),status:='Maj+Min < CN']

print(table(cn.segments.resegmented.dt$status))

cn.segments.resegmented.dt[,nTotal:=round(ifelse(copyNumber<0, 0, copyNumber), 0)]
cn.segments.resegmented.dt[,nMajor:=ifelse(status%in%c('Maj+Min == CN', 'Maj+Min > CN'), round(majorAlleleCopyNumber), ceiling(majorAlleleCopyNumber))]
cn.segments.resegmented.dt[,nMinor:=nTotal-nMajor]

output.dir <- file.path(run.dir, 'result/CN_segments', sep="")

if ( length(list.dirs(output.dir, recursive=FALSE))==0 ){
	system(paste("mkdir -p ", output.dir, sep=""))
}

save( list=c('cn.segments', 'cn.segments.resegmented.dt'),
	file = file.path( output.dir, paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.ver_20220131.Rdata', sep="") ) )

 