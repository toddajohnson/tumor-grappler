#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=350)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

HOME.dir <- '/home/tjohnson'

gpl_prefix <- "GRIDSS-2.12.0"

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})

extra.samples.to.exclude <- c()
tumor.depth.cutoff <- 10;

source( file.path(run.dir, "config_files/common_config.R") )

#install.packages('multidplyr')

#library(purple)
library(tidyr)
library(dplyr)
library(multidplyr)
library(doParallel)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)
library(stringr)

chroms.ls <- paste('chr', c(1:22,'X', 'Y'), sep="")

#ah <- AnnotationHub()
#query(ah, "OrgDb")


# original
#isMaleSexChromosome <- function(gender, chromosome) {
#	return (gender != "FEMALE" & chromosome %in% c('X','Y'))
#}
#
#isLoh <- function(minAllelePloidy, chromosome, gender) {
#	result <- ifelse(isMaleSexChromosome(gender, chromosome), FALSE, minAllelePloidy < 0.5)
#	return (result)
#}
#
#isDel <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
#	result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber < cutoff * ploidy / 2, copyNumber < cutoff * ploidy)
#	return (result)
#}
#
#isAmp <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
#	result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber > cutoff * ploidy / 2,  copyNumber > cutoff * ploidy)
#	return (result)
#}

# change to output absolute CN value based analysis
isMaleSexChromosome <- function(gender, chromosome) {
	return (gender != "FEMALE" & chromosome %in% c('chrX','chrY', 'X', 'Y'))
}

isLoh <- function(minAllelePloidy, chromosome, gender) {
	result <- ifelse(isMaleSexChromosome(gender, chromosome), FALSE, minAllelePloidy < 0.5)
	return (result)
}

isAbsDel <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
	result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber < cutoff * 1, copyNumber < cutoff * 1)
	return (result)
}


isDel <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
	result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber < cutoff * ploidy / 2, copyNumber < cutoff * ploidy)
	return (result)
}

isAmp <- function(copyNumber, chromosome, gender, ploidy, cutoff) {
	result <- ifelse(isMaleSexChromosome(gender, chromosome), copyNumber > cutoff * ploidy / 2,  copyNumber > cutoff * ploidy)
	return (result)
}


load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))
load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.info.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

min.segment.length <- 0
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )

clinical.data <- clinical.data[!tumor.sample.id%in%c(exclude.from.CN.analysis, extra.samples.to.exclude),]

cn.segments.resegmented.dt <- merge(
	x = clinical.data[,list(sample=tumor.sample.id, gender, cancerType=Histology)],
	y = cn.segments.resegmented.dt,
	by = c('sample'))

cn.resegmented.dt <- cn.segments.resegmented.dt[,list(sampleId=sample, ploidy, cancerType, gender, chromosome, start=startpos, end=endpos, actualBaf=BAF, copyNumber)]

cn.resegmented.dt <- cn.resegmented.dt[chromosome!='chrY',]

purple.segment.info <- merge(
		x = clinical.data[,list(tumor.sample.id, purity, ploidy, gender, Histology)],
		y = purple.segment.info,
		by = c('tumor.sample.id'))

chrY.segments <- purple.segment.info[!(germlineStatus=='HET_DELETION' & depthWindowCount<100) & chromosome=='chrY' & !germlineStatus%in%c('NOISE', 'UNKNOWN', 'AMPLIFICATION', 'HOM_DELETION'),
	list(sampleId=tumor.sample.id, ploidy, cancerType=Histology, gender, chromosome, start, end, actualBaf=tumorBAF, copyNumber=refNormalisedCopyNumber)]
#	list(sampleId=tumor.sample.id, ploidy, cancerType=Histology, gender, chromosome, start, end, actualBaf=tumorBAF, copyNumber=tumorCopyNumber)]


cn.resegmented.dt <- rbind(cn.resegmented.dt, chrY.segments)
cn.resegmented.dt <- cn.resegmented.dt[order(sampleId, chromosome, start),]

cn.resegmented <- as.data.frame(cn.resegmented.dt)


processedCopyNumbers = cn.resegmented %>%
		mutate(minAllelePloidy = pmax(0, (1-actualBaf) * copyNumber),
				loh = isLoh(minAllelePloidy, chromosome, gender),
				absDel = isAbsDel(copyNumber, chromosome, gender, 1, 0.5),
				relDel = loh & isDel(copyNumber, chromosome, gender, ploidy=2, 0.6),
				amp1_4 = isAmp(copyNumber, chromosome, gender, ploidy=2, 1.25),
				amp2_0 = isAmp(copyNumber, chromosome, gender, ploidy=2, 1.75),
				amp3_0 =isAmp(copyNumber, chromosome, gender, ploidy=2, 2.25))

processedCopyNumbers$region <- GRanges(processedCopyNumbers$chromosome, ranges = IRanges(start = processedCopyNumbers$start, end = processedCopyNumbers$end))


bins <- tileGenome(seqinfo(Hsapiens), tilewidth=300000, cut.last.tile.in.chrom=TRUE)
bins <- bins[seqnames(bins)%in%chroms.ls]
#bins <- cyto_hg38_bins
length(bins)
# 10306

ol = as.matrix(findOverlaps(bins, processedCopyNumbers$region, type = "any"))
olCopyNumbers = cbind(bin = ol[, 1], processedCopyNumbers[ol[, 2], c("sampleId", "cancerType", "loh", "absDel", "relDel", "amp1_4", "amp2_0", "amp3_0")])

no_cores <- 4
cl <- new_cluster(no_cores)

#binSampleSummary = olCopyNumbers %>% partition(bin, sampleId, cluster = cl) %>%
binSampleSummary = olCopyNumbers %>% group_by(bin, sampleId) %>% partition(cluster = cl) %>%
		summarise(cancerType = dplyr::first(cancerType),
			loh = any(loh), absDel = any(absDel), relDel = any(relDel), amp1_4 = any(amp1_4), amp2_0 = any(amp2_0), amp3_0 = any(amp3_0))  %>%
		collect() %>%
		as_tibble()

binCancerTypeSummary = binSampleSummary %>% group_by( bin, cancerType ) %>% partition(cluster = cl) %>%
	summarise(loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0)) %>%
	collect() %>%
	as_tibble()

binSummary <- binSampleSummary %>% 
	group_by(bin) %>% 
	summarise(cancerType = 'All', loh = sum(loh), absDel = sum(absDel), relDel = sum(relDel), amp1_4 = sum(amp1_4), amp2_0 = sum(amp2_0), amp3_0 = sum(amp3_0))
binSummary <- bind_rows(binSummary, binCancerTypeSummary)

primaryTumorLocationCounts <- clinical.data[,list(nTotal=.N), by=list(cancerType=Histology)] 
primaryTumorLocationMaleCounts <- clinical.data[,list(nMale=sum(ifelse(gender=="MALE", 1, 0))), by=list(cancerType=Histology)]

primaryTumorLocationCounts <- merge(primaryTumorLocationCounts, primaryTumorLocationMaleCounts, by = 'cancerType', all = T)
primaryTumorLocationCounts[,nMale:=ifelse(is.na(nMale), 0, nMale)]

allCounts <- clinical.data[,list(nTotal=.N, cancerType='All')]
allCountsMale <- clinical.data[,list(nMale=sum(ifelse(gender=='MALE', 1, 0)), cancerType='All')]
allCounts <- merge(allCounts, allCountsMale, by = 'cancerType', all = T)

primaryTumorLocationCounts <- as.data.frame(rbind(primaryTumorLocationCounts, allCounts))
rm(primaryTumorLocationMaleCounts, allCounts, allCountsMale)

processedCopyNumberSummary <- as.data.frame(left_join(binSummary, primaryTumorLocationCounts, by = "cancerType"))
processedCopyNumberSummary$region <- bins[processedCopyNumberSummary$bin]
processedCopyNumberSummary$chromosome = substring(as.character(seqnames(processedCopyNumberSummary$region)), 4)
processedCopyNumberSummary$start = start(processedCopyNumberSummary$region)
processedCopyNumberSummary$end = end(processedCopyNumberSummary$region)

processedCopyNumberSummary = processedCopyNumberSummary %>%
		mutate(n = ifelse(chromosome %in% c('Y', 'chrY'), nMale, nTotal)) %>%
		mutate(lohPercentage = loh /n, absDelPercentage = absDel / n, relDelPercentage = relDel / n, amp1_4Percentage = amp1_4 / n,  amp2_0Percentage = amp2_0 / n, amp3_0Percentage = amp3_0 / n) %>%
		ungroup() %>%
		select(chromosome, start, end, cancerType, lohPercentage, absDelPercentage, relDelPercentage, amp1_4Percentage, amp2_0Percentage, amp3_0Percentage)

processedCopyNumberSummary.dt <- as.data.table(processedCopyNumberSummary)


save(list=c('processedCopyNumberSummary', 'primaryTumorLocationCounts', 'processedCopyNumberSummary.dt'),
	file = file.path(run.dir, 'result/CNA_circos/processedCopyNumberSummary.RData'))

# copy gaps.txt from one of the purple circos directories
load(file.path(run.dir, 'result/CNA_circos/processedCopyNumberSummary.RData'))

primaryTumorLocations = unique(processedCopyNumberSummary$cancerType)
for (location in primaryTumorLocations) {
	locationString = gsub(" ", "", location, fixed = TRUE)
	locationString = gsub("/", "", locationString, fixed = TRUE)
	
	primaryTumorLocationSummary = processedCopyNumberSummary %>%
			filter(cancerType == location) %>%
			mutate(chromosome = paste0("hs", chromosome))
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, lohPercentage) %>% mutate(lohPercentage = -lohPercentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".loh.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, absDelPercentage) %>% mutate(absDelPercentage = -absDelPercentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".absDel.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, relDelPercentage) %>% mutate(relDelPercentage = -relDelPercentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".relDel.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp1_4Percentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".amp1_4.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp2_0Percentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".amp2_0.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
	
	write.table(primaryTumorLocationSummary %>% select(chromosome, start, end, amp3_0Percentage),
			file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".amp3_0.circos"), sep = "\t", row.names = F, col.names = F, quote = F)
}


#tidyTargets <- function(copyNumberTargets) {
#	names = names(copyNumberTargets)
#	namesInd = setNames(1:ncol(copyNumberTargets), names)
#	
#	topTargets = copyNumberTargets %>% group_by(target) %>% summarise(N = sum(N)) %>% arrange(-N)
#	tidyCopyNumberTargets = copyNumberTargets %>% 
#			ungroup() %>%
#			filter(target %in% topTargets$target) %>%
#			mutate(gene = target) %>%
#			select(gene, chromosome, start, end, c(namesInd["Biliary"]:namesInd["Uterus"])) %>%
#			gather(cancerType, N, -gene, -chromosome, -start, -end) %>%
#			filter(!is.na(N)) %>%
#			group_by(gene, chromosome, cancerType) %>%
#			summarise(N = sum(N), start = min(start), end = max(end))
#	
#	return (tidyCopyNumberTargets)
#}


label_size <- function(genes) {
	label_size_single <- function(gene) {
		if (nchar(gene) < 7) {
			return ("label_size=40p")
		}
		
		if (nchar(gene) < 10) {
			return ("label_size=30p")
		}
		
		return ("label_size=24p")
	}
	
	sapply(genes,label_size_single)
	
}


label_size_by_driverLH <- function(driverLH) {
	if( driverLH>0.9 ){
		base.label.size <- 40
	}else{
		base.label.size <- 30
	}
	paste('label_size=', base.label.size, sep="");
}

label_size_by_gene_ct <- function(genect, label.sizes=c(40,30)) {
	if( genect>1 ){
		base.label.size <- label.sizes[[1]]
	}else{
		base.label.size <- label.sizes[[2]]
	}
	paste('label_size=', base.label.size, sep="");
}

load(file.path(HOME.dir, 'reference/HMF/38/dbs/linx/ensembl.gene.info.Rdata'))

merged.drivers.dt.filt <- merge(
	x = clinical.data[,list(tumor.sample.id, cancerType=Histology)],
	y = merged.drivers.dt.filt,
	by = c('tumor.sample.id'))

sample.gene.info <- merged.drivers.dt.filt[,list(
	cancerType, tumor.sample.id, gene, category, driverLikelihood, biallelic, AnnotationType)][
order(tumor.sample.id, -driverLikelihood),]


sample.gene.info[,all.gene.ct:=.N, by=list(gene)]
sample.gene.info[,cancerType.gene.ct:=.N, by=list(cancerType, gene)]
sample.gene.info[,gene.ct:=.N, by=list(tumor.sample.id)]
sample.gene.info[,idx:=1:unique(gene.ct), by=list(tumor.sample.id)]

#cancerType.gene.ct <- merged.drivers.dt.filt[,list(gene.ct=.N, driverLikelihood=max(driverLikelihood)),
#	by = list(cancerType, gene, category, biallelic, AnnotationType)][order(-gene.ct),]

#cancerType.gene.ct <- sample.gene.info[ (cancerType=='ccRCC' & (driverLikelihood>0.9 | idx==1)) | 
#	(cancerType!='ccRCC'  & (driverLikelihood>0.5 | idx==1)),
#list(gene.ct=.N, driverLikelihood=max(driverLikelihood)),
#by=list(cancerType, gene, category)][order(-gene.ct, -driverLikelihood),]

cancerType.gene.ct <- sample.gene.info[cancerType.gene.ct>1 | driverLikelihood==1 | idx==1,
		list(gene.ct=.N, driverLikelihood=max(driverLikelihood)),
		by=list(cancerType, gene, category)][order(cancerType, -gene.ct, -driverLikelihood),]

#all.gene.ct <- merged.drivers.dt.filt[,list(gene.ct=.N, cancerType='All'), by = list(gene, category)][order(-gene.ct),]

## all.gene.ct <- sample.gene.info[ (cancerType=='ccRCC' & (driverLikelihood>0.9 | idx==1)) | 
##                 (cancerType!='ccRCC'  & (driverLikelihood>0.5 | idx==1)),
##         list(gene.ct=.N, cancerType='All', driverLikelihood=max(driverLikelihood)),
##         by=list(gene, category)][order(-gene.ct, -driverLikelihood),]


all.gene.ct <- sample.gene.info[ all.gene.ct>1 | driverLikelihood==1 | idx==1,
	list(gene.ct=.N, cancerType='All', driverLikelihood=max(driverLikelihood)),
	by=list(gene, category)][order(-gene.ct, -driverLikelihood),]

gene.ct <- rbind(cancerType.gene.ct, all.gene.ct)
gene.ct[,info:=ifelse(category=='ONCO', "color=vdred", "color=black")]

label.sizes.ls <- c(60,40)

for (locationString in primaryTumorLocations) {
	gene.list.dt <- gene.ct[cancerType==locationString,list(gene, gene.ct, driverLikelihood, info)]
	
	if( nrow(gene.list.dt)>0 ){
		gene.list.dt <- merge(
			x = gene.list.dt,
			y = ensembl.gene.data.dt[GeneName%in%gene.list.dt$gene,list(chromosome=paste('hs', Chromosome, sep=""), gene=GeneName, start=GeneStart, end=GeneEnd)],
			by = c('gene'))

		gene.list.dt <- gene.list.dt[order(chromosome, -gene.ct, start),]
	
		genes <- gene.list.dt[,list(label_size = label_size_by_gene_ct(gene.ct, label.sizes.ls)),
			by=list(chromosome, start, end, gene, info)][,
			list(info=paste(info, label_size, sep = ",")),
			by =list(chromosome, start, end, gene)]
	
	}else{
		genes <- data.table(chromosome=NULL, start=NULL, end=NULL, gene=NULL, info=NULL)
	}

	write.table( genes,
		file = paste0(run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".genes.circos"),
		sep = "\t", row.names = F, col.names = F, quote = F)
}

command.file1 <- paste0(run.dir, '/sh/05_signature_analysis/CNA_analysis/CNA_circos/01_process_conf_cmd.sh', sep="")

# modified template gene r0 from 0.600r to 0.300r
for (locationString in primaryTumorLocations) {
	
	if( locationString==primaryTumorLocations[[1]]){
		curr.append <- FALSE
	}else{
		curr.append <- TRUE
	}
	curr.cmd = paste0("sed 's/CANCER/", locationString,"/g' circos_larger_chr.template > ", run.dir, "/result/CNA_circos/copyNumberSummary/", locationString, ".conf\n")

	fwrite(
		x = data.table(
			cmd = curr.cmd),
		file = command.file1,
		append = curr.append,
		col.names = FALSE,
		quote = FALSE)
}

fwrite(data.table(cancerType=primaryTumorLocations), file=file.path(run.dir, "config_files/CNA_circos_sample_info.csv"))
