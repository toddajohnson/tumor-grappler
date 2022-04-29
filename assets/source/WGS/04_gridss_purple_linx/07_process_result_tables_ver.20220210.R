#!/usr/bin/env Rscript

#module use /usr/local/package/modulefiles/
#module load R/4.0.2
# or
# . ~/.R-4.1.0_setup.sh

### NEED TO RUN SCRIPT problem_sample.info.R in 00_Preparation folder
### to output problematic samples: no matched normal, same CNA profile in tumor and normal (QDNAseq)
###
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
	library(parallel)})

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
tumor.normal.pair.info <- tumor_normal_pair_info.GPL

message( paste("Loading merge results tables for ", study.name, sep=""))

load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

#germline.variant.info.dt2 <- germline.variant.info.dt[is.na(AF_EAS) | (MAF_EAS<maf.cutoff & MAF_TOT<maf.cutoff & AF_EAS<af.cutoff),]
#print(paste('Filter on AF_EAS, MAF_EAS, and MAF_TOT: ', nrow(germline.variant.info.dt2), ' variants', sep=''))
#	"Filter on AF_EAS, MAF_EAS, and MAF_TOT: 21008 variants"
#print(paste('Filter on AF_EAS, MAF_EAS, MAF_TOT, TOMMO, KOREAN, and dbGaP_PopFreq: ', nrow(germline.variant.info.dt1), ' variants', sep=''))
#	"Filter on AF_EAS, MAF_EAS, MAF_TOT, TOMMO, KOREAN, and bGaP_PopFreq: 5883 variants"

load(file = file.path(cg_pipeline_dir, 'common_config_files/DriverGenePanel.38.list.Rdata'))

if (study.name%in%c('BHD', 'BTC', 'KakimiLab', 'KeioOrganoid') ){
	gene_panel.dt <- driver_gene_panel.dt.ls[[study.name]]
}else{
	gene_panel.dt <- driver_gene_panel.dt.ls[['Other']]
}

message( paste("Loading all germline variant info ", study.name, "\n", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/germline.variant.info.Rdata", sep=""))
## Added filter for germline variants not betwee 0.01 and 0.99 AF and not >= 0.99 AF ALFA EAS, ALFA TOT, or dbSNP 155 dbGaP (should be ALFA TOT), TOMMO, or KOREAN
germline.variant.info.dt <- germline.variant.info.dt[
	!((AF_TOT>=maf.cutoff & AF_TOT<af.cutoff) | AF_TOT>=af.cutoff | 
	(AF_EAS>=maf.cutoff & AF_EAS<af.cutoff) | AF_EAS>=af.cutoff | 
	(TOMMO>=maf.cutoff & TOMMO<af.cutoff) | TOMMO>=af.cutoff |
	(KOREAN>=maf.cutoff & KOREAN<af.cutoff) | KOREAN>=af.cutoff |
	(dbGaP_PopFreq>=maf.cutoff & dbGaP_PopFreq<af.cutoff) | dbGaP_PopFreq>=af.cutoff),]

# filter out germline variants that are not in germline gene panel
germline.variant.info.dt <- germline.variant.info.dt[gene%in%gene_panel.dt$gene,]
germline.variant.funcClass2.summary <- germline.variant.info.dt[,list(sample.SEW.funcClass2.ct=.N), by=list(subject.id, tumor.sample.id, SEW.funcClass2)]
germline.variant.info.dt.candidates <- germline.variant.info.dt[SEW.funcClass2%in%c('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE'),]

load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.info.Rdata", sep=""))
somatic.variant.funcClass2.summary <- somatic.variant.info.dt[,list(sample.SEW.funcClass2.ct=.N), by=list(subject.id, tumor.sample.id, SEW.funcClass2)]
somatic.variant.info.dt.candidates <- somatic.variant.info.dt[SEW.funcClass2%in%c('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE'),]

#modify 2021/10/09 to filter germline and somatic additional candidates to those just in gene panel
germline.variant.info.dt.rdcd <- germline.variant.info.dt[,
	list(subject.index, tumor.index, subject.id, tumor.sample.id,
		snp.id, rsid=ID, chr, CHROM, start, end, REF, ALT,
		AF_TOT, AF_EAS, TOMMO, KOREAN,
		AF.N, AF.T, BIALLELIC, PURPLE_AF, PURPLE_CN, PURPLE_VCN,
		PATH, gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]

somatic.variant.info.dt.rdcd <- somatic.variant.info.dt[,
	list(subject.index, tumor.index, subject.id, tumor.sample.id, snp.id, rsid=ID, chr, CHROM, POS=start, REF, ALT,
		AF_TOT, AF_EAS, TOMMO, KOREAN, AF, AF_raw, AF_eas,
		AF.N, AF.T, BIALLELIC, PURPLE_VCN, PURPLE_AF, PURPLE_CN, PURPLE_MACN, PURPLE_GERMLINE, SUBCL,
		gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]

message( paste("Extracting annotations from germline driver variants for ", study.name, "\n", sep=""))

if (length(germline.driver.variant.info.dt)>0){
	germline.driver.variant.info.dt.rdcd <- germline.driver.variant.info.dt[,
		list(subject.index, tumor.index, subject.id, tumor.sample.id,
		snp.id, rsid=ID, chr, CHROM, start, end, REF, ALT,
		AF_TOT, AF_EAS, TOMMO, KOREAN,
		AF.N, AF.T, BIALLELIC, PURPLE_AF, PURPLE_CN, PURPLE_VCN,
		PATH, gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]
	
	germline.driver.variant.info.dt.rdcd <- germline.driver.variant.info.dt.rdcd[order(tumor.index, chr, start),]
	
	germline.driver.genes.ls <- unique(germline.driver.variant.info.dt.rdcd$gene)
}else{
	germline.driver.genes.ls <- c()
}

message( paste("Extracting SEC annotations from somatic driver variants for ", study.name, "\n", sep=""))

somatic.driver.variant.info.dt.rdcd <- somatic.driver.variant.info.dt[,
	list(subject.index, tumor.index, subject.id, tumor.sample.id, snp.id, rsid=ID, chr, CHROM, POS=start, REF, ALT,
	AF_TOT, AF_EAS, TOMMO, KOREAN, AF, AF_raw, AF_eas,
	AF.N, AF.T, BIALLELIC, PURPLE_VCN, PURPLE_AF, PURPLE_CN, PURPLE_MACN, PURPLE_GERMLINE, SUBCL,
	gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]

somatic.driver.variant.info.dt.rdcd <- somatic.driver.variant.info.dt.rdcd[order(tumor.index, chr, POS),]
somatic.variant.ct <- nrow(somatic.driver.variant.info.dt.rdcd)

somatic.driver.genes.ls <- unique(somatic.driver.variant.info.dt.rdcd$gene)

message( paste( "There were ", somatic.variant.ct, " somatic variants across the ", study.name, " samples", sep="") )


message( paste("Filtering merged driver table for ", study.name, "\n", sep=""))

setnames(merged.drivers.dt, old='chromosome', new='CHROM')
merged.drivers.dt[,chr:=chrom.map.ls[CHROM]]

setnames(purple.cnv.gene.info, old='chromosome', new='CHROM')
purple.cnv.gene.info[,chr:=chrom.map.ls[CHROM]]

linx.vis_gene_exon.info.summary <- linx.vis_gene_exon.info.dt[,
	list(linx.start=min(ExonStart), linx.end=max(ExonEnd),
		clusterIds=paste(unique(ClusterId), collapse=",")),
	by = list(subject.id, tumor.sample.id, subject.index, tumor.index, gene=Gene, Chromosome, AnnotationType)]

linx.drivers.dt <- merge(
	x = linx.drivers.dt,
	y = linx.vis_gene_exon.info.summary,
	by = c("subject.id", "tumor.sample.id", "subject.index", "tumor.index", "gene"),
	all.x = TRUE)

non.driver.disruption.info <- linx.vis_gene_exon.info.summary[AnnotationType!='DRIVER',]

merged.drivers.dt <- merge(
	x = merged.drivers.dt,
	y = linx.drivers.dt,
	by = c("subject.id", "tumor.sample.id", "subject.index", "tumor.index", "gene"),
	all = TRUE)

merged.drivers.dt[,eventType:=ifelse(is.na(eventType), 'none', eventType)]
merged.drivers.dt[,clusterId:=ifelse(is.na(clusterId), -99L, clusterId)]

merged.drivers.dt.filt <- merged.drivers.dt[
	gene%in%germline.driver.genes.ls |
	gene%in%somatic.driver.genes.ls |
	eventType!='none' | driver%in%c("AMP", "DEL", "PARTIAL_AMP"),]
merged.drivers.dt.filt <- merged.drivers.dt.filt[order(subject.index, tumor.index, chr),]

message( paste("Adding variant status flags to merged driver table for ", study.name, "\n", sep=""))

merged.drivers.dt.filt[,germline.status:=ifelse(driver=="GERMLINE", 1, 0)]
merged.drivers.dt.filt[,mutation.status:=ifelse(driver=="MUTATION", 1, 0)]
merged.drivers.dt.filt[,deletion.status:=ifelse(driver%in%c("DEL", "HOM_DISRUPTION"), 1, 0)]
#merged.drivers.dt.filt[,deletion.status:=ifelse(driver=="DEL", 1, 0)]
merged.drivers.dt.filt[,amplification.status:=ifelse(driver%in%c("AMP", "PARTIAL_AMP"), 1, 0)]

linx.annotated.cluster.info.dt <- merge(
	x = linx.clusters.dt,
	y = linx.vis_sv_data.info.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, clusterId=ClusterId,
		svId=SvId, Type, resolvedType=ResolvedType, ChrStart, ChrEnd, PosStart, PosEnd, OrientStart, OrientEnd, InfoStart, InfoEnd,
		JunctionCopyNumber, InDoubleMinute)],
	by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'clusterId', 'resolvedType'))

reportable.biotypes <- c('protein_coding')

linx.annotated.cluster.info.dt <- merge(
	x = linx.annotated.cluster.info.dt[category!='ARTIFACT'],
#	y = linx.breakend.dt[reportedDisruption==TRUE,
#	y = linx.breakend.dt[disruptive==TRUE,
	y = linx.breakend.dt[disruptive==TRUE & biotype%in%reportable.biotypes,
		list(subject.id, tumor.sample.id, subject.index, tumor.index, svId, disrupted.gene=gene, reportedDisruption, regionType)],
	by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'svId'),
	all.x = TRUE)
#linx.annotated.cluster.info.dt[,disrupted.gene:=ifelse(is.na(disrupted.gene), 'none', disrupted.gene)]

linx.annotated.cluster.info.dt[,disrupted.gene:=ifelse(is.na(disrupted.gene), 'none', disrupted.gene)]
#linx.annotated.cluster.info.dt[,reportedDisruption:=ifelse(is.na(reportedDisruption), FALSE, reportedDisruption)]

merge_disrupted_gene_info <- function( curr.reportedDisruptions, curr.genes ){
	curr.genes <- unique(curr.genes[which(curr.genes!='none')])
	curr.gene.ct <- length(curr.genes)
	
	if (TRUE%in%curr.reportedDisruptions){
		curr.reportedDisruption.status <- TRUE
	}else{
		curr.reportedDisruption.status <- FALSE
	}
	
	if (curr.gene.ct>0){
		curr.genes.str <- paste(curr.genes[order(curr.genes)], collapse=";")
	}else{
		curr.genes.str <- 'none'
	}
	data.table(
		disrupted.genes = curr.genes.str,
		disrupted.gene.ct = curr.gene.ct,
		reportedDisruption = curr.reportedDisruption.status)
}

linx.SV.CN.info.dt <- linx.annotated.cluster.info.dt[,
	list(
		clusterDesc = ifelse(clusterCount>1, clusterDesc,
			ifelse(ChrStart==ChrEnd,
				ifelse( is.na(disrupted.gene),
					paste(unique(clusterDesc), "(", PosEnd - PosStart + 1, "bp)", ":", unique(paste(ChrStart, ":", PosStart, "-", PosEnd, sep="")), sep=""),
					paste(unique(clusterDesc), "(", disrupted.gene, " ", PosEnd - PosStart + 1, "bp)", ":", unique(paste(ChrStart, ":", PosStart, "-", PosEnd, sep="")), sep="")),
				paste(unique(clusterDesc), ":", unique(paste(ChrStart, ":", PosStart, "-", ChrEnd, ":", PosEnd, sep="")), sep="")))),
		by = list(subject.id, tumor.sample.id, subject.index, tumor.index, clusterId, disrupted.gene, resolvedType, category, reportedDisruption, clusterCount)]

linx.SV.CN.info.dt <- linx.SV.CN.info.dt[,merge_disrupted_gene_info(reportedDisruption, disrupted.gene),
	by=list(subject.id, tumor.sample.id, subject.index, tumor.index, clusterId, resolvedType, category, clusterCount, clusterDesc)]

linx.SV.CN.drivers.info.dt <- merge(
	x = linx.drivers.dt,
	y = linx.SV.CN.info.dt,
	by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'clusterId'),
	all.x = TRUE)

linx.SV.CN.drivers.info.dt <- linx.SV.CN.drivers.info.dt[,list(
	eventType,
	cluster.SV.ct=ifelse(is.na(clusterCount), 0L, clusterCount),
	disrupted.genes=ifelse(is.na(disrupted.genes), 'none', disrupted.genes),
	disrupted.gene.ct=ifelse(is.na(disrupted.gene.ct), 0L, disrupted.gene.ct),
	clusterDesc=ifelse(is.na(clusterDesc), 'ARM_CHR_LEVEL', clusterDesc),
	resolvedType=ifelse(is.na(resolvedType), 'N/A', resolvedType)),
	by = list(subject.id, tumor.sample.id, subject.index, tumor.index, ClusterId=clusterId, gene)]

message( paste("Extracting copy-number drivers from merged driver table for ", study.name, "\n", sep=""))

whole.gene.CN.info <- merged.drivers.dt.filt[eventType!='none' | biallelic==TRUE | driver%in%c('AMP', 'DEL', 'PARTIAL_AMP'),
	list(subject.id, tumor.sample.id, subject.index, tumor.index, chr, CHROM, chromBand=chromosomeBand, gene, driver, category, eventType, ClusterId=clusterId)]
#list(subject.id, tumor.sample.id, subject.index, tumor.index, chr, CHROM, chromBand=chromosomeBand, gene, driver, category, eventType, clusterId)]


message( paste("Merging copy-number drivers with purple CNV gene info for ", study.name, "\n", sep=""))

merged.SV.CN.info.dt <- merge(
	x = whole.gene.CN.info,
	y = linx.SV.CN.drivers.info.dt,
	by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'ClusterId', 'eventType', 'gene'),
	all.x = TRUE)
merged.SV.CN.info.dt[,eventType:=ifelse(eventType=='none', 'biallelic', eventType)]
merged.SV.CN.info.dt[,cluster.SV.ct:=ifelse(is.na(cluster.SV.ct), 0, cluster.SV.ct)]
merged.SV.CN.info.dt[,clusterDesc:=ifelse(is.na(clusterDesc), 'N/A', clusterDesc)]
merged.SV.CN.info.dt[,resolvedType:=ifelse(is.na(resolvedType), 'N/A', resolvedType)]
merged.SV.CN.info.dt[,disrupted.genes:=ifelse(is.na(disrupted.genes), 'none', disrupted.genes)]
merged.SV.CN.info.dt[,disrupted.gene.ct:=ifelse(is.na(disrupted.gene.ct), 0, disrupted.gene.ct)]

driver.purple.cnv.gene.info <- merge(
	x = merged.SV.CN.info.dt,
	y = purple.cnv.gene.info,
	by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'chr', 'CHROM', 'gene'))

SV.CN.variant.genes.ls <- unique(driver.purple.cnv.gene.info$gene)

HRD.genes <- c('BRCA1', 'BRCA2', 'PALB2', 'RAD51C', 'RAD51D', 'ATM', 'BRIP1')


message( paste("Merging germline and somatic driver variant info for ", study.name, "\n", sep=""))

# 2022/01/20: merged driver variant table should included SV variants
if (length(germline.driver.genes.ls)>0){
	merged.driver.variant.info.dt <- merge(
		x = germline.driver.variant.info.dt.rdcd[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, 
			snp.id.gl=snp.id, POS.gl=start, REF.gl=REF, ALT.gl=ALT,
			AF.N.gl=AF.N, AF.T.gl=AF.T,
			BIALLELIC.gl=BIALLELIC, P_AF.gl=PURPLE_AF, MACN.gl=PURPLE_CN, VCN.gl=PURPLE_VCN,
			func1.gl=funcClass1, func2.gl=funcClass2, bp.descr.gl=bp.descr, aa.descr.gl=aa.descr)],
		y = somatic.driver.variant.info.dt.rdcd[PURPLE_AF>0,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, 
			snp.id.som=snp.id, POS.som=POS, REF.som=REF, ALT.som=ALT,
			AF.N.som=AF.N, AF.T.som=AF.T,
			BIALLELIC.som=BIALLELIC, P_AF.som=PURPLE_AF, CN.som=PURPLE_CN, VCN.som=PURPLE_VCN, MACN.som=PURPLE_MACN, GERMLINE.som=PURPLE_GERMLINE, SUBCL=SUBCL,
			func1.som=funcClass1, func2.som=funcClass2, bp.descr.som=bp.descr, aa.descr.som=aa.descr)],
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'chr', 'CHROM', 'gene'),
		all = TRUE)
}else{
	merged.driver.variant.info.dt <- unique(merge(
		x = somatic.driver.variant.info.dt.rdcd[PURPLE_AF>0,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, 
			snp.id.gl='none', POS.gl=as.integer(NA), REF.gl=as.character(NA), ALT.gl=as.character(NA),
			AF.N.gl=-99, AF.T.gl=-99,
			BIALLELIC.gl=as.character(NA), P_AF.gl=-99, MACN.gl=-99, VCN.gl=-99,
			func1.gl=as.character(NA), func2.gl=as.character(NA), bp.descr.gl=as.character(NA), aa.descr.gl=as.character(NA))],
		y = somatic.driver.variant.info.dt.rdcd[PURPLE_AF>0,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, 
			snp.id.som=snp.id, POS.som=POS, REF.som=REF, ALT.som=ALT,
			AF.N.som=AF.N, AF.T.som=AF.T,
			BIALLELIC.som=BIALLELIC, P_AF.som=PURPLE_AF, CN.som=PURPLE_CN, VCN.som=PURPLE_VCN, MACN.som=PURPLE_MACN, GERMLINE.som=PURPLE_GERMLINE, SUBCL=SUBCL,
			func1.som=funcClass1, func2.som=funcClass2, bp.descr.som=bp.descr, aa.descr.som=aa.descr)],
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'chr', 'CHROM', 'gene'),
		all = TRUE))
}

if (length(SV.CN.variant.genes.ls)>0 ){
	# ? merged SV variant into into somatic driver info?
	SV.CN.driver.variant.info.dt <- driver.purple.cnv.gene.info[,
		list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, minRegionStart, minRegionEnd,
			gene, ClusterId, eventType, driver, category, SV_ct=cluster.SV.ct, clusterDesc, resolvedType,
			max_CN=maxCopyNumber, min_minorCN=minMinorAlleleCopyNumber)]
	
	merged.driver.variant.info.dt <- unique(merge(
		x = merged.driver.variant.info.dt,
		y = SV.CN.driver.variant.info.dt,
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'chr', 'CHROM', 'gene'),
		all = TRUE))

	## worst effect for SV/CN variants
	SV_CN.func.summary <- linx.SV.CN.info.dt[,list(sv.ct=sum(clusterCount), cluster.ct=.N, SV.disrupted.gene.ct=sum(disrupted.gene.ct)), 
		by=list(subject.id, tumor.sample.id)]
}else{
	merged.driver.variant.info.dt <- unique(merge(
		x = merged.driver.variant.info.dt,
		y = merged.driver.variant.info.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, minRegionStart=-99L, minRegionEnd=-99L,
			gene, ClusterId=as.integer(NA), eventTypeas.character(NA), driveras.character(NA), as.character(NA), SV_ct=as.integer(NA), clusterDesc=as.character(NA), resolvedType=as.character(NA),
			max_CN=as.numeric(NA), min_minorCN=as.numeric(NA))],
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'chr', 'CHROM', 'gene'),
		all = TRUE))

	SV_CN.func.summary <- merged.driver.variant.info.dt[,
		list(sv.ct=0, cluster.ct=0, SV.disrupted.gene.ct=0),
		by=list(subject.id, tumor.sample.id)]
}
merged.driver.variant.info.dt <- merged.driver.variant.info.dt[order(tumor.index, chr, POS.gl, POS.som, minRegionStart),]


message( paste("Summarizing hit information for germline and somatic driver variant info for ", study.name, "\n", sep=""))

merged.driver.variant.info.dt[,eventType:=ifelse(is.na(eventType), 'N/A', eventType)]
merged.driver.variant.info.dt[,clusterDesc:=ifelse(is.na(clusterDesc), 'N/A', clusterDesc)]
merged.driver.variant.info.dt[,snp.id.gl:=ifelse(is.na(snp.id.gl), 'none', snp.id.gl)]
merged.driver.variant.info.dt[,AF.N.gl:=ifelse(is.na(AF.N.gl), -99, AF.N.gl)]
merged.driver.variant.info.dt[,AF.T.som:=ifelse(is.na(AF.T.som), -99, AF.T.som)]

merged.driver.variant.info.dt[,gene.GL:=0L]
merged.driver.variant.info.dt[,gene.SOM:=0L]
merged.driver.variant.info.dt[,gene.LOH:=0L]
merged.driver.variant.info.dt[,gene.som.LOH:=0L]
merged.driver.variant.info.dt[,gene.som.SV_CN:=0L]
merged.driver.variant.info.dt[,gene.som.DEL:=0L]
merged.driver.variant.info.dt[,gene.som.GAIN:=0L]
merged.driver.variant.info.dt[,gene.GL:=ifelse(AF.N.gl>0, 1L, gene.GL)]
merged.driver.variant.info.dt[,gene.SOM:=ifelse(AF.T.som>0, 1L, gene.SOM)]
merged.driver.variant.info.dt[which(BIALLELIC.gl==TRUE),gene.LOH:=1L]

merged.driver.variant.info.dt[grep('LOH', eventType),gene.som.LOH:=1L]
merged.driver.variant.info.dt[which(eventType=='biallelic'),gene.som.LOH:=1L]

merged.driver.variant.info.dt[which(SV_ct>0 | gene.som.LOH>0 | eventType=='biallelic'),gene.som.SV_CN:=1L]

merged.driver.variant.info.dt[grep('DEL', eventType),gene.som.DEL:=1L]
merged.driver.variant.info.dt[grep('GAIN', eventType),gene.som.GAIN:=1L]


merged.driver.variant.info.dt[,gene.GL.hits:=sum(gene.GL), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.SOM.hits:=sum(gene.SOM), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.LOH.hits:=sum(gene.LOH), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.som.LOH.hits:=sum(gene.som.LOH), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.som.SV_CN.hits:=sum(gene.som.SV_CN), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.som.DEL.hits:=sum(gene.som.DEL), by=list(subject.id, tumor.sample.id, gene)]
merged.driver.variant.info.dt[,gene.som.GAIN.hits:=sum(gene.som.GAIN), by=list(subject.id, tumor.sample.id, gene)]

#merged.driver.variant.info.dt[,gene.hits:=gene.GL.hits + gene.SOM.hits + gene.LOH.hits + gene.som.LOH.hits]
merged.driver.variant.info.dt[,gene.hits:=gene.GL.hits + gene.SOM.hits + gene.LOH.hits + gene.som.SV_CN.hits]
merged.driver.variant.info.dt[,max.gene.hits:=max(gene.hits), by=list(subject.id, tumor.sample.id)]

merged.driver.variant.info.dt[,fs.hit:=ifelse(is.na(func2.som), 0L, ifelse(func2.som=='NONSENSE_OR_FRAMESHIFT', 1L, 0L))]
merged.driver.variant.info.dt[,max.fs.hit:=max(fs.hit), by=list(tumor.sample.id)]

merged.driver.variant.info.dt[,gene.GL.func1.hits:=sum(gene.GL), by=list(subject.id, tumor.sample.id, gene, func1.gl)]
merged.driver.variant.info.dt[,gene.SOM.func1.hits:=sum(gene.SOM), by=list(subject.id, tumor.sample.id, gene, func1.som)]

merged.driver.variant.info.dt[,func1.gl.ann:=ifelse(is.na(func1.gl), NA, paste(func1.gl, "(", gene.GL.func1.hits, ")", sep=""))]
merged.driver.variant.info.dt[,func1.som.ann:=ifelse(is.na(func1.som), NA, paste(func1.som, "(", gene.SOM.func1.hits, ")", sep=""))]

message( paste("Labelling HRD genes in germline and somatic driver variant info for ", study.name, "\n", sep=""))
merged.driver.variant.info.dt[,HRD.gene.hits:=ifelse(gene%in%HRD.genes, gene.hits, 0L)]



message( paste("Exporting summary table of copy-number related driver genes for ", study.name, "\n", sep=""))
fwrite(driver.purple.cnv.gene.info[order(tumor.index),list(subject.id, tumor.sample.id, subject.index, tumor.index,
	chr, CHROM, gene, driver, category, eventType, clusterId=ClusterId, SV_ct=cluster.SV.ct, minRegionStart, minRegionEnd, maxCopyNumber, minMinorAlleleCopyNumber)],
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.somatic.driver.CNV.gene.summary.txt", sep=""),
	col.names = TRUE,
	sep = '\t')


fusions.all.dt <- linx.fusion.dt[order(tumor.index),
	list(subject.id, tumor.sample.id,
	fivePrimeBreakendId, threePrimeBreakendId, name, reported, reportedType, phased, likelihood, chainLength, chainLinks, chainTerminated, domainsKept, domainsLost, skippedExonsUp, skippedExonsDown, fusedExonUp, fusedExonDown, geneStart, geneContextStart, transcriptStart, geneEnd, geneContextEnd, transcriptEnd, junctionCopyNumber)]


fusions.reported.dt <- linx.vis_fusion.info.dt[order(tumor.index), list(subject.id, tumor.sample.id,
	ClusterId, Reportable, GeneNameUp, TranscriptUp, ChrUp, PosUp, StrandUp, RegionTypeUp, FusedExonUp,
	GeneNameDown, TranscriptDown, ChrDown, PosDown, StrandDown, RegionTypeDown, FusedExonDown)]

message( paste("Exporting all fusion table for ", study.name, "\n", sep=""))
fwrite(fusions.all.dt,
		file = paste(run.dir, "/result_summaries/", gpl_prefix, "/fusions.all.txt", sep=""),
		col.names = TRUE,
		sep = '\t')

message( paste("Exporting reported fusion table for ", study.name, "\n", sep=""))
fwrite(fusions.reported.dt,
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/fusions.reported.txt", sep=""),
	col.names = TRUE,
	sep = '\t')

message( paste("Exporting driver gene table for ", study.name, "\n", sep=""))
fwrite(merged.drivers.dt,
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.gene.summary.txt", sep=""),
	col.names = TRUE,
	sep = '\t')

#SV_CN.type.som=eventType, SV_CN.descr.som=clusterDesc, 
merged.driver.variant.info.dt[,snp.id.som:=ifelse(is.na(snp.id.som), 'none', snp.id.som)]


merged.driver.variant.info.dt[,
	SV_CN.type.som:=paste(eventType, collapse=";"),
	by = list( subject.id, tumor.sample.id, CHROM, gene, snp.id.gl, snp.id.som)]

merged.driver.variant.info.dt[,
		SV_CN.descr.som:=paste(clusterDesc, collapse=";"),
	by = list( subject.id, tumor.sample.id, CHROM, gene, snp.id.gl, snp.id.som)]

#SV_CN.descr.som:=paste(paste(clusterDesc, ' (ClusterID=', ClusterId, ')', sep=""), collapse=";"),

message( paste("Exporting merged table of candidate driver variant info for ", study.name, "\n", sep=""))
fwrite(unique(merged.driver.variant.info.dt[order(tumor.index, -gene.hits, chr),list(
	gene.hits, HRD.gene.hits, gene.GL.hits, gene.SOM.hits, gene.LOH.hits, gene.som.LOH.hits,
	subject.id, tumor.sample.id,
	CHROM, gene,
	snp.id.gl,
	AF.N.gl, AF.T.gl,
	BIALLELIC.gl, P_AF.gl, MACN.gl, VCN.gl,
	func1.gl, func2.gl, bp.descr.gl, aa.descr.gl, ph1="",
	snp.id.som,
	AF.N.som, AF.T.som,
	BIALLELIC.som, P_AF.som, CN.som, VCN.som, MACN.som, GERMLINE.som, SUBCL,
	func1.som, func2.som, bp.descr.som, aa.descr.som, SV_CN.type.som, SV_CN.descr.som, min_minorCN, max_CN)]),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.variant.info.txt", sep=""),
	col.names = TRUE,
	sep = '\t')


message( paste("Extracting hit annotated driver variant info for ", study.name, "\n", sep=""))

SV_CN.descr.str.ls <- c(
		'N/A',
		'unann. LOH',
		'amplification by SV', 'amplification of whole arm', 'amplification of whole chromosome',
		'focal homozygous deletion',
		'focal LOH', 'arm level LOH', 'chromosome level LOH', 'LOH from SV to telomere', 'LOH from SV to centromere',
		'homozygous disruption via cross exonic tandem duplication', 'homozygous disruption without homozygous copy number loss')
names(SV_CN.descr.str.ls) <- c(
		'N/A',
		'biallelic',
		'GAIN', 'GAIN_ARM', 'GAIN_CHR',
		'DEL',
		'LOH', 'LOH_ARM', 'LOH_CHR', 'LOH_SV_TELO', 'LOH_SV_CENTRO',
		'HOM_DUP_DISRUPTION', 'HOM_DEL_DISRUPTION')

merged.driver.variant.info.dt[,SV_CN.str:=SV_CN.descr.str.ls[eventType]]

summarized.driver.variant.info <- unique(merged.driver.variant.info.dt[,
	list(
		func1.gl=ifelse(snp.id.gl=='none', 'N/A', paste(unique(func1.gl.ann), collapse=";")),
		snp.descr.gl=ifelse(snp.id.gl=='none', 'No reported germline variant', paste(snp.id.gl, ' (', bp.descr.gl, ';', aa.descr.gl, ')', sep="")),
		func1.som=ifelse(is.na(func1.som.ann),
			"No reported somatic variant",
			paste(unique(func1.som.ann), collapse=";")),
		snp.descr.som=ifelse(is.na(snp.id.som), 'N/A', ifelse(snp.id.som=='none', 'N/A', paste(snp.id.som, ' (', bp.descr.som, ';', aa.descr.som, ')', sep=""))),
		SV_CN.som=ifelse(eventType=='N/A', 'No reported SV/CN', paste(unique(eventType), collapse=";")),
		SV_CN.str.som=ifelse(SV_CN.str=='N/A', 'No reported SV/CN', paste(unique(SV_CN.str), collapse=" & ")),
		SV_CN.descr=ifelse(eventType=='N/A', 'N/A', ifelse(clusterDesc=='N/A', 'No reported SVs', paste(unique(clusterDesc), collapse="\n")))),
	by = list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, max.gene.hits, gene.hits, HRD.gene.hits, gene.GL.hits, gene.SOM.hits, gene.LOH.hits, gene.som.LOH.hits,
			gene.som.SV_CN.hits, gene.som.DEL.hits, gene.som.GAIN.hits)])

merge_Variant.values <- function( curr.values ){
	curr.values <- curr.values[which(curr.values!='N/A')]
	if( length(curr.values)>0 ){
		curr.values <- unique(curr.values)
		paste(curr.values, collapse="\n")
	}else{
		'N/A'
	}
}

summarized.driver.variant.info.summary <- unique(summarized.driver.variant.info[,
	list(
		Description.gl = ifelse(snp.descr.gl=='No reported germline variant', 'N/A', paste( unique(func1.gl), collapse=';')),
		Variant.gl = ifelse(snp.descr.gl=='No reported germline variant', 'No reported germline variant', paste(unique(snp.descr.gl), collapse="\n")),
		Description.som = ifelse(func1.som=='No reported somatic variant',
				ifelse(SV_CN.str.som=='No reported SV/CN', 'No reported somatic small or SV/CN variants', paste(unique(SV_CN.str.som), collapse='\n')),
				ifelse(SV_CN.str.som=='No reported SV/CN', paste(unique(func1.som), collapse=';'), paste(paste(unique(func1.som), collapse=';'), ' with ', unique(SV_CN.str.som), sep=""))),
#		Variant.cmbd.som = ifelse(func1.som=='No reported somatic variant',
#				ifelse(SV_CN.str.som=='No reported SV/CN' | SV_CN.descr=='No reported SVs', 'N/A', paste(unique(SV_CN.descr), collapse='\n')),
#				ifelse(SV_CN.str.som=='No reported SV/CN' | SV_CN.descr=='No reported SVs', paste(unique(snp.descr.som), collapse='\n'), paste( c(unique(snp.descr.som), unique(SV_CN.descr)), collapse='\n'))),
		Variant.som = ifelse(func1.som=='No reported somatic variant', 'N/A', paste(unique(snp.descr.som), collapse='\n')),
		Variant.SV_CN.som = ifelse(SV_CN.str.som=='No reported SV/CN' | SV_CN.descr=='No reported SVs', 'N/A', paste(unique(SV_CN.descr), collapse='\n'))),
	by = list(subject.id, tumor.sample.id, subject.index, tumor.index,
			chr, CHROM, gene, max.gene.hits, gene.hits, HRD.gene.hits, gene.GL.hits, gene.SOM.hits, gene.LOH.hits, gene.som.LOH.hits,
			gene.som.SV_CN.hits, gene.som.DEL.hits, gene.som.GAIN.hits)])

summarized.driver.variant.info.summary[,Description.gl:=ifelse(Description.gl!='N/A', paste('Germline: ', Description.gl, sep=""), Description.gl)]
summarized.driver.variant.info.summary[,Variant.gl:=ifelse(Variant.gl!='No reported germline variant', paste('Germline: ', Variant.gl, sep=""), Variant.gl)]

summarized.driver.variant.info.summary[,Description.som:=ifelse(Description.som!='No reported somatic small or SV/CN variants', paste('Somatic: ', Description.som, sep=""), Description.som)]
summarized.driver.variant.info.summary[,Variant.som:=ifelse(Variant.som!='N/A', paste('Somatic: ', Variant.som, sep=""), Variant.som)]
summarized.driver.variant.info.summary[,Variant.SV_CN.som:=ifelse(Variant.SV_CN.som!='N/A', paste('SV: ', Variant.SV_CN.som, sep=""), Variant.SV_CN.som)]

summarized.driver.variant.info.summary[,Variant.cmbd.som:=ifelse(Variant.som!='N/A', ifelse(Variant.SV_CN.som!='N/A', paste(Variant.som, Variant.SV_CN.som, sep="\n"), Variant.som),
	ifelse(Variant.SV_CN.som!='N/A', Variant.SV_CN.som, 'N/A'))]


message( paste("Exporting summarized table of candidate driver variant info for ", study.name, "\n", sep=""))
fwrite(summarized.driver.variant.info[order(tumor.index, -gene.hits, chr),],
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.variant.info.txt", sep=""),
	col.names = TRUE,
	sep = '\t')


message( paste("Constructing summary of hits for candidate driver genes for ", study.name, "\n", sep=""))
driver.summary <- merge(
	x = summarized.driver.variant.info.summary[Variant.gl!='No reported germline variant',
		list(Description=Description.gl, Variant=Variant.gl),
		by = list(subject.index, tumor.index, subject.id, tumor.sample.id, chr, CHROM, gene, HRD.gene.hits)],
	y = summarized.driver.variant.info.summary[Description.som!='No reported somatic small or SV/CN variants',
		list(Description.som, Variant.som, Variant.SV_CN.som, Variant.cmbd.som),
		by = list(subject.index, tumor.index, subject.id, tumor.sample.id, chr, CHROM, gene, HRD.gene.hits)],
	by = c('subject.index', 'tumor.index', 'subject.id', 'tumor.sample.id', 'chr', 'CHROM', 'gene', 'HRD.gene.hits'),
	all = TRUE)

driver.summary[,Description.som:=ifelse(is.na(Description.som), 'No reported somatic small or SV/CN variants', Description.som)]
driver.summary[,Variant.som:=ifelse(is.na(Variant.som), 'N/A', Variant.som)]
driver.summary[,Variant.SV_CN.som:=ifelse(is.na(Variant.SV_CN.som), 'N/A', Variant.SV_CN.som)]
driver.summary[,Variant.cmbd.som:=ifelse(is.na(Variant.cmbd.som), 'N/A', Variant.cmbd.som)]

driver.summary[,Description.mrg:=ifelse(is.na(Description), Description.som,
	ifelse( Description.som=='No reported somatic small or SV/CN variants', Description, paste( Description, Description.som, sep="\n")))]

driver.summary[,Variant.mrg:=ifelse(is.na(Variant), Variant.som,
	ifelse( Variant.som=='N/A', Variant, paste(Variant, Variant.som, sep="\n")))]

driver.summary[,Variant.cmbd.mrg:=ifelse(is.na(Variant), Variant.cmbd.som,
	ifelse( Variant.cmbd.som=='N/A', Variant, paste(Variant, Variant.cmbd.som, sep="\n")))]

message( paste("Loading CHORD results for ", study.name, "\n", sep=""))
load( file = paste(run.dir, "/result/CHORD/output/chord_pred.Rdata", sep=""))
chord_output.dt <- as.data.table(chord_output)

message( paste("Merging PURPLE QC results with problem sample info for ", study.name, "\n", sep=""))
message( paste("Loading problem sample info for ", study.name, "\n", sep=""))
load( file = paste(run.dir, "/config_files/problem.sample.info.Rdata", sep=""))

qc.int.cols <- c("CopyNumberSegments", "DeletedGenes", "UnsupportedCopyNumberSegments", "AmberMeanDepth")
qc.int.cols <- qc.int.cols[which(qc.int.cols%in%colnames(purple.qc.dt))]
qc.num.cols <- c("Contamination", "Purity")
qc.num.cols <- qc.num.cols[which(qc.num.cols%in%colnames(purple.qc.dt))]

for (col in qc.int.cols)
	set(purple.qc.dt, j = col, value = as.integer(purple.qc.dt[[col]]))

for (col in qc.num.cols)
	set(purple.qc.dt, j = col, value = as.numeric(purple.qc.dt[[col]]))

if( nrow(problem.sample.info.dt)==0 ){
	purple.qc.dt.combined <- copy(purple.qc.dt)
}else{
	purple.qc.dt.combined <- merge(
		x = problem.sample.info.dt,
		y = purple.qc.dt,
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index'),
		all = TRUE)
	
	purple.qc.dt.combined[,QCStatus:=ifelse(is.na(QCStatus), problem, QCStatus)]
}

purple.qc.dt.combined <- merge(
	x = purple.qc.dt.combined,
	y = chord_output.dt[,list(tumor.sample.id=sample, hr_status, hrd_type)],
	by = c('tumor.sample.id'),
	all.x = TRUE)

if ('AmberMeanDepth'%in%colnames(purple.qc.dt.combined)){
	purple.qc.dt.combined <- purple.qc.dt.combined[order(tumor.index),list(
		subject.id, tumor.sample.id, AmberGender, CobaltGender, Contamination, CopyNumberSegments, DeletedGenes,
		GermlineAberrations, Method, Purity, QCStatus, UnsupportedCopyNumberSegments, AmberMeanDepth, hr_status, hrd_type)]
}else{
	purple.qc.dt.combined <- purple.qc.dt.combined[order(tumor.index),list(
		subject.id, tumor.sample.id, AmberGender, CobaltGender, Contamination, CopyNumberSegments, DeletedGenes,
		GermlineAberrations, Method, Purity, QCStatus, UnsupportedCopyNumberSegments, hr_status, hrd_type)]
}

purple.qc.dt.failed.contamination <- purple.qc.dt[grep("FAIL_CONTAMINATION", QCStatus, fixed=TRUE),]

purple.purity.dt <- purple.purity.dt[order(tumor.index),
	list(subject.id, tumor.sample.id,
		purity, normFactor, score, diploidProportion, ploidy, gender, status, polyclonalProportion,
		minPurity, maxPurity, minPloidy, maxPloidy, minDiploidProportion, maxDiploidProportion, version,
		somaticPenalty, wholeGenomeDuplication, msIndelsPerMb, msStatus, tml, tmlStatus, tmbPerMb, tmbStatus, svTumorMutationalBurden)]

message( paste("Annotating driver gene summary with CHORD analysis results for ", study.name, "\n", sep=""))

sample.info.driver.summary <- merge(
	x = tumor.normal.pair.info[,list(subject.index, tumor.index, subject.id, tumor.sample.id=sample.id.T)],
	y = chord_output.dt[,list(tumor.sample.id=sample, hr_status, hrd_type)],
	by = c('tumor.sample.id'),
	all.x = TRUE)

sample.info.driver.summary <- merge(
	x = sample.info.driver.summary[,list(subject.index, tumor.index, subject.id, tumor.sample.id, hr_status, hrd_type)],
	y = driver.summary[,list(subject.index, tumor.index, subject.id, tumor.sample.id, chr, CHROM, gene,
		Description.mrg, Variant.mrg, Variant.SV_CN.som, Variant.cmbd.mrg, HRD.gene.hits)],
	by = c('subject.index', 'tumor.index', 'subject.id', 'tumor.sample.id'),
	all.x = TRUE)

sample.info.driver.summary[,Description.mrg:=ifelse(is.na(Description.mrg), "No driver variant/mutations identified by PURPLE", Description.mrg)]

sample.info.driver.summary[tumor.sample.id%in%purple.qc.dt.failed.contamination$tumor.sample.id, hr_status:="Not tested: contamination failure"]
sample.info.driver.summary[tumor.sample.id%in%tumor.normal.pair.info[sample.id.N=="none",]$sample.id.T, hr_status:="Not tested: no matching normal"]


message( paste("Exporting driver gene variant summary table for ", study.name, "\n", sep=""))

fwrite(sample.info.driver.summary[,list(subject.id, tumor.sample.id,
	MutationDescription=Description.mrg, VariantDescription=Variant.mrg, HRD_status=hr_status, HRD_type=hrd_type, HRD.gene.hits)],
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.variant.description.txt", sep=""),
	col.names = TRUE,
	sep = '\t')

message( paste("Saving Rdata file of summarized and detailed driver gene information for ", study.name, "\n", sep=""))

save(list=c('germline.variant.info.dt.rdcd', 'somatic.variant.info.dt.rdcd',
	'germline.variant.funcClass2.summary', 'germline.variant.info.dt.candidates',
	'somatic.variant.funcClass2.summary', 'somatic.variant.info.dt.candidates',
	'germline.driver.variant.info.dt', 'somatic.driver.variant.info.dt',
	'merged.driver.variant.info.dt', 'merged.drivers.dt.filt', 'driver.purple.cnv.gene.info',
	'summarized.driver.variant.info', 'summarized.driver.variant.info.summary',
	'driver.summary', 'sample.info.driver.summary',
	'purple.qc.dt.combined', 'chord_output.dt', 'purple.purity.dt',
	'fusions.all.dt', 'fusions.reported.dt',
	'merged.SV.CN.info.dt', 'linx.SV.CN.info.dt', 'SV_CN.func.summary'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))
