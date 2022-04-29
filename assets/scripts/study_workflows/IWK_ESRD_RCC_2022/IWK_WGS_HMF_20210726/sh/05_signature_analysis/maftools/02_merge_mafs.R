#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=350)
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

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "maftools", "figures")

dir.create(file.path( run.dir, "result", "sigminer"))
dir.create(file.path( run.dir, "result", "maftools"))
if ( file.exists(fig.dir, recursive=FALSE)==FALSE ){
	dir.create(fig.dir)
}

load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))

load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.info.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

min.segment.length <- 0
#list=c('cn.segments', 'cn.segments.resegmented.dt'),
#load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.ver_20220131.Rdata', sep="") ) )
# Renamed dated file to not have to change downstream scripts
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )

chr3p.start <- 1L
chr3p.end <- 87000000L

cn.segments.resegmented.dt.chr3 <- cn.segments.resegmented.dt[chromosome=='chr3',]
cn.segments.resegmented.dt.chr3[,chr3p.length:=as.integer(ifelse(startpos<chr3p.end, ifelse(endpos<=chr3p.end, endpos-startpos+1L, chr3p.end-startpos+1L), 0))]
chr3p.loss.summary <- cn.segments.resegmented.dt.chr3[,list(chr3p.loss.length=sum(ifelse(nMinor==0, chr3p.length, 0L))), by=list(sample)]
chr3p.loss.summary[,chr3p.pct.loss:=chr3p.loss.length/chr3p.end]
chr3p.loss.summary[,chr3p.loss:=ifelse(chr3p.pct.loss>0.7, 'chr3p loss', 'No chr3p loss')]

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]

if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined <- purple.qc.dt.combined[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined <- rbind(
		purple.qc.dt.combined[QCStatus%in%c('PASS', 'WARN_HIGH_COPY_NUMBER_NOISE'),],
		purple.qc.dt.combined[tumor.depth>=tumor.depth.cutoff,][grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),])

purple.qc.dt.combined <- purple.qc.dt.combined[!tumor.sample.id%in%extra.samples.to.exclude,]

vep.sample.info.dt <- fread( file.path(run.dir, "result_summaries", paste(gpl_prefix, "", sep=""), "vep_annotation_sample_info.tsv"))

vep.tumor.sample.ids.ls <- vep.sample.info.dt[tumor.sample.id%in%purple.qc.dt.combined$tumor.sample.id,]$tumor.sample.id
names(vep.tumor.sample.ids.ls) <- vep.tumor.sample.ids.ls

length(vep.tumor.sample.ids.ls)
# 33

GPL.driver.summary <- merged.drivers.dt.filt[tumor.sample.id%in%vep.tumor.sample.ids.ls,
	list(driver.ct=.N), by=list(gene, driver, biallelic, eventType)]

GPL.driver.summary <- GPL.driver.summary[order(-driver.ct),]

fwrite(GPL.driver.summary, file=file.path(run.dir, "result", "maftools", "GPL_driver_cts.tsv"), sep='\t')

clinical.data <- merge(
	x = purple.purity.dt[,list(Tumor_Sample_Barcode=tumor.sample.id, gender, purity, ploidy)],
	y = patient.info.dt[,list(Tumor_Sample_Barcode=tumor.sample.id,
			histology, histology.short, histology.long=histology.description, age, Gender=gender, years.dialysis, primary.illness)],
	by = c('Tumor_Sample_Barcode'),
	all.x = TRUE)

clinical.data[,clear.cell:=ifelse(histology=='Clear cell', TRUE, FALSE)]
clinical.data[,chromophobe:=ifelse(histology=='Chromophobe', TRUE, FALSE)]
clinical.data[,ACD:=ifelse(histology=='ACD-RCC', TRUE, FALSE)]
clinical.data[,clear.papillary:=ifelse(histology=='Clear cell papillary', TRUE, FALSE)]

clinical.data[,exclude.sample:=FALSE]
clinical.data[!Tumor_Sample_Barcode%in%vep.tumor.sample.ids.ls,exclude.sample:=TRUE]

clinical.data[,MALE:=ifelse(gender=='MALE', 1L, 0L)]
clinical.data[,FEMALE:=ifelse(gender=='FEMALE', 1L, 0L)]


clinical.data <- merge(
	x = clinical.data,
	y = chr3p.loss.summary[,list(Tumor_Sample_Barcode=sample, chr3p.loss)],
	by = c('Tumor_Sample_Barcode'),
	all.x = TRUE)

clinical.summary.excluded <- clinical.data[exclude.sample==TRUE,
	list(ct.excluded=.N,
		males.excluded = sum(MALE),
		females.excluded = sum(FEMALE)),
	by=list(histology, histology.short)]

clinical.summary <- clinical.data[,
	list(ct=.N,
		males = sum(MALE),
		females = sum(FEMALE),
		purity.mean = mean(purity), purity.min = min(purity), purity.max = max(purity),
		years.dialysis.mean=mean(years.dialysis), years.dialysis.min=min(years.dialysis), years.dialysis.max=max(years.dialysis)),
by=list(histology, histology.short)]

clinical.summary <- merge(
	x = clinical.summary,
	y = clinical.summary.excluded,
	by = c('histology', 'histology.short'),
	all = TRUE)

clinical.summary[,ct:=ifelse(is.na(ct), 0L, ct)]
clinical.summary[,males:=ifelse(is.na(males), 0L, males)]
clinical.summary[,females:=ifelse(is.na(females), 0L, females)]


clinical.summary[,ct.excluded:=ifelse(is.na(ct.excluded), 0L, ct.excluded)]
clinical.summary[,males.excluded:=ifelse(is.na(males.excluded), 0L, males.excluded)]
clinical.summary[,females.excluded:=ifelse(is.na(females.excluded), 0L, females.excluded)]

clinical.summary[,total.str:=paste(ct, " (", ct.excluded, ")", sep="")]
clinical.summary[,males.str:=paste(males, " (", males.excluded, ")", sep="")]
clinical.summary[,females.str:=paste(females, " (", females.excluded, ")", sep="")]

clinical.summary[,purity := paste(
	format(purity.mean, digits=2, nsmall=2, justify='none'),
	' [',
	format(purity.min, digits=2, nsmall=2, justify='none'), '-',
	format(purity.max, digits=2, nsmall=2,  justify='none'), ']', sep=""),
	by=1:nrow(clinical.summary)]
clinical.summary[,years.dialysis := paste(format(years.dialysis.mean, digits=2, nsmall=2, justify='none'),
	' [', 
	format(years.dialysis.min, digits=2, nsmall=2, justify='none'), '-', 
	format(years.dialysis.max, digits=2, nsmall=2, justify='none'),
	']', sep=""), by=1:nrow(clinical.summary)]


clinical.summary.excluded.total <- clinical.data[exclude.sample==TRUE,
		list(histology='Total', ct.excluded=.N,
				males.excluded = sum(MALE),
				females.excluded = sum(FEMALE))]

clinical.summary.total <- clinical.data[,
		list(histology='Total', 
			ct=.N,
				males = sum(MALE),
				females = sum(FEMALE),
				purity.mean = mean(purity), purity.min = min(purity), purity.max = max(purity),
				years.dialysis.mean=mean(years.dialysis), years.dialysis.min=min(years.dialysis), years.dialysis.max=max(years.dialysis))]

clinical.summary.total <- merge(
		x = clinical.summary.total,
		y = clinical.summary.excluded.total,
		by = 'histology',
		all = TRUE)

clinical.summary.total[,ct.excluded:=ifelse(is.na(ct.excluded), 0L, ct.excluded)]
clinical.summary.total[,males.excluded:=ifelse(is.na(males.excluded), 0L, males.excluded)]
clinical.summary.total[,females.excluded:=ifelse(is.na(females.excluded), 0L, females.excluded)]


clinical.summary.total[,total.str:=paste(ct, " (", ct.excluded, ")", sep="")]
clinical.summary.total[,males.str:=paste(males, " (", males.excluded, ")", sep="")]
clinical.summary.total[,females.str:=paste(females, " (", females.excluded, ")", sep="")]

clinical.summary.total[,purity := paste(
				format(purity.mean, digits=2, nsmall=2, justify='none'),
				' [',
				format(purity.min, digits=2, nsmall=2, justify='none'), '-',
				format(purity.max, digits=2, nsmall=2,  justify='none'), ']', sep=""),
		by=1:nrow(clinical.summary.total)]
clinical.summary.total[,years.dialysis := paste(format(years.dialysis.mean, digits=2, nsmall=2, justify='none'),
				' [', 
				format(years.dialysis.min, digits=2, nsmall=2, justify='none'), '-', 
				format(years.dialysis.max, digits=2, nsmall=2, justify='none'),
				']', sep=""), by=1:nrow(clinical.summary.total)]

clinical.summary <- clinical.summary[order(-ct),]

clinical.summary.for.export <- rbind(
	clinical.summary[,list(Histology=paste(histology, ' (', histology.short, ')', sep=""), Total=total.str, Males=males.str, Females=females.str, Purity=purity, years.dialysis)],
	clinical.summary.total[,list(Histology=histology, Total=total.str, Males=males.str, Females=females.str, Purity=purity, years.dialysis)])


fwrite( clinical.summary.for.export,
	file = file.path( run.dir, "result", "maftools", "clinical_summary.tsv"), sep='\t')

fwrite( clinical.summary.excluded,
	file = file.path( run.dir, "result", "maftools", "clinical_summary.excluded.tsv"), sep='\t')

save(list=c('clinical.data', 'clinical.summary', 'clinical.summary.excluded', 'clinical.summary.for.export', 'chr3p.loss.summary'),
		file = file.path( run.dir, "result", "maftools", "clinical_data.Rdata" ))

clinical.summary.for.export
##                     Histology  Total  Males Females           Purity       years.dialysis
## 1:            Clear cell (ccRCC) 18 (0) 15 (0)   3 (0) 0.56 [0.34-0.74]    8.51 [0.29-28.10]
## 2:             ACD-RCC (ACD-RCC)  9 (3)  6 (2)   3 (1) 0.43 [0.08-0.82]   19.36 [6.43-39.48]
## 3:              Papillary (pRCC)  5 (0)  2 (0)   3 (0) 0.52 [0.25-0.73] 13.79 [0.0083-27.63]
## 4:           Chromophobe (chRCC)  3 (1)  2 (1)   1 (0) 0.69 [0.08-1.00]   11.62 [0.49-27.50]
## 5: Clear cell papillary (ccpRCC)  2 (0)  2 (0)   0 (0) 0.50 [0.48-0.52]   11.29 [6.05-16.52]
## 6:         Unclassified (uncRCC)  1 (1)  1 (1)   0 (0) 0.19 [0.19-0.19]  14.21 [14.21-14.21]
## 7:                         Total 38 (5) 28 (4)  10 (1) 0.52 [0.08-1.00] 12.31 [0.0083-39.48]


clinical.summary.excluded
##         histology ct males females
## 1:      ACD-RCC  3     2       1
## 2:  Chromophobe  1     1       0
## 3: Unclassified  1     1       0


clinical.data
t.test(x=clinical.data[histology=='Clear cell',]$years.dialysis, y=clinical.data[histology!='Clear cell',]$years.dialysis )
## t.test(x=clinical.data[histology=='Clear cell',]$years.dialysis, y=clinical.data[histology!='Clear cell',]$years.dialysis )
## 
## Welch Two Sample t-test
## 
## data:  clinical.data[histology == "Clear cell", ]$years.dialysis and clinical.data[histology != "Clear cell", ]$years.dialysis
## t = -2.3832, df = 35.714, p-value = 0.0226
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##         -13.398151  -1.076818
## sample estimates:
##         mean of x mean of y 
## 8.50571  15.74319 

t.test(x=clinical.data[histology=='Clear cell',]$years.dialysis, y=clinical.data[histology=='ACD-RCC',]$years.dialysis )
## t.test(x=clinical.data[histology=='Clear cell',]$years.dialysis, y=clinical.data[histology=='ACD-RCC',]$years.dialysis )
## 
## Welch Two Sample t-test
## 
## data:  clinical.data[histology == "Clear cell", ]$years.dialysis and clinical.data[histology == "ACD-RCC", ]$years.dialysis
## t = -2.6382, df = 13.03, p-value = 0.02043
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##         -19.742464  -1.968339
## sample estimates:
##         mean of x mean of y 
## 8.50571  19.36111

clinical.data <- clinical.data[exclude.sample==FALSE,]

purple.purity.dt <- purple.purity.dt[tumor.sample.id%in%vep.tumor.sample.ids.ls,]

purple.cnv.gene.info <- merge(
	x = purple.purity.dt[,list(subject.id, tumor.sample.id,  purity, ploidy, gender)],
	y = purple.cnv.gene.info,
	by = c('subject.id', 'tumor.sample.id'))

purple.cnv.gene.info <- merge(
	x = purple.cnv.gene.info,
	y = clinical.data[,list(tumor.sample.id=Tumor_Sample_Barcode, histology=histology.short, years.dialysis)],
	by = c('tumor.sample.id'))

autosomes.ls <- paste('chr', 1:22, sep="")

## if cn_value < 0.5:
##             status = 'homozygous_deletion'
## 
## elif 0.5 <= cn_value < (0.6 * ploidy):
##         status = 'heterozygous_deletion'
## 
## elif (0.6 * ploidy) <= cn_value < (1.4 * ploidy):
##         status = 'neutral'
## 
## elif (1.4 * ploidy) <= cn_value < (3 * ploidy):
##         status = 'gain'
## 
## else:
##             status = 'gain_high'
## return status`

purple.cnv.gene.info.filt <- purple.cnv.gene.info[chromosome!='chrY',]

gene.cn.info <- purple.cnv.gene.info.filt[,list(chromosome, gene, sample=tumor.sample.id, ploidy,
	histology, years.dialysis,
	minCopyNumber, maxCopyNumber, minMinorAlleleCopyNumber, minRegionMethod,
	cn.status=ifelse(chromosome%in%autosomes.ls | (chromosome=='chrX' & gender=='FEMALE'),
		ifelse(minCopyNumber<0.5, 'DeepDel',
			ifelse( 0.5<=minCopyNumber & minCopyNumber<0.6*ploidy, 'Del',
				ifelse( 0.6*ploidy<=minCopyNumber & minCopyNumber<1.4*ploidy, ifelse(minMinorAlleleCopyNumber*ploidy<0.1, 'Del', 'neutral'),
					ifelse( 1.4*ploidy<=minCopyNumber & minCopyNumber<3*ploidy, 'ShallowAmp',
						ifelse( 3*ploidy<=minCopyNumber, 'Amp', 'Auto. unclass.'))))),
		ifelse(minCopyNumber<0.5*(ploidy/2), 'DeepDel',
			ifelse( 0.5*(ploidy/2)<=minCopyNumber & minCopyNumber<0.6*(ploidy/2), 'Del',
				ifelse( 0.6*(ploidy/2)<=minCopyNumber & minCopyNumber<1.4*(ploidy/2), 'neutral',
					ifelse( 1.4*(ploidy/2)<=minCopyNumber & minCopyNumber<3*(ploidy/2), 'ShallowAmp',
						ifelse( 3*(ploidy/2)<=minCopyNumber, 'Amp', 'chrX unclass.')))))))]

#	cn.status=ifelse(minCopyNumber<0.5, "Del", ifelse( minCopyNumber>3*ploidy, "Amp", "none")))]
table(gene.cn.info$cn.status)
##	Amp    DeepDel        Del    neutral ShallowAmp 
#	7        394      55938     727231      30870
#gene.cn.info <- gene.cn.info[cn.status!='neutral',]
gene.cn.info <- gene.cn.info[cn.status%in%c('Amp', 'DeepDel'),]

SV.samples <- c('IWK020_T', 'IWK022_T')

get_file_path <- function( curr.tumor.id, variant.type ){
	curr.subject.id <- vep.sample.info.dt[tumor.sample.id==curr.tumor.id,]$subject.id
	if ( variant.type=='SVs' ){
		file.path( run.dir, "result", gpl_prefix, curr.subject.id, "purple", curr.tumor.id, paste(curr.tumor.id, ".linx.sv.vep.maf", sep=""))
	}else{
		file.path( run.dir, "result", gpl_prefix, curr.subject.id, "purple", curr.tumor.id, paste(curr.tumor.id, ".purple.", variant.type, ".vep.maf", sep=""))
	#	file.path( run.dir, "result", gpl_prefix, curr.subject.id, "purple", curr.tumor.id, paste(curr.tumor.id, ".purple.", variant.type, ".vep.maf", sep=""))
	}
}
somatic_maf_files <- unlist(lapply(X=vep.tumor.sample.ids.ls, FUN=get_file_path, variant.type='somatic'))
germline_maf_files <- unlist(lapply(X=vep.tumor.sample.ids.ls, FUN=get_file_path, variant.type='germline'))
SV_maf_files <- unlist(lapply(X=vep.tumor.sample.ids.ls[which(vep.tumor.sample.ids.ls%in%SV.samples)], FUN=get_file_path, variant.type='SVs'))

## do not need to run if already know that now germline variants were reported
#merged_germline_mafs <- merge_mafs(
#	mafs = germline_maf_files, verbose = TRUE)

#print(merged_germline_mafs)

#merged_germline_mafs.filt <- subsetMaf(maf = merged_germline_mafs, query="REPORTED == 1")
#no germline variants

clinical.data[,Histology:=histology]
setnames(clinical.data, old='chr3p.loss', new='chr3p_loss')

merged_somatic_mafs <- merge_mafs(
	mafs = c(somatic_maf_files, SV_maf_files),
	clinicalData = clinical.data[exclude.sample==FALSE,],
	verbose = TRUE)

merged_mafs.with.cn <- merge_mafs(
#	mafs = c(merged_germline_mafs.filt, merged_somatic_mafs),
	mafs = c(merged_somatic_mafs),
	cnTable = gene.cn.info[,list(Gene=gene, Sample_name=sample, CN=cn.status)],
	verbose = TRUE)


##GISTIC.AllSamples.dir <- file.path(run.dir, "result/GISTIC2/20211221_0.1_0.1_AllSamples")

GISTIC.AllSamples.dir <- file.path(run.dir, "result/GISTIC2/20220205_none_noChrX_brlen_0.98_0.3_0.3_AllSamples")
#GISTIC.AllSamples.dir <- file.path(run.dir, "result/GISTIC2/20220205_none_noChrX_brlen_0.98_0.3_0.3_non-ccRCC")

all.lesions <- file.path( GISTIC.AllSamples.dir, "all_lesions.conf_75.txt")
amp.genes <- file.path( GISTIC.AllSamples.dir, "amp_genes.conf_75.txt")
del.genes <- file.path( GISTIC.AllSamples.dir, "del_genes.conf_75.txt")
scores.gis <- file.path( GISTIC.AllSamples.dir, "scores.gistic")

merged.gistic =  readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)

merged_mafs.with.GISTIC2 <- merge_mafs(
	mafs = c(merged_somatic_mafs),
	gisticAllLesionsFile = all.lesions,
	gisticAmpGenesFile = amp.genes,
	gisticDelGenesFile = del.genes,
	gisticScoresFile = scores.gis,
	verbose = TRUE)

save( list=c('merged.gistic', 'merged_somatic_mafs', 'merged_mafs.with.GISTIC2', 'merged_mafs.with.cn', 'GPL.driver.summary'),
	file = file.path( run.dir, "result", "maftools", "merged_mafs_GISTIC2_AllSamples.Rdata" ))
#	file = file.path( run.dir, "result", "maftools", "merged_mafs_GISTIC2_non-ccRCC.Rdata" ))

print( merged_mafs.with.GISTIC2)
#ID summary         Mean Median
#1:             NCBI_Build  GRCh38           NA     NA
#2:                 Center       .           NA     NA
#3:                Samples      33           NA     NA
#4:                 nGenes    2732           NA     NA
#5:        Frame_Shift_Del     130   3.93939394      3
#6:        Frame_Shift_Ins      40   1.21212121      1
#7:           In_Frame_Del      33   1.00000000      1
#8:           In_Frame_Ins       1   0.03030303      0
#9:      Missense_Mutation    1092  33.09090909     32
#10:      Nonsense_Mutation      82   2.48484848      2
#11:       Nonstop_Mutation       2   0.06060606      0
#12:            Splice_Site      54   1.63636364      1
#13: Translation_Start_Site       2   0.06060606      0
#14:                  total    1436  43.51515152     42
#15:                    Amp     453  13.72727273      1
#16:                    Del   12458 377.51515152    438
#17:              CNV_total   12911 391.24242424    451
# before adding SVs

## An object of class MAF 
## ======================
##         ID summary         Mean Median
## 1:             NCBI_Build  GRCh38           NA     NA
## 2:                 Center       .           NA     NA
## 3:                Samples      33           NA     NA
## 4:                 nGenes    2908           NA     NA
## 5:        Frame_Shift_Del     130   3.93939394      3
## 6:        Frame_Shift_Ins      40   1.21212121      1
## 7:           In_Frame_Del      33   1.00000000      1
## 8:           In_Frame_Ins       1   0.03030303      0
## 9:      Missense_Mutation    1092  33.09090909     32
## 10:      Nonsense_Mutation      82   2.48484848      2
## 11:       Nonstop_Mutation       2   0.06060606      0
## 12:            Splice_Site      50   1.51515152      1
## 13: Translation_Start_Site       2   0.06060606      0
## 14:                  total    1432  43.39393939     42
## 15:                    Amp    1760  53.33333333      4
## 16:                    Del   14656 444.12121212    284
## 17:              CNV_total   16416 497.45454545    457
# before resegmentation
#An object of class MAF 
#======================
#		ID summary         Mean Median
#1:             NCBI_Build  GRCh38           NA     NA
#2:                 Center       .           NA     NA
#3:                Samples      33           NA     NA
#4:                 nGenes    2681           NA     NA
#5:        Frame_Shift_Del     130   3.93939394      3
#6:        Frame_Shift_Ins      40   1.21212121      1
#7:           In_Frame_Del      33   1.00000000      1
#8:           In_Frame_Ins       1   0.03030303      0
#9:      Missense_Mutation    1092  33.09090909     32
#10:      Nonsense_Mutation      82   2.48484848      2
#11:       Nonstop_Mutation       2   0.06060606      0
#12:            Splice_Site      50   1.51515152      1
#13: Translation_Start_Site       2   0.06060606      0
#14:                  total    1432  43.39393939     42
#15:                    Amp     532  16.12121212      3
#16:                    Del   15690 475.45454545    567
#17:              CNV_total   16222 491.57575758    568

print(merged_mafs.with.cn)
## An object of class MAF 
## ======================
##         ID summary        Mean Median
## 1:             NCBI_Build  GRCh38          NA     NA
## 2:                 Center       .          NA     NA
## 3:                Samples      33          NA     NA
## 4:                 nGenes    1417          NA     NA
## 5:                DeepDel     394 11.93939394     10
## 6:        Frame_Shift_Del     130  3.93939394      3
## 7:        Frame_Shift_Ins      40  1.21212121      1
## 8:           In_Frame_Del      33  1.00000000      1
## 9:           In_Frame_Ins       1  0.03030303      0
## 10:      Missense_Mutation    1092 33.09090909     32
## 11:      Nonsense_Mutation      82  2.48484848      2
## 12:       Nonstop_Mutation       2  0.06060606      0
## 13:            Splice_Site      50  1.51515152      1
## 14: Translation_Start_Site       2  0.06060606      0
## 15:                  total    1826 55.33333333     57
## 16:                    Amp       7  0.21212121      0
## 17:              CNV_total       7  0.21212121      0


histology.cts <- table(clinical.data$histology)
histology.cts <- histology.cts[order(-histology.cts)]
histcolors = RColorBrewer::brewer.pal(n = length(histology.cts),name = 'Spectral')
names(histcolors) = names(histology.cts)
histcolors  = list(histology = histcolors )

histology.names <- names(histology.cts)
names(histology.names) <- histology.names

#load(file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/merged.refit.activities.Rdata", sep="")))

library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

merged_mafs.with.GISTIC2.tnm = trinucleotideMatrix(maf = merged_mafs.with.GISTIC2,
	ref_genome = "BSgenome.Hsapiens.UCSC.hg38")


library('NMF')
merged_mafs.10.sign = estimateSignatures(mat = merged_mafs.with.GISTIC2.tnm, nMin=2, nTry = 10, nrun=10, parallel=12)

pdf( file = file.path(fig.dir, paste(study.name, "_nmf.10trys.cophenetic.pdf", sep="")))
plotCophenetic(res = merged_mafs.10.sign)
dev.off()

merged_mafs.sig = extractSignatures(mat = merged_mafs.with.GISTIC2.tnm, n = 5)

save( list=c('merged_mafs.with.GISTIC2.tnm', 'merged_mafs.sig'),
		file = file.path( run.dir, "result", "maftools", "merged_tnm.Rdata" ))

#Compate against original 30 signatures 
merged_mags.og30.cosm = compareSignatures(nmfRes = merged_mafs.sig, sig_db = "SBS")
#-Comparing against COSMIC signatures
#------------------------------------
#		--Found Signature_1 most similar to SBS40
#Aetiology: Unknown [cosine-similarity: 0.838]
#--Found Signature_2 most similar to SBS5
#Aetiology: Unknown [cosine-similarity: 0.864]
#--Found Signature_3 most similar to SBS40
#Aetiology: Unknown [cosine-similarity: 0.923]
#--Found Signature_4 most similar to SBS12
#Aetiology: Unknown [cosine-similarity: 0.93]
#--Found Signature_5 most similar to SBS5
#Aetiology: Unknown [cosine-similarity: 0.843]
#------------------------------------

top.mutated.genes <- unique(GPL.driver.summary[driver.ct>1,]$gene)
print(histcolors)
## Clear cell              ACD-RCC            Papillary          Chromophobe Clear cell papillary 
## "#D7191C"            "#FDAE61"            "#FFFFBF"            "#ABDDA4"            "#2B83BA" 

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_apobec.pdf", sep="")))
plotApobecDiff(tnm = merged_mafs.with.GISTIC2.tnm , maf = merged_mafs.with.GISTIC2, pVal = 0.2)
dev.off()

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_GISTIC2.chromPlot.pdf", sep="")))
gisticChromPlot(gistic = merged.gistic, ref.build = "hg38", markBands = "all")
dev.off()

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_maf_summary.pdf", sep="")))
plotmafSummary(maf = merged_mafs.with.GISTIC2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_titv.pdf", sep="")))
merged_somatic.titv = titv(maf = merged_somatic_mafs, plot = FALSE, useSyn = TRUE)
plotTiTv(res = merged_somatic.titv)
dev.off()

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_top_mutated_genes_lollipop.pdf", sep="")))

lapply(
	X = top.mutated.genes,
	FUN = function( curr.gene ){
		lollipopPlot(
				maf = merged_mafs.with.GISTIC2,
				gene = curr.gene,
#				AACol = 'Protein_Change',
				showMutationRate = TRUE,
		)
	})

dev.off()



pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_no_CN_maf_oncoplot.pdf", sep="")))
oncoplot(maf=merged_somatic_mafs, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

## oncoplots with default gene selection
pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_oncoplot.pdf", sep="")))
oncoplot(maf=merged_mafs.with.GISTIC2, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

catch.results.ls <- lapply(
	X = histology.names,
	FUN = function( curr.histology ){
		pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_oncoplot_", curr.histology, ".pdf", sep="")))
		oncoplot(maf=subsetMaf(maf=merged_mafs.with.GISTIC2,
			clinQuery=paste("histology == '", curr.histology, "'", sep="")),
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
		dev.off()
	})


pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_histology.pdf", sep="")))
oncoplot(
		maf = merged_mafs.with.GISTIC2,
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_histology.top_altered.pdf", sep="")))
oncoplot(
		maf = merged_mafs.with.GISTIC2,
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2,
		altered = TRUE)
dev.off()

pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_CN_calls.oncoplot_with_histology.pdf", sep="")))
oncoplot(
		maf = merged_mafs.with.cn,
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=7, height=10, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_pathways_auto.pdf", sep="")))
oncoplot(maf=merged_mafs.with.GISTIC2, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, pathways = 'auto', legend_height = 2)
dev.off()

catch.results.ls <- lapply(
	X = histology.names,
	FUN = function( curr.histology ){
		pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_pathways_", curr.histology, ".pdf", sep="")))
		oncoplot(maf=subsetMaf(maf=merged_mafs.with.GISTIC2,
			clinQuery=paste("histology == '", curr.histology, "'", sep="")),
			showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, pathways = 'auto', legend_height = 2)
		dev.off()
	})

pdf( width=10, height=7, file = file.path(fig.dir, paste(study.name, "_Vaf.pdf", sep="")))
plotVaf(maf=merged_mafs.with.GISTIC2)
dev.off()



# Oncoplots with GPL driver genes
pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_GPL_drivers.oncoplot.pdf", sep="")))
oncoplot(maf=merged_mafs.with.GISTIC2, genes = unique(GPL.driver.summary$gene), showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_GPL_drivers.oncoplot_with_histology.pdf", sep="")))
oncoplot(
		maf = merged_mafs.with.GISTIC2,
		genes = unique(GPL.driver.summary$gene), 
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE,
		gene_mar = 6,
		barcode_mar = 5,
		legend_height = 2)
dev.off()

GPL.driver.summary[,gene.mod:=ifelse(gene=='UGT1A8', 'UGT1A1', gene)]

chr3p.loss.genes.ls <- c('VHL', 'PBRM1', 'SETD2', 'BAP1')
GPL.driver.summary2 <- GPL.driver.summary[gene.mod!='OR11H1',list(driver.ct=sum(driver.ct)), by=list(gene, gene.mod)]
GPL.driver.summary2[,chr3p.loss.gene:=ifelse(gene.mod%in%chr3p.loss.genes.ls, 1L, 0L)]
GPL.genes.ls <- GPL.driver.summary2[order(-chr3p.loss.gene, -driver.ct),]$gene.mod

chr3p.loss.colors = c('black', 'lightgrey')
names(chr3p.loss.colors) = c('chr3p loss', 'No chr3p loss')

annot_colors = list(Histology = histcolors$histology, chr3p_loss = chr3p.loss.colors)

pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_GPL_drivers.oncoplot_with_histology_and_chr3ploss.pdf", sep="")))
oncoplot(
		maf = merged_somatic_mafs,
#		genes = unique(GPL.driver.summary$gene), 
		genes = GPL.genes.ls, 
		clinicalFeatures = c("Histology", "chr3p_loss"),
		sortByAnnotation = TRUE,
		annotationColor = annot_colors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2,
		includeColBarCN = FALSE)
dev.off()

pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_GPL_drivers.oncoplot_with_histology_and_chr3ploss.allsamples.pdf", sep="")))
oncoplot(
		maf = merged_somatic_mafs,
#		genes = unique(GPL.driver.summary$gene), 
		genes = GPL.genes.ls, 
		clinicalFeatures = c("Histology", "chr3p_loss"),
		sortByAnnotation = TRUE,
		annotationColor = annot_colors,
		showTumorSampleBarcodes = TRUE,
		keepGeneOrder = TRUE,
		gene_mar = 6,
		barcode_mar = 5,
		fontSize = 0.6,
		legend_height = 2,
		titleFontSize = 1.0,
		annotationFontSize = 1.0,
		legendFontSize = 1.0,
		includeColBarCN = FALSE,
		removeNonMutated = FALSE)
dev.off()


pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_GPL_drivers.CN_calls.oncoplot_with_histology.pdf", sep="")))
oncoplot(
		maf = merged_mafs.with.cn,
		genes = unique(GPL.driver.summary$gene), 
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()