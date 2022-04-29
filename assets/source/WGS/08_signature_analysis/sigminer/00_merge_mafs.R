#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=300)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

gpl_prefix <- "GRIDSS-2.12.0"

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

dir.create(file.path( run.dir, "result", "maftools"))
if ( file.exists(fig.dir, recursive=FALSE)==FALSE ){
	dir.create(fig.dir)
}

load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

vep.sample.info.dt <- fread( file.path(run.dir, "result_summaries", gpl_prefix, "vep_annotation_sample_info.tsv"))

vep.tumor.sample.ids.ls <- vep.sample.info.dt$tumor.sample.id
names(vep.tumor.sample.ids.ls) <- vep.tumor.sample.ids.ls

load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.Rdata", sep=""))

GPL.driver.summary <- merged.drivers.dt.filt[,list(driver.ct=.N), by=list(gene, driver, biallelic, eventType)]

GPL.driver.summary <- GPL.driver.summary[order(-driver.ct),]

fwrite(GPL.driver.summary, file=file.path(run.dir, "result", "maftools", "GPL_driver_cts.tsv"), sep='\t')


sex.dt <- purple.qc.dt.combined[,list(sample=tumor.sample.id, sex=tolower(AmberGender))]

save( sex.dt,
	file = file.path( run.dir, "result", "sigminer", "sample.sex.Rdata" ))

purple.cnv.segments <- purple.cnv.somatic.info[,list(
	tumor.sample.id, chromosome, start, end, seg.length=end-start, copyNumber, minorAlleleCopyNumber, majorAlleleCopyNumber, segmentStartSupport, segmentEndSupport, method, depthWindowCount, bafCount, observedBAF)]
purple.cnv.segments[,minorAlleleCopyNumber.int:=round(ifelse(minorAlleleCopyNumber<0, 0, minorAlleleCopyNumber), 0)]
purple.cnv.segments[,majorAlleleCopyNumber.int:=round(ifelse(majorAlleleCopyNumber<0, 0, majorAlleleCopyNumber), 0)]
purple.cnv.segments[,copyNumber.int:=round(ifelse(copyNumber<0, 0, copyNumber), 0)]

options(sigminer.sex = sex.dt)

merged.cn <- read_copynumber(
	purple.cnv.segments[seg.length>0,list(sample=tumor.sample.id, chromosome, start, end, segVal=log2(copyNumber.int+0.01), minor_cn=minorAlleleCopyNumber.int)],
	seg_cols = c("chromosome", "start", "end", "segVal"),
	genome_build = "hg38",
	genome_measure = "wg",
	join_adj_seg = TRUE,
	add_loh=TRUE,
	loh_min_length=0,
	complement = FALSE, verbose = TRUE
)

save( merged.cn,
	file = file.path( run.dir, "result", "sigminer", "merged_CN_W.Rdata" ))



merged.cn <- read_copynumber(
	purple.cnv.segments[seg.length>0,list(sample=tumor.sample.id, chromosome, start, end, segVal=copyNumber.int, minor_cn=minorAlleleCopyNumber.int)],
	seg_cols = c("chromosome", "start", "end", "segVal"),
	genome_build = "hg38",
	genome_measure = "wg",
	join_adj_seg = FALSE,
	add_loh=TRUE,
	loh_min_length=0,
	complement = FALSE, verbose = TRUE
)

save( merged.cn,
		file = file.path( run.dir, "result", "sigminer", "merged_CN_S.Rdata" ))

## read_sv_as_rs( linx.vis_sv_data.info.dt[,list(sample=tumor.sample.id,
## chr1=ChrStart, start1=PosStart, end1=PosEnd,
## chr2=ChrEnd, start2=, end2, strand1, strand2, svclass)]

## curr.segs <- purple.cnv.segments[tumor.sample.id=='IWK011_T',]
## curr.segs <- curr.segs[,index:=as.integer(1:nrow(curr.segs))]
## 
## curr.linx.sv <- linx.vis_sv_data.info.dt[tumor.sample.id=='IWK011_T',]

purple.cnv.gene.info <- merge(
	x = purple.purity.dt[,list(subject.id, tumor.sample.id,  purity, ploidy)],
	y = purple.cnv.gene.info,
	by = c('subject.id', 'tumor.sample.id'))

gene.cn.info <- purple.cnv.gene.info[,list(gene, sample=tumor.sample.id, cn.status=ifelse(minCopyNumber<0.5, "Del", ifelse( minCopyNumber>3*ploidy, "Amp", "none")))]
table(gene.cn.info$cn.status)
##	Amp    Del   none 
##	 7   5859 934634


clinical.data <- merge(
	x = purple.purity.dt[,list(Tumor_Sample_Barcode=tumor.sample.id, purity, ploidy)],
	y = patient.info.dt[,list(Tumor_Sample_Barcode=paste(subject.id, '_T', sep=""), histology, histology.description)],
	by = c('Tumor_Sample_Barcode'),
	all.x = TRUE)
clinical.data[,histology:=ifelse(is.na(histology), 'Unknown', histology)]
clinical.data[,histology.description:=ifelse(is.na(histology.description), 'Unknown', histology.description)]

get_file_path <- function( curr.tumor.id, variant.type ){
	curr.subject.id <- vep.sample.info.dt[tumor.sample.id==curr.tumor.id,]$subject.id
	file.path( run.dir, "result", gpl_prefix, curr.subject.id, "purple", curr.tumor.id, paste(curr.tumor.id, ".purple.", variant.type, ".vep.maf", sep=""))
}
somatic_maf_files <- unlist(lapply(X=vep.tumor.sample.ids.ls, FUN=get_file_path, variant.type='somatic'))
germline_maf_files <- unlist(lapply(X=vep.tumor.sample.ids.ls, FUN=get_file_path, variant.type='germline'))


merged_germline_mafs <- merge_mafs(
	mafs = germline_maf_files, verbose = TRUE)

## > merged_germline_mafs
## An object of class  MAF 
## ID summary    Mean Median
## 1:        NCBI_Build  GRCh38      NA     NA
## 2:            Center       .      NA     NA
## 3:           Samples      38      NA     NA
## 4:            nGenes     235      NA     NA
## 5:   Frame_Shift_Del       3   0.079    0.0
## 6:   Frame_Shift_Ins       1   0.026    0.0
## 7:      In_Frame_Del      63   1.658    1.5
## 8:      In_Frame_Ins      34   0.895    1.0
## 9: Missense_Mutation    7878 207.316  203.5
## 10: Nonsense_Mutation      44   1.158    1.0
## 11:       Splice_Site      13   0.342    0.0
## 12:             total    8036 211.474  207.5

merged_germline_mafs.filt <- subsetMaf(maf = merged_germline_mafs, query="REPORTED == 1")
#no germline variants

merged_somatic_mafs <- merge_mafs(
	mafs = somatic_maf_files,
	clinicalData = clinical.data,
	verbose = TRUE)

merged_mafs <- merge_mafs(
#	mafs = c(merged_germline_mafs.filt, merged_somatic_mafs),
	mafs = c(merged_somatic_mafs),
	cnTable = gene.cn.info[which(cn.status!='none'),],
	verbose = TRUE)

save( merged_mafs,
	file = file.path( run.dir, "result", "maftools", "merged_mafs.Rdata" ))

## > merged_mafs
## ID summary         Mean Median
## 1:             NCBI_Build  GRCh38           NA     NA
## 2:                 Center       .           NA     NA
## 3:                Samples      38           NA     NA
## 4:                 nGenes    4569           NA     NA
## 5:        Frame_Shift_Del     131   3.44736842      3
## 6:        Frame_Shift_Ins      40   1.05263158      1
## 7:           In_Frame_Del      34   0.89473684      1
## 8:           In_Frame_Ins       1   0.02631579      0
## 9:      Missense_Mutation    1104  29.05263158     30
## 10:      Nonsense_Mutation      82   2.15789474      2
## 11:       Nonstop_Mutation       2   0.05263158      0
## 12:            Splice_Site      50   1.31578947      1
## 13: Translation_Start_Site       2   0.05263158      0
## 14:                  total    1446  38.05263158     38
## 15:                    Amp       7   0.18421053      0
## 16:                    Del    5859 154.18421053     17
## 17:              CNV_total    5866 154.36842105     17

merged_mafs.MT <- subsetMaf(maf=merged_mafs, query="Chromosome == 'chrM'")

merged_mafs.nonMT <- filterMaf(maf=merged_mafs, gene = c(paste('MT-ND', 1:6, sep=""), paste('MT-ATP', 6:8, sep=""), paste('MT-CO', 1:2, sep=""), 'MT-CYB', 'MT-ND4L'))

saveRDS( merged_mafs.nonMT,
	file = file.path( run.dir, "result", "maftools", "merged_mafs_nonMT.RDS" ))

histology.cts <- table(clinical.data$histology)
histology.cts <- histology.cts[order(-histology.cts)]
histcolors = RColorBrewer::brewer.pal(n = length(histology.cts),name = 'Spectral')
names(histcolors) = names(histology.cts)
histcolors  = list(histology = histcolors )

print(fabcolors)
#> $FAB_classification
#>        M0        M1        M2        M3        M4        M5        M6        M7 
#> "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#E6F598" "#ABDDA4" "#66C2A5" "#3288BD"

pdf( width=10, height=8, file = file.path(fig.dir, paste(study.name, "_maf_summary.pdf", sep="")))
plotmafSummary(maf = merged_mafs.nonMT, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_oncoplot_noFilter.pdf", sep="")))
oncoplot(maf=merged_mafs, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=7, height=8, file = file.path(fig.dir, paste(study.name, "_oncoplot.pdf", sep="")))
oncoplot(maf=merged_mafs.nonMT, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=8, height=10, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_histology.pdf", sep="")))
oncoplot(
		maf = merged_mafs.nonMT,
		clinicalFeatures = "histology",
		sortByAnnotation = TRUE,
		annotationColor = histcolors,
		showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, legend_height = 2)
dev.off()

pdf( width=7, height=10, file = file.path(fig.dir, paste(study.name, "_oncoplot_with_pathways_auto.pdf", sep="")))
oncoplot(maf=merged_mafs.nonMT, showTumorSampleBarcodes = TRUE, gene_mar = 9, barcode_mar = 8, pathways = 'auto', legend_height = 2)
dev.off()

pdf( width=10, height=7, file = file.path(fig.dir, paste(study.name, "_Vaf.pdf", sep="")))
plotVaf(maf=merged_mafs.nonMT)
dev.off()


RCC.sig = oncodrive(maf = merged_mafs.nonMT, AACol = 'HGVSp_Short', minMut = 2, pvalMethod = 'zscore')
