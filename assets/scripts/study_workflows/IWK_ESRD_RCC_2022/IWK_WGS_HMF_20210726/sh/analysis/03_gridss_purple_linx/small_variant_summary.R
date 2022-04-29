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
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )


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


## clinical.data <- merge(
##     x = purple.purity.dt[,list(tumor.sample.id, gender, purity, ploidy)],
##     y = patient.info.dt[,list(tumor.sample.id,
##         histology, histology.short, years.dialysis)],
##     by = c('tumor.sample.id'),
##     all.x = TRUE)

load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))
#Histology.colors.ls <- histology.colors.ls

merged.drivers.dt <- merge(
	x = merged.drivers.dt,
	y = clinical.data,
	by = 'tumor.sample.id')

driver.short.map <- c('mut', 'del', 'amp')
names(driver.short.map) <- c('MUTATION', 'DEL', 'AMP')

## work on summarizing GPL driver counts
merged.drivers.dt <- merged.drivers.dt[order(subject.id, tumor.sample.id, chromosome, chromosomeBand, -driverLikelihood),]
merged.drivers.dt[,driver.short:=driver.short.map[driver]]
merged.drivers.dt[,band.driver:=paste(gene, '(', chromosomeBand, ';', driver.short, ')', sep='')]
merged.drivers.dt[,band.driver.status:=1L]

sample.variant.mat <- dcast(
	data = merged.drivers.dt[,list(subject.id, tumor.sample.id, Histology, years.dialysis, gene, gene.driver.ct=1L)],
	formula = subject.id + tumor.sample.id + Histology + years.dialysis ~ gene,
	value.var = 'gene.driver.ct',
	fill = 0)

merged.drivers.dt.merged <- merged.drivers.dt[,list(band.driver=paste(band.driver, collapse=","), band.driver.status=sum(band.driver.status)),
	by=list(subject.id, tumor.sample.id, histology, Histology, chromosome)]
merged.drivers.dt.merged[,chr.drivers:=paste(chromosome, ':', band.driver, sep='')]

#merged.drivers.dt.merged[,list(drivers.list=paste(chr.drivers, sep='|'), driver.ct=sum(band.driver.status)), by=list(subject.id, tumor.sample.id, histology, Histology, chro)]

merged.drivers.histology.ct <- merged.drivers.dt[,list(gene.ct=sum(band.driver.status), LOH.ct=sum(biallelic)),
	by=list(histology, Histology, chromosome, chromosomeBand, gene)]
merged.drivers.histology.ct[,max.gene.ct:=max(gene.ct), by=list(gene)]
merged.drivers.histology.ct[,total.gene.ct:=sum(gene.ct), by=list(gene)]
merged.drivers.histology.ct[,total.gene.LOH.ct:=sum(LOH.ct), by=list(gene)]

merged.drivers.ct <- merged.drivers.dt[,list(gene.ct=sum(band.driver.status)), by=list(chromosome, chromosomeBand, gene)]
length(unique(merged.drivers.dt$tumor.sample.id))

merged.drivers.histology.ct[order(-total.gene.ct, gene, -gene.ct),]
## > merged.drivers.histology.ct[order(-total.gene.ct, gene, -gene.ct),]
##             histology Histology chromosome chromosomeBand    gene gene.ct LOH.ct max.gene.ct total.gene.ct total.gene.LOH.ct
## 1:           Clear cell     ccRCC       chr3          p25.3     VHL      11     11          11            12                12
## 2:            Papillary      pRCC       chr3          p25.3     VHL       1      1          11            12                12
## 3:           Clear cell     ccRCC       chr3          p21.1   PBRM1       8      8           8             8                 8
## 4:            Papillary      pRCC       chr4          q35.2    FAT1       1      0           1             3                 0
## 5:              ACD-RCC   ACD-RCC       chr4          q35.2    FAT1       1      0           1             3                 0
## 6:           Clear cell     ccRCC       chr4          q35.2    FAT1       1      0           1             3                 0
## 7:           Clear cell     ccRCC      chr22         q11.23 SMARCB1       2      0           2             3                 1
## 8:            Papillary      pRCC      chr22         q11.23 SMARCB1       1      1           2             3                 1
## 9:           Clear cell     ccRCC       chrX          q21.1    ATRX       2      1           2             2                 1
## 10:           Clear cell     ccRCC       chrX         p11.22   KDM5C       1      1           1             2                 2
## 11:            Papillary      pRCC       chrX         p11.22   KDM5C       1      1           1             2                 2
## 12:              ACD-RCC   ACD-RCC      chr19          p13.2   KEAP1       2      0           2             2                 0
## 13:           Clear cell     ccRCC       chr2          q31.2  NFE2L2       1      0           1             2                 0
## 14:              ACD-RCC   ACD-RCC       chr2          q31.2  NFE2L2       1      0           1             2                 0
## 15:              ACD-RCC   ACD-RCC       chr2          q33.1   SF3B1       1      0           1             2                 0
## 16:          Chromophobe     chRCC       chr2          q33.1   SF3B1       1      0           1             2                 0
## 17:           Clear cell     ccRCC      chr17          p13.1    TP53       2      1           2             2                 1
## 18:           Clear cell     ccRCC       chr2    p23.1-p23.2     ALK       1      0           1             1                 0
## 19:           Clear cell     ccRCC       chrX            q12      AR       1      1           1             1                 1
## 20:              ACD-RCC   ACD-RCC      chr10          q21.2  ARID5B       1      0           1             1                 0
## 21:           Clear cell     ccRCC      chr11          q22.3     ATM       1      0           1             1                 0
## 22:            Papillary      pRCC       chr1          p13.1  ATP1A1       1      0           1             1                 0
## 23:              ACD-RCC   ACD-RCC       chr3            q23     ATR       1      0           1             1                 0
## 24:           Clear cell     ccRCC       chr3          p21.1    BAP1       1      1           1             1                 1
## 25:            Papillary      pRCC      chr17            q12   CDK12       1      0           1             1                 0
## 26:           Clear cell     ccRCC      chr12         p13.31    CHD4       1      0           1             1                 0
## 27:              ACD-RCC   ACD-RCC      chr16          p13.3  CREBBP       1      0           1             1                 0
## 28:           Clear cell     ccRCC       chr5            q32   CSF1R       1      0           1             1                 0
## 29:           Clear cell     ccRCC       chr1          p34.3   CSF3R       1      1           1             1                 1
## 30:           Clear cell     ccRCC       chrX          p11.4   DDX3X       1      0           1             1                 0
## 31:           Clear cell     ccRCC       chr1          p21.3    DPYD       1      0           1             1                 0
## 32:           Clear cell     ccRCC       chr5          p13.3  DROSHA       1      0           1             1                 0
## 33:          Chromophobe     chRCC       chr6            q13  EEF1A1       1      0           1             1                 0
## 34:           Clear cell     ccRCC      chr12          q13.2   ERBB3       1      0           1             1                 0
## 35:           Clear cell     ccRCC       chr2          p16.1   FANCL       1      0           1             1                 0
## 36:           Clear cell     ccRCC       chr4          p16.3   FGFR3       1      0           1             1                 0
## 37:           Clear cell     ccRCC      chr17          p11.2    FLCN       1      1           1             1                 1
## 38:           Clear cell     ccRCC      chr13          q12.2    FLT3       1      0           1             1                 0
## 39:            Papillary      pRCC      chr12         q13.12   KMT2D       1      0           1             1                 0
## 40:           Clear cell     ccRCC      chr17      p11.2-p12   NCOR1       1      0           1             1                 0
## 41:           Clear cell     ccRCC      chr22          q11.1  OR11H1       1      0           1             1                 0
## 42:           Clear cell     ccRCC       chr5          q13.1  PIK3R1       1      0           1             1                 0
## 43:           Clear cell     ccRCC      chr19         q13.33   POLD1       1      0           1             1                 0
## 44:            Papillary      pRCC       chr8          p21.2 PPP2R2A       1      0           1             1                 0
## 45:            Papillary      pRCC      chr10         q23.31    PTEN       1      1           1             1                 1
## 46:              ACD-RCC   ACD-RCC       chr1          p34.1  RAD54L       1      0           1             1                 0
## 47:              ACD-RCC   ACD-RCC       chr1         p36.31   RPL22       1      0           1             1                 0
## 48:           Clear cell     ccRCC       chr3         p21.31   SETD2       1      1           1             1                 1
## 49:           Clear cell     ccRCC      chr19          p13.2 SMARCA4       1      0           1             1                 0
## 50:           Clear cell     ccRCC      chr15            q14  SPRED1       1      0           1             1                 0
## 51:            Papillary      pRCC      chr12         q24.21    TBX3       1      0           1             1                 0
## 52:           Clear cell     ccRCC       chr8          q22.3    UBR5       1      0           1             1                 0
## 53:           Clear cell     ccRCC       chr2          q37.1  UGT1A8       1      0           1             1                 0
## 54: Clear cell papillary    ccpRCC      chr11            p13     WT1       1      0           1             1                 0
## 55:           Clear cell     ccRCC       chr3         q13.31  ZBTB20       1      0           1             1                 0
## 56:           Clear cell     ccRCC       chr1          q21.3  ZBTB7B       1      0           1             1                 0
##                 histology Histology chromosome chromosomeBand    gene gene.ct LOH.ct max.gene.ct total.gene.ct total.gene.LOH.ct

chr3.gene.drivers.dt <- merged.drivers.dt[gene%in%c('VHL', 'PBRM1', 'SETD2', 'BAP1'),list(subject.id, tumor.sample.id, Histology, chromosome, chromosomeBand, gene)]
chr3.gene.drivers.dt <- chr3.gene.drivers.dt[order(tumor.sample.id, gene),]
chr3.gene.drivers.dt[,list(sample.ct=.N), by=list(gene)]
#> chr3.gene.drivers.dt[,list(sample.ct=.N), by=list(gene)]
## gene sample.ct
## 1:  BAP1         1
## 2:   VHL        12
## 3: PBRM1         8
## 4: SETD2         1



chr3.gene.driver.cmbns <- chr3.gene.drivers.dt[,list(genes=paste(gene, collapse=",")), by=list(subject.id, tumor.sample.id, Histology)]

chr3.gene.driver.cmbns[,list(gene.combo.ct=.N), by=list(genes)][order(-gene.combo.ct),]
## chr3.gene.driver.cmbns
## subject.id tumor.sample.id Histology           genes
## 1:     IWK002        IWK002_T     ccRCC        BAP1,VHL
## 2:     IWK005        IWK005_T     ccRCC           PBRM1
## 3:     IWK006        IWK006_T     ccRCC       PBRM1,VHL
## 4:     IWK009        IWK009_T     ccRCC             VHL
## 5:     IWK010        IWK010_T     ccRCC PBRM1,SETD2,VHL
## 6:     IWK012        IWK012_T     ccRCC           PBRM1
## 7:     IWK015        IWK015_T     ccRCC             VHL
## 8:     IWK020        IWK020_T     ccRCC             VHL
## 9:     IWK021        IWK021_T     ccRCC             VHL
## 10:     IWK022        IWK022_T     ccRCC       PBRM1,VHL
## 11:     IWK024        IWK024_T     ccRCC       PBRM1,VHL
## 12:     IWK025        IWK025_T     ccRCC           PBRM1
## 13:     IWK029        IWK029_T     ccRCC             VHL
## 14:     IWK033        IWK033_T     ccRCC       PBRM1,VHL
## 15:     IWK048        IWK048_T      pRCC             VHL

total.chr3.gene.driver.samples <- sum(chr3.gene.driver.cmbns[,list(gene.combo.ct=.N), by=list(genes)][order(-gene.combo.ct),]$gene.combo.ct)
##[1] 15

chr3.gene.driver.cmbns.ct <- chr3.gene.driver.cmbns[,list(gene.combo.ct=.N), by=list(genes)][order(-gene.combo.ct),]
chr3.gene.driver.cmbns.ct[,gene.combo.propn:=gene.combo.ct/total.chr3.gene.driver.samples]
##                 genes gene.combo.ct gene.combo.propn
## 1:             VHL             6       0.40000000
## 2:       PBRM1,VHL             4       0.26666667
## 3:           PBRM1             3       0.20000000
## 4:        BAP1,VHL             1       0.06666667
## 5: PBRM1,SETD2,VHL             1       0.06666667
fwrite(chr3.gene.driver.cmbns.ct,
	file = file.path( run.dir, "result", "variant_summaries", "chr3.LOH.gene.driver.combo.cts.tsv"), sep='\t')

# 15 samples with genes in chr3p but 16 with LOH?

merge(x=linx.clusters.dt, y=merge(x=linx.breakend.dt[tumor.sample.id=='IWK047_T',], y=linx.svs.dt[,list(tumor.sample.id, svId, clusterId, geneStart, geneEnd)], by=c('tumor.sample.id', 'svId')), by=c('subject.id', 'subject.index', 'tumor.index', 'tumor.sample.id', 'clusterId'))


merged.drivers.histology.ct[,list(gene=paste(unique(gene), collapse='\n')), by=list(Histology)]
#> merged.drivers.histology.ct[,list(gene=paste(unique(gene), collapse=',')), by=list(Histology)]
#Histology                                                                                                                                                                                                  gene
#1:      pRCC                                                                                                                                           TBX3,FAT1,SMARCB1,PTEN,KMT2D,CDK12,ATP1A1,VHL,PPP2R2A,KDM5C
#2:     ccRCC BAP1,VHL,PBRM1,SMARCA4,FANCL,OR11H1,CHD4,NFE2L2,ERBB3,NCOR1,ALK,SETD2,AR,SPRED1,ZBTB20,CSF3R,ATRX,ATM,POLD1,ZBTB7B,FLT3,FGFR3,DROSHA,PIK3R1,DPYD,KDM5C,TP53,SMARCB1,UBR5,FLCN,UGT1A8,FAT1,CSF1R,DDX3X
#3:   ACD-RCC                                                                                                                                                RAD54L,NFE2L2,KEAP1,SF3B1,FAT1,ATR,CREBBP,RPL22,ARID5B
#4:    ccpRCC                                                                                                                                                                                                   WT1
#5:     chRCC                                                                                                                                                                                          SF3B1,EEF1A1


# summary of variant counts
somatic.variant.info.dt <- merge(
	x = somatic.variant.info.dt,
	y = clinical.data,
	by = 'tumor.sample.id')

bases <- c('A', 'C', 'G', 'T')
somatic.variant.info.dt.QCd <- somatic.variant.info.dt[tumor.sample.id%in%vep.tumor.sample.ids.ls,]
somatic.variant.info.dt.QCd[,var.class:=ifelse(REF%in%bases & ALT%in%bases, 'SNV', ifelse(nchar(REF)>nchar(ALT), 'DEL', ifelse(nchar(ALT)>nchar(REF), "INS", "MNV")))]
somatic.variant.info.dt.QCd[,funcClass2:=ifelse(is.na(funcClass2), 'Unannotated', funcClass2)]

somatic.variant.info.dt.QCd[,list(var.class.ct=.N), by=list(var.class)]
##     var.class var.class.ct
## 1:       SNV       158885
## 2:       INS         5818
## 3:       DEL        16256
## 4:       MNV         2245

somatic.variant.info.dt.QCd[,list(var.class.ct=.N), by=list(funcClass2, var.class)][order(funcClass2, var.class)]

var.class.type.map <- c('SNV', 'INDEL', 'INDEL', 'MNV')
names(var.class.type.map) <- c('SNV', 'INS', 'DEL', 'MNV')
somatic.variant.info.dt.QCd[,var.type:=var.class.type.map[var.class]]

somatic.variant.info.dt.QCd.effectors <- somatic.variant.info.dt.QCd[funcClass2%in%c('NONSENSE_OR_FRAMESHIFT', 'MISSENSE', 'SPLICE'),]
somatic.variant.info.dt.QCd.effectors[,list(var.class.ct=.N), by=list(var.class)]
## var.class var.class.ct
## 1:       SNV         1182
## 2:       DEL          170
## 3:       INS           45
## 4:       MNV           21

#load("/Users/tajohnson/Documents/workspace/RemoteSystemsTempFiles/SLOGIN.HGC.JP/yshare2/ZETTAI_path_WA_slash_home_KARA/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/result/CN_segments/IWK_merged_CN_gt_0bp.Rdata")
load( file=paste(run.dir, "/result/CN_segments/IWK_merged_CN_gt_0bp.Rdata", sep=""))
sample.cn.segment.length.by.chrom <- cn.segments.resegmented.dt[,list(total.seg.length=as.numeric(sum(seg.length))), by=list(sample, chromosome)]
sample.cn.segment.length <- sample.cn.segment.length.by.chrom[,list(total.seg.length=sum(total.seg.length)), by=list(tumor.sample.id=sample)]
sample.cn.segment.length[,total.seg.Mb:=total.seg.length/1000000]


somatic.variant.info.dt.QCd <- merge(
	x = somatic.variant.info.dt.QCd,
	y = sample.cn.segment.length,
	by = 'tumor.sample.id')

sample.variant.type.cts <- somatic.variant.info.dt.QCd[,list(ct=.N), by=list(tumor.sample.id, histology, Histology, years.dialysis, total.seg.Mb, var.type)]
sample.variant.type.cts[,ct.per.Mb:=ct/total.seg.Mb]

sample.variant.cts <- somatic.variant.info.dt.QCd[,list(ct=.N), by=list(tumor.sample.id, histology, Histology, years.dialysis, total.seg.Mb)]
sample.variant.cts[,ct.per.Mb:=ct/total.seg.Mb]

sample.variant.summary <- sample.variant.cts[,
	list(sample.ct=.N, ct.mean=mean(ct), ct.sd=sd(ct), ct.per.Mb.mean=mean(ct.per.Mb), ct.per.Mb.sd=sd(ct.per.Mb)),
	by=list(histology, Histology)]
##             histology Histology sample.ct  ct.mean     ct.sd ct.per.Mb.mean ct.per.Mb.sd
## 1:            Papillary      pRCC         5 3287.000 1849.5481      1.1013985   0.61214975
## 2:           Clear cell     ccRCC        18 6836.667 1341.4475      2.2791130   0.44629502
## 3:              ACD-RCC   ACD-RCC         6 6096.167 2339.9968      2.0372630   0.77496573
## 4: Clear cell papillary    ccpRCC         2 1882.000  673.1657      0.6252699   0.22365050
## 5:          Chromophobe     chRCC         2 1684.000  185.2620      0.5644871   0.05447967

sample.variant.type.summary <- sample.variant.type.cts[,
	list(sample.ct=.N, ct.mean=mean(ct), ct.sd=sd(ct), ct.per.Mb.mean=mean(ct.per.Mb), ct.per.Mb.sd=sd(ct.per.Mb)),
	by=list(histology, Histology, var.type)]
##             histology Histology var.type sample.ct    ct.mean       ct.sd ct.per.Mb.mean ct.per.Mb.sd
## 1:            Papillary      pRCC      SNV         5 2882.00000 1602.301002    0.965697625 0.5299705696
## 2:            Papillary      pRCC    INDEL         5  360.40000  228.900197    0.120774917 0.0761527881
## 3:            Papillary      pRCC      MNV         5   44.60000   30.220854    0.014925947 0.0100519881
## 4:           Clear cell     ccRCC      SNV        18 5928.83333 1199.831912    1.976459237 0.3990117924
## 5:           Clear cell     ccRCC    INDEL        18  825.33333  259.763467    0.275150727 0.0864706578
## 6:           Clear cell     ccRCC      MNV        18   82.50000   29.609716    0.027503074 0.0098695404
## 7:              ACD-RCC   ACD-RCC      SNV         6 5256.16667 2125.282044    1.756523689 0.7041608752
## 8:              ACD-RCC   ACD-RCC    INDEL         6  760.33333  276.113503    0.254125126 0.0914113115
## 9:              ACD-RCC   ACD-RCC      MNV         6   79.66667   32.580158    0.026614167 0.0107658124
## 10: Clear cell papillary    ccpRCC      SNV         2 1632.50000  562.149891    0.542376808 0.1867669608
## 11: Clear cell papillary    ccpRCC    INDEL         2  229.00000   91.923882    0.076082260 0.0305405093
## 12: Clear cell papillary    ccpRCC      MNV         2   20.50000   19.091883    0.006810857 0.0063430289
## 13:          Chromophobe     chRCC      SNV         2 1477.00000  137.178716    0.495157061 0.0392923486
## 14:          Chromophobe     chRCC    INDEL         2  198.00000   49.497475    0.066307715 0.0157027095
## 15:          Chromophobe     chRCC      MNV         2    9.00000    1.414214    0.003022329 0.0005153865



cor.test(x=sample.variant.cts$years.dialysis, y=sample.variant.cts$ct)
## Pearson's product-moment correlation
## 
##         data:  sample.variant.cts$years.dialysis and sample.variant.cts$ct
##         t = -1.3319, df = 31, p-value = 0.1926
##         alternative hypothesis: true correlation is not equal to 0
##         95 percent confidence interval:
##         -0.5333565  0.1202683
##         sample estimates:
##         cor 
##         -0.2326465 

cor.test(x=sample.variant.cts[Histology=='ccRCC',]$years.dialysis, y=sample.variant.cts[Histology=='ccRCC',]$ct)
cor.test(x=sample.variant.cts[Histology!='ccRCC',]$years.dialysis, y=sample.variant.cts[Histology!='ccRCC',]$ct)


library(ggpubr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(multcomp)


p.SNV.vs.years <- ggscatter(sample.variant.cts, x = "years.dialysis", y = "ct", 
	color = "Histology",
#	legend.title = "RCC sub-type",
	legend.nrow = 2,
	add = "reg.line",  # Add regressin line
	add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
	conf.int = TRUE, # Add confidence interval
	xlab = "Years of dialysis", ylab = "Somatic variant count") +
	stat_cor(method = "pearson", label.x = 3, label.y = 30) +
	rremove("legend.title") +
	guides(col = guide_legend(nrow = 2)) +
	scale_color_manual(values = Histology.colors.ls)

ggexport(p.SNV.vs.years,  width=5.0, height=3.5,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.variant.cts.vs.years.dialysis.pdf'))




barplot.variant.ct.per.Mb.by.RCC.type <- ggbarplot(sample.variant.type.summary,
				x = "Histology", y = "ct.per.Mb.mean", fill="var.type", 
				color = "var.type",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Mean (#variants/Mb)", xlab = "RCC histology sub-type") +
		rotate_x_text(45) +
#		rremove('legend') +
		rremove('legend.title')

ggexport(barplot.variant.ct.per.Mb.by.RCC.type, width=4.5, height=3.0,
		filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.type.cts.per.Mb.pdf'))

barplot.variant.ct.by.RCC.type <- ggbarplot(sample.variant.type.summary,
				x = "Histology", y = "ct.mean", fill="var.type", 
				color = "var.type",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Mean (somatic variant count)", xlab = "RCC histology sub-type") +
		rotate_x_text(45) +
#		rremove('legend') +
		rremove('legend.title')

ggexport(barplot.variant.ct.by.RCC.type, width=4.0, height=3.0,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.type.cts.pdf'))


boxplot.variant.ct.by.RCC.type <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC histology sub-type") +
		rotate_x_text(45) +
#		rremove('legend') +
		rremove('legend.title') +
		scale_color_manual(values = Histology.colors.ls)

ggexport(boxplot.variant.ct.by.RCC.type, width=4.0, height=3.0,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.boxplot.pdf'))



pdf( width=5.0, height=3.5,
	file=file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.boxplot.pdf'))

my_comparisons <- list(
	c("ccpRCC", "chRCC"),
	c("pRCC", "chRCC"),
	c("pRCC", "ccpRCC"),
	c("ACD-RCC", "chRCC"),
	c("ccRCC", "pRCC"),
	c("ACD-RCC", "ccpRCC"),
	c("ccRCC", "ccpRCC"),
	c("ccRCC", "chRCC"),
	c("ACD-RCC", "pRCC"), 
	c("ccRCC", "ACD-RCC"))

boxplot.variant.ct.by.RCC.type.with.p <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend') +
		rremove('legend.title') +
		stat_compare_means(comparisons = my_comparisons, method="t.test") +
		scale_color_manual(values = Histology.colors.ls)

ggexport(boxplot.variant.ct.by.RCC.type.with.p,  width=5.0, height=3.5,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.boxplot.with.p.pdf'))


boxplot.variant.ct.by.RCC.type.no.p <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend') +
		rremove('legend.title') +
		stat_compare_means(comparisons = my_comparisons, method="t.test",
			symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
		scale_color_manual(values = Histology.colors.ls)

ggexport(boxplot.variant.ct.by.RCC.type.no.p,  width=5.0, height=3.5,
		filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.boxplot.no.p.pdf'))


my_comparisons1 <- list(
#		c("ccpRCC", "chRCC"),
#		c("pRCC", "chRCC"),
#		c("pRCC", "ccpRCC"),
		c("ACD-RCC", "chRCC"),
		c("ccRCC", "pRCC"),
		c("ACD-RCC", "ccpRCC"),
		c("ccRCC", "ccpRCC"),
		c("ccRCC", "chRCC"))
#		c("ACD-RCC", "pRCC")),
#		c("ccRCC", "ACD-RCC"))


barplot.variant.ct.per.Mb.by.RCC.type.for.figure <- ggbarplot(sample.variant.type.summary,
		x = "Histology", y = "ct.per.Mb.mean", fill="var.type", 
		legend.title = "Variant class",
		title = 'A',
		order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
		ylab = "Mean (Somatic variants/Mb)", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend.title')

boxplot.variant.ct.by.RCC.type.with.p.for.figure <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend') +
		rremove('legend.title') +
		stat_compare_means(comparisons = my_comparisons1, method="t.test") +
		scale_color_manual(values = Histology.colors.ls)


boxplot.variant.ct.by.RCC.type.no.p.for.figure <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend') +
		rremove('legend.title') +
		stat_compare_means(comparisons = my_comparisons1, method="t.test",
				symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
		scale_color_manual(values = Histology.colors.ls)

variant.analysis.plots.with.p <- ggarrange(barplot.variant.ct.per.Mb.by.RCC.type.for.figure, boxplot.variant.ct.by.RCC.type.with.p.for.figure, nrow=1, widths=c(1, 1.25))
variant.analysis.plots.no.p <- ggarrange(barplot.variant.ct.per.Mb.by.RCC.type.for.figure, boxplot.variant.ct.by.RCC.type.no.p.for.figure, nrow=1, widths=c(1, 1.25))

ggexport(variant.analysis.plots.with.p,  width=7.0, height=4.0,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.figure.with.p.pdf'))


ggexport(variant.analysis.plots.no.p,  width=7.0, height=4.0,
		filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.cts.figure.no.p.pdf'))


sample.variant.cc.noncc.summary <- sample.variant.cts[,list(sample.ct=.N, ct.mean=mean(ct), ct.sd=sd(ct)), by=list(RCC.group=ifelse(Histology=='ccRCC', 'ccRCC', 'non-ccRCC'))]
sample.variant.cc.noncc.summary
##         RCC.group sample.ct  ct.mean    ct.sd
## 1: non-ccRCC        15 4009.600 2537.606
## 2:     ccRCC        18 6836.667 1341.448


sample.variant.type.cts.wide <- dcast(
	data = sample.variant.type.cts,
	formula = tumor.sample.id + histology + Histology + years.dialysis ~ var.type,
	value.var = c('ct', 'ct.per.Mb'))

## sample.variant.type.cts.wide
## 		tumor.sample.id            histology Histology years.dialysis ct_INDEL ct_MNV ct_SNV ct.per.Mb_INDEL ct.per.Mb_MNV ct.per.Mb_SNV
## 1:        IWK001_T            Papillary      pRCC   15.833333333      217     21   2199      0.07349274   0.007112200     0.7447490
## 2:        IWK002_T           Clear cell     ccRCC    2.216666667      424     41   4116      0.14086847   0.013621715     1.3674873
## 3:        IWK004_T            Papillary      pRCC   12.666666667      327     69   2702      0.10864148   0.022924349     0.8977042
## 4:        IWK005_T           Clear cell     ccRCC    6.619444444      818     84   6271      0.27176982   0.027907903     2.0834579
## 5:        IWK006_T           Clear cell     ccRCC    6.780555556      798     76   4800      0.26526025   0.025262881     1.5955504
## 6:        IWK008_T           Clear cell     ccRCC   28.102777778      931    101   5791      0.31530756   0.034206298     1.9612740
## 7:        IWK009_T           Clear cell     ccRCC    3.083333333      661     65   8445      0.21960862   0.021595401     2.8057410
## 8:        IWK010_T           Clear cell     ccRCC    0.416666667      999     81   6727      0.33190471   0.026911193     2.2349580
## 9:        IWK012_T           Clear cell     ccRCC    8.166666667     1082    118   7638      0.35948146   0.039204078     2.5376335
## 10:        IWK013_T              ACD-RCC   ACD-RCC   24.444444444      734     73   5532      0.24858834   0.024723364     1.8735568
## 11:        IWK015_T           Clear cell     ccRCC   25.647222222      683     80   5379      0.22691783   0.026578955     1.7871025
## 12:        IWK016_T Clear cell papillary    ccpRCC   16.525000000      294     34   2030      0.09767766   0.011296056     0.6744410
## 13:        IWK017_T           Clear cell     ccRCC    3.150000000      838     91   4248      0.27841456   0.030233562     1.4113425
## 14:        IWK018_T              ACD-RCC   ACD-RCC    6.433333333      336     37   2831      0.11163162   0.012292767     0.9405628
## 15:        IWK019_T            Papillary      pRCC   12.813888889      498     59   3177      0.16866075   0.019981897     1.0759743
## 16:        IWK020_T           Clear cell     ccRCC    7.672222222      988     88   5319      0.32939157   0.029338520     1.7733135
## 17:        IWK021_T           Clear cell     ccRCC    5.502777778      544     39   4695      0.18429573   0.013212378     1.5905670
## 18:        IWK022_T           Clear cell     ccRCC    0.294444444      688     64   5151      0.22857906   0.021263168     1.7113528
## 19:        IWK024_T           Clear cell     ccRCC    9.025000000     1540    118   7256      0.51164489   0.039203959     2.4107113
## 20:        IWK025_T           Clear cell     ccRCC    4.666666667      738     55   5523      0.24519086   0.018273032     1.8349446
## 21:        IWK026_T              ACD-RCC   ACD-RCC   18.447222222      677     63   4001      0.22928380   0.021336601     1.3550435
## 22:        IWK028_T              ACD-RCC   ACD-RCC   12.316666667      910    134   9063      0.30233563   0.044519752     3.0110635
## 23:        IWK029_T           Clear cell     ccRCC   22.030555556      585     81   5671      0.19435878   0.026911215     1.8841173
## 24:        IWK030_T          Chromophobe     chRCC   27.502777778      233      8   1574      0.07741121   0.002657896     0.5229409
## 25:        IWK031_T              ACD-RCC   ACD-RCC   10.105555556     1173     94   5544      0.38971393   0.031230273     1.8419216
## 26:        IWK033_T           Clear cell     ccRCC    1.094444444      556     46   6931      0.18472416   0.015282934     2.3027394
## 27:        IWK034_T              ACD-RCC   ACD-RCC   13.638888889      732     77   4566      0.24319744   0.025582245     1.5169939
## 28:        IWK036_T Clear cell papillary    ccpRCC    6.052777778      164      7   1235      0.05448686   0.002325659     0.4103126
## 29:        IWK037_T           Clear cell     ccRCC    5.600000000     1025    154   5724      0.34054287   0.051164489     1.9017243
## 30:        IWK040_T          Chromophobe     chRCC    0.491666667      163     10   1380      0.05520422   0.003386762     0.4673732
## 31:        IWK042_T            Papillary      pRCC   27.630555556       90      4    984      0.03048086   0.001354705     0.3332574
## 32:        IWK047_T           Clear cell     ccRCC   13.033333333      958    103   7034      0.32445189   0.034883659     2.3822491
## 33:        IWK048_T            Papillary      pRCC    0.008333333      670     70   5348      0.22259875   0.023256586     1.7768032
## 		tumor.sample.id            histology Histology years.dialysis ct_INDEL ct_MNV ct_SNV ct.per.Mb_INDEL ct.per.Mb_MNV ct.per.Mb_SNV

sample.variant.type.cts.wide[,list(
	mean.SNV=mean(ct_SNV), mean.INDEL=mean(ct_INDEL), mean.MNV=mean(ct_MNV),mean.SNV.per.Mb=mean(ct.per.Mb_SNV), mean.INDEL.per.Mb=mean(ct.per.Mb_INDEL), mean.MNV.per.Mb=mean(ct.per.Mb_MNV))]
##	mean.SNV mean.INDEL mean.MNV mean.SNV.per.Mb mean.INDEL.per.Mb mean.MNV.per.Mb
##1: 4814.697   668.9091  68.0303        1.606635         0.2232157      0.02269807


sample.variant.type.cts.wide[,list(
	mean.SNV=mean(ct_SNV), mean.INDEL=mean(ct_INDEL), mean.MNV=mean(ct_MNV), mean.SNV.per.Mb=mean(ct.per.Mb_SNV), mean.INDEL.per.Mb=mean(ct.per.Mb_INDEL), mean.MNV.per.Mb=mean(ct.per.Mb_MNV)),
	by = list(Histology)]
##     Histology mean.SNV ct_INDEL mean.MNV mean.SNV.per.Mb mean.INDEL.per.Mb mean.MNV.per.Mb
## 1:      pRCC 2882.000 360.4000 44.60000       0.9656976        0.12077492     0.014925947
## 2:     ccRCC 5928.833 825.3333 82.50000       1.9764592        0.27515073     0.027503074
## 3:   ACD-RCC 5256.167 760.3333 79.66667       1.7565237        0.25412513     0.026614167
## 4:     chRCC 1632.500 229.0000 20.50000       0.5423768        0.07608226     0.006810857
## 5:    ccpRCC 1477.000 198.0000  9.00000       0.4951571        0.06630772     0.003022329


fwrite( sample.variant.type.cts.wide[,list(
	mean.SNV=mean(ct_SNV), mean.INDEL=mean(ct_INDEL), mean.MNV=mean(ct_MNV),mean.SNV.per.Mb=mean(ct.per.Mb_SNV), mean.INDEL.per.Mb=mean(ct.per.Mb_INDEL), mean.MNV.per.Mb=mean(ct.per.Mb_MNV)),
	by = list(Histology)],
	file = file.path( run.dir, "result", "maftools", "sample_variant_count_summary.tsv"), sep='\t')
		
sample.variant.cts <- somatic.variant.info.dt.QCd[,list(ct=.N), by=list(tumor.sample.id, histology, Histology, years.dialysis, total.seg.Mb)]
sample.variant.cts[,ct.per.Mb:=ct/total.seg.Mb]

sample.variant.type.cts.wide[,list(mean.SNV=mean(ct_SNV), mean.INDEL=mean(ct_INDEL), mean.MNV=mean(ct_MNV),
	median.SNV=median(ct_SNV), median.INDEL=median(ct_INDEL), median.MNV=median(ct_MNV))]
#		mean.SNV mean.INDEL mean.MNV median.SNV median.INDEL median.MNV
#	1: 4814.697   668.9091  68.0303       5319          683         70

sample.variant.effectors.summary <- dcast(
		data = somatic.variant.info.dt.QCd.effectors[,list(var.type.ct=.N), by=list(tumor.sample.id, var.type)],
		formula = tumor.sample.id ~ var.type,
		value.var = 'var.type.ct',
		fill = 0)

save(list=c('merged.drivers.dt', 'merged.drivers.histology.ct', 'merged.drivers.ct',
				'somatic.variant.info.dt.QCd', 'sample.variant.cts', 'sample.variant.summary', 'sample.variant.type.cts', 'sample.variant.type.summary',
				'sample.variant.type.cts.wide', 'sample.variant.effectors.summary'),
		file = file.path(run.dir, 'result/variant_summaries/small_variant_summaries.Rdata'))

#> sample.variant.effectors.summary
##     tumor.sample.id INDEL MNV SNV
## 1:        IWK001_T     0   0  15
## 2:        IWK002_T     6   0  23
## 3:        IWK004_T     1   0  24
## 4:        IWK005_T    11   0  31
## 5:        IWK006_T    12   0  46
## 6:        IWK008_T    11   1  40
## 7:        IWK009_T     7   1  59
## 8:        IWK010_T    13   0  65
## 9:        IWK012_T     8   3  53
## 10:        IWK013_T    13   2  51
## 11:        IWK015_T     7   1  39
## 12:        IWK016_T     3   0  17
## 13:        IWK017_T     4   0  39
## 14:        IWK018_T     3   0  33
## 15:        IWK019_T     5   1  29
## 16:        IWK020_T    15   0  35
## 17:        IWK021_T     5   0  31
## 18:        IWK022_T     2   1  34
## 19:        IWK024_T    17   2  59
## 20:        IWK025_T     5   0  39
## 21:        IWK026_T     9   0  28
## 22:        IWK028_T     5   2  71
## 23:        IWK029_T     1   0  25
## 24:        IWK030_T     2   0  14
## 25:        IWK031_T    10   0  43
## 26:        IWK033_T     6   0  32
## 27:        IWK034_T    11   2  36
## 28:        IWK036_T     0   0   7
## 29:        IWK037_T     7   2  57
## 30:        IWK040_T     0   0   9
## 31:        IWK042_T     3   0  10
## 32:        IWK047_T     9   2  55
## 33:        IWK048_T     4   1  33
##     tumor.sample.id INDEL MNV SNV

sample.variant.effectors.summary[,list(mean.SNV=mean(SNV), mean.INDEL=mean(INDEL), mean.MNV=mean(MNV),
				median.SNV=median(SNV), median.INDEL=median(INDEL), median.MNV=median(MNV))]

#	mean.SNV mean.INDEL  mean.MNV median.SNV median.INDEL median.MNV
#	35.81818   6.515152 0.6363636         34            6          0


somatic.variant.info.dt.QCd[funcClass2%in%c('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE'),list(tumor.sample.id, Histology, snp.id, CHROM, REF, ALT, gene, funcClass1, bp.descr, aa.descr, CLNSIG, BIALLELIC)][,list(var.ct=.N), by=list(tumor.sample.id, Histology, gene)][,list(sample.ct=.N, samples=paste(tumor.sample.id, collapse=';'), Histologies=paste(Histology, collapse=';')), by=list(gene)][sample.ct>2,]
