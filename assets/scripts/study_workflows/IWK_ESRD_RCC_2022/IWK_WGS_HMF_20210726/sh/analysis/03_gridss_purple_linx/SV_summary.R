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
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggpubr)
library(openxlsx)

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "SVs", "figures")

Shale.sample.info <- fread( file.path(run.dir, 'result/SVs/Shale2022/Shale et al 2022_TableS2.csv'))
Shale.SV.info <- fread( file.path(run.dir, 'result/SVs/Shale2022/Shale et al 2022_TableS3.csv'))
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

Shale.SV.info <- merge(
	x = Shale.sample.info,
	y = Shale.SV.info,
	by = c('HmfSampleId'))

Shale.SV.info[,purity.level:=cut(Purity, breaks=4, ordered=TRUE)]

Shale.SV.CancerSubType.with.purity.ct.summary <- Shale.SV.info[,list(sample.ct=.N, median.TotalSVs=median(TotalSVs)),
	by=list(PrimaryTumorLocation, CancerType, CancerSubType, purity.level)]

Shale.SV.CancerSubType.ct.summary <- Shale.SV.info[,list(sample.ct=.N, median.TotalSVs=median(TotalSVs)),
		by=list(PrimaryTumorLocation, CancerType, CancerSubType)]

Shale.SV.CancerType.ct.summary <- Shale.SV.info[,list(sample.ct=.N, median.TotalSVs=median(TotalSVs)),
	by=list(CancerType)]
> Shale.SV.CancerType.ct.summary[order(median.TotalSVs),]
##         CancerType sample.ct median.TotalSVs
## 1:          Thyroid        21            45.0
## 2:   Neuroendocrine       146            48.5
## 3:         Lymphoid        24            59.0
## 4:           Kidney       125            65.0
## 5:   Nervous system        79           111.0
## 6: Bone/Soft tissue       252           118.5
## 7:             Skin       326           125.5
## 8:            Liver        56           130.5
## 9:         Pancreas        97           141.0
## 10:          Biliary        81           152.0
## 11:    Head and neck        66           158.0
## 12:            Other       188           168.0
## 13:           Uterus        64           174.0
## 14:     Mesothelioma        39           190.0
## 15:     Colon/Rectum       594           233.5
## 16:             Lung       514           250.0
## 17:           Breast       769           294.0
## 18:    Urinary tract       175           303.0
## 19:         Prostate       394           312.5
## 20:            Ovary       161           344.0
## 21:          Stomach        45           463.0
## 22:        Esophagus       142           605.0
##         CancerType sample.ct median.TotalSVs

#mean(Shale.SV.info[CancerType=='Kidney',]$TotalSVs)
#[1] 115.976

load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

##load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))
load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))

load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.info.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

min.segment.length <- 0
#list=c('cn.segments', 'cn.segments.resegmented.dt'),
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )
#load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.ver_20220131.Rdata', sep="") ) )

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]

if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined.QCd <- rbind(
		purple.qc.dt.combined.QCd[QCStatus%in%c('PASS', 'WARN_HIGH_COPY_NUMBER_NOISE'),],
		purple.qc.dt.combined.QCd[tumor.depth>=tumor.depth.cutoff,][grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),])

purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[!tumor.sample.id%in%extra.samples.to.exclude,]

QC.pass.samples <- purple.qc.dt.combined.QCd$tumor.sample.id
length(QC.pass.samples)
# 33

#merged.SV.CN.info.dt', 'linx.SV.CN.info.dt', 'SV_CN.func.summary

clinical.data <- merge(
	x = purple.purity.dt[,list(tumor.sample.id, gender, purity, ploidy)],
	y = patient.info.dt[,list(tumor.sample.id,
		histology, histology.short, years.dialysis, age)],
	by = c('tumor.sample.id'),
	all.x = TRUE)

## clinical.data[!is.na(age),list(sample.ct=.N, mean.age=mean(age)), by=list(histology.short)]
##     histology.short sample.ct mean.age
## 1:            pRCC         4 44.75000
## 2:           ccRCC        17 65.05882
## 3:         ACD-RCC         8 58.12500
## 4:          ccpRCC         2 50.00000
## 5:           chRCC         3 60.00000
## 6:          uncRCC         1 63.00000

clinical.data <- clinical.data[tumor.sample.id%in%QC.pass.samples,]
clinical.data[,Histology:=factor(histology.short, levels=c('ccRCC', 'ACD-RCC', 'pRCC', 'ccpRCC', 'chRCC'))]

Histology.ls <- unique(clinical.data$Histology)
Histology.ls <- Histology.ls[order(Histology.ls)]
Histology.colors.ls <- get_palette(palette='jco', k=length(Histology.ls))
names(Histology.colors.ls) <- Histology.ls
##> Histology.colors.ls
##	ccRCC     ACD-RCC        pRCC      ccpRCC       chRCC 
##	"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" 
##	pRCC       ccRCC     ACD-RCC      ccpRCC       chRCC 
##	"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" 

## Note that passing Histology here must be as an integer index into the factor...
#clinical.data[,Histology.color:=Histology.colors.ls[as.character(Histology)]]
clinical.data[,Histology.color:=Histology.colors.ls[Histology]]

save(list=c('Histology.colors.ls', 'clinical.data'),
	file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))


linx.annotated.cluster.info.dt <- merge(
		x = linx.clusters.dt,
		y = linx.vis_sv_data.info.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, clusterId=ClusterId,
						svId=SvId, Type, resolvedType=ResolvedType, ChrStart, ChrEnd, PosStart, PosEnd, OrientStart, OrientEnd, InfoStart, InfoEnd,
						JunctionCopyNumber, InDoubleMinute)],
		by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'clusterId', 'resolvedType'))

linx.annotated.cluster.info.dt <- linx.annotated.cluster.info.dt[category!='ARTIFACT' & tumor.sample.id%in%QC.pass.samples,]
linx.annotated.cluster.info.dt <- merge(
	x = clinical.data[,list(tumor.sample.id, Histology, years.dialysis, purity, ploidy)],
	y = linx.annotated.cluster.info.dt,
	by = c('tumor.sample.id'))

linx.annotated.sv.cluster.cts <- linx.annotated.cluster.info.dt[,list(ct=.N), by=list(tumor.sample.id, Histology, category, Type)]


linx.annotated.sample.sv.cts <- linx.annotated.sv.cluster.cts[,list(ct=sum(ct)), by=list(tumor.sample.id, Histology)]
linx.sv.summary <- linx.annotated.sample.sv.cts[,list(mean.ct=mean(ct), median.ct=median(ct))]
#	mean.ct median.ct
#	1: 20.81818        14
linx.Histology.sv.summary <- linx.annotated.sample.sv.cts[,list(mean.ct=mean(ct), median.ct=median(ct)), by=list(Histology)]
#	Histology  mean.ct median.ct
#1:      pRCC 14.80000       9.0
#2:     ccRCC 27.77778      16.5
#3:   ACD-RCC 15.83333      15.5
#4:    ccpRCC  1.50000       1.5
#5:     chRCC  7.50000       7.5

linx.cluster.category.cts <- dcast(
	data = linx.annotated.sv.cluster.cts,
	formula = tumor.sample.id ~ category,
	value.var = 'ct',
	fun.aggregate = sum)

linx.SV.type.cts <- dcast(
		data = linx.annotated.sv.cluster.cts,
		formula = tumor.sample.id ~ Type,
		value.var = 'ct',
		fun.aggregate = sum)

merged.drivers.dt.filt <- merge(
	x = clinical.data,
	y = merged.drivers.dt.filt,
	by = "tumor.sample.id")

linx.SV.CN.info.dt <- merge(
	x = linx.SV.CN.info.dt,
	y = clinical.data,
	by = "tumor.sample.id")

linx.SV.CN.info.dt <- linx.SV.CN.info.dt[tumor.sample.id%in%QC.pass.samples,]

linx.vis_sv_data.info.dt[,SV.status:=1L]

linx.SV.info.dt <- merge(
	x = linx.SV.CN.info.dt[,list(subject.id, tumor.sample.id, histology, Histology, Histology.color, purity, ploidy, years.dialysis,
		ClusterId=clusterId, cluster.status=1L, category)],
	y = linx.vis_sv_data.info.dt,
	by = c('subject.id', 'tumor.sample.id', 'ClusterId'),
	all = TRUE)

linx.SV.info.dt[,SV.status:=ifelse(is.na(SV.status), 0L, SV.status)]
linx.SV.info.dt[,cluster.status:=ifelse(is.na(cluster.status), 0L, cluster.status)]

linx.SV.info.dt <- linx.SV.info.dt[tumor.sample.id%in%QC.pass.samples,]
# should be no clusters without SVs in table
table(linx.SV.info.dt[,list(SV.status, cluster.status)])
##         cluster.status
## SV.status    0    1
##         1 1011  687

disrupted.gene.breakends <- linx.breakend.dt[biotype=='protein_coding' & disruptive==TRUE,]
table(disrupted.gene.breakends$tumor.sample.id)
#IWK002_T IWK006_T IWK008_T IWK009_T IWK010_T IWK011_T IWK012_T IWK018_T IWK019_T IWK020_T IWK022_T IWK028_T IWK031_T IWK033_T IWK040_T IWK047_T IWK048_T 
#	2       25        9        6        2       23        5        1       13        2        2        5        2        9        2        6        6 

table(disrupted.gene.breakends[,list(type, regionType)])
##         regionType
## type  EXONIC INTRONIC
## BND      0       32
## DEL      3       15
## DUP      6        7
## INV      2       55

disrupted.gene.breakends.exonic <- disrupted.gene.breakends[regionType=='EXONIC',]

linx.SV.info.dt.disrupted.genes <- merge(
	x = linx.SV.info.dt[,list(subject.id, tumor.sample.id, Histology, Histology.color, purity, ploidy, years.dialysis, ClusterId, ChainId, svId=SvId, category, Type, ResolvedType,
			ChrStart, ChrEnd, PosStart, PosEnd, OrientStart, OrientEnd, JunctionCopyNumber)],
	y = disrupted.gene.breakends.exonic[,list(subject.id, tumor.sample.id, svId, isStart, chromosome, gene, transcriptId, geneOrientation,
		undisruptedCopyNumber, codingContext)],
	by = c('subject.id', 'tumor.sample.id', 'svId'))

linx.SV.info.dt.disrupted.genes.DEL <- linx.SV.info.dt.disrupted.genes[Type=='DEL' & ResolvedType=='DEL',]
## linx.SV.info.dt.disrupted.genes.DEL
## subject.id tumor.sample.id svId histology.short.f Histology.color purity ploidy years.dialysis ClusterId ChainId Type ResolvedType ChrStart ChrEnd PosStart   PosEnd OrientStart OrientEnd
## 1:     IWK020        IWK020_T   18             ccRCC       #0073C2FF   0.57    3.2      7.6722222         7       1  DEL          DEL    chr11  chr11 66052394 66052521           1        -1
## 2:     IWK022        IWK022_T    1             ccRCC       #0073C2FF   0.50    1.9      0.2944444        10       1  DEL          DEL     chr3   chr3 10149734 10149804           1        -1
## JunctionCopyNumber isStart chromosome  gene    transcriptId geneOrientation undisruptedCopyNumber codingContext
## 1:             1.9453    TRUE      chr11 SF3B2 ENST00000322535        Upstream                 1.070        CODING
## 2:             0.9340   FALSE       chr3   VHL ENST00000256474      Downstream                 0.072        CODING

linx.SV.info.dt.disrupted.genes.DEL[,VAF:=(JunctionCopyNumber-undisruptedCopyNumber)/JunctionCopyNumber]

genome <- BSgenome.Hsapiens.UCSC.hg38
linx.SV.info.dt.disrupted.genes.DEL[,REF:=getSeq(x=genome,
	names = linx.SV.info.dt.disrupted.genes.DEL$ChrStart,
	start = linx.SV.info.dt.disrupted.genes.DEL$PosStart,
	end = linx.SV.info.dt.disrupted.genes.DEL$PosEnd,
	strand="+",
	as.character = TRUE)]

linx.SV.info.dt.disrupted.genes.DEL[,ALT:=paste(substring(REF, 1, 1), substring(REF, nchar(REF), nchar(REF)), sep="")]
linx.SV.info.dt.disrupted.genes.DEL[,IWK020_T.GT:=ifelse(subject.id=='IWK020', ifelse( VAF>0.4 & VAF<0.6, "0/1", ifelse( VAF>0.6, "1/1", "0/0")), "0/0")]
linx.SV.info.dt.disrupted.genes.DEL[,IWK022_T.GT:=ifelse(subject.id=='IWK022', ifelse( VAF>0.4 & VAF<0.6, "0/1", ifelse( VAF>0.6, "1/1", "0/0")), "0/0")]

vcf.dt <- linx.SV.info.dt.disrupted.genes.DEL[,list(CHROM=ChrStart, POS=PosStart, ID='.', REF, ALT, QUAL='.', FILTER='PASS',
	INFO='.', FORMAT='GT', IWK020_T=IWK020_T.GT, IWK020_R='0/0', IWK022_T=IWK022_T.GT, IWK022_R='0/0')]
setnames(vcf.dt, old='CHROM', new='#CHROM')
#	IWK020_T=paste(IWK020_T.GT, ':', formatC(VAF, digits=4, width=4, format='f'), sep=''),
#	IWK022_T=paste(IWK022_T.GT, ':', formatC(VAF, digits=4, width=4, format='f'), sep=''))]

output.vcf <- file.path(run.dir, 'result', gpl_prefix, 'merged.sv.vcf')
write(x="##fileformat=VCFv4.2", file=output.vcf, append=FALSE)
write(x="##fileDate=20220125", file=output.vcf, append=TRUE)
write(x='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=output.vcf, append=TRUE)
fwrite(vcf.dt, file=output.vcf, append=TRUE, col.names=TRUE, sep='\t')
## cd /home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/result/GRIDSS-2.12.0
##	export PATH=~/.local/bin:$PATH
##	bgzip merged.sv.vcf
## tabix -p vcf merged.sv.vcf.gz

artefect.ResolvedTypes <- c('DUP_BE', 'LOW_VAF', 'PAIR_INF', 'SGL_MAPPED_INF')
## artefact SVs appears not to be in cluster, so use cluster status flag
linx.SV.info.dt.filt <- linx.SV.info.dt[cluster.status==1,]
#sum(table(linx.SV.info.dt.filt$tumor.sample.id))
#[1] 687
##sum(table(linx.SV.CN.info.dt$tumor.sample.id))
##[1] 433

linx.SV.CN.info.dt[,event.class:=ifelse(category=='SIMPLE',
	ifelse( resolvedType%in%c('DEL', 'SGL_PAIR_DEL'), 'SIMPLE (DEL)',
		ifelse(resolvedType%in%c('DUP', 'SGL_PAIR_DUP'), 'SIMPLE (DUP)', 'SIMPLE (other)')),
			ifelse(resolvedType=='LINE', 'LINE INSERTION', category))]

linx.SV.info.dt.filt[,event.class:=ifelse(category=='SIMPLE',
		ifelse( ResolvedType%in%c('DEL', 'SGL_PAIR_DEL'), 'SIMPLE (DEL)',
				ifelse(ResolvedType%in%c('DUP', 'SGL_PAIR_DUP'), 'SIMPLE (DUP)', 'SIMPLE (other)')),
		ifelse(ResolvedType=='LINE', 'LINE INSERTION', category))]

linx.sample.sv.cluster.type.ct <- linx.SV.CN.info.dt[,list(event.ct=.N),
	by=list(tumor.sample.id, histology, Histology, Histology.color,
		purity, ploidy, years.dialysis, event.class)]

linx.sample.sv.cluster.type.summary <- linx.sample.sv.cluster.type.ct[,list(event.ct.mean=mean(event.ct)),
	by=list(tumor.sample.id, histology, Histology, purity, ploidy, years.dialysis, event.class)]

linx.sv.sample.type.summary <- linx.SV.info.dt.filt[,list(Cluster.SV.Type.ct=.N),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
		purity, ploidy, years.dialysis, ClusterId, ResolvedType, event.class, Type)]


linx.sv.sample.type.summary[,Type.ct.str:=paste(Type, "(", Cluster.SV.Type.ct, ")", sep="")]

linx.sv.sample.Cluster.Type.summary <- linx.sv.sample.type.summary[,
	list(Cluster.SV.ct=sum(Cluster.SV.Type.ct),
		Cluster.Type.str=paste(unique(ResolvedType), paste(Type.ct.str, collapse=","), sep=":")),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
		purity, ploidy, years.dialysis, ClusterId, ResolvedType, event.class)]

linx.sv.sample.Cluster.summary <- linx.sv.sample.Cluster.Type.summary[,
	list(Cluster.ct=.N, Cluster.SV.ct=sum(Cluster.SV.ct), mean.Cluster.SV.ct=mean(Cluster.SV.ct)),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
		purity, ploidy, years.dialysis, ResolvedType, event.class)]

linx.sv.sample.SV.ResolvedType.summary <- linx.sv.sample.Cluster.summary[,
	list(Cluster.ct=sum(Cluster.ct), Cluster.SV.ct=sum(Cluster.SV.ct)),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
	purity, ploidy, years.dialysis, ResolvedType, event.class)]

linx.sv.sample.SV.summary <- linx.sv.sample.Cluster.summary[,
	list(Cluster.ct=sum(Cluster.ct), Cluster.SV.ct=sum(Cluster.SV.ct)),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
	purity, ploidy, years.dialysis)]

linx.sv.sample.Type.summary <- linx.sv.sample.type.summary[,
	list(Cluster.SV.Type.ct=sum(Cluster.SV.Type.ct)),
	by=list(subject.id, tumor.sample.id, histology, Histology, Histology.color,
		purity, ploidy, years.dialysis, Type)]

linx.sv.subtype.summary <- linx.sv.sample.SV.summary[,list(sample.ct=.N, mean.Cluster.ct=mean(Cluster.ct), mean.Cluster.SV.ct=mean(Cluster.SV.ct)),
	by=list(histology, Histology, Histology.color)]
## linx.sv.subtype.summary
## histology Histology Histology.color sample.ct mean.Cluster.ct mean.Cluster.SV.ct
## 1:            Papillary      pRCC       #868686FF         5         9.40000           14.80000
## 2:           Clear cell     ccRCC       #0073C2FF        18        16.33333           27.77778
## 3:              ACD-RCC   ACD-RCC       #EFC000FF         6        12.33333           15.83333
## 4: Clear cell papillary    ccpRCC       #CD534CFF         2         1.50000            1.50000
## 5:          Chromophobe     chRCC       #7AA6DCFF         2         7.50000            7.50000



fwrite(linx.sv.subtype.summary,
	file = file.path(run.dir, 'result/variant_summaries/SV_counts_by_RCC_subtype.tsv'),
	sep = '\t', col.names = TRUE)

linx.sv.subtype.ResolvedType.summary <- linx.sv.sample.SV.ResolvedType.summary[,
	list(Cluster.ct=sum(Cluster.ct), Cluster.SV.ct=sum(Cluster.SV.ct), mean.Cluster.ct=mean(Cluster.ct), mean.Cluster.SV.ct=mean(Cluster.SV.ct)),
	by=list(histology, Histology, Histology.color, ResolvedType, event.class)]

linx.sv.subtype.ResolvedType.summary[,status:=ifelse(Cluster.ct>0, 1L, 0L)]
linx.sv.subtype.ResolvedType.summary[,total.Cluster.SV.ct:=sum(Cluster.SV.ct), by=list(Histology)]
linx.sv.subtype.ResolvedType.summary[,Cluster.SV.propn:=Cluster.SV.ct/total.Cluster.SV.ct]

linx.sv.subtype.ResolvedType.summary[,total.mean.Cluster.SV.ct:=sum(mean.Cluster.SV.ct), by=list(Histology)]
linx.sv.subtype.ResolvedType.summary[,mean.Cluster.SV.ct.propn:=mean.Cluster.SV.ct/total.mean.Cluster.SV.ct, by=list(Histology)]

#linx.sv.multiCluster.ct <- linx.sv.subtype.ResolvedType.summary[,list(status.ct=sum(status)), by=list(ResolvedType)]
#ResolvedTypes.to.plot <- linx.sv.multiCluster.ct[status.ct>1,][order(-status.ct),]$ResolvedType
#last.ResolvedType <- ResolvedTypes.to.plot[[length(ResolvedTypes.to.plot)]]

linx.sv.multiCluster.ct <- linx.sv.subtype.ResolvedType.summary[,list(status.ct=sum(status)), by=list(event.class)]
event.classes.to.plot <- linx.sv.multiCluster.ct[status.ct>0,][order(-status.ct),]$event.class
last.event.class <- event.classes.to.plot[[length(event.classes.to.plot)]]

## SV.Cluster.ct.plot.ls <- lapply(
##     X = ResolvedTypes.to.plot,
##     FUN = function ( curr.ResolvedType ){
##         p <- ggboxplot(linx.sv.sample.SV.ResolvedType.summary[ResolvedType==curr.ResolvedType,],
##             x = "Histology", y = "Cluster.ct", 
##             outlier.shape = NA,
##             title = curr.ResolvedType,
##             color = "Histology",
##             order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
##             add = "jitter",
##             ylab = "Cluster count", xlab = "RCC sub-type") +
##             theme( plot.margin = margin(0.1,0.1,0.1,0.1, "cm") ) +
##             theme( plot.title = element_text(hjust = 0.5) ) +
##             rremove('legend') +
##             rremove('xlab') +
##             rotate_x_text(45) +
##             scale_color_manual(values = Histology.colors.ls)
##         p;
## })

SV.Cluster.ct.plot.ls <- lapply(
	X = event.classes.to.plot,
	FUN = function ( curr.event.class ){
		p <- ggboxplot(linx.sv.sample.SV.ResolvedType.summary[event.class==curr.event.class,],
						x = "Histology", y = "Cluster.ct", 
						outlier.shape = NA,
						title = curr.event.class,
						color = "Histology",
						order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
						add = "jitter",
						ylab = "Cluster count", xlab = "RCC sub-type") +
				theme( plot.margin = margin(0.1,0.1,0.1,0.1, "cm") ) +
				theme( plot.title = element_text(hjust = 0.5) ) +
				rremove('legend') +
				rremove('xlab') +
				rotate_x_text(45) +
				scale_color_manual(values = Histology.colors.ls)
		p;
	})


SV.Cluster.ct.plots <- ggarrange( plotlist=SV.Cluster.ct.plot.ls, nrow=4 , ncol=2)
ggexport(SV.Cluster.ct.plots, width=7.0, height=2.0*(length(event.classes.to.plot)/2),
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.cluster.summary.figure.pdf'))


linx.sv.subtype.Type.summary <- linx.sv.sample.Type.summary[,
	list(Cluster.SV.Type.ct=sum(Cluster.SV.Type.ct), mean.Cluster.SV.Type.ct=mean(Cluster.SV.Type.ct)),
	by=list(histology, Histology, Histology.color, Type)]

linx.sv.subtype.Type.summary[,status:=ifelse(Cluster.SV.Type.ct>0, 1L, 0L)]

linx.sv.multiSV.ct <- linx.sv.subtype.Type.summary[,list(status.ct=sum(status)), by=list(Type)]
Types.to.plot <- linx.sv.multiSV.ct[status.ct>0,][order(-status.ct),]$Type

SV.plot.ls <- lapply(
	X = Types.to.plot,
	FUN = function ( curr.Type ){
		p <- ggboxplot(linx.sv.sample.Type.summary[Type==curr.Type,],
				x = "Histology", y = "Cluster.SV.Type.ct", 
				outlier.shape = NA,
				title = curr.Type,
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				add = "jitter",
				ylab = "SV count", xlab = "RCC sub-type") +
			theme( plot.margin = margin(0.1,0.1,0.1,0.1, "cm") ) +
			theme( plot.title = element_text(hjust = 0.5) ) +
			rremove('legend') +
			rremove('xlab') +
			rotate_x_text(45) +
			scale_color_manual(values = Histology.colors.ls)
		p;
	})

SV.plots <- ggarrange( plotlist=SV.plot.ls, nrow=4 , ncol=2)
ggexport(SV.plots, width=7.0, height=2.0*(length(Types.to.plot)/2), filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.Type.summary.figure.pdf'))


#c('merged.drivers.dt', 'merged.drivers.histology.ct', 'merged.drivers.ct',
#'somatic.variant.info.dt.QCd', 'sample.variant.cts', 'sample.variant.summary', 'sample.variant.type.cts', 'sample.variant.type.summary'),
load( file = file.path(run.dir, 'result/variant_summaries/small_variant_summaries.Rdata'))
#setnames(sample.variant.cts, old=c('histology.short.f'), new=c('Histology'))

sample.small_variant.SV.cts <- merge(
	x = linx.sv.sample.SV.summary,
	y = sample.variant.cts,
	by = c('tumor.sample.id', 'histology', 'Histology', 'years.dialysis'))

barplot.SV.cluster.ct.by.RCC.type <- ggbarplot(linx.sv.subtype.summary, x = "Histology", y = "mean.Cluster.ct", 
#				legend.title = "Variant class",
				title = 'A',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Mean (Cluster ct.)", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		rremove('legend') +
		rremove('legend.title') +
		scale_color_manual(values = Histology.colors.ls)


## Figure 2
var_ct_comparisons <- list(
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
		theme(legend.text = element_text(size = 6, color = "black", face = "plain"),
				legend.key.size = unit(0.7,"line"),
				legend.position = c(0.5, 0.8)) +
		rotate_x_text(45) +
		rremove('legend.title') +
		rremove('x.text') +
		rremove('xlab') +
		theme(title = element_text(size=10, face='bold', color='black'), axis.title = element_text(size=8), axis.text = element_text(size=8))


boxplot.variant.ct.by.RCC.type.no.p.for.figure <- ggboxplot(sample.variant.cts,
				x = "Histology", y = "ct", 
				title = 'B',
				add = 'jitter',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Somatic variant count", xlab = "RCC sub-type") +
		rotate_x_text(45) +
		stat_compare_means(comparisons = var_ct_comparisons, method="t.test", vjust=0.4,
				symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
		scale_color_manual(values = Histology.colors.ls) +
		scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
		rremove('legend') +
		rremove('legend.title') +
		rremove("xlab") +
		rremove("x.text") +
		theme(plot.margin = margin(0.15,0.1,0.2,0.0, "cm") ) +
		theme(title = element_text(size=10, face='bold', color='black'), axis.title = element_text(size=8), axis.text = element_text(size=8))


barplot.SV.cluster.Type.ct.by.RCC.type <- ggbarplot(linx.sv.subtype.ResolvedType.summary,
				x = "Histology", y = "mean.Cluster.ct", fill="event.class", color="event.class", 
				legend.title = "SV cluster class",
#				legend = 'right',
				title = 'C',
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				ylab = "Mean (Cluster ct.)", xlab = "RCC sub-type") +
		theme(legend.text = element_text(size = 6, color = "black", face = "plain"),
			legend.key.size = unit(0.4,"line"),
			legend.position = c(0.6, 0.9)) +
		rotate_x_text(45) +
		rremove('legend.title') +
		theme(title = element_text(size=10, face='bold', color='black'), axis.title = element_text(size=8), axis.text = element_text(size=8))

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

boxplot.SV.cluster.ct.by.RCC.type <- ggboxplot(linx.sv.sample.SV.summary, x = "Histology", y = "Cluster.ct", 
				outlier.shape = NA,
				title = 'D',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				add = "jitter",
				ylab = "Somatic SV cluster ct.", xlab = "RCC sub-type") +
		stat_compare_means(comparisons = my_comparisons, method="t.test") +
		rremove('legend') +
		rotate_x_text(45) +
		scale_color_manual(values = Histology.colors.ls)

boxplot.SV.ct.by.RCC.type <- ggboxplot(linx.sv.sample.SV.summary, x = "Histology", y = "Cluster.SV.ct", 
				outlier.shape = NA,
				title = 'D',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				add = "jitter",
				ylab = "Somatic SV ct.", xlab = "RCC sub-type") +
		stat_compare_means(comparisons = my_comparisons, method="t.test") +
		rremove('legend') +
		rotate_x_text(45) +
		scale_color_manual(values = Histology.colors.ls)

SV.analysis.plot.with.p <- ggarrange(barplot.SV.cluster.ct.by.RCC.type, boxplot.SV.cluster.ct.by.RCC.type, boxplot.SV.ct.by.RCC.type, nrow=1, widths=c(1, 1.25, 1.25))
ggexport(SV.analysis.plot.with.p, width=7.0, height=4.0, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.cts.figure.with.p.pdf'))



my_comparisons1 <- list(
		c("ccpRCC", "chRCC"),
#		c("pRCC", "chRCC"),
#		c("pRCC", "ccpRCC"),
		c("ACD-RCC", "chRCC"),
#		c("ccRCC", "pRCC"),
		c("ACD-RCC", "ccpRCC"),
		c("ccRCC", "ccpRCC"))
#		c("ccRCC", "chRCC"))
#		c("ACD-RCC", "pRCC"), 
#		c("ccRCC", "ACD-RCC"))

boxplot.SV.cluster.ct.by.RCC.type.no.p <- ggboxplot(linx.sv.sample.SV.summary, x = "Histology", y = "Cluster.ct", 
				outlier.shape = NA,
				title = 'B',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				add = "jitter",
				ylab = "Somatic SV cluster ct.", xlab = "RCC sub-type") +
		stat_compare_means(comparisons = my_comparisons1, method="t.test",
				symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
		rremove('legend') +
		rotate_x_text(45) +
		scale_color_manual(values = Histology.colors.ls)

my_comparisons2 <- list(
		c("ccpRCC", "chRCC"),
#		c("pRCC", "chRCC"),
#		c("pRCC", "ccpRCC"),
		c("ACD-RCC", "chRCC"),
#		c("ccRCC", "pRCC"),
		c("ACD-RCC", "ccpRCC"),
		c("ccRCC", "ccpRCC"),
		c("ccRCC", "chRCC"))
#		c("ACD-RCC", "pRCC"), 
#		c("ccRCC", "ACD-RCC"))


boxplot.SV.ct.by.RCC.type.no.p <- ggboxplot(linx.sv.sample.SV.summary, x = "Histology", y = "Cluster.SV.ct", 
				outlier.shape = NA,
				title = 'D',
				color = "Histology",
				order = c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'),
				add = "jitter",
				ylab = "Somatic SV ct.", xlab = "RCC sub-type") +
		stat_compare_means(comparisons = my_comparisons2, method="t.test", vjust=0.4,
				symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
		rremove('legend') +
		rotate_x_text(45) +
		scale_color_manual(values = Histology.colors.ls) +
		scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
		theme(plot.margin = margin(0.1,0.1,0.3,0.45, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'), axis.title = element_text(size=8), axis.text = element_text(size=8))


SV.analysis.plot.no.p <- ggarrange(barplot.SV.cluster.ct.by.RCC.type, boxplot.SV.cluster.ct.by.RCC.type.no.p, boxplot.SV.ct.by.RCC.type.no.p, nrow=1, widths=c(1, 1.25, 1.25))
ggexport(SV.analysis.plot.no.p, width=7.0, height=4.0, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.cts.figure.no.p.pdf'))

SV.analysis.with.Cluster.type.plot.no.p <- ggarrange(barplot.SV.cluster.Type.ct.by.RCC.type, boxplot.SV.ct.by.RCC.type.no.p, ncol=1, nrow=2, widths=c(1, 1), heights=c(1, 1))
ggexport(SV.analysis.with.Cluster.type.plot.no.p, width=3.5, height=7.0, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.Cluster.type.SV.cts.figure.no.p.pdf'))



variant.analysis.plot.no.p <- ggarrange(
	barplot.variant.ct.per.Mb.by.RCC.type.for.figure,
	boxplot.variant.ct.by.RCC.type.no.p.for.figure,
	barplot.SV.cluster.Type.ct.by.RCC.type,
	boxplot.SV.ct.by.RCC.type.no.p,
	ncol=2, nrow=2, widths=c(1, 1), heights=c(0.8, 1))

ggexport(variant.analysis.plot.no.p, width=4.5, height=4.5, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variant.SV.cts.figure.no.p.pdf'))


p.SNVs.vs.years <- ggscatter(sample.small_variant.SV.cts, x = "years.dialysis", y = "ct", 
				color = "Histology",
				title = 'A',
				legend.nrow = 2,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Years of dialysis", ylab = "Somatic small variant ct.") +
		stat_cor(method = "pearson", label.x = 9, label.y = 1500, size=2.5) +
		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.2, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

p.SVs.vs.years <- ggscatter(sample.small_variant.SV.cts, x = "years.dialysis", y = "Cluster.SV.ct", 
				color = "Histology",
				title = 'B',
				legend.nrow = 2,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Years of dialysis", ylab = "Somatic SV ct.") +
		stat_cor(method = "pearson", label.x = 5, label.y = 50, size=2.5) +
		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.2, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

p.SVs.vs.years.no.outlier <- ggscatter(sample.small_variant.SV.cts[tumor.sample.id!='IWK006_T',], x = "years.dialysis", y = "Cluster.SV.ct", 
				color = "Histology",
				title = 'B',
				legend.nrow = 2,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Years of dialysis", ylab = "Somatic SV ct.") +
		stat_cor(method = "pearson", label.x = 5, label.y =35, size=2.5) +
		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


variants.years.dialysis.plot <- ggarrange(p.SNVs.vs.years, p.SVs.vs.years, ncol=2, nrow=1, widths=c(1, 1), common.legend = TRUE, legend="bottom")
ggexport(variants.years.dialysis.plot, width=4.5, height=2.5,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variants.years.dialysis.figure.pdf'))

variants.years.dialysis.plot.no.outlier <- ggarrange(p.SNVs.vs.years, p.SVs.vs.years.no.outlier, ncol=2, nrow=1, widths=c(1, 1), common.legend = TRUE, legend="bottom")
ggexport(variants.years.dialysis.plot.no.outlier, width=4.5, height=2.5,
	filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variants.years.dialysis.no_outliers.figure.pdf'))

p.SV.vs.SNV <- ggscatter(sample.small_variant.SV.cts, x = "ct", y = "Cluster.SV.ct", 
				color = "Histology",
				title = 'A',
				legend.nrow = 2,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic small variant ct.", ylab = "Somatic SV ct.") +
		stat_cor(method = "pearson", label.x = 3, label.y = 60, size=2.5) +
#		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
			legend.key.size = unit(0.6,"line"),
			legend.text = element_text(size=8), legend.title = element_text(size=8),
			axis.title = element_text(size=8), axis.text = element_text(size=8))

p.SV.vs.SNV.no.outlier <- ggscatter(sample.small_variant.SV.cts[tumor.sample.id!='IWK006_T',], x = "ct", y = "Cluster.SV.ct", 
				color = "Histology",
				title = 'B',
				legend.nrow = 2,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic small variant ct.", ylab = "Somatic SV ct.") +
				stat_cor(method = "pearson", label.x = 3, label.y = 40, size=2.5) +
		rremove("ylab") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.1, 0.2, 0.1, 0.0, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8), legend.title = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


p.SV.vs.SNV.cmbd.plot <- ggarrange(p.SV.vs.SNV, p.SV.vs.SNV.no.outlier, ncol=2, nrow=1, widths=c(1, 0.90), heights=c(1, 1), common.legend = TRUE, legend="bottom")
ggexport(p.SV.vs.SNV.cmbd.plot, width=4.5, height=2.5, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.vs.small_variant.cts.cmbd.figure.pdf'))


ggexport(p.SV.vs.SNV, width=2.5, height=2.5, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.vs.small_variant.cts.figure.pdf'))
ggexport(p.SV.vs.SNV.no.outlier, width=2.5, height=2.5, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.SV.vs.small_variant.no_outliers.cts.figure.pdf'))



p.ploidy.vs.SNV <- ggscatter(sample.small_variant.SV.cts, x = "ct", y = "ploidy", 
				color = "Histology",
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic small variant ct.", ylab = "Tumor ploidy") +
		stat_cor(method = "pearson", label.x = 3, label.y = 1.5) +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls)

p.ploidy.vs.SV <- ggscatter(sample.small_variant.SV.cts, x = "Cluster.SV.ct", y = "ploidy", 
				color = "Histology",
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic SV ct.", ylab = "Tumor ploidy") +
		stat_cor(method = "pearson", label.x = 3, label.y = 1.5) +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls)

p.purity.vs.SNV <- ggscatter(sample.small_variant.SV.cts, x = "ct", y = "purity", 
				color = "Histology",
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic small variant ct.", ylab = "Tumor ploidy") +
		stat_cor(method = "pearson", label.x = 3, label.y = 0.8) +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls)

p.purity.vs.SV <- ggscatter(sample.small_variant.SV.cts, x = "Cluster.SV.ct", y = "purity", 
				color = "Histology",
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Somatic SV ct.", ylab = "Tumor purity") +
		stat_cor(method = "pearson", label.x = 3, label.y = 0.8) +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls)




variants.purity.ploidy.plot <- ggarrange(p.ploidy.vs.SNV, p.ploidy.vs.SV, p.purity.vs.SNV, p.purity.vs.SV, ncol=2, nrow=2, widths=c(1, 1, 1, 1), heights=c(1, 1, 1, 1), common.legend = TRUE, legend="bottom")
ggexport(variants.purity.ploidy.plot, filename = file.path(run.dir, 'result/variant_summaries/RCC.sub.type.variants.purity.ploidy.figure.pdf'))


save(list=c('linx.SV.info.dt.filt', 'disrupted.gene.breakends', 'linx.SV.info.dt.disrupted.genes', 'linx.sv.subtype.summary', 
	'linx.sv.sample.Type.summary', 'linx.sv.sample.SV.summary', 'linx.sv.sample.SV.ResolvedType.summary',
	'linx.sv.sample.Cluster.summary', 'linx.sv.sample.Cluster.Type.summary', 'linx.sv.sample.type.summary', 'linx.sample.sv.cluster.type.summary'),
		file = file.path(run.dir, 'result/variant_summaries/SV_summaries.Rdata'))

cn.segments.resegmented.dt <- merge(
	x = clinical.data[,list(sample=tumor.sample.id, Histology, years.dialysis)],
	y = cn.segments.resegmented.dt,
	by = c('sample'))

check[,list(seg.length=sum(seg.length)), by=list(nTotal)]/sum(check$seg.length)
cn.segments.resegmented.dt[,chrom.seg.length:=sum(seg.length), by=list(sample, chromosome)]
cn.segments.resegmented.dt[,seg.chrom.propn:=seg.length/chrom.seg.length]
cn.segments.summary <- cn.segments.resegmented.dt[,list(seg.length=sum(seg.length)), by=list(sample, chromosome, nTotal)]

cn.segments.resegmented.dt[,CN.mult:=copyNumber*seg.length]



disrupted.exon.gene.samples <- linx.vis_gene_exon.info.dt[,list(row.ct=.N), by=list(tumor.sample.id, ClusterId, Chromosome, Gene, AnnotationType)][AnnotationType=='DISRUPTION',]

disrupted.gene.SV.info <- merge(x = linx.clusters.dt, y=linx.vis_gene_exon.info.dt[AnnotationType=='DISRUPTION',list(row.ct=.N), by=list(tumor.sample.id, ClusterId, Chromosome, Gene, AnnotationType)][,list(tumor.sample.id, clusterId=ClusterId, Chromosome, Gene)], by=c('tumor.sample.id', 'clusterId'))[category!='ARTIFACT',]
#tumor.sample.id clusterId category synthetic resolvedType clusterCount clusterDesc subject.id subject.index tumor.index Chromosome Gene
#1:        IWK008_T        22  COMPLEX     FALSE      COMPLEX            3 DEL=2_DUP=1     IWK008             6           6       chr4 FAT1

disrupted.gene.SV.info.with.vcfId <- merge(
	x = disrupted.gene.SV.info,
	y = linx.svs.dt[,list(tumor.sample.id, clusterId, vcfId, svId, geneStart, geneEnd)],
	by = c('tumor.sample.id', 'clusterId'))

disrupted.gene.SV.info.with.vcfId[geneStart==Gene | geneEnd==Gene,]
disrupted.gene.SV.info.with.vcfId.matches <- disrupted.gene.SV.info.with.vcfId[geneStart==Gene | geneEnd==Gene,]
#	tumor.sample.id clusterId category synthetic resolvedType clusterCount clusterDesc subject.id subject.index tumor.index Chromosome Gene            vcfId svId geneStart geneEnd
#1:        IWK008_T        22  COMPLEX     FALSE      COMPLEX            3 DEL=2_DUP=1     IWK008             6           6       chr4 FAT1 gridss305fb_909o   15              FAT1
setnames(disrupted.gene.SV.info.with.vcfId.matches, old='vcfId', new='ID')
purple.sv.variant.info.dt <- purple.sv.variant.info.dt[,unique(colnames(purple.sv.variant.info.dt)), with=FALSE]

purple.sv.variant.info.dt.disrupted.gene.SV <- merge(
	x = disrupted.gene.SV.info.with.vcfId.matches,
	y = purple.sv.variant.info.dt,
	by = c('subject.id', 'tumor.sample.id', 'ID'))



SV.vcf.dt <- purple.sv.variant.info.dt.disrupted.gene.SV[,list(CHROM, POS, ID, REF, ALT, QUAL='.', INFO='', FORMAT='GT', IWK008_T='0/1', IWK008_N='0/0')]
setnames(SV.vcf.dt, old='CHROM', new='#CHROM')

linx.SV.info.dt.disrupted.genes.DEL[,list(CHROM=ChrStart, POS=PosStart, ID='.', REF, ALT, QUAL='.', FILTER='PASS',
				INFO='.', FORMAT='GT', IWK020_T=IWK020_T.GT, IWK020_R='0/0', IWK022_T=IWK022_T.GT, IWK022_R='0/0')]
setnames(vcf.dt, old='CHROM', new='#CHROM')
