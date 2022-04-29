#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=325)
options(width=215)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )
#run.dir <- '/Users/tajohnson/Documents/RIKEN/NakagawaLab/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'

library(data.table)
library(parallel)
library(ggpubr)
library(ggsignif)

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result_summaries", "adhoc")


#TCMA.base.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC'
#TCMA.base.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC'
TCMA.base.dir <- file.path(run.dir, 'result/mutations/somatic/')
load( file = file.path(TCMA.base.dir, 'MitochondrialMutations/TCMA/TCMA.processed.Rdata'))

TCMA.CN.RCC <- TCMA.CN[cancer_type%in%c('Kidney-ChRCC', 'Kidney-RCC'),]
TCMA.SNV.INDEL.with.Consequence.RCC <- TCMA.SNV.INDEL.with.Consequence[cancer_type%in%c('Kidney-ChRCC', 'Kidney-RCC'),]
#	unique(TCMA.CN.RCC$dcc_project_code)
#	"KICH-US" "KIRC-US" "KIRP-US" "RECA-EU"
TCMA.Histology.ls <- c('KIRC-US', 'RECA-EU', 'KIRP-US', 'KICH-US')
names(TCMA.Histology.ls) <- TCMA.Histology.ls;

load(file = paste(run.dir, '/result_summaries/adhoc/germline.variant.coverage.Rdata', sep=""))

load(file = file.path(run.dir, 'result/variant_summaries/mtDNA.variant.summaries.Rdata'))

load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))
#Histology.colors.ls <- histology.colors.ls
Histology.ls <- names(Histology.colors.ls)
names(Histology.ls) <- Histology.ls


Histology.ls <- c(Histology.ls, TCMA.Histology.ls)
Histology.colors.ls <- get_palette(palette='jco', k=length(Histology.ls))
names(Histology.colors.ls) <- Histology.ls
##	> Histology.colors.ls
##	ccRCC     ACD-RCC        pRCC      ccpRCC       chRCC     KIRC-US     RECA-EU     KIRP-US     KICH-US 
##	"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" "#003C67FF" "#8F7700FF" "#3B3B3BFF" "#A73030FF" 



TCMA.SNV.INDEL.age.summary <- TCMA.SNV.INDEL.with.Consequence.RCC[,list(var.ct=.N),
		by=list(sample_id, Histology=dcc_project_code, age=donor_age_at_diagnosis)]

excluded.samples.ls <- c('IWK011_T', 'IWK023_T', 'IWK039_T', 'IWK041_T', 'IWK049_T')

IWK.germline.variant.coverage.dt <- germline.variant.coverage.dt[study.group%in%c('IWK') & !tumor.sample.id%in%excluded.samples.ls,
		list(study.group, subject.id, tumor.sample.id, histology, mean.DP.N_chrM, mean.DP.T_chrM, median.DP.N_chrM, median.DP.T_chrM, mtDNA.copy.est.N, mtDNA.copy.est.T)]

IWK.germline.variant.coverage.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, sex=ifelse(gender=='MALE', 'M', ifelse(gender=='FEMALE', 'F', 'UNK')), age, years.dialysis, Histology, Histology.color)],
		y = IWK.germline.variant.coverage.dt,
		by = c('tumor.sample.id'))

IWK.germline.variant.mtCN.dt.long <- melt(
		data = IWK.germline.variant.coverage.dt[,list(study.group, tumor.sample.id, Histology, Histology.color,
			sex, age, years.dialysis,
			Normal = mtDNA.copy.est.N, Tumor=mtDNA.copy.est.T)],
		id.vars = c('study.group', 'tumor.sample.id', 'Histology', 'Histology.color', 'sex', 'age', 'years.dialysis'),
		value.vars = c('Normal', 'Tumor'),
		variable.name = 'Sample.type',
		value.name = 'Mito.CN')

IWK.germline.variant.mtCN.dt.long[,Histology:=factor(as.character(Histology), levels=Histology.ls, ordered=TRUE)]
IWK.germline.variant.mtCN.dt.long[,Histology.color:=Histology.colors.ls[as.character(Histology)]]

TCMA.CN.RCC[,Histology:=factor(dcc_project_code, levels=Histology.ls, ordered=TRUE)]
TCMA.CN.RCC[,Histology.color:=Histology.colors.ls[as.character(Histology)]]

IWK.germline.variant.mtCN.dt.long <- rbind(
	IWK.germline.variant.mtCN.dt.long[,list(study.group, tumor.sample.id, Histology, Histology.color, Sample.type,
		sex, age, years.dialysis, Mito.CN)],
	TCMA.CN.RCC[,list(study.group='TCMA', tumor.sample.id=sample_id, Histology, Histology.color, Sample.type='Tumor',
		sex = ifelse(donor_sex=='male', 'M', ifelse(donor_sex=='female', 'F', 'UNK')),
		age = donor_age_at_diagnosis, years.dialysis=as.numeric(NA), Mito.CN=tumor_copy_number)])

IWK.germline.variant.mtCN.dt <- dcast(
		IWK.germline.variant.mtCN.dt.long,
		formula = study.group + tumor.sample.id + Histology + Histology.color + sex + age + years.dialysis ~ Sample.type,
		value.var = 'Mito.CN',
		fill = 0)
setnames(IWK.germline.variant.mtCN.dt, old=c('Normal', 'Tumor'), new=c('mtDNA.copy.est.N', 'mtDNA.copy.est.T'))

IWK.germline.variant.coverage.dt[,list(sample.ct=.N), by=list(study.group, Histology)]
#1:         IWK      pRCC         5
#2:         IWK     ccRCC        18
#3:         IWK   ACD-RCC         6
#4:         IWK    ccpRCC         2
#5:         IWK     chRCC         2

IWK.germline.variant.mtCN.dt[,list(sample.ct=.N), by=list(study.group, Histology)]
## study.group Histology sample.ct
## 1:         IWK      pRCC         5
## 2:         IWK     ccRCC        18
## 3:         IWK   ACD-RCC         6
## 4:         IWK    ccpRCC         2
## 5:         IWK     chRCC         2
## 6:        TCMA   KICH-US        42
## 7:        TCMA   KIRC-US        32
## 8:        TCMA   KIRP-US        31
## 9:        TCMA   RECA-EU        58

#mean depth of coverage at chrM for normals and tumor tissues
#germline.variant.coverage.dt[study.name=='IWK',list(mean.DP.N_chrM=mean(mean.DP.N_chrM), mean.DP.T_chrM=mean(mean.DP.T_chrM))]
#mean.DP.N_chrM mean.DP.T_chrM
#1:       6701.951       18640.39

IWK.germline.variant.coverage.dt[,list(mean.DP.N_chrM=mean(mean.DP.N_chrM), mean.DP.T_chrM=mean(mean.DP.T_chrM),
	median.DP.N_chrM=median(median.DP.N_chrM), median.DP.T_chrM=median(median.DP.T_chrM)), by=list(study.group, Histology)]
##         Histology mean.DP.N_chrM mean.DP.T_chrM median.DP.N_chrM median.DP.T_chrM
## 1:      pRCC      11808.888       24488.02          7154.00         17043.00
## 2:     ccRCC       4138.500       11047.57          2147.75          9446.50
## 3:   ACD-RCC       4307.611       23320.54          2767.00         22788.50
## 4:    ccpRCC       4408.418       12503.85          4450.25         12675.75
## 5:     chRCC      11909.500       68905.01         11965.50         72226.50

#	study.group Histology median.mtCN.N median.mtCN.T mean.mtCN.N mean.mtCN.T
#1:         IWK      pRCC      474.8389      851.7313    474.8389   1220.4137
#2:         IWK     ccRCC      397.2651      469.7088    397.2651    491.5200
#3:         IWK   ACD-RCC      116.5452     1049.1144    116.5452   1083.3385
#4:         IWK    ccpRCC      142.9048      654.8807    142.9048    654.8807
#5:         IWK     chRCC      558.4279     3132.7198    558.4279   3132.7198
#6:        TCMA   KICH-US            NA      491.3415         NaN    624.7170
#7:        TCMA   KIRC-US            NA      340.1201         NaN    337.0147
#8:        TCMA   KIRP-US            NA      242.3362         NaN    353.4045
#9:        TCMA   RECA-EU            NA      323.4332         NaN    356.1543
#

IWK.germline.variant.mtCN.dt[,list(
	median.mtCN.N=median(ifelse(study.group=='TCMA', as.numeric(NA), mtDNA.copy.est.N), na.rm=TRUE),
	median.mtCN.T=median(mtDNA.copy.est.T),
	mean.mtCN.N=mean(ifelse(study.group=='TCMA', as.numeric(NA), mtDNA.copy.est.N), na.rm=TRUE),
	mean.mtCN.T=mean(mtDNA.copy.est.T)),
by = list(study.group, Histology)]
#   study.group Histology median.mtCN.N median.mtCN.T mean.mtCN.N mean.mtCN.T
#1:         IWK      pRCC      474.8389      851.7313    474.8389   1220.4137
#2:         IWK     ccRCC      397.2651      469.7088    397.2651    491.5200
#3:         IWK   ACD-RCC      116.5452     1049.1144    116.5452   1083.3385
#4:         IWK    ccpRCC      142.9048      654.8807    142.9048    654.8807
#5:         IWK     chRCC      558.4279     3132.7198    558.4279   3132.7198
#6:        TCMA   KICH-US            NA      491.3415         NaN    624.7170
#7:        TCMA   KIRC-US            NA      340.1201         NaN    337.0147
#8:        TCMA   KIRP-US            NA      242.3362         NaN    353.4045
#9:        TCMA   RECA-EU            NA      323.4332         NaN    356.1543


ESRD.TCMA.var.ct.summary <- rbind(
		IWK.age.max.AF.T.summary[,list(sample_id=tumor.sample.id, Histology, age, var.ct)],
		TCMA.SNV.INDEL.age.summary)

ESRD.TCMA.var.ct.summary[,list(sample.ct=.N, mean.var.ct=mean(var.ct), median.var.ct=median(var.ct)), by=list(Histology)]
##     Histology sample.ct mean.var.ct median.var.ct
## 1:      pRCC         5    4.000000           3.0
## 2:     ccRCC        18    4.500000           4.0
## 3:   ACD-RCC         6    5.666667           6.0
## 4:    ccpRCC         2    1.500000           1.5
## 5:     chRCC         2    7.000000           7.0
## 6:   KICH-US        44    6.227273           5.5
## 7:   KIRC-US        31    4.419355           4.0
## 8:   KIRP-US        31    7.387097           7.0
## 9:   RECA-EU        72    4.430556           4.0

comparison.list <- list(
		c("ccRCC", "chRCC"),
		c("ccRCC", "ACD-RCC"),
		c("ccRCC", "pRCC"),
		c("ccRCC", "ccpRCC"),
		c("ccRCC", "chRCC"),
		c("ACD-RCC", "pRCC"),
		c("ACD-RCC", "ccpRCC"),
		c("ACD-RCC", "chRCC"),
		c("pRCC", "ccpRCC"),
		c("pRCC", "chRCC"),
		c("ccpRCC", "chRCC"),
		c('ccRCC', 'KIRC-US'),
		c('ccRCC', 'RECA-EU'),
		c("ACD-RCC", "KIRC-US"),
		c("ACD-RCC", "RECA-EU"),
		c("ACD-RCC", 'KIRP-US'),
		c("ACD-RCC", 'KICH-US'),
		c('pRCC', 'KIRP-US'),
		c('chRCC', 'KICH-US'))

Mito.CN.group.comparison.p.values.ls <- lapply(
		X = comparison.list,
		FUN = function( curr.comparison ){
			curr.group1 <- curr.comparison[[1]]
			curr.group2 <- curr.comparison[[2]]
			
			curr.wilcox.exact <- wilcox.test(
					x = IWK.germline.variant.mtCN.dt.long[Sample.type=='Tumor' & Histology==curr.group1,]$Mito.CN,
					y = IWK.germline.variant.mtCN.dt.long[Sample.type=='Tumor' & Histology==curr.group2,]$Mito.CN,
					exact=TRUE)
			
			curr.wilcox.not.exact <- wilcox.test(
					x = IWK.germline.variant.mtCN.dt.long[Sample.type=='Tumor' & Histology==curr.group1,]$Mito.CN,
					y = IWK.germline.variant.mtCN.dt.long[Sample.type=='Tumor' & Histology==curr.group2,]$Mito.CN,
					exact=FALSE)
			
			data.table(
					group1 = curr.group1,
					group2 = curr.group2,
					p.value.exact = curr.wilcox.exact$p.value,
					p.value.not.exact = curr.wilcox.not.exact$p.value)
		})
Mito.CN.group.comparison.p.values <- rbindlist(Mito.CN.group.comparison.p.values.ls)

#> Mito.CN.group.comparison.p.values[,list(group1, group2, p.value=formatC(p.value.exact, format='f', digits=6))]
#	group1  group2  p.value
#1:   ccRCC   chRCC 0.010526
#2:   ccRCC ACD-RCC 0.000951
#3:   ccRCC    pRCC 0.024250
#4:   ccRCC  ccpRCC 0.210526
#5:   ccRCC   chRCC 0.010526
#6: ACD-RCC    pRCC 0.536797
#7: ACD-RCC  ccpRCC 0.142857
#8: ACD-RCC   chRCC 0.071429
#9:    pRCC  ccpRCC 0.571429
#10:    pRCC   chRCC 0.380952
#11:  ccpRCC   chRCC 0.333333
#12:   ccRCC KIRC-US 0.067049
#13:   ccRCC RECA-EU 0.075310
#14: ACD-RCC KIRC-US 0.000005
#15: ACD-RCC RECA-EU 0.000000
#16: ACD-RCC KIRP-US 0.000256
#17: ACD-RCC KICH-US 0.011558
#18:    pRCC KIRP-US 0.008494
#19:   chRCC KICH-US 0.002114

Mito.CN.group.comparison.p.values[p.value.exact<0.05,list(group1, group2, p.value=formatC(p.value.exact, format='f', digits=6))]
##     group1  group2  p.value
## 1:   ccRCC   chRCC 0.010526
## 2:   ccRCC ACD-RCC 0.000951
## 3:   ccRCC    pRCC 0.024250
## 4:   ccRCC   chRCC 0.010526
## 5: ACD-RCC KIRC-US 0.000005
## 6: ACD-RCC RECA-EU 0.000000
## 7: ACD-RCC KIRP-US 0.000256
## 8: ACD-RCC KICH-US 0.011558
## 9:    pRCC KIRP-US 0.008494
## 10:   chRCC KICH-US 0.002114

var.ct.group.comparison.p.values.ls <- lapply(
	X = comparison.list,
	FUN = function( curr.comparison ){
		curr.group1 <- curr.comparison[[1]]
		curr.group2 <- curr.comparison[[2]]
		
		curr.wilcox.exact <- wilcox.test(
				x = ESRD.TCMA.var.ct.summary[Histology==curr.group1,]$var.ct,
				y = ESRD.TCMA.var.ct.summary[Histology==curr.group2,]$var.ct,
				exact=TRUE)
		
		curr.wilcox.not.exact <- wilcox.test(
				x = ESRD.TCMA.var.ct.summary[Histology==curr.group1,]$var.ct,
				y = ESRD.TCMA.var.ct.summary[Histology==curr.group2,]$var.ct,
				exact=FALSE)
		
		data.table(
				group1 = curr.group1,
				group2 = curr.group2,
				p.value.exact = curr.wilcox.exact$p.value,
				p.value.not.exact = curr.wilcox.not.exact$p.value)
	})
var.ct.group.comparison.p.values <- rbindlist(var.ct.group.comparison.p.values.ls)

#> var.ct.group.comparison.p.values[p.value<0.05,]
#group1  group2    p.value
#1:  ccRCC  ccpRCC 0.03436913
#2:   pRCC KIRP-US 0.04075758

#> var.ct.group.comparison.p.values
#	group1  group2     p.value
#1:  ccRCC  ccpRCC 0.034369132
#2:   pRCC KIRP-US 0.002844281
#:  chRCC KICH-US 0.935143317

Histology.color.ls <- Histology.colors.ls[1:5]


boxplot.Mito.CN.by.IWK.type1 <- ggboxplot(
				droplevels(IWK.germline.variant.mtCN.dt.long[
						study.group=='IWK',
						list(tumor.sample.id, Sample.type, Histology, Mito.CN)]),
				x = "Sample.type", y = "Mito.CN", 
				outlier.shape = NA,
				title = 'A',
				color = "Histology",
#				order = c(
#					paste('N-', c('ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'ACD-RCC'), sep=""),
#					paste('T-', c('ccRCC', 'KIRC-US', 'RECA-EU', 'pRCC', 'KIRP-US', 'ccpRCC', 'chRCC', 'KICH-US', 'ACD-RCC'), sep="")),
				add = "jitter",
				ylab = "mtCN") +
#		stat_compare_means(method="t.test", paired=TRUE, label.y=2500, label.x=1, size=3) +
		stat_compare_means(method="wilcox.test", paired=TRUE, label.y=2500, label.x=1, size=3) +
#		stat_compare_means(comparisons = my_comparisons, method="t.test") +
#		stat_pvalue_manual(type1.aov.dt, size=2.5, vjust= -0.4) +
#		rremove('legend') +
#		rremove('legend.title') +
#		rotate_x_text(45) +
		scale_colour_manual(values = Histology.colors.ls) +
#		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.5,"line"),
				legend.text = element_text(size=6),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

Mito.CN.comparison.idxs <- which(Mito.CN.group.comparison.p.values$p.value.exact<0.05)

Mito.CN.comparisons <- lapply(
		X = Mito.CN.comparison.idxs,
		FUN = function(curr.row){
			unlist(Mito.CN.group.comparison.p.values[curr.row,list(group1, group2)])
		})

boxplot.Mito.CN.by.IWK.type2 <- ggboxplot(
				IWK.germline.variant.mtCN.dt.long[
					Sample.type=='Tumor',
					list(tumor.sample.id, Sample.type, Histology, Mito.CN)],
				x = "Histology", y = "Mito.CN", 
				outlier.shape = NA,
				title = 'B',
				color = "Histology",
				order = c('ccRCC', 'KIRC-US', 'RECA-EU', 'ACD-RCC', 'pRCC', 'KIRP-US', 'ccpRCC', 'chRCC', 'KICH-US'),
				add = "jitter",
				ylab = "mtCN", xlab = "RCC sub-type") +
#		stat_compare_means(method="wilcox.test", paired=TRUE, label.y=2500, label.x=1, size=3) +
		stat_compare_means(comparisons = Mito.CN.comparisons, method="wilcox.test", exact=TRUE, size=3, vjust=0.4,
			symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
#		stat_compare_means(comparisons = mito.CN.type2.comparisons, method="wilcox.test", exact=FALSE, size=3) +
		rremove('legend') +
#		rremove('legend.title') +
		rotate_x_text(45) +
		scale_colour_manual(values = Histology.colors.ls) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


#IWK.mtDNA_copies.vs.years.dialysis.N <- ggscatter(
#				droplevels(IWK.germline.variant.mtCN.dt[study.group=='IWK' & tumor.sample.id!='IWK031_T',]),
#				x = "years.dialysis", y = "mtDNA.copy.est.N",
#				title = 'C',
#				color = "Histology",
#				legend.nrow = 2,
#				add = "reg.line",  # Add regression line
#				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
##				conf.int = TRUE, # Add confidence interval
#				xlab = "Dialysis period", ylab = "Normal mtCN") +
#		stat_cor(method = "pearson", label.x = 5, label.y = 600, size=2.5) +
#		rremove("legend") +
#		rremove("legend.title") +
#		theme( plot.title = element_text(hjust = 0, size=10) ) +
#		guides(col = guide_legend(nrow = 2)) +
#		scale_color_manual(values = Histology.colors.ls[1:5]) +
#		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
#		theme(title = element_text(size=10, face='bold', color='black'),
#				legend.key.size = unit(0.6,"line"),
#				legend.text = element_text(size=8),
#				axis.title = element_text(size=8), axis.text = element_text(size=8))


IWK.mtDNA_copies.vs.years.dialysis.T <- ggscatter(
				droplevels(IWK.germline.variant.mtCN.dt[study.group=='IWK',]),
				x = "years.dialysis", y = "mtDNA.copy.est.T",
				title = 'C',
				color = "Histology",
				add = "reg.line",  # Add regression line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#				conf.int = TRUE, # Add confidence interval
				xlab = "Dialysis period (years)", ylab = "Tumor mtCN") +
		stat_cor(method = "pearson", label.x = 2, label.y = 1400, size=2.5) +
		rremove("legend") +
		rremove("legend.title") +
		scale_color_manual(values = Histology.colors.ls[1:5]) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=10, face='bold', color='black'),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


var.ct.comparison.idxs <- which(var.ct.group.comparison.p.values$p.value.not.exact<0.05)

var.ct.comparisons <- lapply(
	X = var.ct.comparison.idxs,
	FUN = function(curr.row){
		unlist(var.ct.group.comparison.p.values[curr.row,list(group1, group2)])
	})

ESRD.TCMA.var.ct.RCC.type <- ggboxplot(ESRD.TCMA.var.ct.summary,
				x = "Histology", y = "var.ct", 
				outlier.shape = NA,
				title = 'D',
				color = "Histology",
				order = c('ccRCC', 'KIRC-US', 'RECA-EU', 'ACD-RCC', 'pRCC', 'KIRP-US', 'ccpRCC', 'chRCC', 'KICH-US'),
				add = "jitter",
#				ylim = c(0,13.0),
				ylab = "mtDNA var. count", xlab = "RCC sub-type") +
		stat_compare_means(comparisons = var.ct.comparisons, method="wilcox.test", exact=FALSE, size=3, vjust=0.4,
				symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
#		stat_compare_means(comparisons = var.ct.comparisons, method="wilcox.test", exact=FALSE, size=3, vjust=0.4) +
		rremove('legend') +
#		rremove('legend.title') +
		rotate_x_text(45) +
		scale_colour_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

mtDNA.p <- ggarrange(
		ggarrange(
			boxplot.Mito.CN.by.IWK.type1 + 
					rremove("xlab") + 
					theme(plot.title = element_text(vjust = -15)) +
					rremove('legend.title') +
					theme(legend.margin = margin(0, 0.1, 0.3, 0.1, "cm")) +
					theme(plot.margin = margin(0.1, 0.1, 0.2, 0.15, "cm")) +
					theme(axis.title.y = element_text(margin = margin(r = 10))) +
					theme(axis.title.x = element_text(margin = margin(r = 10))) +
					guides(col = guide_legend(nrow = 3)),
			IWK.mtDNA_copies.vs.years.dialysis.T +
				theme(plot.margin = margin(1.5, 0.1, 0.2, 0.15, "cm")) +
				theme(axis.title.x = element_text(margin = margin(r = 10))),
			boxplot.Mito.CN.by.IWK.type2 + rremove('legend') + 
					rremove('legend') + 
					theme(plot.margin = margin(0.2, 0.1, 0.1, 0.15, "cm")) +
					theme(axis.title.y = element_text(margin = margin(r = 10))) + 
					theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
							panel.background=element_blank(), panel.border=element_blank()),
			ESRD.TCMA.var.ct.RCC.type + 
					rremove('legend') + 
					theme(plot.margin = margin(0.2, 0.1, 0.1, 0.15, "cm")) +
					theme(axis.title.y = element_text(margin = margin(r = 10))),
			ncol = 2, nrow=2,
			heights = c(1.1,1),
			widths = c(1,1)))

fig.1col.width <- 85/25.4
fig.1.5col.width <- 114/25.4
fig.2col.width <- 174/25.4

fig.height.to.width <- 0.85

ggexport(mtDNA.p, width=fig.1.5col.width, height=fig.height.to.width*fig.1.5col.width, filename=file.path(run.dir, 'result/mtDNA/figures/mtDNA.Figure4.1.5col.pdf'))
ggexport(mtDNA.p, width=fig.2col.width, height=fig.height.to.width*fig.2col.width, filename=file.path(run.dir, 'result/mtDNA/figures/mtDNA.Figure4.2col.pdf'))
#ggexport(mtDNA.p, width=6.4, height=5.5, filename=file.path(run.dir, 'result/variant_summaries/mtDNA.Figure4.pdf'))

save(list=c('ESRD.TCMA.var.ct.summary', 'IWK.germline.variant.mtCN.dt.long', 'IWK.germline.variant.mtCN.dt',
	'Mito.CN.group.comparison.p.values', 'var.ct.group.comparison.p.values'),
	file = file.path(run.dir, 'result/mtDNA/mtDNA_analysis_summary.Rdata'))
