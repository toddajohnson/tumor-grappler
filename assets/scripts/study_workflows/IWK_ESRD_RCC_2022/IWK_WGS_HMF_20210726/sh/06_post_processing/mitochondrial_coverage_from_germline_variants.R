library(data.table)
library(parallel)

load("/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/config_files/fastq_file_info.Rdata")
patient.info.dt.IWK <- copy(patient.info.dt)

load("/home/tjohnson/workspace/runs/BHD_WGS_HMF_20210201/config_files/fastq_file_info.Rdata")
patient.info.dt <- rbind(
	patient.info.dt.IWK[,list(study.name='IWK', subject.id, tumor.sample.id, histology, age, gender, years.dialysis)],
	tumor_normal_pair_info.GPL[,list(study.name='BHD', subject.id, tumor.sample.id=sample.id.T,
	histology=ifelse(tumor.class=='ccRCC', 'Clear cell', ifelse(tumor.class=='chRCC', 'Chromophobe', tumor.class)), age, gender, years.dialysis=NA)])

base.run.dir <- "/home/tjohnson/workspace/runs"

run.directories.ls <- c(
	IWK = "IWK_WGS_HMF_20210726",
	BHD = "BHD_WGS_HMF_20210201",
	BTC1 = "BTC_WGS_GRIDSS_PURPLE_LINX_CHORD_20210315",
	BTC2 = "BTC_WGS_GRIDSS_PURPLE_LINX_CHORD_20210702",
	Esophageal_cancer = "Esophageal_cancer_WGS_GPL_20211117",
	KakimiLab = "KakimiLab_WGS_HMF_20210426",
	KeioOrganoid1 = "KeioOrganoid_WGS_GRIDSS_PURPLE_LINX_CHORD_20210315",
	KeioOrganoid2 = "KeioOrganoid_WGS_GRIDSS_PURPLE_LINX_CHORD_20210607",
	KeioOrganoid3 = "KeioOrganoid_WGS_GRIDSS_PURPLE_LINX_CHORD_20211108")

study.names.ls <- names(run.directories.ls)
names(study.names.ls) <- study.names.ls

study.groups <- c('IWK', 'BHD', 'BTC', 'BTC', 'Esophageal_cancer', 'KakimiLab', 'KeioOrganoid', 'KeioOrganoid', 'KeioOrganoid')
names(study.groups) <- study.names.ls

cancer.types.ls <- c('RCC', 'RCC', 'BTC', 'BTC', 'EC', 'UNK', 'UNK', 'UNK', 'UNK')
names(cancer.types.ls) <- study.names.ls

study.colors.ls <- c('red', 'light green', 'blue', 'blue', 'orange', 'brown', 'yellow', 'yellow', 'yellow')
names(study.colors.ls) <- study.names.ls

target.chroms.ls <- paste('chr', c(as.character(1:22), 'M'), sep="")

germline.variant.coverage.ls <- mclapply(
	X = study.names.ls,
	FUN = function( curr.study.name ){
		curr.run.dir <- run.directories.ls[[curr.study.name]]
		curr.full.run.dir <- file.path( base.run.dir, curr.run.dir)
		print(paste('Running ', curr.study.name, sep=""))
		source( file.path( curr.full.run.dir, 'config_files/common_config.R'), local=TRUE )
		load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/merged.germline.chrM.autosomes.depth_summary.Rdata", sep=""))
		
		germline.variant.depth.summary.dt[,study.name:=curr.study.name]
		germline.variant.depth.summary.dt
	}, mc.cores=3)

germline.variant.coverage.dt <- rbindlist(germline.variant.coverage.ls)

germline.variant.coverage.dt[,cancer.type:=cancer.types.ls[study.name]]
germline.variant.coverage.dt[,study.group:=study.groups[study.name]]

germline.variant.coverage.dt <- dcast(
 data = germline.variant.coverage.dt,
 formula = cancer.type + study.group + study.name + subject.id + tumor.sample.id + purity + ploidy ~ chrom.class,
 value.var = c('mean.DP.N', 'mean.DP.T', 'median.DP.N', 'median.DP.T', 'variant.ct'))
 
#germline.variant.coverage.dt[,mtDNA.copy.est.N:=mean.DP.N_chrM/(mean.DP.N_autosome/2)]
#germline.variant.coverage.dt[,mtDNA.copy.est.T:=mean.DP.T_chrM/(mean.DP.T_autosome/ploidy)]
#Scaling formula for mt CN from Yuan et al.
germline.variant.coverage.dt[,mtDNA.copy.est.N:=(mean.DP.N_chrM/mean.DP.N_autosome)*((0*2) + ((1-0)*2))]
germline.variant.coverage.dt[,mtDNA.copy.est.T:=(mean.DP.T_chrM/mean.DP.T_autosome)*((purity*ploidy) + ((1-purity)*2))]

germline.variant.coverage.dt <- merge(
	x = germline.variant.coverage.dt,
	y = patient.info.dt[,list(study.name, subject.id, tumor.sample.id, histology, age, gender, years.dialysis)],
	by = c('study.name', 'subject.id', 'tumor.sample.id'),
	all.x = TRUE)

# ploidy = 3
# AUTO.T.DP = 40
# AUTO.T.DP.1 = 40/3 = 13.3
# MT.T.DP.obs = 3000 
# MT.T.DP.est = 13.3*MT.T.CN


germline.variant.coverage.dt[,histology:=ifelse(is.na(histology), study.name, histology)]

histology.ls <- c('Clear cell', 'ACD-RCC', 'Chromophobe', 'Papillary', 'Clear cell papillary', 'HOCT', 'Esoph. cancer',
	'Unclass.', 'Unclass.', 'Unclass.', 'Unclass.', 'Unclass.', 'Unclass.', 'Unclass.', 'Unclass.')
	
names(histology.ls) <- c('Clear cell', 'ACD-RCC', 'Chromophobe', 'Papillary', 'Clear cell papillary', 'HOCT', 'Esophageal_cancer',
	'Unclassified',
	'BTC1', 'BTC2', 'KakimiLab', 'KeioOrganoid1', 'KeioOrganoid2', 'KeioOrganoid3')

germline.variant.coverage.dt[,histology:=histology.ls[histology]]

histology.shapes <- c(25, 23, 22, 24, 24, 14, 7, rep(21, length=8))
names(histology.shapes) <- histology.ls

germline.variant.coverage.dt[,study.color:=study.colors.ls[study.name]]
germline.variant.coverage.dt[,histology.shape:=histology.shapes[histology]]

germline.mtDNA.copy.summary <- germline.variant.coverage.dt[,
	list(mean.mtDNA.copy.est.N=mean(mtDNA.copy.est.N), mean.mtDNA.copy.est.T=mean(mtDNA.copy.est.T), sample.ct=.N),
	by = list(cancer.type, study.group, study.name, histology, histology.shape, study.color)]

germline.mtDNA.copy.summary <- germline.mtDNA.copy.summary[order(cancer.type, -mean.mtDNA.copy.est.T)]

save( list=c('germline.variant.coverage.dt', 'germline.mtDNA.copy.summary'),
	file = paste(base.run.dir, '/IWK_WGS_HMF_20210726/result_summaries/adhoc/germline.variant.coverage.Rdata', sep=""))

#	cancer.type       study.group        study.name            histology histology.shape study.color mean.mtDNA.copy.est.N mean.mtDNA.copy.est.T sample.ct
#1:         BTC               BTC              BTC2             Unclass.              21        blue              341.6226              616.5480        12
#2:         BTC               BTC              BTC1             Unclass.              21        blue              274.8707              602.6313        43
#3:          EC Esophageal_cancer Esophageal_cancer        Esoph. cancer               7      orange              594.8937              540.6654        33
#4:         RCC               BHD               BHD             Unclass.              21 light green             1958.4480             4462.7771         2
#5:         RCC               BHD               BHD          Chromophobe              22 light green             1248.5552             2870.2510         9
#6:         RCC               IWK               IWK          Chromophobe              22         red              717.4633             2427.4341         3
#7:         RCC               BHD               BHD                 HOCT              14 light green             1691.5302             1925.5202         4
#8:         RCC               IWK               IWK             Unclass.              21         red             1013.4038             1295.1978         1
#9:         RCC               BHD               BHD           Clear cell              25 light green             1418.9169             1246.3087         3
#10:         RCC               IWK               IWK            Papillary              24         red              638.4667             1220.4137         5
#11:         RCC               IWK               IWK              ACD-RCC              23         red              255.4548              896.3804         9
#12:         RCC               IWK               IWK Clear cell papillary              24         red              280.1942              654.8807         2
#13:         RCC               IWK               IWK           Clear cell              25         red              236.0333              491.5200        18
#14:         UNK      KeioOrganoid     KeioOrganoid1             Unclass.              21      yellow              187.1975             1575.4783        39
#15:         UNK      KeioOrganoid     KeioOrganoid3             Unclass.              21      yellow              122.0607             1231.7082         3
#16:         UNK      KeioOrganoid     KeioOrganoid2             Unclass.              21      yellow              320.2522              941.8473        23
#17:         UNK         KakimiLab         KakimiLab             Unclass.              21       brown              511.6537              464.0912        10

# Using PURPLE germline effect annotated variants
#    cancer.type  study.group    study.name                   histology histology.shape study.color mean.mtDNA.copy.est.N mean.mtDNA.copy.est.T sample.ct
# 1:         BTC          BTC          BTC1                    Unclass.              21        blue              309.8748              682.1997        43
# 2:         BTC          BTC          BTC2                    Unclass.              21        blue              309.2813              554.4058        12
# 3:         RCC          BHD           BHD                    Unclass.              21 light green             1850.7901             4476.8019         2
# 4:         RCC          BHD           BHD                 Chromophobe              22 light green             1150.3075             2922.7188         9
# 5:         RCC          IWK           IWK                 Chromophobe              22         red              771.5751             2653.3771         3
# 6:         RCC          BHD           BHD                        HOCT              21 light green             1613.1011             1827.9847         4
# 7:         RCC          IWK           IWK                    Unclass.              21         red             1094.7370             1383.1298         1
# 8:         RCC          IWK           IWK                   Papillary              24         red              631.6944             1301.8292         5
# 9:         RCC          BHD           BHD                  Clear cell              25 light green             1409.7434             1247.5616         3
#10:         RCC          IWK           IWK                     ACD-RCC              23         red              278.0039              997.5111         9
#11:         RCC          IWK           IWK Papillary (Clear papillary)              24         red              306.8673              727.9998         2
#12:         RCC          IWK           IWK                  Clear cell              25         red              255.5955              546.3332        18
#13:         UNK KeioOrganoid KeioOrganoid1                    Unclass.              21      yellow              189.5987             1637.9356        39
#14:         UNK KeioOrganoid KeioOrganoid2                    Unclass.              21      yellow              314.4171              912.8380        23
#15:         UNK    KakimiLab     KakimiLab                    Unclass.              21       brown              482.8851              446.1990        10

# Previous divided by 2 or ploidy
germline.mtDNA.copy.summary
## cancer.type  study.group    study.name                   histology histology.shape study.color mean.mtDNA.copy.est.N mean.mtDNA.copy.est.T sample.ct
## 1:         BTC          BTC          BTC1                    Unclass.              21        blue           309.8748           731.6481        43
## 2:         BTC          BTC          BTC2                    Unclass.              21        blue           309.2813           703.3533        12
## 3:         RCC          BHD           BHD                    Unclass.              21 light green          1850.7901          5469.7010         2
## 4:         RCC          BHD           BHD                 Chromophobe              22 light green          1150.3075          3805.8690         9
## 5:         RCC          IWK           IWK                 Chromophobe              22         red           771.5751          2653.3771         3
## 6:         RCC          BHD           BHD                        HOCT              21 light green          1613.1011          2200.3135         4
## 7:         RCC          BHD           BHD                  Clear cell              25 light green          1409.7434          1632.8608         3
## 8:         RCC          IWK           IWK                   Papillary              24         red           631.6944          1308.6979         5
## 9:         RCC          IWK           IWK                    Unclass.              21         red          1094.7370          1292.1197         1
## 10:         RCC          IWK           IWK                     ACD-RCC              23         red           278.0039          1059.5518         9
## 11:         RCC          IWK           IWK Papillary (Clear papillary)              24         red           306.8673           725.7250         2
## 12:         RCC          IWK           IWK                  Clear cell              25         red           255.5955           543.8028        18
## 13:         UNK KeioOrganoid KeioOrganoid1                    Unclass.              21      yellow           189.5987          1645.7332        39
## 14:         UNK KeioOrganoid KeioOrganoid2                    Unclass.              21      yellow           314.4171           919.2103        23
## 15:         UNK    KakimiLab     KakimiLab                    Unclass.              21       brown           482.8851           517.5104        10

library(ggpubr)
fwrite(germline.mtDNA.copy.summary, file=paste(base.run.dir, '/IWK_WGS_HMF_20210726/result_summaries/adhoc/germline.variant.mtDNA.copy.summary.csv', sep=""))

pdf( file=paste(base.run.dir, '/IWK_WGS_HMF_20210726/result_summaries/adhoc/germline.variant.mtDNA.copy.summary.pdf', sep=""))

plot(mtDNA.copy.est.T ~ mtDNA.copy.est.N,
	data=germline.variant.coverage.dt,
	pch=germline.variant.coverage.dt$histology.shape,
	col = 'black',
	bg=germline.variant.coverage.dt$study.color,
#	xlim=c(0,800), ylim=c(0,1250),
	xlab = 'mtDNA copies in normal', ylab='mtDNA copies in tumor')

legend.dt <- unique(germline.mtDNA.copy.summary[,list(legend.text=paste(study.group, histology, sep=': '), study.color, histology.shape, cancer.type, study.group, histology)])
legend.dt <- legend.dt[order(cancer.type, histology, study.group),]

legend(x = 0, y = 6100,
	legend = legend.dt$legend.text,
	col = 'black',
	pt.bg = legend.dt$study.color,
	cex = 0.7,
    pch = legend.dt$histology.shape, bty='n',
	xpd = NA,
	ncol = 3)

dev.off()


