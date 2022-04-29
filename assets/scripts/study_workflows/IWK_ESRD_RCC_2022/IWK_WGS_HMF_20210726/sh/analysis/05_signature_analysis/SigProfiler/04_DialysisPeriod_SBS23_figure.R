#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=350)
options(width=215)
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
			library(sigminer)
			library(ggplot2)
			library(ggpubr)})

library(stringr)

## arguments that could be added to common_config.R
remove.copy.number.noise.samples <- TRUE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

HOME.dir <- "/home/tjohnson"


sig.type.ls <- c("SBS", "DBS", "ID")
names(sig.type.ls) <- sig.type.ls

sig.class.ls <- c("SBS96", "DBS78", "ID83")
names(sig.class.ls) <- sig.type.ls

curr.categorical.annotation.col <- 'histology'
curr.numerical.annotation.col <- 'years.dialysis'
curr.numerical.annotation.label <- 'Years of dialysis'

base.cols <- c('tumor.sample.id', 'gender', 'purity', 'ploidy', 'cat.annotation', 'num.annotation')

fig.dir <- file.path(run.dir, 'result/SigProfiler/figures')

if (file.exists(fig.dir)!=TRUE){
	dir.create(fig.dir, recursive=TRUE)
}

load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

load(file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/merged.refit.activities.Rdata", sep="")))

purity_legend <- function(x.propn, y.propn){
	x.min = par()$usr[[1]]
	x.max = par()$usr[[2]]
	y.min = par()$usr[[3]]
	y.max = par()$usr[[4]]
	
	x.width <- x.max - x.min
	y.height <- y.max - y.min
	
	x.pos <- x.propn*x.width + x.min
	y.pos <- y.propn*y.height + y.min
	
	legend(
			x = x.pos,
			y = y.pos,
			legend = round(seq(0.1, 1, 0.1),1),
			pch = 16,
			pt.cex = 3*seq(0.1, 1, 0.1),
			horiz = TRUE,
			title = 'Purity',
			xpd = NA)
}


cor_text <- function(x.propn, y.propn, curr.text){
	x.min = par()$usr[[1]]
	x.max = par()$usr[[2]]
	y.min = par()$usr[[3]]
	y.max = par()$usr[[4]]
	
	x.width <- x.max - x.min
	y.height <- y.max - y.min
	
	x.pos <- x.propn*x.width + x.min
	y.pos <- y.propn*y.height + y.min
	
	text(
			x = x.pos,
			y = y.pos,
			labels = curr.text,
			xpd = NA)
}


curr.x.propn <- 0.2
curr.y.propn <- -0.2

curr.type <- 'SBS'
curr.class <- sig.class.ls[[curr.type]]

dt.long <- sig.dt.long.filt.ls[[curr.type]][['dt.long']]
dt.long.filt <- sig.dt.long.filt.ls[[curr.type]][['dt.long.group.sig.class.filt']]



sig.classes.ls <- unique(dt.long.filt$sig.class)

dt.long.filt[,Purity:=purity]
dt.long.filt[,Signature:=factor(sig.class, ordered=FALSE)]

curr.sig.class <- 'SBS23'


p.SBS23.vs.years <- ggscatter(dt.long.filt[sig.class==curr.sig.class,list(Histology=cat.annotation, num.annotation, sig.mut.ct, sig.mut.propn)], x = "num.annotation", y = "sig.mut.propn", 
				color = "Histology",	
				title = "C",
#				legend.nrow = 2,
				add = "reg.line",  # Add regression line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#				conf.int = TRUE, # Add confidence interval
				xlab = "Years of dialysis", ylab = paste("Propn. of mutations in ", curr.sig.class, sep="")) +
#		stat_cor(method = "pearson", label.x = 1, label.y = 0.045, size=2.5) +
		stat_cor(method = "spearman", label.x = 1, label.y = 0.045, size=2.5) +
		rremove("legend") +
		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=9, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

#ggexport(p.SBS23.vs.years,  width=5.0, height=3.5,
#		filename = file.path(run.dir, 'result/SigProfiler/figures/RCC.propn.SBS23.vs.years.dialysis.pdf'))


ct.SBS23.vs.years <- ggscatter(dt.long.filt[sig.class==curr.sig.class,list(Histology=cat.annotation, num.annotation, sig.mut.ct, sig.mut.propn)], x = "num.annotation", y = "sig.mut.ct", 
				color = "Histology",	
				title = "B",
#				legend.nrow = 2,
				add = "reg.line",  # Add regression line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#				conf.int = TRUE, # Add confidence interval
				xlab = "Years of dialysis", ylab = paste("Number of mutations in ", curr.sig.class, sep="")) +
#		stat_cor(method = "pearson", label.x = 1, label.y = 250, size=2.5) +
		stat_cor(method = "spearman", label.x = 1, label.y = 250, size=2.5) +
		rremove("legend") +
		rremove("legend.title") +
		guides(col = guide_legend(nrow = 2)) +
		scale_color_manual(values = Histology.colors.ls) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(title = element_text(size=9, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

#ggexport(ct.SBS23.vs.years,  width=5.0, height=3.5,
#		filename = file.path(run.dir, 'result/SigProfiler/figures/RCC.cts.SBS23.vs.years.dialysis.pdf'))

SBS23.vs.years <- ggarrange(ct.SBS23.vs.years, p.SBS23.vs.years, nrow=1, widths=c(1, 1), legend='bottom', common.legend=TRUE)

ggexport(SBS23.vs.years, width=4.488, height=2.5,
#	filename = file.path(run.dir, 'result/SigProfiler/figures/RCC.SBS23.vs.years.dialysis.pdf'))
	filename = file.path(run.dir, 'result/SigProfiler/figures/RCC.SBS23.vs.years.dialysis.Spearman.pdf'))

clinical.data[,dialysis.period:=cut(years.dialysis, breaks=2, labels=c('Short (<14 years)', 'Long (>14 years)'), ordered=TRUE)]

dt.long.filt <- merge(
	x = clinical.data[,list(tumor.sample.id, years.dialysis, dialysis.period)],
	y = dt.long.filt,
	by = c('tumor.sample.id'))

dt.long.filt.SBS23 <- dt.long.filt[sig.class==curr.sig.class,]
dt.long.filt.SBS23[,SBS23.status:=ifelse(sig.mut.propn>0, 1, 0)]
SBS23.summary.dt <- dt.long.filt.SBS23[,list(SBS23.mut.ct=sum(SBS23.status), total.ct=.N), by=list(dialysis.period)]

#	dialysis.period SBS23.mut.ct total.ct
#1:  Long (>14 years)            5        9
#2: Short (<14 years)            0       24

prop.test(
	x = SBS23.summary.dt$SBS23.mut.ct,
	n = SBS23.summary.dt$total.ct,
	correct = FALSE)

#2-sample test for equality of proportions without continuity correction
#
#data:  SBS23.summary.dt$SBS23.mut.ct out of SBS23.summary.dt$total.ct
#X-squared = 15.714, df = 1, p-value = 7.367e-05
#alternative hypothesis: two.sided
#95 percent confidence interval:
#		0.2309176 0.8801935
#sample estimates:
#		prop 1    prop 2 
#0.5555556 0.0000000 
#
#Warning message:
#In prop.test(x = SBS23.summary.dt$SBS23.mut.ct, n = SBS23.summary.dt$total.ct,  :
#	Chi-squared approximation may be incorrect
