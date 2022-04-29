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

sig.denovo.stats.ls <- lapply(
		X = sig.type.ls,
		FUN = function(curr.type){
			curr.class <- sig.class.ls[[curr.type]]
			dt.long <- sig.denovo.dt.long.filt.ls[[curr.type]][['dt.long']]
			dt.long.filt <- sig.denovo.dt.long.filt.ls[[curr.type]][['dt.long.group.sig.class.filt']]
			
			pdf( onefile=TRUE, width=10, height=8, file = file.path( fig.dir, paste(curr.type, '.denovo.signature.mutation.summary.AllSamples.pdf', sep="")))
			
			op <- par(mai=c(2, 1, 0.3, 0.3))
			bp1 <- boxplot( sig.mut.ct ~ cat.annotation + sig.class, data=dt.long, outline=FALSE, las=3, xlab='', ylab='Signature mutation (counts)', main=paste(curr.type, ' denovo signatures', sep="") )
			

			op <- par(mai=c(2, 1, 0.3, 0.3))
			bp2 <- boxplot( sig.mut.propn ~ cat.annotation + sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutations (Propn. of total)', main=paste(curr.type, ' denovo signatures (median ct. > 0)', sep="") )

			op <- par(mai=c(2, 1, 0.3, 0.3))
			bp3 <- boxplot( sig.mut.propn ~ sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutations (Propn. of total)', main=paste(curr.type, ' denovo signatures (median ct. > 0)', sep="") )

			sig.classes.ls <- unique(dt.long.filt$sig.class)
			
			dt.long.filt[,Purity:=3*purity]
			dt.long.filt[,Signature:=factor(sig.class, ordered=FALSE)]
			
			ggp1 <- ggplot( dt.long.filt , aes(x=num.annotation, y=sig.mut.ct, color= Signature  )) + 
					geom_point( aes(size = Purity) ) +  
					facet_wrap(~sig.class , dir="v", scales = "free_y") +
					stat_cor(method = "pearson", label.x.npc = 0.65, label.y.npc = 0.85, colour='black') +
#					theme(legend.position="none") +
					labs(x = "Years of dialysis") +
					labs(y = "Signature mutations (count)") +
					labs(title = paste("De novo ", curr.type, " signature", sep=""))
			print(ggp1)

			ggp2 <- ggplot( dt.long.filt , aes(x=num.annotation, y=sig.mut.propn, color=Signature )) + 
					geom_point( aes(size = Purity) ) + 
					facet_wrap(~sig.class , dir="v", scales = "free_y")  +
					stat_cor(method = "pearson", label.x.npc = 0.65, label.y.npc = 0.85, colour='black') +
#					theme(legend.position="none") +
					labs(x = "Years of dialysis") +
					labs(y = "Signature mutations (Propn. of total mutations)") +
					labs(title = paste("De novo ", curr.type, " signature", sep=""))
			print(ggp2)

			sig.mut.ct.stats.ls <- lapply(
				X = sig.classes.ls,
				FUN = function(curr.sig.class){
					curr.dt <- dt.long.filt[sig.class==curr.sig.class,]
					curr.cor <- cor.test(x=curr.dt$num.annotation, y=curr.dt$sig.mut.ct)
					
					op <- par(mai=c(2, 1, 0.3, 0.3))
					plot( sig.mut.ct ~ num.annotation,
							data = curr.dt,
							pch = 16,
							cex = 3*curr.dt$purity,
							xlab=curr.numerical.annotation.label,
							ylab='Signature mutations (count)',
							main=paste("De novo ", curr.type, " signature: ", curr.sig.class, ' (median ct. > 0)', sep="") )
					cor_text(0.7, 0.8, paste("r = ", formatC(curr.cor$estimate, digits=5, format='f'), ", P-value = ", formatC(curr.cor$p.value, digits=5, format='f'), sep=""))
					purity_legend(curr.x.propn, curr.y.propn)
					
					data.table(
							sig.class=curr.sig.class,
							p.value = curr.cor$p.value,
							r = curr.cor$estimate)
				})
		
			sig.mut.propn.stats.ls <- lapply(
				X = sig.classes.ls,
				FUN = function(curr.sig.class){
					curr.dt <- dt.long.filt[sig.class==curr.sig.class,]
					curr.cor <- cor.test(x=curr.dt$num.annotation, y=curr.dt$sig.mut.propn)
					
					op <- par(mai=c(2, 1, 0.3, 0.3))
					plot( sig.mut.propn ~ num.annotation,
							data = curr.dt,
							pch = 16,
							cex = 3*curr.dt$purity,
							xlab=curr.numerical.annotation.label,
							ylab='Signature mutations (Propn. of total)',
							main=paste("De novo ", curr.type, " signature: ", curr.sig.class, ' (median ct. > 0)', sep="") )
					cor_text(0.7, 0.8, paste("r = ", formatC(curr.cor$estimate, digits=5, format='f'), ", P-value = ", formatC(curr.cor$p.value, digits=5, format='f'), sep=""))
					purity_legend(curr.x.propn, curr.y.propn)
					
					data.table(
							sig.class=curr.sig.class,
							p.value = curr.cor$p.value,
							r = curr.cor$estimate)
				})
			dev.off()
			
			list(
				sig.mut.ct.stats = rbindlist(sig.mut.ct.stats.ls),
				sig.mut.propn.stats = rbindlist(sig.mut.propn.stats.ls))
		})


sig.COSMIC.stats.ls <- lapply(
	X = sig.type.ls,
	FUN = function(curr.type){
		curr.class <- sig.class.ls[[curr.type]]
		
		dt.long <- sig.dt.long.filt.ls[[curr.type]][['dt.long']]
		dt.long.filt <- sig.dt.long.filt.ls[[curr.type]][['dt.long.group.sig.class.filt']]
		
		pdf( onefile=TRUE, width=10, height=8, file = file.path( fig.dir, paste(curr.type, '.signature.mutation.summary.AllSamples.pdf', sep="")))
		
		op <- par(mai=c(3.0, 1, 0.3, 0.3))
		bp1 <- boxplot( sig.mut.ct ~ cat.annotation + sig.class, data=dt.long, outline=FALSE, las=3, xlab='', ylab='Signature mutation (counts)', main=paste(curr.type, ' decomposed COSMIC signatures', sep="") )
		
		op <- par(mai=c(3.0, 1, 0.3, 0.3))
		bp1 <- boxplot( sig.mut.ct ~ cat.annotation + sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutation (counts)', main=paste(curr.type, ' decomposed COSMIC signatures (median ct. > 0)', sep="") )
		
		op <- par(mai=c(3.0, 1, 0.3, 0.3))
		bp2 <- boxplot( sig.mut.propn ~ cat.annotation + sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutations (Propn. of total)', main=paste(curr.type, ' decomposed COSMIC signatures (median ct. > 0)', sep="") )

		op <- par(mai=c(3.0, 1, 0.3, 0.3))
		bp3 <- boxplot( sig.mut.ct ~ sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutation (counts)', main=paste(curr.type, ' decomposed COSMIC signatures (median ct. > 0)', sep="") )
		
		op <- par(mai=c(3.0, 1, 0.3, 0.3))
		bp4 <- boxplot( sig.mut.propn ~ sig.class, data=dt.long.filt, outline=FALSE, las=3, xlab='', ylab='Signature mutations (Propn. of total)', main=paste(curr.type, ' decomposed COSMIC signatures (median ct. > 0)', sep="") )
		
		
		sig.classes.ls <- unique(dt.long.filt$sig.class)

		dt.long.filt[,Purity:=purity]
		dt.long.filt[,Signature:=factor(sig.class, ordered=FALSE)]

		ggp1 <- ggplot( dt.long.filt, aes(x=num.annotation, y=sig.mut.ct, color=Signature )) + 
				geom_point( aes(size = Purity) ) +  
				facet_wrap(~sig.class , dir="v", scales = "free_y") +
				stat_cor(method = "pearson", label.x.npc = 0.70, label.y.npc = 0.80, colour='black') +
				theme(legend.position=c(0.75,0.2), legend.direction='horizontal') +
				labs(x = "Years of dialysis") +
				labs(y = "Signature mutations (count)") +
				labs(title = paste("COSMIC 3 signature", sep=""))
		
		print(ggp1)

		ggp2 <- ggplot( dt.long.filt, aes(x=num.annotation, y=sig.mut.propn, color=Signature )) + 
				geom_point( aes(size = Purity) ) + 
				facet_wrap(~sig.class , dir="v", scales = "free_y") +
				stat_cor(method = "pearson", label.x.npc = 0.65, label.y.npc = 0.70, colour='black') +
				theme(legend.position=c(0.75,0.2), legend.direction='horizontal') +
				labs(x = "Years of dialysis") +
				labs(y = "Signature mutations (Propn. of total mutations)") +
				labs(title = paste("COSMIC 3 signature", sep=""))
		
		print(ggp2)

		ggp3 <- ggplot( dt.long.filt, aes(x=num.annotation, y=sig.mut.propn, color=Signature )) + 
				geom_point( size = 2) + 
				facet_wrap(~sig.class + cat.annotation , dir="v", scales = "free_y") +
				theme(legend.position='none') +
				labs(x = "Years of dialysis") +
				labs(y = "Signature mutations (Propn. of total mutations)") +
				labs(title = paste("COSMIC 3 signature", sep=""))
		
		print(ggp3)
		
		sig.mut.ct.stats.ls <- lapply(
				X = sig.classes.ls,
				FUN = function(curr.sig.class){
					curr.dt <- dt.long.filt[sig.class==curr.sig.class,]
					curr.cor <- cor.test(x=curr.dt$num.annotation, y=curr.dt$sig.mut.ct)
					
					op <- par(mai=c(2, 1, 0.3, 0.3))
					plot( sig.mut.ct ~ num.annotation,
							data = curr.dt,
							pch = 16,
							cex = 3*curr.dt$purity,
							xlab=curr.numerical.annotation.label,
							ylab='Signature mutations (count)',
							main=paste("De novo ", curr.type, " signature: ", curr.sig.class, ' (median ct. > 0)', sep="") )
					cor_text(0.7, 0.8, paste("r = ", formatC(curr.cor$estimate, digits=5, format='f'), ", P-value = ", formatC(curr.cor$p.value, digits=5, format='f'), sep=""))
					purity_legend(curr.x.propn, curr.y.propn)
					
					data.table(
							sig.class=curr.sig.class,
							p.value = curr.cor$p.value,
							r = curr.cor$estimate)
				})
		
		sig.mut.propn.stats.ls <- lapply(
				X = sig.classes.ls,
				FUN = function(curr.sig.class){
					curr.dt <- dt.long.filt[sig.class==curr.sig.class,]
					curr.cor <- cor.test(x=curr.dt$num.annotation, y=curr.dt$sig.mut.propn)

					op <- par(mai=c(2, 1, 0.3, 0.3))
					plot( sig.mut.propn ~ num.annotation,
							data = curr.dt,
							pch = 16,
							cex = 3*curr.dt$purity,
							xlab=curr.numerical.annotation.label,
							ylab='Signature mutations (Propn. of total)',
							main=paste("De novo ", curr.type, " signature: ", curr.sig.class, ' (median ct. > 0)', sep="") )
					cor_text(0.7, 0.8, paste("r = ", formatC(curr.cor$estimate, digits=5, format='f'), ", P-value = ", formatC(curr.cor$p.value, digits=5, format='f'), sep=""))
					purity_legend(curr.x.propn, curr.y.propn)
					
					data.table(
							sig.class=curr.sig.class,
							p.value = curr.cor$p.value,
							r = curr.cor$estimate)
				})

		dev.off()

		list(
			sig.mut.ct.stats = rbindlist(sig.mut.ct.stats.ls),
			sig.mut.propn.stats = rbindlist(sig.mut.propn.stats.ls))
	})

save(list=c('sig.denovo.dt.long.filt.ls', 'sig.dt.long.filt.ls', 'sig.denovo.stats.ls', 'sig.COSMIC.stats.ls'),
	file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/modeled.activities.Rdata", sep="")))

fwrite(sig.COSMIC.stats.ls[['SBS']][['sig.mut.propn.stats']],
	file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/mutation_activities.propn.tsv", sep="")), sep='\t')

fwrite(sig.COSMIC.stats.ls[['SBS']][['sig.mut.ct.stats']],
		file = file.path(run.dir,paste("result/SigProfiler/", study.name, "_mutation_signature_analysis/results/mutation_activities.ct.tsv", sep="")), sep='\t')

#> sig.denovo.stats.ls
#$SBS
#$SBS$sig.mut.ct.stats
#sig.class    p.value          r2
#1:    SBS96A 0.64528823 -0.08319989
#2:    SBS96B 0.61951220  0.08972168
#3:    SBS96C 0.06301386 -0.32726897
#
#$SBS$sig.mut.propn.stats
#sig.class    p.value         r2
#1:    SBS96A 0.69723467 -0.0703547
#2:    SBS96B 0.07335685  0.3158485
#3:    SBS96C 0.09342829 -0.2968515
#
#
#$DBS
#$DBS$sig.mut.ct.stats
#sig.class   p.value        r2
#1:    DBS78A 0.5742102 -0.101467
#
#$DBS$sig.mut.propn.stats
#sig.class p.value r2
#1:    DBS78A      NA NA
#
#
#$ID
#$ID$sig.mut.ct.stats
#sig.class   p.value         r2
#1:     ID83A 0.7515950 -0.0572648
#2:     ID83B 0.2418072 -0.2095629
#
#$ID$sig.mut.propn.stats
#sig.class   p.value          r2
#1:     ID83A 0.9523934 -0.01080920
#2:     ID83B 0.9371606  0.01427382



## > sig.COSMIC.stats.ls
## $SBS
## $SBS$sig.mut.ct.stats
## sig.class     p.value          r2
## 1:      SBS1 0.159672764 -0.25051930
## 2:     SBS12 0.147531808 -0.25777336
## 3:     SBS23 0.002274235  0.51286253
## 4:     SBS40 0.368357679 -0.16179554
## 5:      SBS5 0.716951107 -0.06556979
## 
## $SBS$sig.mut.propn.stats
## sig.class      p.value          r2
## 1:      SBS1 3.884591e-01  0.15519876
## 2:     SBS12 2.494971e-01 -0.20625112
## 3:     SBS23 5.107695e-05  0.64481678
## 4:     SBS40 8.593133e-01  0.03208419
## 5:      SBS5 9.532359e-01  0.01061769
## 
## 
## $DBS
## $DBS$sig.mut.ct.stats
## sig.class   p.value           r2
## 1:      DBS2 0.9669071 -0.007511637
## 2:      DBS4 0.6635781 -0.078634610
## 3:      DBS9 0.6257270 -0.088139294
## 
## $DBS$sig.mut.propn.stats
## sig.class   p.value          r2
## 1:      DBS2 0.1132917  0.28091194
## 2:      DBS4 0.2566658 -0.20322438
## 3:      DBS9 0.7502623 -0.05758205
## 
## 
## $ID
## $ID$sig.mut.ct.stats
## sig.class   p.value         r2
## 1:       ID1 0.1457867 -0.2588512
## 2:       ID3 0.8701400 -0.0295932
## 3:       ID5 0.5301308 -0.1133030
## 4:       ID8 0.3943752 -0.1532953
## 
## $ID$sig.mut.propn.stats
## sig.class   p.value           r2
## 1:       ID1 0.8624013  0.031373142
## 2:       ID3 0.7625920  0.054653196
## 3:       ID5 0.8006009 -0.045706100
## 4:       ID8 0.9687958  0.007082695
## 



dt.long.filt <- sig.dt.long.filt.ls[['SBS']][['dt.long.group.sig.class.filt']]
dt.long.filt[,Signature:=factor(sig.class, ordered=FALSE)]

unique(dt.long.filt$Signature)

check.lm <- lm(sig.mut.propn ~ num.annotation + cat.annotation, data=dt.long.filt[sig.class=='SBS23',])
summary(check.lm)

check.lm <- lm(sig.mut.propn ~ num.annotation + cat.annotation, data=dt.long.filt)
summary(check.lm)


