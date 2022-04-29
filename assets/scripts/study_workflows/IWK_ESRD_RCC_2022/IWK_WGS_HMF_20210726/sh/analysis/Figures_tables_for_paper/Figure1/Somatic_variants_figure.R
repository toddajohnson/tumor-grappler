#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=210)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

HOME.dir <- '/home/tjohnson'

gpl_prefix <- "GRIDSS-2.12.0"

extra.samples.to.exclude <- c()
tumor.depth.cutoff <- 10;

source( file.path(run.dir, "config_files/common_config.R") )

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(magick)

## Extract data to process for making CNA figure
load( file = file.path(macbook.run.dir , 'result/sample_summaries/clinical_data_with_colors.Rdata'))


# Produce Figure 2 from panels
library(magick)

macbook.run.dir <- '/Users/toddjohnson/HGC_mounts/HGC/workspace/runs/IWK_WGS_HMF_20210726'




panelA.file <- file.path(macbook.run.dir, "result/maftools/figures/IWK_GPL_drivers.oncoplot_with_histology_and_chr3ploss.allsamples.pdf")
panelB.file <- file.path(macbook.run.dir, 'result/SVs/Linx/plots/manually_modified/chromosomes/IWK047_T.chr3.015.png')




fig.width <- 19/2.54
fig.height <- 13/2.54
fig.res <- 600

panelA <- image_read_pdf(panelA.file) %>%
		image_annotate("A", size = 14, color = "black", gravity = 'northwest', location="+50+0")

panelB <- image_read(panelB.file) %>%
		image_annotate("B", size = 110, color = "black", gravity = 'northwest', location="+250+100")

blank.img <- image_blank(width=fig.width*fig.res, height=fig.height*fig.res, color = "none", pseudo_image = "", defines = NULL)

#Figure1.file <- file.path(macbook.run.dir, "figures/CellReports_final_figures/JohnsonTA_Figure_1_somatic_variants.pdf")
#quartz( width=fig.width, height=fig.height, type='pdf', file=Figure1.file);

Figure1.file <- file.path(macbook.run.dir, "figures/CellReports_final_figures/JohnsonTA_Figure_1_somatic_variants_composited.pdf")
pdf( width=fig.width, height=fig.height, file=Figure1.file);
#quartz( width=fig.width, height=fig.height)
par(mai=c(0,0,0,0));
merged.panels <- image_composite(image_composite(blank.img, image_scale(panelA, "x3100"), operator='over', offset="+0+0"), image_scale(panelB, "x2240"), operator='over', offset="+2260+820")
plot(merged.panels)
dev.off()


fig.width <- 19/2.54
fig.height <- 11/2.54

panelA <- image_read_pdf(panelA.file) %>%
		image_annotate("A", size = 16, color = "black", gravity = 'northwest', location="+50+0")

panelB <- image_read(panelB.file) %>%
		image_annotate("B", size = 110, color = "black", gravity = 'northwest', location="+250+100")

blank.img <- image_blank(width=fig.width*fig.res, height=fig.height*fig.res, color = "none", pseudo_image = "", defines = NULL)

Figure1.file <- file.path(macbook.run.dir, "figures/CellReports_final_figures/JohnsonTA_Figure_1_somatic_variants_layout.pdf")
pdf( width=fig.width, height=fig.height, file=Figure1.file);
#quartz(width=fig.width, height=fig.height)
#quartz( width=fig.width, height=fig.height)
par(mai=c(0,0,0,0));
layout(matrix(c(1,2,1,2), ncol=2, byrow=TRUE), widths=c(0.50, 0.50), heights=c(0.15, 0.85))
par(mai=c(0,0,0,0));
plot(panelA)
plot(panelB)
dev.off()