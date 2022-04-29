#!/usr/bin/env Rscript

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
	library(parallel)
	library(openxlsx)})

message( paste("Loading processed results for export to Excel for ", study.name, "\n", sep=""))

#run.dirs.ls <- run.dir
#names(run.dirs.ls) <- basename(run.dirs.ls)
study.run.names <- names(run.dirs.ls)

study.stored.object.names <- c('germline.variant.funcClass2.summary', 'somatic.variant.funcClass2.summary', 'purple.qc.dt.combined',
	'purple.purity.dt', 'merged.driver.variant.info.dt', 'sample.info.driver.summary', 'driver.purple.cnv.gene.info',
	'fusions.all.dt', 'fusions.reported.dt', 'germline.variant.info.dt.candidates', 'somatic.variant.info.dt.candidates',
	'summarized.driver.variant.info',
	'linx.SV.CN.info.dt', 'SV_CN.func.summary')
names(study.stored.object.names) <- study.stored.object.names

script.gpl_prefix <- copy(gpl_prefix)
rm( gpl_prefix )

study.data.ls <- lapply(
	X = run.dirs.ls,
	FUN = function( curr.run.dir ){
		source( file.path( curr.run.dir, "config_files/common_config.R" ), local=TRUE )

		load( file = paste(curr.run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_ver.20220210.Rdata", sep=""))
		mget(x=study.stored.object.names)
	})
#load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))

gpl_prefix <- copy(script.gpl_prefix)

merged.study.data.ls <- lapply(
	X = study.stored.object.names,
	FUN = function( curr.object.name ){
		rbindlist(lapply(
			X = study.run.names,
			FUN = function( curr.study.run.name ){
				curr.dt <- study.data.ls[[curr.study.run.name]][[curr.object.name]]
#				curr.dt[,study.run.name:= curr.study.run.name]
				curr.dt
			}))
})

germline.variant.funcClass2.summary <- merged.study.data.ls[['germline.variant.funcClass2.summary']]
somatic.variant.funcClass2.summary <- merged.study.data.ls[['somatic.variant.funcClass2.summary']]
purple.qc.dt.combined <- merged.study.data.ls[['purple.qc.dt.combined']]
purple.purity.dt <- merged.study.data.ls[['purple.purity.dt']]
merged.driver.variant.info.dt <- merged.study.data.ls[['merged.driver.variant.info.dt']]
sample.info.driver.summary <- merged.study.data.ls[['sample.info.driver.summary']]
driver.purple.cnv.gene.info <- merged.study.data.ls[['driver.purple.cnv.gene.info']]
fusions.all.dt <- merged.study.data.ls[['fusions.all.dt']]
fusions.reported.dt <- merged.study.data.ls[['fusions.reported.dt']]
germline.variant.info.dt.candidates <- merged.study.data.ls[['germline.variant.info.dt.candidates']]
somatic.variant.info.dt.candidates <- merged.study.data.ls[['somatic.variant.info.dt.candidates']]
linx.SV.CN.info.dt <- merged.study.data.ls[['linx.SV.CN.info.dt']]
SV_CN.func.summary <- merged.study.data.ls[['SV_CN.func.summary']]
summarized.driver.variant.info  <- merged.study.data.ls[['summarized.driver.variant.info']]

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]

#modified QC categories 3/11/2022
if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

if (remove.deleted.genes.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_DELETED_GENES', QCStatus, invert=TRUE),]
}

if (remove.gender.mismatch.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_GENDER_MISMATCH', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined.QCd[,low.purity:=0L]
purple.qc.dt.combined.QCd[grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),low.purity:=1L]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[which(!(low.purity==1L & tumor.depth<tumor.depth.cutoff)),]

#purple.qc.dt.combined.QCd <- rbind(
#	purple.qc.dt.combined.QCd[QCStatus%in%c('PASS'),],
#	purple.qc.dt.combined.QCd[tumor.depth>=tumor.depth.cutoff,][grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),])

purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[!tumor.sample.id%in%extra.samples.to.exclude,]

QC.pass.samples <- purple.qc.dt.combined.QCd$tumor.sample.id

message( paste("Preparing tables for Excel workbook insertion for ", study.name, "\n", sep=""))

message( paste("Preparing PURPLE variable key for ", study.name, "\n", sep=""))
purple.key <- data.table(
	PURPLE_variable_key = c(
		"Somatic enrichment",
		"PURPLE_CN: Purity adjusted copy number surrounding variant location",
		"PURPLE_MACN: Purity adjusted minor allele copy number surrounding variant location",
		"PURPLE_AF: Purity adjusted allelic frequency of variant",
		"PURPLE_VCN: Purity adjusted copy number of variant",
		"SUBCL: Subclonal likelihood between 0 and 1",
		"BIALLELIC: Flag to indicate variant is biallelic",
		"",
		"Germline enrichment",
		"PURPLE_CN: Purity adjusted copy number surrounding variant location",
		"PURPLE_MACN: Purity adjusted minor allele copy number surrounding variant location",
		"PURPLE_AF: Purity adjusted allelic frequency of variant",
		"PURPLE_VCN: Purity adjusted copy number of variant",
		"BIALLELIC: Flag to indicate variant is biallelic") )

setnames(purple.key, old="PURPLE_variable_key", new="PURPLE variable key" )

germline.variant.funcClass2.summary[,var.type:='gl']
somatic.variant.funcClass2.summary[,var.type:='som']

SV_CN.func.summary[,var.type:='SV_CN']

SV_CN.func.summary.long <- melt(
	data = SV_CN.func.summary[,list(subject.id, tumor.sample.id, var.type, SV=sv.ct, SV_CLUSTER=cluster.ct, DISRUPTED_GENE=SV.disrupted.gene.ct)],
	id.vars = c('subject.id', 'tumor.sample.id', 'var.type'),
	measure.vars = c('SV', 'SV_CLUSTER', 'DISRUPTED_GENE'),
	variable.name = 'funcClass',
	value.name = 'sample.funcClass.ct',
	variable.factor = FALSE,
	value.factor = FALSE)

variant.worst.effect.info <- rbind(
	germline.variant.funcClass2.summary[,list(subject.id, tumor.sample.id, var.type, funcClass=SEW.funcClass2, sample.funcClass.ct=sample.SEW.funcClass2.ct)],
	somatic.variant.funcClass2.summary[,list(subject.id, tumor.sample.id, var.type, funcClass=SEW.funcClass2, sample.funcClass.ct=sample.SEW.funcClass2.ct)],
	SV_CN.func.summary.long[,list(subject.id, tumor.sample.id, var.type, funcClass, sample.funcClass.ct)])

variant.worst.effect.info[,funcClass:=ifelse(is.na(funcClass), "NOT_ANNOTATED", funcClass)]

## 20211216 note: SEW.funcClass2 may be NA, so those will disappear in below table without update
variant.worst.effect.summary <- dcast(
	formula = subject.id + tumor.sample.id ~ funcClass + var.type,
	data = variant.worst.effect.info[,],
	value.var = "sample.funcClass.ct",
	sep = '.')

# germline variants are prefiltered for those in gene panel, so all have annotation
purple.qc.dt.combined <- merge(
	x = purple.qc.dt.combined,
	y = variant.worst.effect.summary[,list(subject.id, tumor.sample.id,
		NONSENSE_OR_FRAMESHIFT.gl, SPLICE.gl, MISSENSE.gl, SYNONYMOUS.gl, NONE.gl,
		NONSENSE_OR_FRAMESHIFT.som, SPLICE.som, MISSENSE.som, SYNONYMOUS.som, NONE.som, NOT_ANNOTATED.som,
		SV.SV_CN, SV_CLUSTER.SV_CN, DISRUPTED_GENE.SV_CN)],
	by = c('subject.id', 'tumor.sample.id'),
	all.x = TRUE)

message( paste("Renaming PURPLE QC table columns for ", study.name, "\n", sep=""))
if ('AmberMeanDepth'%in%colnames(purple.qc.dt.combined)){
	setnames(purple.qc.dt.combined,
		old = c("subject.id", "tumor.sample.id", "AmberGender", "CobaltGender", "Contamination", "CopyNumberSegments", "DeletedGenes",
			"GermlineAberrations", "Method", "Purity", "QCStatus", "UnsupportedCopyNumberSegments", "AmberMeanDepth", "hr_status", "hrd_type",
			"NONSENSE_OR_FRAMESHIFT.gl", "SPLICE.gl", "MISSENSE.gl", "SYNONYMOUS.gl", "NONE.gl",
			"NONSENSE_OR_FRAMESHIFT.som", "SPLICE.som", "MISSENSE.som", "SYNONYMOUS.som", "NONE.som", "NOT_ANNOTATED.som",
			"SV.SV_CN", "SV_CLUSTER.SV_CN", "DISRUPTED_GENE.SV_CN"),
		new = c("Subject ID", "Tumor sample ID", "Amber Gender", "Cobalt Gender", "Contamination", "Copy Number Segments", "Deleted Genes",
			"Germline Aberrations", "Method", "Purity", "QCStatus", "Unsupported Copy Number Segments", "Amber Mean Depth", "HRD status", "HRD type",
			paste( "Germline:\n", c("nonsense or frameshift", "splice", "missense", "synonymous", "no effect"), sep=" "),
			paste( "Somatic:\n", c("nonsense or frameshift", "splice", "missense", "synonymous", "no effect", "not annotated"), sep=" "),
			paste( "SVs:\n", c("SVs", "SV clusters", "Disrupted genes"), sep=" ")))
}else{
	setnames(purple.qc.dt.combined,
		old = c("subject.id", "tumor.sample.id", "AmberGender", "CobaltGender", "Contamination", "CopyNumberSegments", "DeletedGenes",
			"GermlineAberrations", "Method", "Purity", "QCStatus", "UnsupportedCopyNumberSegments", "hr_status", "hrd_type"),
		new = c("Subject ID", "Tumor sample ID", "Amber Gender", "Cobalt Gender", "Contamination", "Copy Number Segments", "Deleted Genes",
			"Germline Aberrations", "Method", "Purity", "QCStatus", "Unsupported Copy Number Segments", "HRD status", "HRD type"))
}

setnames(purple.purity.dt,
	old = c("subject.id", "tumor.sample.id"),
	new = c("Subject ID", "Tumor sample ID"))

message( paste("Preparing driver variant info table for ", study.name, "\n", sep=""))
driver.variant.info <- unique(merged.driver.variant.info.dt[order(tumor.index, -gene.hits, chr),list(
	gene.hits, HRD.gene.hits, gene.GL.hits, gene.SOM.hits, gene.LOH.hits, gene.som.LOH.hits,
	subject.id, tumor.sample.id,
	CHROM, gene,
	snp.id.gl,
	AF.N.gl, AF.T.gl,
	BIALLELIC.gl, P_AF.gl, MACN.gl, VCN.gl,
	func1.gl, func2.gl, bp.descr.gl, aa.descr.gl, ph1="",
	SV_CN.type.som, SV_CN.descr.som, min_minorCN, max_CN,
	ph2="",
	snp.id.som,
	AF.N.som, AF.T.som,
	BIALLELIC.som, P_AF.som, CN.som, VCN.som, MACN.som, GERMLINE.som, SUBCL,
	func1.som, func2.som, bp.descr.som, aa.descr.som)])

driver.variant.info[,AF.N.gl:=ifelse(AF.N.gl==-99, NA, AF.N.gl)]
driver.variant.info[,AF.T.gl:=ifelse(AF.T.gl==-99, NA, AF.T.gl)]
driver.variant.info[,P_AF.gl:=ifelse(P_AF.gl==-99, NA, P_AF.gl)]

driver.variant.info[,AF.N.som:=ifelse(AF.N.som==-99, NA, AF.N.som)]
driver.variant.info[,AF.T.som:=ifelse(AF.T.som==-99, NA, AF.T.som)]
driver.variant.info[,P_AF.som:=ifelse(P_AF.som==-99, NA, P_AF.som)]

driver.variant.info <- driver.variant.info[tumor.sample.id%in%QC.pass.samples,]

message( paste("Preparing candidate driver summary for ", study.name, "\n", sep=""))
candidate.driver.summary <- sample.info.driver.summary[tumor.sample.id%in%QC.pass.samples,list(subject.id, tumor.sample.id,
		CHROM, gene, Description.mrg, Variant.mrg, Variant.SV_CN.som, hr_status, hrd_type, HRD.gene.hits)]
#	CHROM, gene, Description.mrg, Variant.mrg, hr_status, hrd_type, HRD.gene.hits)]

setnames(candidate.driver.summary,
	old = c("subject.id", "tumor.sample.id", "CHROM", "gene", "Description.mrg", "Variant.mrg", "Variant.SV_CN.som", "hr_status", "hrd_type", "HRD.gene.hits"),
	new = c("Subject ID", "Tumor sample ID", "Chrom.", "Gene", "Mutation status description", "Small variant(s) summary", "SV summary", "HRD status", "HRD type", "Hit ct. in BRCA1/BRCA2/PALB2/RAD51"))

## message( paste("Preparing candidate copy-number driver summary for ", study.name, "\n", sep=""))
## CN.drivers <- driver.purple.cnv.gene.info[order(tumor.index),list(subject.id, tumor.sample.id,
##                 CHROM, gene, driver, category, eventType, clusterId, start, end, maxCopyNumber, minMinorAlleleCopyNumber)]
## setnames(CN.drivers,
##         old = c("subject.id", "tumor.sample.id",
##                 "CHROM", "gene", "driver", "category", "eventType","clusterId", "start", "end", "maxCopyNumber", "minMinorAlleleCopyNumber"),
##         new = c("Subject ID", "Tumor sample ID",
##                 "Chrom", "Gene", "Driver type", "Category", "Event type", "Cluster ID", "start", "end", "maxCopyNumber", "minMinorAlleleCopyNumber"))


fusions.all.dt <- fusions.all.dt[tumor.sample.id%in%QC.pass.samples,]
fusions.reported.dt <- fusions.reported.dt[tumor.sample.id%in%QC.pass.samples,]

setnames(fusions.all.dt,
	old = c("subject.id", "tumor.sample.id"),
	new = c("Subject ID", "Tumor sample ID"))

setnames(fusions.reported.dt,
	old = c("subject.id", "tumor.sample.id"),
	new = c("Subject ID", "Tumor sample ID"))

message( paste("Creating Excel workbook of summarized and detailed driver gene information for ", study.name, "\n", sep=""))
require(ggplot2)

wb <- createWorkbook(title = paste(study.name, " study: GRIDSS-PURPLE-LINX-CHORD (", gpl_prefix, "), QC, driver summary, and variant details", sep=""))

modifyBaseFont(wb, fontSize = 12, fontName = "Calibri")

hs <- createStyle(fontColour = "black", fgFill="grey",
	halign = "center", valign = "bottom", textDecoration = "Bold",
	border = "TopBottomLeftRight", borderStyle = "thin", borderColour = "black",
	textRotation = 0, wrapText=TRUE)

forDatStl <- createStyle(fontColour = "#000000", halign="center")

num.2.s <- createStyle(numFmt = "0.00", halign="center")
num.3.s <- createStyle(numFmt = "0.000", halign="center")
left.s <- createStyle(halign="left")

message( paste("Adding worksheets to Excel workbook for ", study.name, "\n", sep=""))

addWorksheet(wb, sheetName = "PURPLE variable key", gridLines = FALSE)
addWorksheet(wb, sheetName = "PURPLE QC summary", gridLines = TRUE)
addWorksheet(wb, sheetName = "PURPLE purity summary", gridLines = TRUE)
addWorksheet(wb, sheetName = "Candidate driver summary", gridLines = TRUE)
addWorksheet(wb, sheetName = "Candidate driver gene variants", gridLines = TRUE)
#addWorksheet(wb, sheetName = "Candidate copy-number drivers", gridLines = TRUE)
addWorksheet(wb, sheetName = "Reported gene fusions", gridLines = TRUE)
addWorksheet(wb, sheetName = "Candidate germline variants", gridLines = TRUE)
addWorksheet(wb, sheetName = "Candidate somatic variants", gridLines = TRUE)
addWorksheet(wb, sheetName = "All gene fusions", gridLines = TRUE)

message( paste("Adding data to worksheets to Excel workbook for ", study.name, "\n", sep=""))

writeDataTable(wb,
	sheet = "PURPLE variable key",
	startCol = 1, startRow = 1,
	x = purple.key,
	colNames = TRUE, rowNames = FALSE,
	withFilter = FALSE,
	headerStyle = hs,
	tableStyle = 'none')
setColWidths(wb, sheet = "PURPLE variable key", cols = "A", widths = 72)
addStyle(wb, sheet = "PURPLE variable key", style = createStyle(fontColour = "#000000", textDecoration = "bold", halign="left"), cols = "A", rows = 1 + which(purple.key[["PURPLE variable key"]]=="Somatic enrichment") )
addStyle(wb, sheet = "PURPLE variable key", style = createStyle(fontColour = "#000000", textDecoration = "bold", halign="left"), cols = "A", rows = 1 + which(purple.key[["PURPLE variable key"]]=="Germline enrichment") )
addStyle(wb, sheet = "PURPLE variable key", style = createStyle(fontColour = "#000000", halign="left"), cols = "A", rows = 1 + which(!purple.key[["PURPLE variable key"]]%in%c("Somatic enrichment", "Germline enrichment")))


freezePane(wb, sheet = "PURPLE QC summary", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "PURPLE QC summary",
	startCol = 1, startRow = 1,
	x = purple.qc.dt.combined,
	colNames = TRUE, rowNames = FALSE,
	withFilter = FALSE,
	headerStyle = hs,
	bandedRows=FALSE,
	tableStyle = 'none')

addStyle(wb, sheet = "PURPLE QC summary", style = num.3.s, cols = c("E", "J"), rows = 2:nrow(purple.qc.dt.combined), gridExpand = T )

setRowHeights(wb, sheet = "PURPLE QC summary", rows = 1, heights = 35)
letters2 <- c(LETTERS, paste('A', LETTERS, sep=""))

addStyle(wb, sheet = "PURPLE QC summary", style = createStyle(fontColour = "#000000", halign="center"),
	cols = letters2[c(1,2,3,4,5,6,7,8,9,10,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)], rows = 1 + 1:nrow(purple.qc.dt.combined), gridExpand = T )

setColWidths(wb, sheet = "PURPLE QC summary",
	cols = c("A", "B", "C", "D", "E",
		"F", "G", "H", "I", "J",
		"K", "L", "M", "N", "O",
		"P", "Q", "R", "S", "T",
		"U", "V", "W", "X", "Y", "Z",
		"AA", "AB", "AC"),
	widths = c(15, 18, 11, 11, 12.33,
		11.5, 6.83, 17.67, 12, 7.17,
		59, 17.83, 13, 15, 13,
		20, 9, 9, 12, 9,
		20, 9, 9, 12, 9, 9,
		12, 12, 12))

freezePane(wb, sheet = "PURPLE purity summary", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "PURPLE purity summary",
		startCol = 1, startRow = 1,
		x = purple.purity.dt,
		colNames = TRUE, rowNames = FALSE,
		withFilter = FALSE,
		headerStyle = hs,
		bandedRows=FALSE,
		tableStyle = 'none')

addStyle(wb, sheet = "PURPLE purity summary", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(purple.purity.dt), rows = 1 + 1:nrow(purple.purity.dt), gridExpand = T )
addStyle(wb, sheet = "PURPLE purity summary", style = num.2.s, cols = c(3,7,11,12,13,14), rows = 2:nrow(purple.purity.dt), gridExpand = T )
addStyle(wb, sheet = "PURPLE purity summary", style = num.3.s, cols = c(4,5,6,10,15,16,18,20,24), rows = 2:nrow(purple.purity.dt), gridExpand = T )

setRowHeights(wb, sheet = "PURPLE purity summary", rows = 1, heights = 20)

setColWidths(wb, sheet = "PURPLE purity summary",
	cols = c(1:ncol(purple.purity.dt)),
	widths = c(15, 18, rep(15, ncol(purple.purity.dt) - 2)))


freezePane(wb, sheet = "Candidate driver summary", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "Candidate driver summary",
	startCol = 1, startRow = 1,
	x = candidate.driver.summary,
	colNames = TRUE, rowNames = FALSE,
	withFilter = FALSE,
	headerStyle = hs,
	bandedRows=FALSE,
	tableStyle = 'none')

setColWidths(wb, sheet = "Candidate driver summary",
		cols = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
	widths = c(15, 18, 15, 15, 72, 60, 60, 15, 15, 25))

setRowHeights(wb, sheet = "Candidate driver summary", rows = 1, heights = 35)


addStyle(wb, sheet = "Candidate driver summary", style = createStyle(halign="center"), cols = "F", rows = 1 + 1:nrow(candidate.driver.summary) )


addStyle(wb, sheet = "Candidate driver summary", style = createStyle(wrapText=TRUE, valign='top'), cols = c("E", "F", "G"), rows = 2:nrow(candidate.driver.summary), gridExpand = T )

variant.header.row.1 <- transpose(data.frame(
	c(rep("", 10), "Germline variants", rep("", 10), "", "SV/CN event", rep("", 3), "", "Somatic variants", rep("", 13))))

variant.header.row.2 <- transpose(data.frame(
	c("Variant/mutation status in gene", rep("", 5),
	rep("", 8),
	"PURPLE Statistics", rep("", 2),
	rep("", 4), "",
	rep("",4), "",
	"", "", "",
	"PURPLE Statistics", rep("", 6),
	rep("", 4))))

variant.header.row.3 <- transpose(data.frame(
	c("Hit ct. (GL+SOM+LOH)", "HRD gene hits", "Germline", "Somatic", "Germline variant biallelic in tumor (LOH)", "Somatic variant biallelic in tumor (LOH)",
		"Subject ID", "Tumor sample ID",
		"Chr", "Gene", "SNP ID (Chr#-POS-REF-ALT) X=23, Y=24, MT=25", "AF.N", "AF.T", "BIALLELIC", "PAF", "MACN", "VCN",
		"Function class 1", "Function class 2", "DNA Tx. description", "Amino acid Tx. description",
		"",
		"SV/CN event-type", "SV/CN descr.", "Min. minor CN", "Max. CN",
		"",
		"SNP ID (Chr#-POS-REF-ALT) X=23, Y=24, MT=25", "AF.N", "AF.T",
		"BIALLELIC ", "PAF", "CN", "VCN", "MACN", "GERMLINE", "SUBCL",
		"Function class 1", "Function class 2", "DNA Tx. description", "Amino acid Tx. Description")))

writeData(wb, sheet = "Candidate driver gene variants",
	startCol = 1, startRow = 1,
	x = variant.header.row.1,
	colNames = FALSE, rowNames = FALSE,
	withFilter = FALSE)

writeData(wb, sheet = "Candidate driver gene variants",
	startCol = 1, startRow = 2,
	x = variant.header.row.2,
	colNames = FALSE, rowNames = FALSE,
	withFilter = FALSE)

writeData(wb, sheet = "Candidate driver gene variants",
	startCol = 1, startRow = 3,
	x = variant.header.row.3,
	colNames = FALSE, rowNames = FALSE,
	withFilter = FALSE)

writeData(wb, sheet = "Candidate driver gene variants",
	startCol = 1, startRow = 4,
	x = driver.variant.info,
	colNames = FALSE, rowNames = FALSE,
	withFilter = FALSE)


addStyle(wb, sheet = "Candidate driver gene variants", style = hs, cols = 1:ncol(driver.variant.info), rows = 1:3, gridExpand = T )
addStyle(wb, sheet = "Candidate driver gene variants", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(driver.variant.info), rows = 3 + 1:nrow(driver.variant.info), gridExpand = T )

addStyle(wb, sheet = "Candidate driver gene variants", style = num.3.s, cols = c("M", "N", "O", "P", "Q", "Y", "Z", "AD", "AE", "AF", "AG", "AH", "AI"), rows = 3 + 1:nrow(driver.variant.info), gridExpand = T )
addStyle(wb, sheet = "Candidate driver gene variants", style = left.s, cols = c("K", "AB"), rows = 3 + 1:nrow(driver.variant.info), gridExpand = T )

setRowHeights(wb, sheet = "Candidate driver gene variants", rows = 1, heights = 20)
setRowHeights(wb, sheet = "Candidate driver gene variants", rows = 2, heights = 20)
setRowHeights(wb, sheet = "Candidate driver gene variants", rows = 3, heights = 53)

setColWidths(wb, sheet = "Candidate driver gene variants",
	cols = c("A", "B", "C", "D", "E", "F",
		"G", "H", "I", "J",
		"K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
		"V",
		"W", "X", "Y", "Z",
		"AA",
		"AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK",
		"AL", "AM", "AN", "AO"),
	widths = c(15, 18, 11, 11, 13.33, 13.33,
		12, 17.67, 8, 12,
		27, 9, 7, 9, 7, 7, 7, 32, 32, 32, 32,
		2,
		20, 32, 8, 8,
		2,
		59, 9, 7, 7, 7, 7, 7, 7, 9, 8,
		32, 32, 32, 32))

mergeCells(wb, sheet = "Candidate driver gene variants", cols = 11:21, rows = 1)
mergeCells(wb, sheet = "Candidate driver gene variants", cols = 23:26, rows = 1)
mergeCells(wb, sheet = "Candidate driver gene variants", cols = 28:41, rows = 1)
mergeCells(wb, sheet = "Candidate driver gene variants", cols = 1:6, rows = 2)
mergeCells(wb, sheet = "Candidate driver gene variants", cols = 15:17, rows = 2)
mergeCells(wb, sheet = "Candidate driver gene variants", cols = 31:36, rows = 2)

freezePane(wb, sheet = "Candidate driver gene variants", firstActiveRow = 4, firstActiveCol = 10)

## freezePane(wb, sheet = "Candidate copy-number drivers", firstRow = TRUE, firstCol = FALSE)
## writeDataTable(wb, sheet = "Candidate copy-number drivers",
##     x = CN.drivers,
##     colNames = TRUE, rowNames = FALSE,
##     withFilter = FALSE,
##     headerStyle = hs,
##     tableStyle = 'none')
## 
## addStyle(wb, sheet = "Candidate copy-number drivers", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(CN.drivers), rows = 1 + 1:nrow(CN.drivers), gridExpand = T )
## 
## setRowHeights(wb, sheet = "Candidate copy-number drivers", rows = 1, heights = 20)
## 
## setColWidths(wb, sheet = "Candidate copy-number drivers",
##     cols = c("A", "B", "C", "D", "E",
##     "F", "G", "H", "I", "J", "K", "L"),
##     widths = c(15, 18, 11, 11, 12.33,
##     11.5, 11, 11, 11.5, 11.5, 30, 40))


freezePane(wb, sheet = "Reported gene fusions", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "Reported gene fusions",
	x = fusions.reported.dt,
	colNames = TRUE, rowNames = FALSE,
	withFilter = FALSE,
	headerStyle = hs,
	tableStyle = 'none')

addStyle(wb, sheet = "Reported gene fusions", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(fusions.reported.dt), rows = 1 + 1:nrow(fusions.reported.dt), gridExpand = T )

setRowHeights(wb, sheet = "Reported gene fusions", rows = 1, heights = 20)

setColWidths(wb, sheet = "Reported gene fusions",
	cols = c(1:18),
	widths = c(15, 18, rep(15, ncol(fusions.reported.dt) - 2)))

germline.variants.dt <- germline.variant.info.dt.candidates[tumor.sample.id%in%QC.pass.samples,
	list(subject.id, tumor.sample.id,
	snp.id, CHROM, start, end, REF, ALT,
	AF.N, AF.T, BIALLELIC, PURPLE_AF, PURPLE_CN, PURPLE_VCN,
	PATH, gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]

freezePane(wb, sheet = "Candidate germline variants", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "Candidate germline variants",
		x = germline.variants.dt,
		colNames = TRUE, rowNames = FALSE,
		withFilter = FALSE,
		headerStyle = hs,
		tableStyle = 'none')

addStyle(wb, sheet = "Candidate germline variants", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(germline.variants.dt), rows = 1 + 1:nrow(germline.variants.dt), gridExpand = T )

setRowHeights(wb, sheet = "Candidate germline variants", rows = 1, heights = 20)

#setColWidths(wb, sheet = "Candidate germline variants",
#		cols = c(1:ncol(germline.variants.dt)),
#		widths = c(15, 18, rep(15, ncol(fgermline.variants.dt) - 2)))


somatic.variants.dt <- somatic.variant.info.dt.candidates[tumor.sample.id%in%QC.pass.samples,
	list(subject.id, tumor.sample.id, snp.id, CHROM, POS=start, REF, ALT,
	AF.N, AF.T, BIALLELIC, PURPLE_VCN, PURPLE_AF, PURPLE_CN, PURPLE_MACN, PURPLE_GERMLINE, SUBCL,
	gene, ENST, funcClass1, funcClass2, bp.descr, aa.descr)]

freezePane(wb, sheet = "Candidate somatic variants", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "Candidate somatic variants",
		x = somatic.variants.dt,
		colNames = TRUE, rowNames = FALSE,
		withFilter = FALSE,
		headerStyle = hs,
		tableStyle = 'none')

addStyle(wb, sheet = "Candidate somatic variants", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(somatic.variants.dt), rows = 1 + 1:nrow(somatic.variants.dt), gridExpand = T )

setRowHeights(wb, sheet = "Candidate somatic variants", rows = 1, heights = 20)

#setColWidths(wb, sheet = "Candidate somatic variants",
#		cols = c(1:ncol(somatic.variants.dt)),
#		widths = c(15, 18, rep(15, ncol(somatic.variants.dt) - 2)))


freezePane(wb, sheet = "All gene fusions", firstRow = TRUE, firstCol = FALSE)
writeDataTable(wb, sheet = "All gene fusions",
		x = fusions.all.dt,
		colNames = TRUE, rowNames = FALSE,
		withFilter = FALSE,
		headerStyle = hs,
		tableStyle = 'none')

addStyle(wb, sheet = "All gene fusions", style = createStyle(fontColour = "#000000", halign="center"), cols = 1:ncol(fusions.all.dt), rows = 1 + 1:nrow(fusions.all.dt), gridExpand = T )

setRowHeights(wb, sheet = "All gene fusions", rows = 1, heights = 20)

setColWidths(wb, sheet = "All gene fusions",
		cols = c(1:ncol(fusions.reported.dt)),
		widths = c(15, 18, rep(15, ncol(fusions.reported.dt) - 2)))


saveWorkbook(wb,
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/", study.name, "_variant_information_ver.20220210_", Sys.Date(), ".xlsx", sep=""),
	overwrite = TRUE)

message( paste("Finished creating Excel workbook for ", study.name, "\n", sep=""))
