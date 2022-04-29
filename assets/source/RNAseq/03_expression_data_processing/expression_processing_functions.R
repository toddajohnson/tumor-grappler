tximport_processing <- function(input.sample.info, input.read.length = 126){
	##
	# curr.sample.info: each sample ID should have SAMPLEID.isf.transcript_data.csv and SAMPLEID.isf.gene_data.csv
	# under SAMPLEID directoryin isofox results directory 
	##
	require(tximport)
	
	isofox.dir <- file.path(run.dir, 'result', 'isofox')
	output.dir <- file.path(run.dir, 'result/processed_expression_data')

	print( paste("Checking if output directory exists", sep=""))
	if ( file.exists(output.dir, recursive=FALSE)==FALSE ){
		print( paste("Creating output directory: ", output.dir, sep=""))
	    dir.create(output.dir)
	}
	
	sample.ids.ls <- input.sample.info$sample.id
	
	print( paste("Adding file names to sample info table ", sep=""))
	input.sample.info[,tx.expr.file:=paste(isofox.dir, '/', sample.id, '/', sample.id, '.isf.transcript_data.csv', sep="")]
	input.sample.info[,gene.expr.file:=paste(isofox.dir, '/', sample.id, '/', sample.id, '.isf.gene_data.csv', sep="")]
	
	print( paste("Extracting transcript expression file names ", sep=""))
	tx.expr.files <- input.sample.info$tx.expr.file
	names(tx.expr.files) <-  sample.ids.ls
	
	print( paste("Extracting gene expression file names ", sep=""))
	gene.expr.files <- input.sample.info$gene.expr.file
	names(gene.expr.files) <- sample.ids.ls
	
	print( paste("Reading first transcript file to access transcript and gene names ", sep=""))
	one.tx.file.dt <- fread(tx.expr.files[[1]])
	tx2gene.dt <- one.tx.file.dt[,list(TransName, GeneName)]

	one.gene.file.dt <- fread(gene.expr.files[[1]])
	geneid2gene.dt <- one.gene.file.dt[,list(GeneId, GeneName)]

	imported.tx.expr.data <- tximport(
			files = tx.expr.files,
			type = 'none',
			txIn = TRUE,
			txOut = FALSE,
			countsFromAbundance = "no",
			tx2gene = as.data.frame(tx2gene.dt),
			geneIdCol = 'GeneName',
			txIdCol = 'TransName',
			abundanceCol = 'AdjTPM',
			countsCol = 'FittedFragments',
			lengthCol = 'EffectiveLength',
			importer = fread,
			readLength = input.read.length)
	
	imported.gene.expr.data <- tximport(
			files = gene.expr.files,
			type = 'none',
			txIn = FALSE,
			txOut = FALSE,
			countsFromAbundance = "no",
			tx2gene = as.data.frame(geneid2gene.dt),
			geneIdCol = 'GeneName',
			abundanceCol = 'AdjTPM',
			countsCol = 'SplicedFragments',
			lengthCol = 'GeneLength',
			importer = fread,
			readLength = input.read.length)

	list(
		sample.ids.ls = sample.ids.ls,
		sample.info = input.sample.info,
		imported.tx.expr.data = imported.tx.expr.data,
		imported.gene.expr.data = imported.gene.expr.data)
}

normalize_and_filter_counts <- function( input.processed.tximport.data.ls, input.grouping.column = 'CancerType', input.samples.to.remove = c() ){
	require(edgeR)
	options(width=350)

	txi <- copy(input.processed.tximport.data.ls[['imported.gene.expr.data']])
	
	adjTPM <- txi$abundance
	adjTPM.column.sums <- colSums(adjTPM)
	
	print( adjTPM.column.sums )

	cts <- txi$counts
	normMat <- txi$length
	
	gene.expr.ct <- nrow(cts)
	sample.ct <- ncol(cts)
	
	print( paste('There were ', gene.expr.ct, ' genes and ', sample.ct, ' samples for ', study.name, sep=""))

	cols.to.keep <- colnames(cts)
	cols.to.keep <- cols.to.keep[which(!cols.to.keep%in%input.samples.to.remove)]
	
	sample.ids.filtered.ls <- cols.to.keep
	
	cts <- cts[,cols.to.keep]
	normMat <- normMat[,cols.to.keep]

	# Obtaining per-observation scaling factors for length, adjusted to avoid
	# changing the magnitude of the counts.
	normMat <- normMat/exp(rowMeans(log(normMat)))
	normCts <- cts/normMat
	
	# Computing effective library sizes from scaled counts, to account for
	# composition biases between samples.
	eff.lib <- calcNormFactors(normCts) * colSums(normCts)
	
	# Combining effective library sizes with the length factors, and calculating
	# offsets for a log-link GLM.
	normMat <- sweep(normMat, 2, eff.lib, "*")
	normMat <- log(normMat)

	# Creating a DGEList object for use in edgeR.
	
	groups <- input.processed.tximport.data.ls[['sample.info']][[input.grouping.column]]
	names(groups) <- input.processed.tximport.data.ls[['sample.info']]$sample.id
	groups <- groups[colnames(cts)]
	groups.f <- factor(x = groups, ordered = FALSE)
	
	y <- DGEList(cts, group=groups.f)
	
	y <- scaleOffset(y, normMat)
	
	# filtering on expression, reduce library sizes to 
	keep <- filterByExpr(y)
	y <- y[keep,,keep.lib.sizes=FALSE]
	
# TMM normalization and log transformation
	mat.all.samples.natural.cpm <- cpm(y, log=FALSE)
	
	mat.all.samples <- cpm(y, log=TRUE)

	list(
			sample.ids.filtered = sample.ids.filtered.ls,
			expr_filtered_DGEList = y,
			mat.all.samples.natural.cpm = mat.all.samples.natural.cpm,
			mat.all.samples = mat.all.samples)
}