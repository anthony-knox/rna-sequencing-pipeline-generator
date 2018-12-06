#Instructions for DESeq2

You will need to move the abundance.h5 file for each sample to their respective directory within the DESeq2 directory. Enter these commands, one by one in a terminal or command prompt window.

	cd @@@DESEQ2_SHELLSCRIPTS_ABSOLUTEPATH@@@

	vi cp_abundance_to_DESeq2.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

	#!/bin/bash
	while IFS= read -r line || [[ -n "$line" ]]; do
		cp @@@KALLISTO_ABSOLUTEPATH@@@/$line/output/abundance.h5 @@@DESEQ2_KALLISTO_SUBDIR_ABSOLUTEPATH@@@/$line
	done < "$1"

Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	chmod u+x cp_abundance_to_DESeq2.sh

	./cp_abundance_to_DESeq2.sh @@@KALLISTO_SAMPLES_ABSOLUTEPATH@@@

You will now need to download R and RStudio if you have not already. To do so, visit [the RStudio download page](https://www.rstudio.com/products/rstudio/download/)

Next, you will have to download the list of transcript stable id versions from [BiomaRt](https://www.ensembl.org/biomart). Choose your dataset (for example Human genes (GRCh38.p12)). Make sure this matches the same dataset that you chose to download for the kallisto index. If you need a previous version, Ensembl provides a BiomaRt page for any previous versions, just google the version name + ensembl + BiomaRt. Once you choose your dataset, navigate to Attributes and make sure that only "Transcript stable ID version" is selected. Click on Results in the top left and choose to export all results to a File CSV. Now rename this file to `ensembl_transcript_id_version.csv` (note the csv extension) and place it in `@@@DESEQ2_ABSOLUTEPATH@@@`

Now you are ready to perform differential expression analysis using DESeq2. Open RStudio, click File -> New File -> R Script and copy and paste the code below. You can execute commands by navigating to a specific line, the pressing Ctrl+Enter (or Command+Enter). Follow instructions and voila! You now have the tools to perform RNA-Sequencing analysis. I have included links to the vignettes that relate to this specific pipeline. I highly suggest reading through these vignettes (especially the one for DESeq2) to better understand the performance of each function and the parameters that are associated with each function. You can play around with these parameters and shape them to your own needs. I also suggest exploring different data visualization tools, such as these below after you have obtained your differential expression dataset:

[GenePattern Notebook](http://genepattern-notebook.org/)  
[GENE-E and Morpheus](https://software.broadinstitute.org/GENE-E/)  
[Multiplot](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/multiplot/2)  

Thank you for following along and I hope this tool was useful!

	#Installations only need to be run once on your computer
	if (!requireNamespace("BiocManager", quietly = TRUE))
	  install.packages("BiocManager")
	BiocManager::install("DESeq2", version = "3.8")
	BiocManager::install("biomaRt", version = "3.8")
	BiocManager::install("tximport", version = "3.8")
	BiocManager::install("rhdf5", version = "3.8")
	install.packages("ggplot2")
	BiocManager::install("vsn", version = "3.8")
	install.packages("pheatmap")
	BiocManager::install("pcaExplorer", version = "3.8")
	BiocManager::install("dexus", version = "3.8")
	if (!require("RColorBrewer")) {
	  install.packages("RColorBrewer")
	}
	install.packages("SIBERG")
	install.packages("doParallel")
	install.packages("doSNOW")
	install.packages("hexbin")
	
	#Load libraries - this needs to be performed any time you create a new workspace in RStudio
	library("DESeq2")
	library("biomaRt")
	library("tximport")
	library("rhdf5")
	library("ggplot2")
	library("vsn")
	library("pheatmap")
	library("pcaExplorer")
	library("dexus")
	library("RColorBrewer")
	library("SIBERG")
	library("doParallel")
	library("doSNOW")
	library("hexbin")
	
	#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
	#NOTE: I would highly suggest reading through this vignette, as the majority of this pipeline is based on this DESeq2 vignette
	#For Help on functions - ?functionName, or for a full list:
	help(package="DESeq2")
	
	#To remove scientific notation in csv files
	options(scipen=999)
	
	#set DESeq2 directory
	dir <- "@@@DESEQ2_ABSOLUTEPATH@@@"
	
	#to save files to the results directory, which will be your working directory
	setwd("@@@DESEQ2_RESULTS_ABSOLUTEPATH@@@")
	
	#read in samples
	samples <- read.table(file.path(dir,"samples_deseq2.txt"), header=TRUE)
	
	#specify path to files
	files <- file.path(dir,"kallisto", samples$run, "abundance.h5")
	names(files) <- samples$run
	
	#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
	#convert transcripts to genes with biomaRt - this is the preferred method
	#You need to download the list of transcript id's (use transcript stable id version) from BiomaRt - https://www.ensembl.org/biomart
	#Download and save it as ensembl_transcript_id_version.csv and place it in the DESeq2 directory
	#You can choose to convert to gene name (Eg. GAPDH) or geneID (Eg. ENSG##############.#), or any other attributes found in the listFilters function
	#To do so, replace 'external gene name' with any other attribute
	listMarts()
	mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
	listFilters(mart)
	listDatasets(mart)
	values <- read.csv(file.path(dir, "ensembl_transcript_id_version.csv"), header=TRUE, sep = ",", stringsAsFactors=FALSE)
	values <- values[,1]
	head(values)
	tx2geneName <- getBM(attributes=c('ensembl_transcript_id_version', 'external_gene_name'), filters='ensembl_transcript_id_version', values=values, mart=mart)
	head(tx2geneName)
	write.csv(as.data.frame(tx2geneName), file="tx2geneName.csv")
	tx2geneID <- getBM(attributes=c('ensembl_transcript_id_version', 'ensembl_gene_id'), filters='ensembl_transcript_id_version', values=values, mart=mart)
	head(tx2geneID)
	write.csv(as.data.frame(tx2geneID), file="tx2geneID.csv")
	
	#OR read in transcripts to genes file if you've defined your own  - tx2gene.csv should reside in your DESeq2 directory
	#The first column should have header TXNAME and contain a list of transcript id's with version id (Eg. ENST###############.#)
	#The second column should have header GENEID and contain a list of gene names (Eg. GAPDH)
	tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
	
	#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
	#import quantification data (transcript level estimates) - FOR TRANSCRIPT TO GENE NAME (EX A1BG)
	txi <- tximport(files, type="kallisto", tx2gene=tx2geneName)
	head(txi)
	write.csv(as.data.frame(txi), file="tx_level_estimates_geneName.csv")
	#import quantification data (transcript level estimates) - FOR TRANSCRIPT TO GENE ID (EX ENSG00000000003)
	txi <- tximport(files, type="kallisto", tx2gene=tx2geneID)
	head(txi)
	write.csv(as.data.frame(txi), file="tx_level_estimates_geneID.csv")
	
	#construct DESeqDataSet object from txi object
	dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
	
	#To get normalized counts matrix (for DEXUS, SIBER, and prcomp analysis later)
	dds <- estimateSizeFactors(dds)
	normalizedCountsMatrix <- counts(dds, normalized=TRUE)
	head(normalizedCountsMatrix)
	write.csv(as.data.frame(normalizedCountsMatrix), file="normalizedCountsMatrix.csv")
	
	#Differential expression analysis steps
	#In case you got normalized counts matrix previously - re-construct
	dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
	#Pre-filter for low-count genes (optional but recommended)
	keep <- rowSums(counts(dds)) >= 3
	dds <- dds[keep,]
	#this will use the Wald test - if you want to use LRT, add in parameters: test="LRT", reduced=~1
	dds <- DESeq(dds)
	#obtain results
	#this will have independent filtering with alpha = 0.1, which gives adjusted p value cutoff = 0.1
	#to change it to 0.05, add in parameter alpha=0.05
	#can specify contrast - which conditions to compare - it will give results of log(A/B) - contrast=c("condition","B","A")
	res <- results(dds)
	res
	resultsNames(dds)
	
	#Only run this if you want unfiltered results
	dds <- DESeq(dds, minReplicatesForReplace=Inf)
	res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
	
	#see more information about which variables and tests were used
	mcols(res)$description
	
	#perform lfcshrink (shrinkage of log2 fold change)
	resLFC <- lfcShrink(dds, coef=2, res=res)
	resLFC
	
	#order results by lowest p value
	resOrdered <- res[order(res$padj),]
	resOrdered
	#or for lfcshrink data
	resOrdered <- resLFC[order(res$padj),]
	resOrdered
	
	#Output the results to a csv file with certain parameters such as ordered data with padj<0.1
	resSigUpregulated <- subset(resOrdered, padj < 0.1, padj > 0.0)
	resSigUpregulated <- subset(resSigUpregulated, log2FoldChange > 0)
	resSigDownregulated <- subset(resOrdered, padj < 0.1, padj > 0.0)
	resSigDownregulated <- subset(resSigDownregulated, log2FoldChange < 0)
	write.csv(as.data.frame(resSigUpregulated), file="results_Upregulated.csv")
	write.csv(as.data.frame(resSigDownregulated), file="/results_Downregulated.csv")
	
	#summarize results (change res to whatever object you want summarized)
	summary(res)
	
	#How many adjusted p values were less than 0.1?
	sum(res$padj < 0.1, na.rm=TRUE)
	#Example: plotting mean normalized counts vs log10(pvalue) - shows outliers
	plot(res$baseMean+1, -log10(res$pvalue),
	     log="x", xlab="mean of normalized counts",
	     ylab=expression(-log[10](pvalue)),
	     ylim=c(0,30),
	     cex=.4, col=rgb(0,0,0,.3))
	
	#the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the 	DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
	plotMA(res, ylim=c(-2,2))
	#it is more useful to visualize the lfcshrink data, change parameters (-2,2) to fit data
	plotMA(resLFC, ylim=c(-2,2))
	
	#view plot counts in each condition for a specific gene - here the lowest padj value
	#use gene="geneName" or a function such as lowest padj value
	plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
	#for a customizable plot, use ggplot2
	d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
	ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))
	#Transform the data for visualization and clustering - removes the dependency of variance on the mean (vst = variance stabilizing transformation)
	#blind=false is used when we expect a large portion of genes to have differences in counts due to experimental design
	#a second way to do this is to use the rlog function, but it takes very long for large # of samples, but similar output
	vsd <- vst(dds, blind=FALSE)
	head(assay(vsd), 3)
	
	#To show the effect of transforming the data, view plot of mean vs standard deviation
	ntd <- normTransform(dds)
	meanSdPlot(assay(ntd))
	meanSdPlot(assay(vsd))
	#Now view heatmaps of non-transformed data (ntd)
	select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:500]
	log2.norm.counts <- assay(ntd)[select,]
	pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=TRUE, 
         fontsize_row = 6, fontsize_col = 6, main = "Heatmap_ntd_topCounts", width = 30, height = 60, 
         filename = "Heatmap_ntd_topCounts.pdf", treeheight_row = 100, treeheight_col = 100)
	#and transformed data (vsd)
	log2.norm.counts <- assay(vsd)[select,]
	pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=TRUE, 
         fontsize_row = 6, fontsize_col = 6, main = "Heatmap_vst_topCounts", width = 30, height = 60, 
         filename = "Heatmap_vst_topCounts.pdf", treeheight_row = 100, treeheight_col = 100)
	
	#to get a heatmap of transformed data - first X# genes with lowest padj
	resOrdered <- res[order(res$padj),]
	#Check how many have padj < 0.05
	resSig <- subset(resOrdered, padj < 0.05)
	write.csv(as.data.frame(resSig), file="results_padj_lessThan_5pct.csv")
	#open csv file to check # padj < 0.05, continue with heatmap and change value to # of values padj < 0.05 (1:# at the end)
	select <- rownames(resOrdered[resOrdered$padj<0.05 & !is.na(resOrdered$padj), ])[1:24]
	log2.norm.counts <- assay(vsd)[select,]
	pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, show_colnames=TRUE, 
         fontsize_row = 6, fontsize_col = 20, main = "Heatmap_vst_padj_lessThan_5pct", width = 30, height = 60, 
         filename = "Heatmap_vst_padj_lessThan_5pct.pdf", treeheight_row = 100, treeheight_col = 100)
	
	#to get a heatmap of transformed data - first X# genes with highest absolute value of log2 fold change and padj<0.05
	resOrdered <- res[order(abs(res$log2FoldChange), decreasing=TRUE),]
	#Check how many have padj < 0.05
	resSig <- subset(resOrdered, padj < 0.05)
	write.csv(as.data.frame(resSig), file="results_padj_lessThan_5pct_largestLog2FC.csv")
	#open csv file to check # padj < 0.05, continue with heatmap and change value to # of values padj < 0.05 (1:# at the end)
	select <- rownames(resSig[resSig$padj<0.05 & !is.na(resSig$padj), ])[1:24]
	log2.norm.counts <- assay(vsd)[select,]
	pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=TRUE, 
         fontsize_row = 6, fontsize_col = 20, main = "Heatmap_vst_padj_lessThan_5pct_largestLog2FC.pdf", width = 30, height = 60, 
         filename = "Heatmap_vst_padj_lessThan_5pct_largestLog2FC.pdf", treeheight_row = 100, treeheight_col = 100)
	
	
	#PCA Explorer
	pcaExplorer(dds = dds, dst = vsd)
	
	#Or use prcomp - using normalizedCountsMatrix object
	#https://www.biostars.org/p/282685/ and https://www.biostars.org/p/280615/
	project.pca <- prcomp(t(normalizedCountsMatrix))
	summary(project.pca)
	project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
	
	#bar plot of proportion of variances
	pdf("prcomp_proportionOfVariances.pdf", width=8.5, height=11)
	barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of 	variation (%)", main="Scree plot", ylim=c(0,100))
	dev.off()
	
	#bi-plot for PC 1-10
	par(cex=1.0, cex.axis=0.8, cex.main=0.8)
	pdf("prcomp_bi-plot_1-10.pdf", width=9, height=9)
	pairs(project.pca$x[,1:10], col="black", main="Principal components analysis bi-plot\nPCs 1-10", pch=16)
	dev.off()
	#Output coordinates for bi-plot
	write.csv(project.pca$x, file = "prcomp_Bi-plot_coordinates.csv")
	
	#bi-plot for PC 1 and 2
	par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
	pdf("prcomp_bi-plot_1vs2.pdf", width=5, height=5)
	plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), 	"%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
	points(project.pca$x, col="black", pch=16, cex=1)
	dev.off()
	
	#View PCA components
	head(project.pca$rotation)
	#Output PCA components to data file
	write.csv(project.pca$rotation, file = "prcomp_PCA_components.csv")
	#DEXUS - bimodality analysis
	#https://bioconductor.org/packages/release/bioc/vignettes/dexus/inst/doc/dexus.pdf
	
	
	#DEXUS results object - normalization performed previously for dds object
	dexus_results <- dexus(normalizedCountsMatrix, normalization = "none")
	dexus_results
	
	#Genes shown are the top 15 based on I/NI values - top is highest I/NI
	pdf("DEXUS_heatmap_15.pdf", width = 100, height = 8)
	plot(sort(dexus_results))
	dev.off()
	#plot the top 100 genes
	pdf("DEXUS_heatmap_100.pdf", width = 100, height = 50)
	plot(sort(dexus_results), idx=1:100)
	dev.off()
	
	#export the entire results
	dexus_results_data_frame <- as.data.frame(sort(dexus_results))
	write.csv(dexus_results_data_frame, file="DEXUS_INI_results.csv")
	
	
	#SIBER - bimodality analysis
	#https://cran.r-project.org/web/packages/SIBERG/vignettes/SIBER.pdf
	#remove genes with 0 count
	keep <- rowSums(normalizedCountsMatrix) > 0
	normalizedCountsMatrixNonZero <- normalizedCountsMatrix[keep,]
	#Parallel processing
	cl <- makeCluster(3, type = "SOCK")
	registerDoParallel(cl)
	#Perform function on each gene
	func <- function(i) {
	  SIBER(y=normalizedCountsMatrixNonZero[i, ], model='LN')
	}
	BIinfo_LN <- foreach(i=1:nrow(normalizedCountsMatrixNonZero),
	                     .combine='rbind',
	                     .packages='SIBERG') %dopar% {
	                       func(i)
	                     }
	#results from SIBER are processed in the order of the normalizedCountsMatrixNonZero
	#open the SIBER_results_order.csv file and your first column will be the names of your SIBER results
	write.csv(as.data.frame(normalizedCountsMatrixNonZero), file="SIBER_results_order.csv")
	write.csv(file="SIBER_results.csv", BIinfo_LN)

