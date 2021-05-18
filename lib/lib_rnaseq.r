library(tximport)
library(readr)
library(tximportData)
library("S4Vectors")
library(GenomicFeatures)
library(stringi)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(openxlsx)
library(EnrichmentBrowser)
library(npGSEA)
library(GSEABase)
library(annotate)
library(XML)
library(dbplyr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(uwot)
library(ReportingTools)
library(EnsDb.Hsapiens.v86) #BiocManager::install("EnsDb.Hsapiens.v86")
library(ggrepel)

source(file.path(Sys.getenv('R_UTIL'),'lib_workflow.r'))

combine_read_cnts_to_mat <- function(rcnts,comp_info_dt) {
	# browser()
	rcnts2 <- list()
	
	sgroups <- split(comp_info_dt,by="group")
	for (gname in names(sgroups)) {
		#sgroup_name <- "ctrl"
		for (sample in sgroups[[gname]]$sample) {
			if (!is.na(sample)) {
				#sample <- "D1"
				message(sprintf("gname[%s],sample[%s]",gname,sample))
				rcnts2[[sample]] <- rcnts[[sample]][,!"gene_names"]
				colnames(rcnts2[[sample]]) <- sprintf("%s_%s_%d",gname,sample,1:dim(rcnts2[[sample]])[2])
			}
		}
	}
	
	rcnt2 <- rcnts2[[1]]
	if (length(rcnts2)>1) {
		for (i in 2:length(rcnts2)) {
			rcnt2 <- cbind(rcnt2,rcnts2[[i]])
		}
		rm(rcnts2)
	}
	
	rcnt2[is.na(rcnt2)] <- 0
	cnt_mat <- as.matrix(rcnt2)
	rownames(cnt_mat) <- rcnts[[comp_info_dt[1,sample]]]$gene_names
	colname_orginal <- colnames(cnt_mat)
	
	parsed <- tstrsplit(colname_orginal,"_")
	
	meta_row <- rownames(cnt_mat)
	
	meta_col <- data.table(orig=colname_orginal,
												 comp=parsed[[1]],
												 sample=parsed[[2]],
												 sindex=parsed[[3]])
	
	return(list(mtx=cnt_mat,
							meta_row=meta_row,
							meta_col=meta_col))
}

run_deseq2 <- function(deseq2_input,out_dir,design_f_str="~ comp",shrink_lfc=FALSE,alpha=0.1,padj_co=1,exptag="testing",debug2=0) {
	
	
	dds <- DESeqDataSetFromMatrix(countData = deseq2_input$cts,
																colData = deseq2_input$coldata,
																design= as.formula(design_f_str))
	dds <- DESeq(dds)
	
	if (shrink_lfc){ #apply fc shrinkage with apeglm
		res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
	} else {
		res <- results(dds,alpha=alpha)
	}
	
	# plotMA(res)
	
	res <- res[res$baseMean>0.,]
	exp_description <- mcols(res)$description
	
	if (debug2==1){browser()}
	res_dt <- as.data.table(res[order(res$padj),],keep.rownames=TRUE)
	
	wb <- createWorkbook("deseq2_report")
	sheet_name <- 'readme'
	addWorksheet(wb,sheetName = sheet_name)
	exp_description <- c(exp_description,sprintf('BH alpha=%g;design=%s;shrink_lfc=%s',alpha,design_f_str,shrink_lfc))
	writeData(wb,sheet_name,
						data.table(desc=exp_description))
	
	sheet_name <- 'comp_ctrl_vs_expr'
	addWorksheet(wb,sheetName = sheet_name)
	
	writeData(wb,sheet_name,
						res_dt)
	
	# --------------
	read_cnt <- as.data.table(counts(dds),keep.rownames = T)
	sheet_name <- 'raw_read_count'
	addWorksheet(wb,sheetName = sheet_name)
	writeData(wb,sheet_name,
						read_cnt)
	
	# --------------
	read_cnt <- as.data.table(counts(dds,normalized=TRUE),keep.rownames = T)
	sheet_name <- 'normalized_read_count'
	addWorksheet(wb,sheetName = sheet_name)
	writeData(wb,sheet_name,
						read_cnt)
	
	xlsx_fn <- file.path(out_dir,sprintf('deseq2_result_%s.xlsx',exptag))
	message(xlsx_fn)
	saveWorkbook(wb,file=xlsx_fn,overwrite = TRUE)
	
	return(res_dt)
}


dge_analysis_v2 <- function(cnt_mat,gtf_gene_dt,outd,plot_pca=1,exp_tag='dge12',design_str="~ comp",reuse=0,debug=0) {
	if (debug==1){browser()}
	
	out_rd <- file.path(outd,'deseq2.rd')
	
	if (FALSE & file.exists(out_rd)) {
		return(NA)
	}
	
	message('storing Ensembl Gene ID and Gene Symbol mapping table ...')
	
	#TOD: the following step to get tx2gene should b moved to lib_salmon.r
	if (is.na(gtf_gene_dt)) {
		feature_map <- data.table(orig=cnt_mat$meta_row,
															ensembl=cnt_mat$meta_row)
	} else {
		cnt_mat.ensembl <- gtf_gene_dt[match(cnt_mat$meta_row,gtf_gene_dt$symbol),ensembl]
		
		feature_map <- data.table(orig=cnt_mat$meta_row,
															ensembl=cnt_mat.ensembl)
		keep <- !is.na(feature_map$ensembl)
		
		cnt_mat$count <- cnt_mat$count[keep,]
		cnt_mat$meta_row <- feature_map[keep,ensembl]
		rownames(cnt_mat$count) <- cnt_mat$meta_row
	}
	
	cnt_mat$meta_col$comp <- factor(cnt_mat$meta_col$comp)
	cnt_mat$meta_col$Timepoint <- factor(cnt_mat$meta_col$Timepoint)
	if (debug==1){browser()}
	message(design_str)
	
	# dds <- DESeqDataSetFromMatrix(cnt_mat$count,
	# 															colData=cnt_mat$meta_col,
	# 															as.formula(design_str))
	
	dds <- DESeqDataSetFromMatrix(cnt_mat$count,
																colData=cnt_mat$meta_col,
																~ comp)
	
	# -------------
	# prefiltering
	message('keep gene features that has 10+ reads total')
	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep,]
	
	# -------------
	message('read count in log10() across samples ...')
	pdf_file <- file.path(outd,sprintf('01_read_count_log10.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	boxplot(log10(counts(dds)+1),ylab='log10(read_count)+1')
	dev.off()
	
	if (plot_pca==1) {
		message('do some PCA analysis ...')
		
		if (dim(dds)[1]<1000){
			vsd <- varianceStabilizingTransformation(dds, 
																							 blind = TRUE,
																							 fitType = "parametric")
		} else {
			vsd <- vst(dds)
		}
		# browser()
		pcaData <- plotPCA(vsd, intgroup=c("comp","Timepoint"),returnData=TRUE)
		graphics.off()
		
		pdf_file <- file.path(outd,sprintf('02_pca_from_vst.pdf'))
		pdf(pdf_file, width=8.5, height=11)
		
		percentVar <- round(100 * attr(pcaData, "percentVar"))
		
		p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = comp, shape=Timepoint)) +
			geom_point(size = 3) +
			xlab(paste0("PC1: ", percentVar[1], "% variance")) +
			ylab(paste0("PC2: ", percentVar[2], "% variance")) +
			coord_fixed()
		
		plot(p)
		dev.off()
		
		# --------------
		message('heatmap of the sample-to-sample distance')
		
		sampleDists <- dist(t(assay(vsd)))
		
		sampleDistMatrix <- as.matrix(sampleDists)
		rownames(sampleDistMatrix) <- vsd$rep_sample
		colnames(sampleDistMatrix) <- NULL
		
		colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
		pdf_file <- file.path(outd,sprintf('03_sample_to_sample_dist_vst.pdf'))
		pdf(pdf_file, width=8.5, height=11)
		
		pheatmap(sampleDistMatrix,
						 clustering_distance_rows=sampleDists,
						 clustering_distance_cols=sampleDists,
						 col=colors)
		dev.off()
	}
	# =============
	message('performing DESeq analysis ...')
	dds <- DESeq(dds)
	# =============
	
	message('check result name ...')
	resultsNames(dds)
	
	res <- results(dds,alpha=0.05,contrast = c('comp','expr','ctrl'))
	
	message('adding original feature name (gene_symbols) ...')
	if (!is.na(gtf_gene_dt)) {
		res$gene_symbol <- feature_map[match(rownames(res),feature_map$ensembl),orig]
	}
	# --------------
	dds <- estimateSizeFactors(dds)
	norm_cnt <- counts(dds, normalized=TRUE)
	
	M <- dim(norm_cnt)[2]
	n_neighbors <- 15
	if (M < 15) {n_neighbors <- round(M * 0.3)}
	
	message(sprintf("target n_neighbors[%d]",n_neighbors))
	umap_out <- umap(t(norm_cnt), n_neighbors=n_neighbors, init="spca")
	umap_out <- data.frame(umap_out)
	colnames(umap_out) <- c('umap1','umap2')
	rownames(umap_out) <- colnames(norm_cnt)
	
	umap_out <- cbind(umap_out,cnt_mat$meta_col[,c('Patient','Timepoint','comp')])
	
	pdf_file <- file.path(outd,sprintf('02b_umap_from_vst.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	
	p <- ggplot(data=umap_out,aes(x=umap1,y=umap2,shape=Timepoint,color=comp)) + geom_point(alpha=0.4,size=4)
	p <- p + ggtitle(exp_tag)
	plot(p)
	dev.off()
	
	# ---------------------

	message('save diff testing result...')
	
	exp_description <- mcols(res)$description
	
	wb <- createWorkbook("deseq2_report")
	sheet_name <- 'readme'
	addWorksheet(wb,sheetName = sheet_name)
	exp_description <- c(exp_description,'min 10 reads;BH alpha=0.05')
	writeData(wb,sheet_name,
						data.table(desc=exp_description))
	
	sheet_name <- 'comp_ctrl_vs_expr'
	addWorksheet(wb,sheetName = sheet_name)
	
	res_dt <- as.data.table(res)
	res_dt$gene_id <- rownames(res)
	
	res_dt <- res_dt[order(padj)]
	
	writeData(wb,sheet_name,
						res_dt)
	
	# --------------
	genes <- rownames(norm_cnt)
	norm_cnt.dt <- as.data.table(norm_cnt)
	norm_cnt.dt$gene_id <- genes
	if (debug==1){browser()}
	norm_cnt.dt$gene_symbol <- feature_map[match(genes,ensembl),orig]
	
	sheet_name <- 'normalized_read_count'
	addWorksheet(wb,sheetName = sheet_name)
	writeData(wb,sheet_name,
						norm_cnt.dt)
	
	saveWorkbook(wb,file=file.path(outd,sprintf('%s_deseq2_results.xlsx',exp_tag)),overwrite=T)
	
	# ----------------------
	
	hsa.grn <- NA
	comb.res <- NA
	if (!is.na(gtf_gene_dt)) {
		message('reporting ...')
		
		rep <- HTMLReport(shortName="deseq2", 
											title=exp_tag,
											basePath=outd,
											reportDirectory=sprintf("%s_dge_table",exp_tag))
		
		if(debug==1) {browser()}
		res2 <- res
		rownames(res2) <- res$gene_symbol
		dds2 <- dds
		rownames(dds2) <- feature_map[match(rownames(dds),feature_map$ensembl),orig]
		
		publish(res2, rep, dds2, n=20, make.plots=TRUE, factor=dds$comp, pvalueCutoff = 1.)
		finish(rep)
		rm(dds2)
		rm(res2)
		# ---------------
		
		#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment
		
		midx <- match(rownames(dds),gtf_gene_dt$ensembl)
		
		rowRanges <- GRanges(seqnames = gtf_gene_dt[midx,chr],
												 strand = gtf_gene_dt[midx,strand],
												 ranges = IRanges(start=gtf_gene_dt[midx,start],
												 								 end=gtf_gene_dt[midx,end]),
												 feature_id = rownames(dds))
		
		colData <- colData(dds)
		
		colData$GROUP <- ifelse(colData$comp == 'expr',1,0)
		# colData$BLOCK <-  colData$subgroup
		
		se <- SummarizedExperiment(assays=list(counts=norm_cnt),
															 rowRanges=rowRanges,
															 colData=colData)
		
		row_data <- res[,c('log2FoldChange','stat','pvalue','padj')]
		colnames(row_data) <- c('FC','DESeq2.STAT','PVAL','ADJ.PVAL')
		rowData(se) <- row_data
		
		#Set-based enrichment analy
		
		se2 <- idMap(se, org="hsa", from="ENSEMBL", to="ENTREZID")
		
		# message('get.go.genesets() ...')
		# go.gs <- getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db",cache=TRUE)
		
		query_db <- c("go","kegg")
		
		#Network-based enrichment analy
		hsa.grn <- compileGRN(org="hsa",db="kegg")
		#head(hsa.grn)
		
		if (debug==1) {browser()}
		
		for (qdb in query_db) {
			message('getGenesets() ...')
			
			if (qdb == "kegg") {
				hsa.gs <- getGenesets(org="hsa",db=qdb, cache=TRUE)
			} else {
				hsa.gs <- getGenesets(org="hsa", db=qdb, go.onto="BP", go.mode="GO.db", cache=TRUE)
			}
			
			#sbeaMethods()
			message('run sbea(x,se2) ...')
			sbea.res <- sbea(method="gsea", se=se2, gs=hsa.gs, perm=0, alpha=0.1)
			#gsRanking(sbea.res)
			
			
			#nbeaMethods()
			message('run nbea(x,se2) ...')
			nbea.res <- nbea(method="ggea", se=se2, gs=hsa.gs, grn=hsa.grn)
			#gsRanking(nbea.res)
			
			res.list <- list('gsea'=sbea.res, 'ggea'=nbea.res)
			comb.res <- combResults(res.list)
			
			graphics.off()
			if (debug==1) {browser()}
			eaBrowse(comb.res, graph.view=hsa.grn, nr.show=10, out.dir=file.path(outd,sprintf('eabrowse_%s',qdb)),html.only=TRUE,report.name=sprintf("%s_%s",exp_tag,qdb))
		}
	}
	# --------------
	message('MA plot ...')
	pdf_file <- file.path(outd,sprintf('04_ma.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotMA(res,colNonSig='blue')
	dev.off()
	
	# --------------
	message('dispersion ...')
	pdf_file <- file.path(outd,sprintf('05_dispersion.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotDispEsts(dds)
	dev.off()
	
	# --------------
	message('data quality assessment using heatmap of the count matrix ...')
	
	select <- order(rowMeans(norm_cnt),
									decreasing=TRUE)[1:20]
	
	df <- as.data.frame(colData(dds)[,c("comp","sample")])
	if (plot_pca==1) {
		pdf_file <- file.path(outd,sprintf('06_sample_to_sample_dist_norm_topk.pdf'))
		pdf(pdf_file, width=8.5, height=11)
		pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
						 cluster_cols=FALSE, annotation_col=df)
		dev.off()
	}
	# --------------
	message('generating p-value distribution ...')
	pdf_file <- file.path(outd,sprintf('07_padj_hist.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	
	hist(res$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 (adj) P-value", ylab="Number of genes")
	dev.off()
	
	# --------------
	message('the most diff gene read count ...')
	pdf_file <- file.path(outd,sprintf('08_the_most_diff_gene_read_count.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotCounts(dds, gene=which.min(res$padj), intgroup=c('comp'))
	dev.off()
	
	# --------------
	if (debug==1){browser()}
	message('volcano plot ...')
	alpha <- 0.05 # Threshold on the adjusted p-value
	cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
	
	pdf_file <- file.path(outd,sprintf('09_volcano.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	
	plot(res$log2FoldChange, 
			 -log10(res$padj), 
			 col=cols,
			 panel.first=grid(),
			 main="Volcano plot", 
			 xlab="Effect size: log2(fold-change)", 
			 ylab="-log10(adjusted p-value)",
			 pch=20, cex=0.6)
	
	abline(v=0)
	abline(v=c(-1,1), col="brown")
	abline(h=-log10(alpha), col="brown")
	
	gn.selected <- abs(res$log2FoldChange) > 2 & res$padj < alpha
	if (any(gn.selected)) {
		text(res$log2FoldChange[gn.selected],
				 -log10(res$padj)[gn.selected],
				 lab=res$gene_symbol[gn.selected], cex=0.4)
	}
	dev.off()
	graphics.off()
	
	
	
	save(outd,feature_map,dds,comb.res,hsa.grn,file=out_rd,compress = T)
	
	if (debug==1){browser()}
	
	return(out_rd)
	
}

salmons_to_tximport <- function(salmon_dt2,outd,reuse=0,debug2=0) {
	
	if (debug2==1) {browser()}
	
	salmon_reformatted<-file.path(outd,'salmon_reformatted')
	if (!file.exists(salmon_reformatted)) {
		dir.create(salmon_reformatted,showWarnings = FALSE, recursive = TRUE)
	}
	
	salmon_dt2$parsed_out1 <- sapply(1:dim(salmon_dt2)[1],function(i) {
		fpath <- salmon_dt2$out1[i]
		fname <- decompose_fname(fpath)
		parsed_annot_fn <- file.path(salmon_reformatted,sprintf('%s_quant.sf',salmon_dt2$sample[i]))
		
		if (reuse==1 & file.exists(parsed_annot_fn)) {
			message(sprintf('reuse [%s]...',parsed_annot_fn))
		} else {
			message(sprintf('splitting tx/geneID in quant.sf[%s] ...',fpath))
			
			dt <- fread(fpath)
			if (debug2==1){browser()}
			
			annots <- tstrsplit(dt$Name,'\\|')
			dt$Name <- annots[[1]] #ensemb txid
			
			fwrite(dt,file=parsed_annot_fn,sep='\t')
			message(sprintf('Check[%s] ...',parsed_annot_fn))
		}
		return(parsed_annot_fn)
	},simplify=TRUE)
	
	message('storing Ensembl Gene ID and Gene Symbol mapping table ...')
	ignore_gene_ver <- TRUE
	fpath <- salmon_dt2$out1[1]
	dt <- fread(fpath)
	annots <- tstrsplit(dt$Name,'\\|')
	dt$Name <- annots[[1]] #ensemb txid
	
	if (ignore_gene_ver) {
		dt$gene_id <- str_extract(annots[[2]],"(?<=^)[:alnum:]+")
		
	} else {
		dt$gene_id <- annots[[2]]
	}
	
	dt$geneSymbol <- annots[[6]]

	tx2gene <- unique(dt[,list(Name,gene_id,geneSymbol)])
	colnames(tx2gene) <- c('TXNAME','GENEID','GENESYMBOL')
	tx2gene$Length <- dt[match(tx2gene$TXNAME,Name),Length]
	
	
	if (F) {
		geneid_symbol_map2 <- ensembldb::select(EnsDb.Hsapiens.v79,
																						 keys= unique(tx2gene$GENEID), 
																						 keytype = "GENEID", 
																						 columns = c("GENEID","SYMBOL"))
		
		geneid_symbol_map2 <- as.data.table(geneid_symbol_map2)
		
		tx2gene$GENESYMBOL2 <- geneid_symbol_map2[match(tx2gene$GENEID,GENEID),SYMBOL]
		tx2gene[is.na(GENESYMBOL2),GENESYMBOL2:=GENESYMBOL]
	}
	#load salmon data
	files <- salmon_dt2$parsed_out1
	names(files) <- salmon_dt2$sample
	if (debug2==1) {browser()}
	txi <- tximport(files,
									type='salmon',
									tx2gene = tx2gene[,c('TXNAME','GENEID')])
	
	if (debug2==1) {browser()}
	
	return(list(txi=txi,tx2gene=tx2gene))
}

get_canonical_geneid_from_salmon_tx2gene <- function(tx2gene) {
	#tx2gene has the required columns, (GENEID,GENESYMBOL,Length)
	ens2symb <- unique(tx2gene[,c('GENEID','GENESYMBOL')])
	ens2symb$Length <- tx2gene[match(ens2symb$GENEID,GENEID),Length]
	
	maxLen_GeneSymbol <- ens2symb[,.(maxLength=max(Length)),by="GENESYMBOL"]
	
	ens2symb2 <- merge(maxLen_GeneSymbol,
										 ens2symb,
										 by.x=c("GENESYMBOL","maxLength"),
										 by.y=c("GENESYMBOL","Length"))
	ens2symb2
}

dge_analysis <- function(salmon_dt2,outd,exp_tag='dge12',reuse=T,debug2=0) {
	if (debug2==1){browser()}
	
	out_rd <- file.path(outd,'deseq2.rd')
	
	if (FALSE & file.exists(out_rd)) {
		return(NA)
	}
	
	salmon_reformatted<-file.path(outd,'salmon_reformatted')
	
	dir.create(salmon_reformatted,showWarnings = FALSE, recursive = TRUE)
	
	salmon_dt2$parsed_out1 <- sapply(1:dim(salmon_dt2)[1],function(i) {
		#fpath <- salmon_dt2$out_g4[i]
		fpath <- salmon_dt2$out1[i]
		fname <- decompose_fname(fpath)
		parsed_annot_fn <- file.path(salmon_reformatted,sprintf('%s_quant.sf',salmon_dt2$sample[i]))
		
		# parsed_annot_fn <- sprintf('%s_parsed_annot.txt',fpath)
		if (reuse & file.exists(parsed_annot_fn)) {
			message('reuse ...')
		} else {
			message(sprintf('splitting tx/geneID in quant.sf[%s] ...',fpath))
			
			dt <- fread(fpath)
			annots <- tstrsplit(dt$Name,'\\|')
			dt$Name <- annots[[1]] #ensemb txid
			
			fwrite(dt,file=parsed_annot_fn,sep='\t')
			message(sprintf('Check[%s] ...',parsed_annot_fn))
		}
		return(parsed_annot_fn)
	},simplify=TRUE)
	
	message('storing Ensembl Gene ID and Gene Symbol mapping table ...')
	
	#TOD: the following step to get tx2gene should b moved to lib_salmon.r
	if (F) {
		txdb <- makeTxDbFromGFF('~/hwangt-share/Datasets/Reference_seqs/GRCh38_P01/GTF/human_GENCODEv29/gencode.v29.annotation.gtf.gz')
		
		k <- keys(txdb, keytype = "TXNAME")
		tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
	} else {
		ignore_gene_ver <- TRUE
		fpath <- salmon_dt2$out1[1]
		dt <- fread(fpath)
		annots <- tstrsplit(dt$Name,'\\|')
		dt$Name <- annots[[1]] #ensemb txid
		
		if (ignore_gene_ver) {
			dt$gene_id <- str_extract(annots[[2]],"(?<=^)[:alnum:]+")
			
		} else {
			dt$gene_id <- annots[[2]]
		}
		
		dt$geneSymbol <- annots[[6]] 
		
		tx2gene <- unique(dt[,list(Name,gene_id,geneSymbol)])
		colnames(tx2gene) <- c('TXNAME','GENEID','GENESYMBOL')
	}
	
	#load salmon data
	
	files <- salmon_dt2$parsed_out1
	names(files) <- salmon_dt2$sample
	
	txi <- tximport(files,
									type='salmon',
									geneIdCol = 'gene_id',
									tx2gene = tx2gene[,c('TXNAME','GENEID')])
	
	stopifnot(all(colnames(txi$counts) %in% salmon_dt2$sample))
	
	#create DESeq object
	message("[TODO] to include the other variables in the design formula from user input")
	dds <- DESeqDataSetFromTximport(txi,
																	colData = salmon_dt2,
																	design = ~ comp)
	
	# -------------
	# prefiltering
	message('keep gene features that has 10+ reads total')
	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep,]
	
	# -------------
	
	message('read count in log10() across samples ...')
	pdf_file <- file.path(outd,sprintf('01_read_count_log10.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	boxplot(log10(counts(dds)+1),ylab='log10(read_count)+1')
	dev.off()
	
	message('do some PCA analysis ...')
	vsd <- vst(dds)
	
	pcaData <- plotPCA(vsd, intgroup=c("cond1","cond2"),returnData=TRUE)
	graphics.off()
	
	pdf_file <- file.path(outd,sprintf('02_pca_from_vst.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	
	p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = cond1, shape = cond2)) +
		geom_point(size =3) +
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		coord_fixed()
	
	plot(p)
	dev.off()
	# --------------
	message('heatmap of the sample-to-sample distance')
	
	sampleDists <- dist(t(assay(vsd)))
	
	sampleDistMatrix <- as.matrix(sampleDists)
	rowAnnot<- as.data.frame(salmon_dt2[match(rownames(sampleDistMatrix),salmon_dt2$sample),list(cond1,cond2)])
	rownames(rowAnnot) <- rownames(sampleDistMatrix)
	
	# colAnnot<- salmon_dt2[match(colnames(sampleDistMatrix),salmon_dt2$sample),list(sample,cond1,cond2)]
	
	colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
	pdf_file <- file.path(outd,sprintf('03_sample_to_sample_dist_vst.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	
	pheatmap(sampleDistMatrix,
					 clustering_distance_rows=sampleDists,
					 clustering_distance_cols=sampleDists,
					 col=colors,
					 annotation_row=rowAnnot)
	
	dev.off()
	
	# =============
	message('performing DESeq analysis ...')
	dds <- DESeq(dds)
	# =============
	
	# add interaction to see if the result is dependent on subgroup 
	if (F) {
		message('performing DESeq analysis (interact atrributes)...')
		dds.ia <- dds
		dds.ia$group <- factor(paste0(dds.ia$cond2,'.',dds.ia$comp))
		design(dds.ia) <- ~ group
		
		dds.ia <- DESeq(dds.ia)
		resultsNames(dds.ia)
		res.ia <- results(dds.ia)
	}
	
	message('check result name ...')
	resultsNames(dds)
	
	res <- results(dds,alpha=0.05,contrast = c('comp','expr','ctrl'))
	if (debug2==1){browser()}
	message('gene name conversion ...')
	ens2symb <- unique(tx2gene[,c('GENEID','GENESYMBOL')])
	res$gene_symbol <- ens2symb[match(rownames(res),ens2symb$GENEID),GENESYMBOL]
	
	
	# ---------------
	
	#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment
	
	G <- dim(dds)[1]
	
	norm_cnt <- counts(dds,normalized=TRUE)
	
	#TODO: use a grange for the actual gene model instead of this dummy variable
	rowRanges <- GRanges(rep('chr1',G),
											 IRanges(floor(runif(G,1e4,1e7)),width=100),
											 strand=sample(c('+','-'),G,TRUE),
											 feature_id=rownames(dds))
	
	colData <- colData(dds)
	
	colData$GROUP <- ifelse(colData$comp == 'expr',1,0)
	# colData$BLOCK <-  colData$subgroup
	
	se <- SummarizedExperiment(assays=list(counts=norm_cnt),
														 rowRanges=rowRanges,
														 colData=colData)
	
	row_data <- res[,c('log2FoldChange','stat','pvalue','padj')]
	colnames(row_data) <- c('FC','DESeq2.STAT','PVAL','ADJ.PVAL')
	rowData(se) <- row_data
	
	message("Set-based enrichment analy ...")
	
	se2 <- idMap(se, org="hsa", from="ENSEMBL", to="ENTREZID")
	
	
	# message('get.go.genesets() ...')
	# go.gs <- getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db",cache=TRUE)
	
	query_db <- c("go","kegg")
	
	#Network-based enrichment analy
	hsa.grn <- compileGRN(org="hsa",db="kegg")
	#head(hsa.grn)
	
	if (debug2==1) {browser()}
	
	for (qdb in query_db) {
		message('getGenesets() ...')
		
		if (qdb == "kegg") {
			hsa.gs <- getGenesets(org="hsa",db=qdb, cache=TRUE)
		} else {
			hsa.gs <- getGenesets(org="hsa", db=qdb, go.onto="BP", go.mode="GO.db", cache=TRUE)
		}
		
		#sbeaMethods()
		message('run sbea(x,se2) ...')
		sbea.res <- sbea(method="gsea", se=se2, gs=hsa.gs, perm=0, alpha=0.1)
		#gsRanking(sbea.res)
		
		
		#nbeaMethods()
		message('run nbea(x,se2) ...')
		nbea.res <- nbea(method="ggea", se=se2, gs=hsa.gs, grn=hsa.grn)
		#gsRanking(nbea.res)
		
		if (debug2==1) {browser()}
		
		res.list <- list('gsea'=sbea.res, 'ggea'=nbea.res)
		comb.res <- combResults(res.list)
		
		graphics.off()
		
		eaBrowse(comb.res, graph.view=hsa.grn, nr.show=10, out.dir=file.path(outd,sprintf('eabrowse_%s',qdb)),html.only=TRUE,report.name=sprintf("%s_%s",exp_tag,qdb))
		
	}
	
	# --------------
	message('MA plot ...')
	pdf_file <- file.path(outd,sprintf('04_ma.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotMA(res,colNonSig='blue')
	dev.off()
	
	# --------------
	message('dispersion ...')
	pdf_file <- file.path(outd,sprintf('05_dispersion.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotDispEsts(dds)
	dev.off()
	# --------------
	message('data quality assessment using heatmap of the count matrix ...')
	
	select <- order(rowMeans(norm_cnt),
									decreasing=TRUE)[1:20]
	
	df <- as.data.frame(colData(dds)[,c("comp","cond2")])
	
	pdf_file <- file.path(outd,sprintf('06_sample_to_sample_dist_norm_topk.pdf'))
	pdf(pdf_file, width=8.5, height=11)
	pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
					 cluster_cols=FALSE, annotation_col=df)
	dev.off()
	
	# --------------
	message('generating p-value distribution ...')
	pdf_file <- file.path(outd,sprintf('07_padj_hist.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	
	hist(res$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 (adj) P-value", ylab="Number of genes")
	dev.off()
	
	# --------------
	message('the most diff gene read count ...')
	pdf_file <- file.path(outd,sprintf('08_the_most_diff_gene_read_count.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	plotCounts(dds, gene=which.min(res$padj), intgroup=c('comp'))
	dev.off()
	
	# --------------
	message('volcano plot ...')
	alpha <- 0.05 # Threshold on the adjusted p-value
	cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
	
	pdf_file <- file.path(outd,sprintf('09_volcano.pdf'))
	pdf(pdf_file, width=11, height=8.5)
	
	plot(res$log2FoldChange, 
			 -log10(res$padj), 
			 col=cols,
			 panel.first=grid(),
			 main="Volcano plot", 
			 xlab="Effect size: log2(fold-change)", 
			 ylab="-log10(adjusted p-value)",
			 pch=20, cex=0.6)
	
	abline(v=0)
	abline(v=c(-1,1), col="brown")
	abline(h=-log10(alpha), col="brown")
	
	gn.selected <- abs(res$log2FoldChange) > 2 & res$padj < alpha 
	text(res$log2FoldChange[gn.selected],
			 -log10(res$padj)[gn.selected],
			 lab=res$gene_symbol[gn.selected], cex=0.4)
	
	dev.off()
	graphics.off()
	#----------
	message('reporting ...')
	
	tmp <- tempdir(outd)
	rep <- HTMLReport(shortName="deseq2", 
										title="Ovarian Cancer",
										basePath=tmp,
										reportDirectory="html_report")
	
	publish(res, rep, dds, n=20, make.plots=TRUE, factor=dds$comp)
	finish(rep)
	
	# -----------
	message('save diff testing result...')
	
	exp_description <- mcols(res)$description
	
	wb <- createWorkbook("deseq2_report")
	sheet_name <- 'readme'
	addWorksheet(wb,sheetName = sheet_name)
	exp_description <- c(exp_description,'min 10 reads;BH alpha=0.05')
	writeData(wb,sheet_name,
						data.table(desc=exp_description))
	
	sheet_name <- 'comp_expr_vs_ctrl'
	addWorksheet(wb,sheetName = sheet_name)
	
	res_dt <- as.data.table(res)
	res_dt$gene_id <- rownames(res)
	
	res_dt <- res_dt[order(padj)]
	
	writeData(wb,sheet_name,
						res_dt)
	
	# --------------
	genes <- rownames(norm_cnt)
	norm_cnt <- as.data.table(norm_cnt)
	norm_cnt$gene_id <- genes
	
	sheet_name <- 'normalized_read_count'
	addWorksheet(wb,sheetName = sheet_name)
	writeData(wb,sheet_name,
						norm_cnt)
	
	saveWorkbook(wb,file=file.path(outd,'deseq2_results.xlsx'),overwrite=T)
	
	save(outd,txi,dds,comb.res,hsa.grn,file=out_rd,compress = T)
	
	if (debug2==1){browser()}
	
	return(res)
	
}

counts_to_cpm <- function(counts,n=1e6,debug2=0) {
	if (debug2==1){browser()}
	# browser()
	# Process one column at a time.
	cpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
		if (debug2==1){browser()}
		rate = log(counts[,i])
		denom = log(sum(exp(rate)))
		ret <- round(exp(rate - denom + log(n)))
		if (all(is.nan(ret))) {
			ret <- round(n/length(ret))
		}
		return(ret)
	}))
	
	# Copy the row and column names from the original matrix.
	colnames(cpm) <- colnames(counts)
	rownames(cpm) <- rownames(counts)
	# browser()
	return(cpm)
}

#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
	# browser()
	# Ensure valid arguments.
	stopifnot(length(featureLength) == nrow(counts))
	stopifnot(length(meanFragmentLength) == ncol(counts))
	
	# Compute effective lengths of features in each library.
	effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
		featureLength - meanFragmentLength[i] + 1
	}))
	
	# Exclude genes with length less than the mean fragment length.
	idx <- apply(effLen, 1, function(x) min(x) > 1)
	counts <- counts[idx,]
	effLen <- effLen[idx,]
	featureLength <- featureLength[idx]
	
	# Process one column at a time.
	tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
		rate = log(counts[,i]) - log(effLen[,i])
		denom = log(sum(exp(rate)))
		exp(rate - denom + log(1e6))
	}))
	
	# Copy the row and column names from the original matrix.
	colnames(tpm) <- colnames(counts)
	rownames(tpm) <- rownames(counts)
	# browser()
	return(tpm)
}


get_gene_len_from_gtf <- function(GTFfile,reuse=0) {
	#Load the annotation and reduce it
	message(sprintf("reading %s ...",GTFfile))
	rd_file <- sprintf("%s.rd",GTFfile)
	if (reuse==1 & file.exists(rd_file)) {
		load(rd_file)
	} else {
		anno <- import.gff(GTFfile, format="gtf", genome="GRCh38.p12", feature.type="exon")
		
		exons <- anno[anno@elementMetadata$type %in% c("exon","five_prime_UTR","three_prime_UTR"),]
		
		# splitting up isoforms as preparation for the next step
		tmp <- split(exons,elementMetadata(anno)$transcript_id)
		# for each isoform, calculate the sum of all reduced exons
		tx_length <- sum(width(reduce(tmp)))
		
		tx_len_dt <- data.table(transcript_id=names(tx_length),
														len=tx_length)
		
		tx_to_gene_map_dt <- unique(data.table(gene_name=elementMetadata(anno)$gene_name,
																					 transcript_id=elementMetadata(anno)$transcript_id,
																					 ensembl_gene_id=elementMetadata(anno)$gene_id))
		
		tx_len_dt$gene_name <- tx_to_gene_map_dt[match(tx_len_dt$transcript_id,tx_to_gene_map_dt$transcript_id),gene_name]
		tx_len_dt$ensembl_gene_id <- tx_to_gene_map_dt[match(tx_len_dt$transcript_id,tx_to_gene_map_dt$transcript_id),ensembl_gene_id]
		
		gtf_tx_list <- list(tx_len=tx_len_dt,
												gene_len=tx_len_dt[,mean(len),by=gene_name])
		save(gtf_tx_list,file=rd_file,compress=T)
	}
	return(gtf_tx_list)
}



staralign <- function(args,wks,action="script",debug2=1) {
	if (debug2==1){browser()}
	
	if (args$outd=="") {
		outd <- "starout"
	} else {
		outd <- normalizePath(args$outd)
	}
	
	if (!file.exists(outd)) {
		message(sprintf("creating [%s]...",outd))
		dir.create(outd,recursive = T)
	}
	
	message('STAR ...')
	cmethod <- get_prog_info_yaml('star',debug2=debug2)
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		# message(i)
		cmd <- cmethod$bin_path
		cmd <- sprintf("%s --genomeDir %s",cmd,cmethod$resource$hg19_path)
		cmd <- sprintf("%s --runThreadN %d",cmd,cmethod$ncpu)
		
		if (args$inputd!="") {
			r1fpath <- file.path(args$inputd,wks$read1[i])
			r2fpath <- file.path(args$inputd,wks$read2[i])
		} else {
			r1fpath <- wks$read1[i]
			r2fpath <- wks$read2[i]
		}
		
		cmd <- sprintf("%s --readFilesIn %s %s",cmd,r1fpath,r2fpath)
		# message(cmd)
		outdi <- file.path(outd,wks$sample[i])
		if (!file.exists(outdi)) {
			dir.create(outdi,recursive = T)
		}
		out_prefix <- file.path(normalizePath(outdi),"star_")
		
		cmd <- sprintf("%s --outFileNamePrefix %s",cmd,out_prefix)
		cmd <- sprintf("%s --outSAMtype BAM SortedByCoordinate",cmd)
		cmd <- sprintf("%s --outSAMunmapped Within",cmd)
		cmd <- sprintf("%s --outSAMattributes Standard",cmd)
		cmd <- sprintf("%s --twopassMode Basic",cmd)
		
		if (endsWith(r1fpath,".gz")) {
			cmd <- sprintf("%s --readFilesCommand zcat",cmd)
		}
		
		cmd <- sprintf("%s --outSAMattrRGline ID:%s.1 PL:Illumina PU:pu1 SM:%s",cmd, wks$sample[i], wks$sample[i])
		
		# message(cmd)
		
		# message(str(cmethod)) #debug
		
		if (debug2==1){browser()}
		
		if ("default_opt" %in% names(cmethod)) {
			if (!is.null(cmethod$default_opt)) {
				cmd <- sprintf("%s %s",cmd,cmethod$default_opt)
			}
		}
		
		star_bam <- sprintf("%sAligned.sortedByCoord.out.bam",out_prefix)
		
		# message(cmd)
		if (file.exists(star_bam)) {cmd <- sprintf("## reuse [%s]",star_bam)}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		
		return(list(sample=wks$sample[i],read1=r1fpath,read2=r2fpath,cmd=cmd,out1=star_bam))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														out1=sapply(1:M,function(m) {ret[[m]]$out1}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

GSEA_with_tryCatch <- function(geneList, 
															 GO_file, 
															 pval, 
															 desc,
															 use.gage,
															 debug2) {
	gsea_res <- tryCatch(
		{
			gsea_res <- GSEA_fgsea(geneList, 
											 GO_file, 
											 desc=desc,
											 use.gage=use.gage,
											 debug2=debug2)
		},
		error=function(cond){
			message(sprintf('%s,%s',GO_file,desc))
			return(NA)
		},
		finally={
			message("Done.")
		}
	)
	return(gsea_res)
}

GSEA_fgsea <- function(gene_list, gmt_file, pval=0.05, use.gage=0, desc="NA", debug2=0) {
	if (debug2==1){browser()}
	
	if ( any( duplicated(names(gene_list)) )  ) {
		warning("Duplicates in gene names")
		gene_list = gene_list[!duplicated(names(gene_list))]
	}
	
	stopifnot(any(is.finite(gene_list)))
	max_logfc <- max(gene_list[is.finite(gene_list)])
	min_logfc <- min(gene_list[is.finite(gene_list)])
	gene_list[is.infinite(gene_list) & gene_list>0] <- max_logfc + 1
	gene_list[is.infinite(gene_list) & gene_list<0] <- min_logfc - 1
	if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
		warning("Gene list not sorted")
		gene_list = sort(gene_list, decreasing = TRUE)
	}
	
	# message(gmt_file)
	myGmt = fgsea::gmtPathways(gmt_file)
	
	fgRes <- fgsea::fgsea(pathways = myGmt,
												stats = gene_list,
												minSize=3,
												maxSize=600,
												nperm=10000) %>%
		as.data.frame() %>%
		dplyr::filter(padj < !!pval)
	# #print(dim(fgRes))
	
	# fgRes <- fgsea::fgseaMultilevel(pathways = myGmt,
	# 																stats = gene_list,
	# 																minSize=3,
	# 																maxSize=600) %>%
	# 	as.data.frame() %>%
	# 	dplyr::filter(padj < !!pval)
	
	if (use.gage==1) {
		## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
		message("Filter FGSEA by using gage results. Must be significant and in same direction to keep ...")
		gaRes = gage::gage(gene_list, gsets=myGmt, same.dir=TRUE, set.size =c(3,600))
		
		ups = as.data.frame(gaRes$greater) %>% 
			tibble::rownames_to_column("Pathway") %>% 
			dplyr::filter(!is.na(p.geomean) & p.val < pval ) %>%
			dplyr::select("Pathway")
		
		downs = as.data.frame(gaRes$less) %>% 
			tibble::rownames_to_column("Pathway") %>% 
			dplyr::filter(!is.na(p.geomean) & p.val < pval ) %>%
			dplyr::select("Pathway")
		
		#print(dim(rbind(ups,downs)))
		keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
		keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
		
		### Collapse redundant pathways
		Up = fgsea::collapsePathways(keepups, pathways = myGmt, stats = gene_list,  nperm = 500, pval.threshold = pval)
		Down = fgsea::collapsePathways(keepdowns, myGmt, gene_list,  nperm = 500, pval.threshold = pval) 
		
		fgRes = fgRes[ !is.na(match(fgRes$pathway, 
																c( Up$mainPathways, Down$mainPathways))), ] %>% 
			arrange(desc(NES))
	}
	
	if (debug2==1){browser()}
	
	fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
	
	if (F) {
		filtRes = rbind(head(fgRes, n = 10),
										tail(fgRes, n = 10 ))
		
		g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
			geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
			geom_point( size=5, aes( fill = Enrichment),
									shape=21, stroke=2) +
			scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
																	 "Up-regulated" = "firebrick") ) +
			coord_flip() +
			labs(x="GO", y="Normalized Enrichment Score",
					 title=plot_title) + 
			theme_minimal()
		
		output = list("Results" = fgRes, "Plot" = g)
	}
	
	if (nrow(fgRes)>0) {
		fgRes <- as.data.table(fgRes)
		fgRes$desc <- desc
	} else {
		fgRes <- NA
	}
	# message("Done.")
	return(fgRes)
}


## draw volcanoplot : from https://www.biostars.org/p/309962/
EnhancedVolcanoDESeq2 <- function(toptable,column2disp,AdjustedCutoff=0.05, LabellingCutoff=0.05, LFCCutoff=2.0, main="DoubleKD VolcanoPlot",debug2=0)
{
	if (debug2==1){browser()}
	toptable$Significance <- "NS"
	toptable$Significance[(abs(toptable$logFC) > LFCCutoff)] <- "FC"
	toptable$Significance[(toptable$adj.P.Val<AdjustedCutoff)] <- "FDR"
	toptable$Significance[(toptable$adj.P.Val<=AdjustedCutoff) & (abs(toptable$logFC)>=LFCCutoff)] <- "FC_FDR"
	table(toptable$Significance)
	toptable$Significance <- factor(toptable$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
	
	xlim_min <- min(toptable$logFC)*1.2
	xlim_max <- max(toptable$logFC)*1.2
	
	subind <- (toptable$adj.P.Val<=LabellingCutoff & abs(toptable$logFC)>=LFCCutoff)

	p <- ggplot(toptable, aes(x=logFC, y=-log10(adj.P.Val))) +
		
		#Add points:
		#   Colour based on factors set a few lines up
		#   'alpha' provides gradual shading of colour
		#   Set size of points
		geom_point(aes(color=factor(Significance)), alpha=1/2, size=0.8) +
		
		#Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
		scale_color_manual(values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
											 labels=c(NS="NS",
											 				 FC=sprintf("LogFC>|%4.3f|",LFCCutoff),
											 				 FDR=sprintf("FDR adj.Pval<%4.3f",AdjustedCutoff),
											 				 FC_FDR=sprintf("LogFC>=|%4.3f| & FDR adj.Pval<=%4.3f",LFCCutoff,AdjustedCutoff))) +

		#Set the size of the plotting window
		theme_bw(base_size=24) +
		
		#Modify various aspects of the plot text and legend
		theme(legend.background=element_rect(),
					plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
					
					panel.grid.major=element_blank(), #Remove gridlines
					panel.grid.minor=element_blank(), #Remove gridlines
					
					axis.text.x=element_text(angle=0, size=12, vjust=1),
					axis.text.y=element_text(angle=0, size=12, vjust=1),
					axis.title=element_text(size=12),
					
					#Legend
					legend.position="top",            #Moves the legend to the top of the plot
					legend.key=element_blank(),       #removes the border
					legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
					legend.text=element_text(size=8), #Text size
					title=element_text(size=8),       #Title text size
					legend.title=element_blank()) +       #Remove the title
		
		#Change the size of the icons/symbols in the legend
		guides(colour = guide_legend(override.aes=list(size=2.5))) +
		
		#Set x- and y-axes labels
		xlab(bquote(~Log[2]~ "fold change")) +
		ylab(bquote(~-Log[10]~adjusted~italic(P))) +

		ggtitle(main)
	
	p <- p +
		#Tidy the text labels for a subset of genes
		geom_text_repel(data=toptable[subind,],
										aes(label=toptable[subind,get(column2disp)])) +
		
	# p <- p +
	# 	#Tidy the text labels for a subset of genes
	# 	geom_text(data=toptable[subind,],
	# 						aes(label=toptable[subind,get(column2disp)]),
	# 						size=3,
	# 						#size=2.25,
	# 						#segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
	# 						#segment.size=0.01,
	# 						check_overlap=TRUE,
	# 						vjust=1.0) +
		
		#Add a vertical line for fold change cut-offs
		geom_vline(xintercept=c(-LFCCutoff, LFCCutoff), linetype="longdash", colour="black", size=0.4) +
		
		#Add a horizontal line for P-value cut-off
		geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)
	
	p <- p +
		xlim(c(xlim_min,xlim_max))
	
	return(p)
}


select_highly_variable_genes <- function(log2expr_mat,cutoff=0.25) {
	cv <- apply(log2expr_mat, 1, function(x){sd(x) / mean(x)})
	summary(cv)
	
	summary(cv[rank(cv) / length(cv) > 1 - cutoff])
	hvgs_mat <- log2expr_mat[rank(cv) / length(cv) > 1 - cutoff, ]
	dim(hvgs_mat)
	hvgs_mat
}


get_gc_len_on_gtf <- function(gtf_fbase="/media/sammy/apps/cm/resource/human_GENCODEv19/gencode.v19.annotation") {
	library(BSgenome.Hsapiens.UCSC.hg19)

	rds_file <- sprintf("%s_gclen.rds",gtf_fbase)
	if (file.exists(rds_file)) {
		ebg.red <- readRDS(rds_file)
	} else {
		gtf_fname <- sprintf("%s.gtf",gtf_fbase)
		txdb <- makeTxDbFromGFF(gtf_fname)
		ebg <- exonsBy(txdb, by="gene")
		head(names(ebg))
		
		ebg.red <- reduce(ebg)
		
		dna <- extractTranscriptSeqs(Hsapiens, ebg.red)
		all(sum(width(ebg.red)) == width(dna))
		
		mcols(ebg.red)$gc <- as.numeric(letterFrequency(dna, "GC", as.prob=TRUE))
		mcols(ebg.red)$len <- sum(width(ebg.red))
		
		saveRDS(ebg.red,file = rds_file)
	}
	ebg.red
}


get_correlated_variable_genes = function(mat, n = nrow(mat), cor_cutoff = 0, pct_cutoff = 0.1,debug2=0) {
	ind = order(apply(mat, 1, function(x) {
		q = quantile(x, c(0.1, 0.9))
		x = x[x < q[1] & x > q[2]]
		var(x)/mean(x)
	}), decreasing = TRUE)[1:n]
	mat2 = mat[ind, , drop = FALSE]
	dt = cor(t(mat2), method = "spearman")
	diag(dt) = 0
	dt[abs(dt) < cor_cutoff] = 0
	dt[dt < 0] = -1
	dt[dt > 0] = 1
	if(debug2==1){browser()}
	
	n_cutoff = n * pct_cutoff
	message(sprintf("n_cutoff[%g],cor_cutoff[%g]",n_cutoff,cor_cutoff))
	i = colSums(abs(dt)) > n_cutoff
	
	mat3 = mat2[i, ,drop = FALSE]
	return(mat3)
}

scale_and_clip <- function(expriv,scale.center=TRUE) {
	expriv.scaled = t(apply(expriv, 1, function(x) {
		q10 = quantile(x, 0.1)
		q90 = quantile(x, 0.9)
		x[x < q10] = q10
		x[x > q90] = q90
		scale(x,center = scale.center)
	}))
	dim(expriv.scaled)
	colnames(expriv.scaled) = colnames(expriv)
	
	#get rid of any NaN
	j <- rowSums(is.nan(expriv.scaled))==0
	expriv.scaled[j,]
}

get_effective_gene_length <- function() {
	
	rds_fn = file.path(Sys.getenv('REFD'),"human_GENCODEv19","EnsDb.Hsapiens.v86.gene_len.rds")
	
	if (file.exists(rds_fn)) {
		eGeneLen = readRDS(rds_fn)
	} else {
		library(EnsDb.Hsapiens.v86)
		edb <- EnsDb.Hsapiens.v86
		
		eg2symbol = as.data.table(genes(edb,columns=c("gene_id","symbol"),return.type="DataFrame"))
		
		EGene <- exonsBy(edb, by="gene",columns=c("seq_name","symbol","gene_id","exon_id","exon_id","seq_strand"))
		
		effectiveLens <- mclapply(EGene,function(EGenej){
			sum(width(reduce(EGenej)))
		}
		,mc.cores = 12
		)
		
		eGeneLen <- data.table(gene_id=names(effectiveLens),
													 exonBp=unlist(effectiveLens))
		
		eGeneLen$geneSymbol = eg2symbol[match(eGeneLen$gene_id,gene_id),symbol]
		saveRDS(eGeneLen,file=rds_fn)
	}
	eGeneLen
}
