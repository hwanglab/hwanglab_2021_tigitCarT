library(data.table)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(patchwork)
library(ggplotify)
library(dplyr)
library(openxlsx)
library(grid)
library(gridExtra)
library(reshape2)
library(circlize)
library('ComplexHeatmap')
library(ggbeeswarm)
library(AUCell)
library(GSEABase)
library(MASS)

source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_workflow.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_rnaseq.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_seurat3.r'))

# library(randomcoloR)
# source(file.path(Sys.getenv('R_UTIL'),'lib_seurat3.r'))
 
get_projd0 <- function(subd=NA) {
	outd='out'
	if (!is.na(subd)) {
		outd=file.path(outd,subd)
	}
	outd
}

locate_sample_sheet <- function() {
	file.path(get_projd0(),'doc','seqrun_info','cellranger_samplesheet.xlsx')
}

load_sample_worksheet <- function(sheet_name="library",debug2=0) {
	if (debug2==1){browser()}
	# 1. read library and generate library.csv file
	swks <- as.data.table(read.xlsx(locate_sample_sheet(),sheet=sheet_name))
	if (debug2==1){browser()}
	swks
}

load_inndb_genes <- function() {
	inndb_rds_fn <- file.path(get_projd0(),"cart/00a_resource/InnateDB_ImmPort_genes.rds")
	if (file.exists(inndb_rds_fn)) {
		gois <- readRDS_w_msg(inndb_rds_fn)
	} else {
		inndb_imm <- fread(file.path(get_projd0(),'cart/00a_resource/InnateDB_ImmPort_genes.csv'),sep=",")
		gois <- sort(union(get_immreg_genes(extended=1),inndb_imm$name))
		ga <- rownames(seus$P12D14.c2$CD8)
		gois <- ga[ga %in% gois]
		
		saveRDS(gois,file=inndb_rds_fn)
	}
	gois
}

get_immreg_genes <- function(extended=0) {
	# goi <- "CCR7,CREM,FOS,FOSB,FOSL1,FOSL2,GZMM,HAVCR2,IL21,IL21R,JUN,JUNB,JUND,LAG3,NFAT5,NFATC1,NFATC2,NFATC2IP,NFATC3,NFKBIA,NR4A1,NR4A2,NR4A3,PDCD1,PDCD10,PDCD11,PDCD1LG2,SELL,TCF7,HAVCR2,TOX,TOX2,CD226,CD96,CTLA4,FGL1,FGL2,FOXP3,PVR,PVRIG,NECTIN3,VSIR,LRIG1,CXCR5,BTLA,SIGIRR,CD244,CD160,KLRG1,PRDM1,TGFB1,CD274,TIGIT,EOMES,TBX21,GATA3,FOXO1,FOXP1,BATF,IRF4,SPRY2,GZMA,GZMB,GZMK,IFNG,TNF,PRF1,FASLG,IL10,IL21,IL21R,TGFB2,IDO1,CCL5,"
	
	goi <- "CCR7,CREM,FOS,FOSB,FOSL1,FOSL2,GZMM,HAVCR2,IL21,IL21R,JUN,JUNB,JUND,LAG3,NFAT5,NFATC1,NFATC2,NFATC2IP,NFATC3,NFKBIA,NR4A1,NR4A2,NR4A3,PDCD1,PDCD10,PDCD11,PDCD1LG2,SELL,TCF7,HAVCR2,TOX,TOX2,CD226,CD96,CTLA4,FGL1,FGL2,FOXP3,PVR,PVRIG,NECTIN3,VSIR,LRIG1,CXCR5,BTLA,SIGIRR,CD244,CD160,KLRG1,PRDM1,TGFB1,CD274,TIGIT,EOMES,TBX21,GATA3,FOXO1,FOXP1,BATF,IRF4,SPRY2,GZMA,GZMB,GZMK,IFNG,TNF,PRF1,FASLG,IL10,IL21,IL21R,TGFB2,IDO1,CCL5,CD27,CD276,CMKLR1,CXCL9,NKG7,PSMB10,PDL2,STAT1,AKT1,ARG1,B2M,BATF3,BCL2,CCND1,CD3E,CD40,CD47,VEGFA,HIF1A,PTEN,CXCL10,DKK2,EPCAM,FAS,CTNNB1,CSF1R,ICAM1,ITGB2,ITGB8,ITGAM,ITGAV,PECAM1,STAT2,STAT3,KI67"
	
	gois <- comma_string_to_list(goi,sort_flag="ascending")
	
	if (extended==1) {
		gois <- sort(unique(c(gois,c('KLRG1','TBX21','ICOS','TRAF1','TRAF2','TRAF3','TRAF4','TRAIP','NFKB1','TANK','SOCS1','ELOB','RASA2'))))
	}
	gois
}

extract_cell_type <- function(seu,ctype,stype="cart",debug2=0) {
	if (debug2==1){browser()}
	message(sprintf("extracting %s,%s ...",ctype,seu$orig.ident[1]))
	
	if (stype=="pbmc") {
		j <- which(seu@meta.data$pmid_30726743 %in% c('T',sprintf('%s+ T',ctype)))
	} else {
		if (ctype=="CD4") {
			j<-grep("CD8",seu@meta.data$pmid_30726743,invert = TRUE)
		} else {
			j<-grep("CD8",seu@meta.data$pmid_30726743)
		}
	}
	subset(seu,cells=rownames(seu@meta.data)[j])
}

get_CD4_CD8_ctypes <- function(norm.pat_by_cellid,debug2=0) {
	if (debug2==1) {browser()}
	ctype <- colnames(norm.pat_by_cellid)
	if ("B" %in% ctype) {
		ctype[grep("CD4",ctype)] <- "CD4"
		ctype[grep("CD8",ctype)] <- "CD8"
	} else {
		ctype[grep("CD8",ctype)] <- "CD8"
		ctype[grep("CD8",ctype,invert = T)] <- "CD4"
	}
	# ctype[!(ctype %in% c("CD4","CD8"))] <- NA
	return(ctype)
}

get_cellid_seus_rds <- function(stype="cart") {
	if (stype=="nr") {
		rds_fpath = "~/projects/20201123.cart/m01_3seqruns/analyze_NR/01a_get_NR/01a_seus_nr.rds"
	} else {
		rds_fpath = sprintf("~/projects/20201123.cart/m01_3seqruns/a01_combine_samples_all/a01_%s_seus.rds",stype)
	}
	message(rds_fpath)
	rds_fpath
}

get_20210202_seus_rds <- function(stype="cart") {
	sprintf('~/projects/20201123.cart/m01_3seqruns/samples_20210202/%s/seus_20210202.rds',stype)
}

get_meta_rds <- function(stype="cart") {
	if (stype=="nr") {
		rds_fpath = "~/projects/20201123.cart/m01_3seqruns/analyze_NR/01a_get_NR/01a_meta_nr.rds"
	} else {
		rds_fpath = sprintf("~/projects/20201123.cart/m01_3seqruns/a01_combine_samples_all/a01_%s_meta.rds",stype)
	}
	message(rds_fpath)
	rds_fpath
}

get_T_seus_rd <- function(wo_c3=1) {
	if (wo_c3==1) {
 		rds_fpath = "/home/hongc2/projects/20201123.cart/m01_3seqruns/filter_non_CART/a02_extractT/a02_T_cells_wo_c3.rds"
	} else {
		stop("not yet impleted!")
	}
	rds_fpath
}

get_genes_given_chr <- function(gtf_file="",gene_field="symbol",chrs=c("X","Y")) {
	if (gtf_file==""){
		gtf_file <- "~/hwangt-share/Datasets/Reference_seqs/GRCh38_P01/10x_genomics/GRCh38/genes/genes.gtf"
	}
	stopifnot(file.exists(gtf_file))
	gtf_gene_dt <- read_gtf_in_dt(gtf_file)
	return(sort(unique(gtf_gene_dt[chr %in% chrs,get(gene_field)])))
}

readRDS_w_msg <- function(rds_fpath) {
	message(sprintf("reading %s ...",rds_fpath))
	readRDS(file=rds_fpath)
}

saveRDS_w_msg <- function(var2,rds_fpath) {
	message(sprintf("writing %s ...",rds_fpath))
	saveRDS(var2,file=rds_fpath)
}

check_cols <- function(req_fields,dt2) {
	stopifnot(all(req_fields %in% colnames(dt2)))
}

dimplot_by_meta <- function(seui,pdf_outd,debug2=0) {
	if (debug2==1) {browser()}
	reds <- c("tsne")
	
	DefaultAssay(seui) <- "Protein"
	# =====================
	# protein expression on GEX umap
	lapply(reds,function(red2) {
		pdf_file <- file.path(pdf_outd,sprintf("adt_on_gex_cluster_by_tpoint_%s.pdf",red2))
		if (debug2==0){pdf(file=pdf_file,width=12,height = 4)}
		
		ret <- lapply(rownames(seui),function(feat) {
			message(feat)
			#Fig2.F
			p_list <- FeaturePlotSingle(seui, feature= `feat`, metadata_column = "tpoint", pt.size = 0.05, order =TRUE, assay2="Protein", red2=red2)
			
			p <- wrap_plots(p_list ,cols=3,,guides = 'collect')
			p <- p + plot_annotation(title = gsub("-TotalSeqB","",feat))
			if (debug2==1){browser()}
			plot(p)
		})
		if (debug2==0){dev.off()}
	})

	# =====================
	# protein avg expression on GEX cluster
	seu_w_adt <- subset(seui,seqrun %in% c(2))
	DefaultAssay(seu_w_adt) <- "Protein"
	seu_w_adt <- NormalizeData(seu_w_adt,verbose = TRUE)
	adt_avg <- AverageExpression(seu_w_adt,assay="Protein",slot = "data", return.seurat = TRUE)
	
	avgExprAdt <- GetAssayData(adt_avg,assay="Protein",slot="data")
	
	avgExprAdt <- avgExprAdt[!(rownames(avgExprAdt) %in% c("adt-CD4","adt-CD8")),]
	
	
	rownames(avgExprAdt) <- gsub("-TotalSeqB","",rownames(avgExprAdt))
	avgExprAdt2 <- avgExprAdt[,colnames(avgExprAdt) %in% c(0:12)]
	mtitle <- sprintf("ADT expression on GEX-based clusters\n[set mean to 0 across clusters]")
	
	pdf_file <- file.path(pdf_outd,sprintf("avgAdt_by_gex_cluster.pdf"))
	pdf(file=pdf_file)
	heatmap_feat_by_cluster_avgexprscaled(avgExprAdt2,title2=mtitle)
	dev.off()
	
	# =====================
	# low dim w/ patient meta
	pdf_files <- lapply(reds,function(red2) {
		
		pdf_file <- file.path(pdf_outd,sprintf("%s_with_cd4cd8.pdf",red2))
		
		pdf(pdf_file,width = 8,height = 4)
		p <- DimPlot(seui,
								 reduction = red2, 
								 pt.size = 0.07, 
								 group.by = "tpoint",
								 split.by = "cd4_cd8")
		
		p$layers[[1]]$aes_params$alpha = .5
		plot(p)
		dev.off()
		
		pdf_file <- file.path(pdf_outd,sprintf("%s_pmid30726743_tpoint.pdf",red2))
		pdf(pdf_file,width = 12,height = 4)
		lapply(c("CD4","CD8"),function(ctype) {
			p <- DimPlot(subset(seui,cd4_cd8==`ctype`),
							reduction = "tsne", 
							pt.size = 0.1, 
							group.by = "pmid_30726743",
							raster = FALSE,
							split.by = "tpoint")
			
			p$layers[[1]]$aes_params$alpha = .5
			p <- p + ggtitle(ctype)
			plot(p)
		})
		dev.off()
		
		# --------
		pdf_file <- file.path(pdf_outd,sprintf("%s_with_meta.pdf",red2))
		
		if (debug2==0) {pdf(pdf_file,width = 10,height = 8)}
		
		tplot <- DimPlot(seui, reduction = red2, pt.size = 0.07, group.by = "tpoint")
		tplot[[1]]$layers[[1]]$aes_params$alpha = .5
		plist <- list()
		plist[[1]] <- tplot + ggtitle("timepoint")
		
		tplot <- DimPlot(seui, reduction = red2, pt.size = 0.07, group.by = "resp")
		tplot[[1]]$layers[[1]]$aes_params$alpha = .5
		plist[[2]] <- tplot + ggtitle("treatment response")
		
		tplot <- DimPlot(seui, reduction = red2, pt.size = 0.07, group.by = "Phase")
		tplot[[1]]$layers[[1]]$aes_params$alpha = .5
		plist[[3]] <- tplot + ggtitle("cell cycle phase")

		tplot <- DimPlot(seui, reduction = red2, pt.size = 0.07, group.by = "patid")
		tplot[[1]]$layers[[1]]$aes_params$alpha = .5
		plist[[4]] <- tplot + ggtitle("patient id")
		
		p <- wrap_plots(plist,ncol=2)
		plot(p)
		dev.off()
		
		message(pdf_file)
		pdf_file
	})
	pdf_files
}

dimplot_by_cluster <- function(seui,pdf_outd,cluster_cols,debug2=0) {
	if (debug2==1) {browser()}
	reds <- c("tsne")
	
	# =====================
	# low dim w/ patient meta
	pdf_files <- lapply(reds,function(red2) {
		# --------
		pdf_file <- file.path(pdf_outd,sprintf("%s_with_meta_by_cluster.pdf",red2))
		
		if (debug2==0) {pdf(pdf_file,width = 10,height = 8)}

		plist <- lapply(cluster_cols,function(cluster_col){
			tplot <- DimPlot(seui, reduction = red2, pt.size = 0.07, group.by = cluster_col, label=T)
			tplot$layers[[1]]$aes_params$alpha = .5
			tplot + ggtitle(cluster_col)
		})
		p <- wrap_plots(plist,ncol=2,nrow=2)
		plot(p)
		dev.off()
		
		message(pdf_file)
		pdf_file
	})
	pdf_files
}

select_top_deg_feats <- function(deg_dt,out_tsv=NA,n1=10,n2=3,excl_pat="^HLA-|^IG[HJKL]|^RNA|^MT|^RP[SL]",incl=c("TIGIT"),ctype="both",debug2=0) {
	
	D <- length(unique(deg_dt$cluster))
	
	D1 <- round(D/2)
	D2 <- D1 + 1
	
	if (debug2==1){browser()}
	
	top1 <- deg_dt %>%
		dplyr::filter(cluster %in% c(0:D1)) %>%
		dplyr::group_by(cluster) %>% top_n(n = n1, wt = sort_by)
	
	top2 <- deg_dt %>% 
		dplyr::filter(cluster %in% D2:(D-1)) %>%
		dplyr::group_by(cluster) %>% top_n(n = n2, wt = sort_by)
	
	if (!is.na(out_tsv)) {
		message(sprintf("save the topK deg table into tsv file format [%s] ...",out_tsv))
		fwrite(bind_rows(top1,top2),file=out_tsv,sep="\t")
	}
	
	genes = union(top1$gene,top2$gene)
	genes = union(genes,incl)
	
	genes = genes[grep(excl_pat,genes,invert = T)]
	message(sprintf("total num of genes [%d]",length(genes)))
	message(sprintf("%s",paste0(genes,collapse=",")))
	genes[!is.na(genes)]
}

heatmap_annotated_clusters <- function(plot_title,cls_marker_deg,smeta,pdf_pref,majcn=10,mincn=3,apvalco=0.05,zscored=1,incl_genes=c("CD4","CD40LG","TNFRSF4","CD8A","CD8B","TIGIT"),ctypes=c("CD4","CD8"),row_split_cnt=3,col_split_cnt=3,width1=6,height1=10,abs_avg_logFC=0,top_genes=NA,debug2=0) {

	cls_marker_deg <- as.data.table(cls_marker_deg)
	colnames(cls_marker_deg) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")
	
	if (debug2==1){browser()}
	
	if (nrow(cls_marker_deg)>1) {
		message(sprintf("print umap along with topk DGE heatmap[%s]",plot_title))
		
		logfc_mat <- acast(cls_marker_deg, 
											 gene ~ cluster, 
											 value.var='avg_logFC')
		
		logfc_mat[is.na(logfc_mat)] <- 0
		
		if (dim(logfc_mat)[1]==1 | dim(logfc_mat)[2]==1) {
			message(sprintf("either 1 feature or 1 cluster is observed and ignore heatmap plotting!"))
		} else {
			cls_marker_deg = cls_marker_deg[gene %in% rownames(logfc_mat),]
			if (abs_avg_logFC==1) {
				cls_marker_deg[,sort_by:=abs(avg_logFC)]
			} else {
				cls_marker_deg[,sort_by:=avg_logFC]
			}
			
			if (dim(logfc_mat)[1]>100) {
				if (is.na(top_genes)){
					top_genes = select_top_deg_feats(cls_marker_deg,
																					 incl=incl_genes,
																					 n1=majcn,
																					 n2=mincn,
																					 debug2=debug2)
					cls_marker_deg$sort_by <- NULL
					logfc_mat <- logfc_mat[top_genes[top_genes %in% rownames(logfc_mat)],]
				} else {
					logfc_mat <- logfc_mat[top_genes,]
				}
			}
			
			if (debug2==1){browser()}
			
			ret1  = Heatmap_cluster_prof(smeta,
																	 cluster_ids = colnames(logfc_mat), 
																	 ctype = ctypes,
																	 debug2=0)
			
			ret = Heatmap_annot_gen(ret1,ctype=ctypes,debug2=0)
			
			hetmap_legend_tt <- "logFC"
			if (zscored==1) {
				logfc_mat <- get_zscore(logfc_mat)
				hetmap_legend_tt = "zscored_logFC"
			}
			
			myPal <- get_red_white_blue_for_Heatmap(logfc_mat)
			logfc_mat[is.na(logfc_mat)] <- 0.

			htm <- Heatmap(logfc_mat,
										 row_km = row_split_cnt,
										 column_km = col_split_cnt,
										 border = TRUE,
										 name=hetmap_legend_tt,
										 col=myPal,
										 #row_names_side="left",
										 # right_annotation = ha1,
										 top_annotation = ret$annot,
										 show_row_names = TRUE,
										 cluster_rows = TRUE,
										 column_title=plot_title)
			
			
			if (debug2==1){browser()}
			else {
				pdf_file <- sprintf("%s_z%d_deg_heatmap.pdf",pdf_pref,zscored)
				message(pdf_file)
				pdf(file = pdf_file,width = width1, height = height1)
			}
			
			ComplexHeatmap::draw(htm,annotation_legend_side="right",heatmap_legend_list = ret$lgd)
			
			if (debug2==0) {
				dev.off()
				graphics.off()
			}
			
		}
	} else {
		message(sprintf("No DGE observed in [%s]",plot_title))
	}
}

get_harmony_seui_smeta <- function(ctype="CART") {
	csize=readRDS('~/projects/20201123.cart/m01_3seqruns/cart/03a_cluster_size_comp/03a_cluster_size_pcts.rds')
	
	csize$smeta[is.na(pmid_30726743),pmid_30726743:="unk"]
	csize
}

Heatmap_cluster_prof <- function(smeta,cluster_ids,ctypes=c("CD4","CD8"), debug2=0) {
	if (debug2==1){browser()}
	# -------------
	#generate cell counts of each cluster 
	cellcnt.dt <- smeta[,.N,by=c("seurat_clusters")][order(seurat_clusters),]
	cellcnt.dt[,log10c:=log10(N+1)]
	
	cellcnt.dt <- data.frame(y=cellcnt.dt$log10c,
													 row.names = cellcnt.dt$seurat_clusters)
	
	# -------------
	# generate CD4 and CD8 fraction from the norm pct of each sample
	smry <- smeta[,.(cnt=.N),by=c("orig.ident","seurat_clusters","patid","tpoint","pmid_30726743","resp")]
	smry$tgroup <- "post_inf"
	smry[tpoint=="D0",tgroup:="pre_inf"]
	smry[,tgroup:=factor(tgroup,level=c("pre_inf","post_inf"))]
	smry[,resp:=factor(resp,level=c("R","NR"))]
	smry$cd4_cd8 <- "CD4"
	smry[grep("CD8",pmid_30726743),cd4_cd8:="CD8"]
	
	message(sprintf("total # of cells [%d]",sum(smry$cnt)))

	#total number of cells for each sample
	tcell <-smry[,.(denom=sum(cnt)),by=c("orig.ident")]
	smry$pct <- 100. * (smry$cnt / tcell[match(smry$orig.ident,orig.ident),denom])
	
	cgroup_frac <- acast(smry,
											 seurat_clusters ~ cd4_cd8,
											 fun.aggregate = sum,
											 value.var="pct")
	
	cgroupc = sweep(cgroup_frac,1,rowSums(cgroup_frac),"/")
	
	# ---------------
	if (debug2==1){browser()}
	#generate CD4 subtype composition stack (100%)
	cdj_profs = lapply(ctypes,function(ctype) {
		#ctype = "CD8"

		#filter by ctype and sum the subtype abundance
		abund.prof = smry[cd4_cd8==ctype & !is.na(pmid_30726743),.(abund=sum(pct)),by=c("seurat_clusters","pmid_30726743")]
		
		cdxprof <- acast(abund.prof,
												seurat_clusters ~ pmid_30726743, 
												value.var = "abund")
		
		cdxprof[is.na(cdxprof)] <- 0
		
		cdxprof <- sweep(cdxprof,1,rowSums(cdxprof),"/")
		
		#prep template
		cdxprof.sync = matrix(data=0,nrow=nrow(cgroupc),ncol=ncol(cdxprof))
		rownames(cdxprof.sync) <- rownames(cgroupc)
		colnames(cdxprof.sync) <- colnames(cdxprof)
		
		l <- !is.na(match(rownames(cdxprof.sync),rownames(cdxprof)))
		cdxprof.sync[l,] <- cdxprof
		cdxprof.sync[cluster_ids,]
	})
	names(cdj_profs) <- ctypes
	
	ret = list(cnt=cellcnt.dt[cluster_ids,],
			 group_cnt=cgroupc[cluster_ids,],
			 profs=cdj_profs)
	ret
}

Heatmap_annot_gen <- function(ret,ctypes=c("CD4","CD8"),debug2=0) {
	if (debug2==1) {browser()}
	
	if (all(ctypes=="CD4")) {
		# <=======
		
		cluster_ids <- rownames(ret$profs$CD4)
		annot_text_ci = anno_text(cluster_ids,gp=gpar(fontsize=rep(12,length(cluster_ids))))

		N1=dim(ret$profs$CD4)[1]
		cgroup_col = distinctColorPalette(N1)
		
		N2=dim(ret$profs$CD4)[2]
		cd4_col = distinctColorPalette(N2)

		cluster.prof = HeatmapAnnotation(log10_cell_count = anno_barplot(ret$cnt),

																	 cluster_id = annot_text_ci,
																	 
																	 CD4.T.subtypes = anno_barplot(ret$profs[[1]], 
																	 															gp = gpar(fill = cd4_col),
																	 															bar_width = 1, 
																	 															height = unit(2, "cm"),
																	 															show_row_names=FALSE),
																	 show_annotation_name = TRUE)
	
		lgd_list <- list(Legend(labels = colnames(ret$profs$CD4), 
										 			 title = "CD4.T.subtypes",
										 			 legend_gp = gpar(fill = cd4_col)))
	
	} else if (all(ctypes=="CD8")) {
		# <=======
		cluster_ids <- rownames(ret$profs$CD8)
		annot_text_ci = anno_text(cluster_ids,gp=gpar(fontsize=rep(12,length(cluster_ids))))
		
		N1=dim(ret$profs$CD8)[1]
		cgroup_col = distinctColorPalette(N1)

		N3=dim(ret$profs$CD8)[2]
		cd8_col = distinctColorPalette(N3)
		cluster.prof = HeatmapAnnotation(log10_cell_count = anno_barplot(ret$cnt),

																		 cluster_id = annot_text_ci,
																		 
																		 CD8.T.subtypes = anno_barplot(ret$profs$CD8,
																		 															gp = gpar(fill = cd8_col),
																		 															bar_width = 1,
																		 															height = unit(2, "cm"),
																		 															show_row_names=FALSE),
																		 show_annotation_name = TRUE)
		
		lgd_list <- list(Legend(labels = colnames(ret$profs$CD8), 
														title = "CD8.T.subtypes",
														legend_gp = gpar(fill = cd8_col)))
	} else {
		# <=======
		cluster_ids <- rownames(ret$profs$CD8)
		annot_text_ci = anno_text(cluster_ids,gp=gpar(fontsize=rep(12,length(cluster_ids))))
		
		N1=dim(ret$profs$CD8)[1]
		cgroup_col = distinctColorPalette(N1)
		
		N2=dim(ret$profs$CD4)[2]
		cd4_col = distinctColorPalette(N2)
		
		N3=dim(ret$profs$CD8)[2]
		cd8_col = distinctColorPalette(N3)
		cluster.prof = HeatmapAnnotation(log10_cell_count = anno_barplot(ret$cnt),
																		 
																		 CD4_CD8_frac = anno_barplot(ret$group_cnt, 
																		 														gp = gpar(fill = cgroup_col),
																		 														bar_width = 1, 
																		 														height = unit(1, "cm"),
																		 														show_row_names=FALSE),
																		 
																		 cluster_id = annot_text_ci,
																		 
																		 CD4.T.subtypes = anno_barplot(ret$profs$CD4, 
																		 															gp = gpar(fill = cd4_col),
																		 															bar_width = 1, 
																		 															height = unit(2, "cm"),
																		 															show_row_names=FALSE),
																		 CD8.T.subtypes = anno_barplot(ret$profs$CD8,
																		 															gp = gpar(fill = cd8_col),
																		 															bar_width = 1,
																		 															height = unit(2, "cm"),
																		 															show_row_names=FALSE),
																		 show_annotation_name = TRUE)
		
		lgd_list <- list(Legend(labels = colnames(ret$group_cnt), 
														title = "CD4_CD8_frac",
														legend_gp = gpar(fill = cgroup_col)),
										 Legend(labels = colnames(ret$profs$CD4), 
										 			 title = "CD4.T.subtypes",
										 			 legend_gp = gpar(fill = cd4_col)),
										 Legend(labels = colnames(ret$profs$CD8), 
										 			 title = "CD8.T.subtypes",
										 			 legend_gp = gpar(fill = cd8_col)))
	}
	list(annot=cluster.prof,lgd=lgd_list)
}

cluster_with_heatmap <- function(plot_title,cls_marker_deg,seui,cluster_col,pdf_pref,apvalco=0.05,zscored=1,incl_genes=c("CD4","CD40LG","TNFRSF4","CD8A","CD8B","TIGIT"),debug2=0) {
	
	Idents(seui) <- cluster_col
	
	pum <- DimPlot(seui, reduction = "tsne", label = TRUE, combine=T,label.size = 6,raster=FALSE)
	pum[[1]]$layers[[1]]$aes_params$alpha = .5
	
	pum <- pum + ggtitle(plot_title)
	
	pdf_file <- sprintf("%s_cluster.pdf",pdf_pref)
	if (debug2==1){browser()}
	else {
		pdf(pdf_file,width = 10, height = 8)
	}
	plot(pum)
	
	marker_deg_dt <- as.data.table(cls_marker_deg)
	colnames(marker_deg_dt) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")
	
	if (nrow(marker_deg_dt)>1) {
		message(sprintf("print umap along with topk DGE heatmap[%s]",plot_title))
		
		logfc_mat <- acast(marker_deg_dt, 
											 gene ~ cluster, 
											 value.var='avg_logFC')
		
		logfc_mat[is.na(logfc_mat)] <- 0
		
		if (dim(logfc_mat)[1]==1 | dim(logfc_mat)[2]==1) {
			message(sprintf("either 1 feature or 1 cluster is observed and ignore heatmap plotting!"))
		} else {
			a <- as.matrix(table(Idents(seui)))
			cluster_names <- rownames(a)
			cluster_info <- data.frame(frac=a/sum(a),
																 row.names=cluster_names)
			
			cls_marker_deg[,sort_by:=abs(avg_logFC)]
			
			top_genes = select_top_deg_feats(cls_marker_deg,
																			 incl=incl_genes,
																			 n1=10,
																			 n2=3,
																			 debug2=0)
			cls_marker_deg$sort_by <- NULL
			
			logfc_mat <- logfc_mat[top_genes,]
			
			if (debug2==1){browser()}
			
			hetmap_legend_tt <- "logFC"
			if (zscored==1) {
				logfc_mat <- get_zscore(logfc_mat)
				hetmap_legend_tt = "zscored_logFC"
			}
			
			myPal <- get_red_white_blue_for_Heatmap(logfc_mat)
			logfc_mat[is.na(logfc_mat)] <- 0.
			
			cluster_info <- cluster_info[colnames(logfc_mat),]
			
			ha1 <- rowAnnotation(frac=anno_barplot(cluster_info),
													 show_annotation_name = TRUE)
			
			phm <- as.ggplot(Heatmap(t(logfc_mat),
															 name=hetmap_legend_tt,
															 col=myPal,
															 #row_names_side="left",
															 right_annotation = ha1,
															 show_row_names = TRUE,
															 cluster_rows = TRUE,
															 column_title=sprintf("avg_logFC;padj<%g;zscored[%d]",apvalco,zscored)))
			
			pdf_file <- sprintf("%s_z%d_dge_heatmap.pdf",pdf_pref,zscored)
			if (debug2==1){browser()}
			else {
				pdf(pdf_file,width = 24, height = 8)
			}
			plot(phm)
			message(pdf_file)
			if (debug2==0){dev.off()}
			
		}
	} else {
		message(sprintf("No DGE observed in [%s]",plot_title))
	}
	
	if (debug2==0){dev.off()}
	
}

load_20210202_seus <- function(stype="cart",
															 reuse=1,
															 debug2=0) {
	if (debug2==1){browser()}
	suffix = "20210202"
	wkd<-file.path(get_projd0(),sprintf("samples_%s",suffix),stype)
	if (!file.exists(wkd)){dir.create(wkd,showWarnings = F,recursive = T)}
	rds_fpath <- file.path(wkd,sprintf("seus_%s.rds",suffix))
	if (reuse==1 & file.exists(rds_fpath)) {
		seus <- readRDS_w_msg(rds_fpath)
	} else {
		message("loading sample metadata ...")
		seus_rds=get_cellid_seus_rds(stype = stype)
		sample_meta=get_meta_rds(stype = stype)
		sample_meta=readRDS_w_msg(sample_meta)
		sample_meta = sample_meta[batch %in% c(1,2) | sample == "P3D14.c3",]
		seus <- readRDS_w_msg(seus_rds)
		seus <- seus[sample_meta$sample]
		
		seus <- lapply(seus,function(seu) {
			if (stype=="cart") {
				seu$cd4_cd8 <- "CD4"
				seu$cd4_cd8[grep("CD8",seu$pmid_30726743)] <- "CD8"
			} else {
				seu$cd_cd8 <- "unk"
			}
			seu
		})
		saveRDS(seus,file=rds_fpath)
	}
	seus
}

select_samples_from_sdt <- function(sample_meta,suffix="20210202") {
	if (suffix=="20210202") {
		sample_meta = sample_meta[(batch %in% c(1,2)) | (sample == "P3D14.c3"),]
	}
	sample_meta
}

#checked 05/14/2021
harmony_findmarkers <- function(seus,
																		stype="cart",
																		stepd="harmony_findmarkers",
																		reuse=1) {
	
	wkd<-file.path(get_projd0(),stype,stepd)
	if (!file.exists(wkd)){dir.create(wkd,showWarnings = F,recursive = T)}
	# -------------------------
	message(sprintf("aligning samples ..."))
	
	rds_fpath <- file.path(wkd,sprintf("seui.rds"))
	if (reuse==1 & file.exists(rds_fpath)) {
		seui = readRDS_w_msg(rds_fpath)
	} else {
		seui = aligned_by_harmony_w_cc(seus,
																	 D=50,
																	 harmony=1,
																	 do_cluster=1,
																	 min.cell=30,
																	 debug2=0)
		
		saveRDS_w_msg(seui,rds_fpath)
	}
	rm(seus)
	
	# -------------------------
	cls_marker_rds <- file.path(wkd,sprintf("cls_deg.rds")) # <==== out_file
	if (reuse==1 & file.exists(cls_marker_rds)) {
		cls_marker_deg <- readRDS_w_msg(cls_marker_rds)
	} else {
		markder_dt <- batch_solo_FindAllMarkers(list(harmony=seui),
																						ident2="seurat_clusters",
																						assay2="RNA",
																						debug2=0)
		cls_marker_deg <- markder_dt[[1]]
		saveRDS_w_msg(cls_marker_deg,cls_marker_rds)
		rm(seui)
	}
	
	cls_marker_tsv <- file.path(wkd,"cls_marker.tsv") # <==== out_file
	fwrite(cls_marker_deg,file=cls_marker_tsv,sep="\t")
}

locate_harmony_rds <- function(stype="cart",
															stepd="harmony_findmarkers") {
	wkd<-file.path(get_projd0(),stype,stepd)
	list(seu=file.path(wkd,sprintf("seui.rds")),
			 deg=file.path(wkd,sprintf("cls_deg.rds")))
}

vcluster_heatmap <- function(args) {
	if (args$debug==1){browser()}
	# ---------------------
	message("checking input/output ...")
	if (!file.exists(args$outd)) {
		dir.create(args$outd,showWarnings = F,recursive = T)
	}
	
	check_fpath(args$marker_rds)
	cls_marker_deg <- as.data.table(readRDS_w_msg(args$marker_rds))
	
	check_fpath(args$harmony_rds)
	seui <- readRDS_w_msg(args$harmony_rds)
	
	#Fig1.B, Fig2.E, Fig3.B
	dimplot_pdf_fns = dimplot_by_meta(sui,pdf_outd=args$outd)
	
	# -------------------
	message("drawing scaled heatmap w/ cluster markers by only avgLogFC ...")
	tsv_fpath <- get_out_fpath(args,sprintf("doheatmap_topk%d_%d_pos.tsv",8,4))
	
	cls_marker_deg[,sort_by:=abs(avg_logFC)]
	topk_genes = select_top_deg_feats(cls_marker_deg,tsv_fpath,n1=8,n2=4,debug2=args$debug)
	cls_marker_deg$sort_by <- NULL
	
	incl_genes = comma_string_to_list(args$req_genes)
	
	# -------------------
	message("drawing cluster/heatmap w/ significant cluster markers ...")
	pdf_pref2 <- get_out_fpath(args,"cls_markers") # <==== out_file
	
	#Fig1.C
	heatmap_annotated_clusters(args$heatmap_title,
											 cls_marker_deg[cls_marker_deg$p_val_adj<args$pval,],
											 smeta = as.data.table(seui@meta.data),
											 pdf_pref = pdf_pref2,
											 apvalco=args$pval,
											 incl_genes=incl_genes,
											 zscored=0,
											 row_split_cnt=5,
											 col_split_cnt=4,
											 width1=6,
											 height1=14,
											 debug2=args$debug)
	
	if (args$debug==1){browser()}
}

reduce_by_cluster_size = function(dt2_w_seurat_clusters,dominant_clusters=0:10) {
	dt2_w_seurat_clusters[seurat_clusters %in% dominant_clusters,]
}


split_seurat_by_tgroup <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$tgroup)
	seus <- lapply(items, function(item){
		subset(seu,tgroup==`item`)
	})
	names(seus) <- items
	seus
}

split_seurat_by_resp <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$resp)
	seus <- lapply(items, function(item){
		subset(seu,resp==`item`)
	})
	names(seus) <- items
	seus
}

add_tcgroup_meta <- function(seu,stype="cart") {
	seu$tgroup <- "post_inf"
	seu$tgroup[seu$tpoint=="D0"] <- "pre_inf"

	seu$cgroup <- ""
	for (ctype in c("CD4","CD8")) {
		if (stype=="pbmc") {
			j <- which(seu@meta.data$pmid_30726743 %in% c('T',sprintf('%s+ T',ctype)))
			
		} else {
			if (ctype=="CD4") {
				j<-grep("CD8",seu@meta.data$pmid_30726743,invert = TRUE)
			} else {
				j<-grep("CD8",seu@meta.data$pmid_30726743)
			}
		}
		seu$cgroup[j] <- ctype
	}
	seu
}

add_meta2 <- function(seu,cd4_cd8T=1) {
	seu$tgroup <- "post_inf"
	seu$tgroup[seu$tpoint=="D0"] <- "pre_inf"
	
	if (cd4_cd8T==1) {
		seu$cd4_cd8 <- "CD4"
		seu$cd4_cd8[grep("CD8",seu$pmid_30726743)] <- "CD8"
	} else {
		seu$cd4_cd8 = "pbmc"
	}
	seu
}


# ------------
split_seurat_by_cellgroup <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$cd4_cd8)
	seus <- lapply(items, function(item){
		subset(seu,cd4_cd8==`item`)
	})
	names(seus) <- items
	seus
}

# ------------
split_seurat_by_patient <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$patid)
	seus <- lapply(items, function(item){
		subset(seu,patid==`item`)
	})
	names(seus) <- items
	seus
}
# ------------
split_seurat_by_sample <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$orig.ident)
	seus <- lapply(items, function(item){
		subset(seu,orig.ident==`item`)
	})
	names(seus) <- items
	seus
}
# ------------
split_seurat_by_cluster <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$seurat_clusters)
	seus <- lapply(items, function(item){
		subset(seu,seurat_clusters==`item`)
	})
	names(seus) <- items
	seus
}

# ------------

split_seurat_by_tgroup <- function(seu,debug2=0) {
	if (debug2==1){browser()}
	items <- unique(seu@meta.data$tgroup)
	seus <- lapply(items, function(item){
		subset(seu,tgroup==`item`)
	})
	names(seus) <- items
	seus
}

get_cart_markers_by_ctype_tpoint <- function() {
	readRDS('~/projects/20200327.cart/ms_tigit0/CART/02a_singleR/s01_compare_pbmc_tcells_T/a05b02_cart_markers_ctype_tpoint.rds')
}

get_non_responders <- function() {
	c("P3","P8","P20","P23")
}

get_cluster_size_rd <- function() {
	rds_fpath = "~/projects/20201123.cart/m01_3seqruns/cart/03a_cluster_size_comp/03a_cluster_size_pcts.rds"
	readRDS_w_msg(rds_fpath)
}


#SFig1.C
cluster_by_pmid_30726743_composition <- function(smeta,debug2=0) {
	
	if (debug2==1){browser()}
	pmeta <- readRDS(get_meta_rds())
	
	cellcnt_by_cluster <- smeta[,.N,by=c("orig.ident","seurat_clusters","pmid_30726743")]
	cellcnt_by_sample <- smeta[,.(denom=.N),by=c("orig.ident")]
	cellcnt_by_cluster = merge(cellcnt_by_cluster,cellcnt_by_sample,by="orig.ident")
	cellcnt_by_cluster[,pct:=100.*N/denom]
	cellcnt_by_cluster <- merge(cellcnt_by_cluster,pmeta,by.x="orig.ident",by.y="sample")
	
	m <- acast(cellcnt_by_cluster,
						 seurat_clusters ~ pmid_30726743, 
						 fun.aggregate = sum,
						 value.var = "pct")
	
	agg_pct <- as.data.table(melt(m))
	colnames(agg_pct) <- c("seurat_clusters","pmid_30726743","sum.pct")
	agg_pct[,seurat_clusters:=factor(seurat_clusters)]
	
	agg_pct[is.na(pmid_30726743),pmid_30726743:="unk"]
	
	agg_pct$cd4_cd8 <- "CD4"
	agg_pct$cd4_cd8[grep("CD8",agg_pct$pmid_30726743)] <- "CD8"
	
	agg_denom = agg_pct[,.(denom=sum(sum.pct)),by=c("cd4_cd8")]
	agg_pct = merge(agg_pct,agg_denom,by="cd4_cd8")
	agg_pct[,norm.pct:=100*(sum.pct/denom)]
	
	delim <- list(s1=c(0,11),e2=c(10,19))
	plist <- lapply(1:2,function(j) {
		p <- ggplot(agg_pct[seurat_clusters %in% c(delim$s1[j]:delim$e2[j]),], 
								aes(fill=pmid_30726743, y=norm.pct, x=seurat_clusters)) + 
			geom_bar(position="stack", stat="identity") +
			facet_grid(rows = vars(cd4_cd8)) +
			labs(title = "Abundance of the T cell subtype across the samples",
					 y="Normalized cell fraction [pct]")
		p
	})
	list(plt=wrap_plots(plist,nrow=2),
			 cellid_comp=agg_pct)
}

deg_heatmap_w_clusters <- function(deg_dt,
																	 cluster_id="pmid_30726743",
																	 apval.co=0.01,
																	 logFc.co=0.5,
																	 topk1=25,
																	 topk2=15,
																	 width1=4,
																	 height1=12,
																	 pdf_file="out.pdf",
																	 mtitle="heatmap",
																	 debug2=0) {
	
	if (debug2==1){browser()}
	deg_dt <- deg_dt[grep('^HLA-|^IG[HJKL]|^RNA|^MT|^RP',gene,invert = T),]
	
	min_adj_pval <- deg_dt[p_val_adj>0,min(p_val_adj)]
	deg_dt[p_val_adj==0,p_val_adj:=min_adj_pval*0.9]
	# deg_dt2 <- deg_dt[p_val_adj<apval.co,]
	
	high_genes<-unique(unlist(lapply(split(deg_dt,by="cluster_id"),function(degci){
		degci[order(avg_logFC,decreasing = T),gene][1:topk1]
	})))
	
	low_genes<-unique(unlist(lapply(split(deg_dt,by="cluster_id"),function(degci){
		degci[order(avg_logFC,decreasing = F),gene][1:topk2]
	})))
	
	m <- acast(deg_dt[p_val_adj<apval.co & abs(avg_logFC)>=logFc.co,],
						 gene ~ `cluster_id`,
						 value.var = "avg_logFC")
	
	m[is.na(m)] <- 0.
	
	findex2disp <- which(rownames(m) %in% c(low_genes,high_genes))
	ha <- rowAnnotation(deg=anno_mark(at=findex2disp,
																		labels=rownames(m)[findex2disp],
																		labels_gp = gpar(fontsize = 8),
																		extend = unit(2,'mm'),
																		padding = unit(1, "mm"))
	)
	# if (debug2==1) {browser()}
	
	if (debug2==0){
		pdf(file=pdf_file,width=width1,height=height1)
	}
	
	v<-as.vector(m)
	col_fun = colorRamp2(c(min(v), 0, max(v)), c("blue", "white", "red"))
	
	if (debug2==0){dev.off()}
	
	ht_list <- Heatmap(m,
										 column_title = sprintf("%s\n(FDR<%g & |logFC|>=%g)",mtitle,apval.co,logFc.co),
										 cluster_rows = TRUE, 
										 name="logFC",
										 col=col_fun,
										 right_annotation = ha,
										 show_row_names=FALSE,
										 heatmap_legend_param = list(direction = "horizontal"),
										 row_names_gp = gpar(fontsize = 4))
	
	
	ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom")
	# p <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,annotation_legend_side = "bottom"))
	# plot(p)
	if (debug2==0){dev.off()}
}

dysfunction_genes_list <- function(topK=100) {
	
	long2015 <- c('KIAA1324','FN1','ITGA9','SMPDL3A','ITPRIPL2','PRSS23','CLECL1','GRB10','KIF21A','CTLA4','ATP10D','NUDT16','PLD1','PRKCSH','LGMN','SASH1','ZNF704','PLEKHH2','CTHRC1','ID1','FAXC','EGR1','TNFSF11','PTGIS','KLRC1','KLRC2','CD38','RBM47','SPINK2','AKAP5','ENTPD1')
	
	sadefeldman2018.CD8T_exhaust <- c('LAG3','PDCD1','HAVCR2','TIGIT','CD38','ENTPD1')
	leun2020.CD8T_dysfunctional <- c('LAYN','ITGAE','PDCD1','CTLA4','HAVCR2','LAG3','TIGIT','CXCL13','CD38','ENTPD1','CDK1','HSPH1','CCNB1','HSPB1','MKI67','CDK4','GZMB','TOX','IFNG','MIR155HG','TNFRSF9','RB1')
	
	long2015.all <- as.data.table(openxlsx::read.xlsx('~/projects/20201123.cart/m01_3seqruns/cart_r1r2/05a_aucell/ref_dysfunctional_markers/Long_2015_Exh_Gene_Set_GD2-CD28_vs_CD19-CD28.xlsx',sheet=1))
	
	j <- which(long2015.all$`Fold-Change(GD2.28z.vs..CD19.28z).(Description)`=="GD2 28z up vs CD19 28z")
	
	a <- long2015.all[`Fold-Change(GD2.28z.vs..CD19.28z).(Description)`=="GD2 28z up vs CD19 28z",][order(`Fold-Change(GD2.28z.vs..CD19.28z)`,decreasing = T),][1:topK]
	
	b <- unlist(sapply(a$Gene.Symbol,function(gene_str) {
		tstrsplit(gene_str,"\\/\\/\\/")
	}))
	
	long2015.all <- unique(sapply(b,trimws))
	long2015.all <- long2015.all[long2015.all != "---"]
	names(long2015.all) <- NULL
	
	gene_list <- list(long2015.lite=long2015,
										long2015.full=long2015.all,
										sadefeldman2018.CD8T_exhaust=sadefeldman2018.CD8T_exhaust,
										leun2020.CD8T_dysfunctional=leun2020.CD8T_dysfunctional)
	gene_list
}


deg_wo_clusters2 <- function(deg_dt,
														 domain1,
														 domain2,
														 ctrl1,
														 expr2,
														 apval.co=0.01,
														 min_logfc2disp=1.5,
														 logFc.co=1.,
														 genes_incl=NA,
														 topk1=25,
														 topk2=15,
														 pdf_file="out.pdf",
														 mtitle="heatmap",
														 debug2=0) {
	
	deg_dt <- deg_dt[grep('^HLA-|^IG[HJKL]|^RNA|^MT|^RP',gene,invert = T),]
	setnames(deg_dt,'p_val_adj','adj.P.Val')
	if ('avg_logFC' %in% colnames(deg_dt)) {
		setnames(deg_dt,'avg_logFC','logFC')
	} else if ('avg_log2FC' %in% colnames(deg_dt)) {
		setnames(deg_dt,'avg_log2FC','logFC')
	}
	
	min_adj_pval <- deg_dt[adj.P.Val>0,min(adj.P.Val)]
	deg_dt[adj.P.Val==0,adj.P.Val:=min_adj_pval]
	
	if (F) {
		deg_dts <- split(deg_dt,by="cluster_id")
		plist <- list()
		plist[[1]] <- EnhancedVolcanoDESeq2(deg_dts[[domain1]],column2disp="gene",AdjustedCutoff=0.05, LabellingCutoff=0.001, LFCCutoff=1, main=sprintf("Differentially experssed genes between %s and %s in %s",ctrl1,expr2,domain2),debug2=0)
		
		plist[[2]] <- EnhancedVolcanoDESeq2(deg_dts[[domain2]],column2disp="gene",AdjustedCutoff=0.05, LabellingCutoff=0.001, LFCCutoff=1, main=sprintf("Differentially experssed genes between %s and %s in %s",ctrl1,expr2,domain2),debug2=0)
		
		p <- wrap_plots(plist,nrow=2)
		plot(p)
	}
	
	m <- acast(deg_dt[adj.P.Val<apval.co & abs(logFC)>=logFc.co,],
						 gene ~ cluster_id, 
						 value.var = "logFC")
	
	m[is.na(m)] <- 0.
	
	m.dt <- as.data.table(m,keep.rownames = T)
	m.dt$todisp <- 0
	
	if (debug2==1) {browser()}
	
	g0a <- NA
	if (!is.na(genes_incl)) {
		m.dt[rn %in% genes_incl & (abs(get(domain1))>min_logfc2disp | abs(get(domain2))>min_logfc2disp),todisp:=3]
		g0a <- m.dt[todisp==3,rn]
	}
	
	m.dt[abs(get(domain1))>min_logfc2disp & abs(get(domain2))>min_logfc2disp,todisp:=2]
	g1a=m.dt[todisp==2,][order(abs(get(domain1)),decreasing = T),rn][1:topk2]
	g2a=m.dt[todisp==2,][order(abs(get(domain2)),decreasing = T),rn][1:topk2]
	g1b=m.dt[todisp==0,][order(abs(get(domain1)),decreasing = T),rn][1:topk1]
	g2b=m.dt[todisp==0,][order(abs(get(domain2)),decreasing = T),rn][1:topk1]
	
	m.dt$todisp <- 0
	m.dt[rn %in% unique(c(g0a,g1a,g2a,g1b,g2b,"TIGIT")),todisp:=1]
	
	j <- which(m.dt$todisp==1)
	
	ha <- rowAnnotation(deg=anno_mark(at=j,
																		labels=rownames(m)[j],
																		labels_gp = gpar(fontsize = 8),
																		extend = unit(2,'mm'),
																		padding = unit(1, "mm"))
	)
	if (debug2==1) {browser()}
	
	pdf(file=pdf_file,width=4,height=12)
	
	ht_list <- Heatmap(m, 
										 column_title = sprintf("%s;logFC of %s to %s\n(FDR<%g & |logFC|>=%g)",mtitle,expr2,ctrl1,apval.co,logFc.co),
										 cluster_rows = TRUE, 
										 name="logFC",
										 right_annotation = ha,
										 show_row_names=FALSE,
										 heatmap_legend_param = list(direction = "horizontal"),
										 row_names_gp = gpar(fontsize = 4))
	
	ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom")
	# p <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,annotation_legend_side = "bottom"))
	# return(p)
	dev.off()
}


deg_wo_clusters1 <- function(deg_dtj,
														 ctrl1,
														 expr2,
														 apval.co=0.01,
														 min_logfc2disp=1.5,
														 logFc.co=1.,
														 genes_incl=NA,
														 topk1=25,
														 topk2=15,
														 pdf_file="out.pdf",
														 mtitle="heatmap",
														 debug2=0) {
	
	deg_dtj <- deg_dtj[grep('^HLA-|^IG[HJKL]|^RNA|^MT|^RP',gene,invert = T),]
	
	min_adj_pval <- deg_dtj[p_val_adj>0,min(p_val_adj)]
	deg_dtj[p_val_adj==0,p_val_adj:=min_adj_pval]
	
	m <- acast(deg_dtj[p_val_adj<apval.co & abs(avg_logFC)>=logFc.co,],
						 gene ~ cluster_id, 
						 value.var = "avg_logFC")
	
	m[is.na(m)] <- 0.
	
	if (debug2==1) {browser()
	} else {pdf(file=pdf_file,width=4,height=12)}
	
	ht_list <- Heatmap(m, 
										 column_title = sprintf("%s;logFC of %s to %s\n(FDR<%g & |logFC|>=%g)",mtitle,expr2,ctrl1,apval.co,logFc.co),
										 cluster_rows = TRUE, 
										 name="logFC",
										 show_row_names=TRUE,
										 row_names_gp = gpar(fontsize = 12))
	
	ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom")
	# p <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,annotation_legend_side = "bottom"))
	# return(p)
	if (debug2==0){dev.off()}
	
}

load_pseudobulk <- function(seu,uassay="RNA",uq3=1,debug2=0) { #harmony seurat obj
	
	by_cgroups <- split_seurat_by_cellgroup(seu)
	
	bulk.ms <- imap(by_cgroups,function(by_cgroup,cgroup){
		if (debug2==1) {browser()}
		counts.m <- as.matrix(GetAssayData(by_cgroup,assay=uassay,slot = "counts"))
		counts.m[is.na(counts.m)] = 0
		if (uq3==1) {
			uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
			counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
			counts.m[is.nan(counts.m)] = 0
		}
		t(Matrix.utils::aggregate.Matrix(t(counts.m),
																		 groupings = by_cgroup@meta.data$orig.ident, fun = "sum"))
	})
	bulk.ms
}

pseudobulk_with_replicate2 <- function(by_sample,sname,uassay="RNA",R=3,uq3=0,debug2=0) {
	if (debug2==1) {browser()}
	counts.m <- as.matrix(GetAssayData(by_sample,assay=uassay,slot = "counts"))
	counts.m[is.na(counts.m)] = 0
	cbc_on = colnames(counts.m)
	
	if (length(cbc_on) >= 2*R) {
		pseudo_repi <- (sample(1:R, length(cbc_on), replace=T))
	} else {
		pseudo_repi = rep(1:R,length(cbc_on))
		pseudo_repi = pseudo_repi[1:length(cbc_on)]
	}
	
	pbk = NA
	if (length(cbc_on)>0) {
		if (uq3==1) {
			uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
			counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
			counts.m[is.nan(counts.m)] = 0
		}
		
		pbk = t(Matrix.utils::aggregate.Matrix(t(counts.m),
																					 groupings = pseudo_repi, fun = "sum"))
		
		colnames(pbk) = sprintf("%s.%s",sname,colnames(pbk))
	}
	pbk
}


load_pseudobulk_TIGIT <- function(seu,uq3=1,reuse=1) { #harmony seurat obj
	rds_fn = "out/cart/harmony_findmarkers/pseudo_bulk_TIGIT.rds"
	if (reuse==1 & file.exists(rds_fn)){
		bulk.lst = readRDS(rds_fn)
	} else {
		by_cgroups <- split_seurat_by_cellgroup(seu)
		
		bulk.lst <- imap(by_cgroups,function(by_cgroup,cgroup) {
			
			tigit_expr <- GetAssayData(by_cgroup,assay="RNA",slot = "data")['TIGIT',]
			by_cgroup <- AddMetaData(by_cgroup,tigit_expr,col.name = "TIGIT")
			indic <- rep("OFF",length(tigit_expr))
			indic[tigit_expr>0] <- "ON"
			by_cgroup <- AddMetaData(by_cgroup,indic,col.name = "TIGIT_status")
			
			by_cgroup@meta.data$orig.ident_w_TIGIT_status <- paste0(by_cgroup@meta.data$orig.ident,'_',by_cgroup@meta.data$TIGIT_status)
			
			counts.m <- GetAssayData(by_cgroup,slot = "counts")
			
			if (uq3==1) {
				uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
				counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
			}
			bulkm=t(Matrix.utils::aggregate.Matrix(t(counts.m),
																			 groupings = by_cgroup@meta.data$orig.ident_w_TIGIT_status, fun = "sum"))
			list(mtx=bulkm,smeta=by_cgroup@meta.data)
		})
		saveRDS(bulk.lst,file=rds_fn)
	}
	bulk.lst
}

map_geneSymbol_ensemblGenes <- function(genes,convert_to="GENEID",debug2=0) {
	library(EnsDb.Hsapiens.v86)
	if (debug2==1){browser()}
	
	if (convert_to=="GENEID") {
		gene_map <- ensembldb::select(EnsDb.Hsapiens.v86, 
																	keys= genes, 
																	keytype = "SYMBOL", 
																	columns = c("SYMBOL","GENEID"))
	} else {
		gene_map <- ensembldb::select(EnsDb.Hsapiens.v86, 
																	keys= genes, 
																	keytype = "GENEID", 
																	columns = c("SYMBOL","GENEID"))
	}
	
	gene_map2 <- as.data.table(gene_map)
	gene_map2 <- gene_map2[!startsWith(GENEID,'LRG'),]
	dupcnt <- gene_map2[,.N,by="SYMBOL"]
	gene_map2[SYMBOL %in% dupcnt[N>1,SYMBOL],]
	message("[WARNING]pick up the first mapped GENEID for duplicated!!!")
	gene_map2 = gene_map2[!duplicated(SYMBOL),]
	gene_map2[!is.na(GENEID) & !is.na(SYMBOL),]
}
