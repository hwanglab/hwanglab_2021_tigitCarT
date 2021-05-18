#0. load shared library
source('lib/lib_project.r')

heatmap_annotated_clusters_cellid <- function(plot_title,cls_marker_deg,smeta,pdf_pref,majcn=10,mincn=3,apvalco=0.05,zscored=1,incl_genes=c("CD4","CD40LG","TNFRSF4","CD8A","CD8B","TIGIT"),ctypes=c("CD4","CD8"),row_split_cnt=3,col_split_cnt=3,width1=6,height1=10,debug2=0) {
	
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
			cls_marker_deg[,sort_by:=abs(avg_logFC)]
			
			if (dim(logfc_mat)[1]>100) {
				top_genes = select_top_deg_feats(cls_marker_deg,
																				 incl=incl_genes,
																				 n1=majcn,
																				 n2=mincn,
																				 debug2=debug2)
				cls_marker_deg$sort_by <- NULL
				logfc_mat <- logfc_mat[top_genes %in% rownames(logfc_mat),]
			}
			
			if (debug2==1){browser()}
			
			cellcnt.dt <- smeta[,.N,by=c("pmid_30726743")][order(pmid_30726743),]
			cellcnt.dt[,log10c:=log10(N+1)]
			
			cavail = colnames(logfc_mat)
			
			cellcnt.dt = cellcnt.dt[pmid_30726743 %in% cavail,]
			
			cellcnt.dt <- data.frame(y=cellcnt.dt$log10c,
															 row.names = cellcnt.dt$pmid_30726743)
			
			# -------------
			# generate CD4 and CD8 fraction from the norm pct of each sample
			smry <- smeta[,.(cnt=.N),by=c("orig.ident","seurat_clusters","patid","tpoint","pmid_30726743","resp","TIGIT_status")]
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
													 pmid_30726743 ~ TIGIT_status,
													 fun.aggregate = sum,
													 value.var="pct")
			
			cgroupc = sweep(cgroup_frac,1,rowSums(cgroup_frac),"/")
			cgroupc = cgroupc[cavail,]
			
			N3=dim(cgroupc)[2]
			TIGIT_col = distinctColorPalette(N3)
			cluster.prof = HeatmapAnnotation(log10_cell_count = anno_barplot(cellcnt.dt[rownames(cgroupc),]),
																			 TIGIT = anno_barplot(cgroupc,
																			 										 gp = gpar(fill = TIGIT_col),
																			 										 bar_width = 1,
																			 										 height = unit(2, "cm"),
																			 										 show_row_names=FALSE),
																			 show_annotation_name = TRUE)
			
			lgd_list <- list(Legend(labels = colnames(cgroupc), 
															title = "TIGIT_expr",
															legend_gp = gpar(fill = TIGIT_col)))
			
			
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
										 top_annotation = cluster.prof,
										 show_row_names = TRUE,
										 cluster_rows = TRUE,
										 column_title=plot_title)
			
			
			if (debug2==1){browser()}
			else {
				pdf_file <- sprintf("%s_z%d_deg_heatmap.pdf",pdf_pref,zscored)
				message(pdf_file)
				pdf(file = pdf_file,width = width1, height = height1)
			}
			
			ComplexHeatmap::draw(htm,annotation_legend_side="right",heatmap_legend_list = lgd_list)
			
			if (debug2==0) {
				dev.off()
				graphics.off()
			}
			
		}
	} else {
		message(sprintf("No DGE observed in [%s]",plot_title))
	}
}

#####################
fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd=get_wkd(get_projd0(fpath_dt$parentd),fpath_dt$fpath),
									 dge_method="MAST",
									 ncpu=2,
									 reuse=0,
									 debug=0)

message("loading seus ...")

harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)

query <- "TIGIT"

query_expr <- GetAssayData(seu,assay="RNA",slot = "data")[query,]
seu <- AddMetaData(seu,query_expr,col.name = `query`)

seu$comp <- "ctrl"
seu$comp[seu[[query]]>0] <- "expr"

region_col <- "cd4_cd8"

deg_dt <- diff_from_integration(seu,
																cluster_col=region_col,
																min_cell_cnts=3,
																method=args$dge_method,
																ncpu=2,
																debug=0)

deg_dt = deg_dt[p_val_adj<0.05,]
deg_dt$imm_reg_genes <- "N"
deg_dt[gene %in% get_immreg_genes(extended = 1),imm_reg_genes:="Y"]

tsv_fn <- file.path(args$outd,"deg_tigit_on_off_wo_cluster.tsv")
fwrite(deg_dt,file=tsv_fn,sep="\t")

# ===========
# pseudobulk
pbulk.lst = load_pseudobulk_TIGIT(seu)
pmeta <- readRDS(get_meta_rds())
pmeta1 <- pmeta
pmeta1$TIGIT_status <- "ON"
pmeta1[,sample_TIGIT_status:=paste0(sample,'_',TIGIT_status)]

pmeta2 <- pmeta
pmeta2$TIGIT_status <- "OFF"
pmeta2[,sample_TIGIT_status:=paste0(sample,'_',TIGIT_status)]
pmeta3 <- rbind(pmeta1,pmeta2)

sdt = pmeta3[match(colnames(pbulk.lst$CD8$mtx),sample_TIGIT_status),]
sdt$tgroup="post_inf"
sdt[tpoint=="D0",tgroup:="pre_inf"]

deg_immreg.dts = imap(split(sdt,by="tgroup"),function(by_tgroup,tgroup) {
	
	by_tgroup$comp <- "ctrl"
	by_tgroup[TIGIT_status=="ON",comp:="expr"]
	
	de2in <- list()
	de2in$coldata = as.data.frame(by_tgroup,row.names = by_tgroup$sample)
	de2in$cts <- as.matrix(pbulk.lst$CD8$mtx[,by_tgroup$sample_TIGIT_status])
	
	design_f_str = "~ patient + comp"
	deg.dt <- run_deseq2(de2in,args$outd,design_f_str,shrink_lfc=FALSE,alpha=0.1,exptag=sprintf("CD8_TIGIT_%s",tgroup),debug2=0)
	
	deg.dt = deg.dt[!is.na(padj),]
	
	setnames(deg.dt,'rn','gene')
	setnames(deg.dt,'log2FoldChange','logFC')
	setnames(deg.dt,'padj','adj.P.Val')
	
	deg_immreg.dt = deg.dt[gene %in% get_immreg_genes(extended=1) & gene != "TIGIT",]
	# browser()
	# Fig4.A
	p=EnhancedVolcanoDESeq2(deg_immreg.dt,
													column2disp="gene",
													AdjustedCutoff=0.05, 
													LabellingCutoff=0.05, 
													LFCCutoff=0.5, 
													main=sprintf("DEG LogFC(%s/%s);%s;%s in pseudo bulk mRNA","TIGIT_on","TIGIT_off","CD8",tgroup),
													debug2=0)
	
	pdf_fn <- file.path(args$outd,sprintf("deseq2_deg_by_tigit_%s.pdf",tgroup))
	pdf(file=pdf_fn,width=8,height=10)
	plot(p)
	dev.off()
	deg_immreg.dt
})

# ==================
# same analysis but per cell types
# region_cols <- c("seurat_clusters","pmid_30726743")
region_cols <- c("pmid_30726743")
by_cgroups <- split_seurat_by_cellgroup(seu)
csize=get_harmony_seui_smeta()

region_col <- region_cols[[1]]

# region_col <- "seurat_clusters"
by_cgroup = by_cgroups$CD8
cgroup = "CD8"
by_tgroups <- split_seurat_by_tgroup(by_cgroup,debug2=0)
tgroup = "post_inf"
by_tgroup = by_tgroups$post_inf

# tgroup="pre_inf"
message(sprintf("%s;%s;%s",region_col,cgroup,tgroup))
# by_tgroup = by_tgroups[[tgroup]]
by_tgroup = by_cgroup

tigit_expr <- GetAssayData(by_tgroup,assay="RNA",slot = "data")['TIGIT',]
by_tgroup <- AddMetaData(by_tgroup,tigit_expr,col.name = "TIGIT")
indic <- rep("OFF",length(tigit_expr))
indic[tigit_expr>0] <- "ON"
by_tgroup <- AddMetaData(by_tgroup,indic,col.name = "TIGIT_status")

by_tgroup$comp <- "ctrl"
by_tgroup$comp[by_tgroup$TIGIT_status=="ON"] <- "expr"

deg_rds <- file.path(args$outd,sprintf("%s_%s.rds",region_col,cgroup))
if (args$reuse==0 & file.exists(deg_rds)){
	deg_dt <- readRDS(deg_rds)
} else {
	deg_dt <- diff_from_integration(by_tgroup,
																	cluster_col=region_col,
																	min_cell_cnts=3,
																	method=args$dge_method,
																	ncpu=2,
																	debug=0)
	saveRDS(deg_dt,file=deg_rds)
}

deg_dtj <- deg_dt[gene %in% get_immreg_genes(extended = 1),]
deg_dtj <- deg_dtj[p_val_adj<0.05,]
deg_dtj <- deg_dtj[!(gene %in% 'TIGIT'),]
setnames(deg_dtj,'cluster_id','cluster')

mtitle=sprintf("%s;\nDEG(logFC=TIGIT_on/TIGIT_off)\n%s",cgroup,region_col)

pdf_pref2=file.path(args$outd,sprintf("%s_%s_deg_heatmap_immreg_x_%s",cgroup,tgroup,region_col))

heatmap_annotated_clusters_cellid(mtitle,
													 deg_dtj,
													 smeta = as.data.table(by_tgroup@meta.data),
													 pdf_pref = pdf_pref2,
													 apvalco=0.05,
													 incl_genes=NA,
													 zscored=0,
													 ctypes=cgroup,
													 row_split_cnt=3,
													 col_split_cnt=2,
													 width1=5,
													 height1=10,
													 debug2=0)