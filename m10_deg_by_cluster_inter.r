#0. load shared library
source('lib/lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd=get_wkd(get_projd0(fpath_dt$parentd),fpath_dt$fpath),
									 cluster_col="seurat_clusters",
									 dge_method="MAST",
									 ncpu=2,
									 reuse=1,
									 debug=0)

deg_clusters_heatmap2 <- function(deg_dt,
																	ctype="CD8",
														 apval.co=0.05,
														 logFc.co=0.3,
														 pdf_file="out.pdf",
														 mtitle="post_inf.CD8",
														 title2="any",
														 debug2=0) {
	
	if (debug2==1) {browser()}
	
	M <- 40
	deg_dt <- deg_dt[grep('^HLA-|^IG[HJKL]|^RNA|^MT|^RP',gene,invert = T),]
	setnames(deg_dt,'p_val_adj','adj.P.Val')
	setnames(deg_dt,'avg_logFC','logFC')
	
	min_adj_pval <- deg_dt[adj.P.Val>0,min(adj.P.Val)]
	deg_dt[adj.P.Val==0,adj.P.Val:=min_adj_pval]

	m <- acast(deg_dt[adj.P.Val<apval.co & abs(logFC)>=logFc.co,],
						 gene ~ cluster_id, 
						 value.var = "logFC")
	
	m[is.na(m)] <- 0.
	
	if (debug2==1) {browser()}
	
	N <- dim(m)[1]
	
	m1 <- select_topk_from_each_sample(abs(m),
																		 num_feats=10)
	
	m2 <- m[rownames(m1),]
	
	m2.dt <- as.data.table(m2,keep.rownames = TRUE)
	
	#m2.dt$todisp <- FALSE
	#j <- which(m2.dt$todisp==TRUE)
	
	# -----------
	ret1  = Heatmap_cluster_prof(smeta,
															 cluster_ids=colnames(m2),
															 ctype=ctypes)
	
	ret = Heatmap_annot_gen(ret1,ctype=ctypes)

	if (debug2==1) {browser()}
	
	pdf(file=pdf_file,width=8,height=10)
	
	myPal <- get_red_white_blue_for_Heatmap(m2)
	
	ht_list <- Heatmap(m2,
										 row_km = 5,
										 column_km = 4,
										 border = TRUE,
										 column_title = sprintf("%s\nFDR<%g;|logFC(NR/R)|>%g",title2,apval.co,logFc.co),
										 cluster_rows = TRUE, 
										 name="logFC",
										 # right_annotation = ha,
										 top_annotation = ret$annot,
										 show_row_names=TRUE,
										 col = myPal)
	
	ComplexHeatmap::draw(ht_list, heatmap_legend_list=ret$lgd,heatmap_legend_side = "right")
	
	# p <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,annotation_legend_side = "bottom"))
	# return(p)
	dev.off()
}

# ------------
message("highlight immreg genes ...")

rds_file <- get_out_fpath(args,sprintf("cluster_deg_on_treatment.rds"))
if (args$reuse==1 & file.exists(rds_file)) {
	deg_dts <- readRDS(rds_file)
} else {
	harmony_rds=locate_harmony_rds()
	seu <- readRDS(harmony_rds$seu)
	seu <- add_meta2(seu)
	
	seus <- split_seurat_by_tgroup(seu,debug2=0)
	
	deg_dts <- lapply(names(seus),function(tgroup) {
		# tgroup <- "post_inf"
		seu <- seus[[tgroup]]
		seus2 <- split_seurat_by_cellgroup(seu,debug2=0)
		deg_dts2 <- lapply(names(seus2),function(cgroup) {
			# cgroup <- "CD4"
			message(sprintf("%s,%s",tgroup,cgroup))
			seu2 <- seus2[[cgroup]]
			seu2$comp <- "ctrl"
			seu2$comp[seu2$resp=="NR"] <- "expr"
			
			# browser()
			region_col <- "seurat_clusters"
			deg_dt <- diff_from_integration(seu2,
																			cluster_col=region_col,
																			min_cell_cnts=3,
																			method=args$dge_method,
																			ncpu=8,
																			debug=0)
			deg_dt
			
		})
		names(deg_dts2) <- names(seus2)
		deg_dts2
	})
	names(deg_dts) <- names(seus)
	saveRDS(deg_dts,file=rds_file)
}

lapply(names(deg_dts),function(tgroup) {
	lapply(names(deg_dts[[tgroup]]),function(cgroup) {
		if (cgroup!="unk") {
			mtitle2 <- paste0(tgroup,'.',cgroup)
			message(mtitle2)
			deg_dt2<-deg_dts[[tgroup]][[cgroup]]
			
			deg_clusters_heatmap2(deg_dt2,
														ctype=cgroup,
														pdf_file=get_out_fpath(args,sprintf("%s_heatmap.pdf",mtitle2)),
														logFc.co=0.5,
														mtitle=mtitle2,
														title2=sprintf("%s;%s",mtitle2,'any'),
														debug2=1)
		}
	})
})

harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)
seus <- split_seurat_by_tgroup(seu,debug2=0)
rm(seu)

lapply(names(deg_dts),function(tgroup) {
	seus2 <- split_seurat_by_cellgroup(seus[[tgroup]])
	lapply(names(deg_dts[[tgroup]]),function(cgroup) {
		seu2 <- seus2[[cgroup]]
		if (cgroup!="unk") {
			mtitle2 <- paste0(tgroup,'.',cgroup)
			
			message(mtitle2)
			# browser()
			deg_dt2<-deg_dts[[tgroup]][[cgroup]]

			#Fig3.D
			heatmap_annotated_clusters(sprintf("%s;%s",mtitle2,'immreg'),
																 deg_dt2[gene %in% get_immreg_genes(extended=1) & abs(avg_logFC)>=0.25,],
																 smeta = as.data.table(seu2@meta.data),
																 ctypes=list(cgroup),
																 pdf_pref = get_out_fpath(args,sprintf("%s_immreg_heatmap",mtitle2)),
																 majcn=10,
																 mincn=5,
																 apvalco=0.05,
																 incl_genes="TIGIT",
																 zscored=0,
																 row_split_cnt=3,
																 col_split_cnt=3,
																 width1=6,
																 height1=14,
																 debug2=0)

		}
	})
})
