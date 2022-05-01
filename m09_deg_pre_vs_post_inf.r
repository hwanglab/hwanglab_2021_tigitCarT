#0. load shared library
source('lib/lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd="out/cart/deg_pre_vs_post_inf",
									 inrds1="out/cart/vcluster_heatmap/pmid_30726743_pct_by_seurat_clusters.rds",
									 cluster_col="seurat_clusters",
									 dge_method="MAST",
									 ncpu=2,
									 reuse=1,
									 debug=0)

deg_wo_clusters <- function(deg_dt,
														domain1,
														domain2,
														ctrl1,
														expr2,
														apval.co=0.01,
														logFc.co=1.,
														topk1=25,
														topk2=15,
														pdf_file="out.pdf",
														mtitle="heatmap",
														debug2=0) {

	deg_dt <- deg_dt[grep('^HLA-|^IG[HJKL]|^RNA|^MT|^RP',gene,invert = T),]
	setnames(deg_dt,'p_val_adj','adj.P.Val')
	setnames(deg_dt,'avg_logFC','logFC')
	
	min_adj_pval <- deg_dt[adj.P.Val>0,min(adj.P.Val)]
	deg_dt[adj.P.Val==0,adj.P.Val:=min_adj_pval*0.5]
	
	if (debug2==1) {browser()}

	m <- acast(deg_dt[adj.P.Val<apval.co & abs(logFC)>=logFc.co,],
						 gene ~ cluster_id, 
						 value.var = "logFC")
	
	m[is.na(m)] <- 0.
	
	v<-as.vector(m)
	col_fun = colorRamp2(c(min(v), 0, max(v)), c("blue", "white", "red"))
	
	min_logfc2disp <- 1
	
	m.dt <- as.data.table(m,keep.rownames = T)
	m.dt$todisp <- FALSE
	
	if (debug2==1) {browser()}
	
	m.dt[abs(get(domain1))>min_logfc2disp & abs(get(domain2))>min_logfc2disp,todisp:=TRUE]
	g1a=m.dt[todisp==TRUE,][order(abs(get(domain1)),decreasing = T),rn][1:topk2]
	g2a=m.dt[todisp==TRUE,][order(abs(get(domain2)),decreasing = T),rn][1:topk2]
	g1b=m.dt[todisp==FALSE,][order(abs(get(domain1)),decreasing = T),rn][1:topk1]
	g2b=m.dt[todisp==FALSE,][order(abs(get(domain2)),decreasing = T),rn][1:topk1]
	m.dt$todisp <- FALSE
	
	m.dt[rn %in% unique(c(g1a,g2a,g1b,g2b,"TIGIT")),todisp:=TRUE]
	
	j <- which(m.dt$todisp==TRUE)
	
	ha <- rowAnnotation(deg=anno_mark(at=j,
																		labels=rownames(m)[j],
																		labels_gp = gpar(fontsize = 8),
																		extend = unit(2,'mm'),
																		padding = unit(1, "mm"))
	)
	if (debug2==1) { browser() 
	} else {
		pdf(file=pdf_file,width=4,height=12)
	}
	
	ht_list <- Heatmap(m, 
										 column_title = sprintf("%s;logFC of %s to %s\n(FDR<%g & |logFC|>=%g)",mtitle,expr2,ctrl1,apval.co,logFc.co),
					cluster_rows = TRUE, 
					name="logFC",
					# right_annotation = ha,
					col=col_fun,
					show_row_names=TRUE
					# row_names_gp = gpar(fontsize = 4),
					# heatmap_legend_param = list(direction = "horizontal")
					)
	
	ComplexHeatmap::draw(ht_list, heatmap_legend_side = "right")
	# p <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,annotation_legend_side = "bottom"))
	# return(p)
	if (debug2==0){dev.off()}
}


# ------------
harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)
seu$whole <- "unite"

# split by cell types
by_cgroups <- split_seurat_by_cellgroup(seu)

region_col <- "whole"

cellid_frac.dt <- readRDS(args$inrds1)

# pdf_fn <- file.path(args$outd,"deg_tgroup_per_cluster.pdf")
# pdf(file=pdf_fn,width=6,height=12)
# ====================
imap(by_cgroups,function(seu,ctype) {

	seu$comp <- "ctrl"
	seu$comp[seu$tgroup=="post_inf"] <- "expr"
	region_col <- "seurat_clusters"
	
	rds_fn <- file.path(args$outd,sprintf("%s_tgroup_per_cluster.rds",ctype))
	message(rds_fn)
	if (file.exists(rds_fn)) {
		deg_dt <- readRDS(rds_fn)
	} else {
		deg_dt <- diff_from_integration(seu,
																		cluster_col=region_col,
																		min_cell_cnts=3,
																		method=args$dge_method,
																		ncpu=3,
																		debug=args$debug)
		saveRDS(deg_dt,file=rds_fn)
	}
	
	tsv_fn = file.path(args$outd,sprintf("%s_btn_tgroup.tsv",ctype))
	fwrite(deg_dt,file=tsv_fn,sep="\t")
	
	title2 <- sprintf("DEG(Post/Pre);%s;FDR<0.05",ctype)
	
	pdf_pref2 <- file.path(args$outd,sprintf("tgroup_per_cluster_%s",ctype))
	
	#removed from ms_fig (was in Fig2C, 2021.04.26 version)
	heatmap_annotated_clusters(title2,
														 deg_dt[(gene %in% get_immreg_genes(extended=1) & p_val_adj<0.05),],
														 smeta = as.data.table(seu@meta.data),
														 pdf_pref = pdf_pref2,
														 apvalco=args$pval,
														 incl_genes="TIGIT",
														 zscored=0,
														 ctypes=ctype,
														 majcn=20,
														 mincn=10,
														 row_split_cnt=5,
														 col_split_cnt=4,
														 width1=6,
														 height1=10,
														 debug2=0)
	
	#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
	#grid.text()
	NA
})
# dev.off()

# ------------------------
message("(ADT)DEG between timepoint samples at each cell group")
region_col <- "whole"
plist <- imap(by_cgroups,function(seu,ctype) {
	# ctype="CD8"
	# ------->
	# DEG(D14/D0)
	seu$comp <- "unk"
	seu$comp[seu$tpoint=="D0"] <- "ctrl"
	seu$comp[seu$tpoint=="D14"] <- "expr"
	
	deg_d0d14_dt <- diff_from_integration(seu,
																				cluster_col=region_col,
																				assay2="Protein",
																				min_cell_cnts=3,
																				method=args$dge_method,
																				ncpu=1,
																				debug=args$debug)
	
	deg_d0d14_dt$diff<-"D0_D14"
	
	# ------->
	# DEG(D30/D14)
	seu$comp <- "unk"
	seu$comp[seu$tpoint=="D14"] <- "ctrl"
	seu$comp[seu$tpoint=="D30"] <- "expr"
	
	deg_d14d30_dt <- diff_from_integration(seu,
																				 cluster_col=region_col,
																				 assay2="Protein",
																				 min_cell_cnts=3,
																				 method=args$dge_method,
																				 ncpu=1,
																				 debug=args$debug)
	
	deg_d14d30_dt$diff<-"D14_D30"
	
	deg_dt=rbind(deg_d0d14_dt,deg_d14d30_dt)

	deg_dt[,protein:=gsub("-TotalSeqB","",gene)]
	deg_dt <- deg_dt[(!(protein %in% c("CD4","CD8A")) & p_val_adj<0.05),]
	deg_dt$ctype = ctype
	
	
	by_ctypes <- split(deg_dt,by="cluster_id")
	
	m1<-acast(deg_dt,protein ~ diff,value.var = "avg_logFC")
	
	m1[is.na(m1)]=0
	
	myPal <- get_red_white_blue_for_heatmap(m1,
																					debug2=0)
	#Fig2.G
	as.ggplot(pheatmap(m1,
																color=myPal$color,
																breaks=myPal$breaks,
																show_colnames = T,
																main=sprintf("ADT(DEG(D14/D0),DEG(D30/D14))\n%s;FDR<0.05",ctype),
																cluster_cols = F,
																silent=T))
})
pdf_fn = file.path(args$outd,"adt_deg_per_tpoint_ctype.pdf")
pdf(file=pdf_fn,width=8,height=6)
wrap_plots(plist,ncol=2)
dev.off()

# ===========
# pseudobulk
pbulk.m = load_pseudobulk(seu)
pmeta <- readRDS(get_meta_rds())

de2in <- list()
de2in$coldata = pmeta[match(colnames(pbulk.m$CD8),sample),]

de2in$coldata$comp <- "ctrl"
de2in$coldata[tpoint!="D0",comp:="expr"]
de2in$coldata = as.data.frame(de2in$coldata,row.names = de2in$coldata$sample)

de2in$cts <- as.matrix(pbulk.m$CD8)
design_f_str = "~ patient + comp"
deg.dt <- run_deseq2(de2in,args$outd,design_f_str,shrink_lfc=FALSE,alpha=0.1,exptag="CD8_tgroup",debug2=0)

deg.dt = deg.dt[!is.na(padj),]

setnames(deg.dt,'rn','gene')
setnames(deg.dt,'log2FoldChange','logFC')
setnames(deg.dt,'padj','adj.P.Val')

#Fig2.B
p=EnhancedVolcanoDESeq2(deg.dt[gene %in% get_immreg_genes(extended=1),],
											 column2disp="gene",
											 AdjustedCutoff=0.05, 
											 LabellingCutoff=0.05, 
											 LFCCutoff=0.5, 
											 main=sprintf("DEG LogFC(%s/%s);%s in pseudo bulk mRNA","post_inf","pre_inf","CD8"),
											 debug2=0)

pdf_fn <- file.path(args$outd,"deseq2_deg_by_tgroup_CD8.pdf")
pdf(file=pdf_fn,width=8,height=10)
plot(p)
dev.off()
# ========================
message("split cells into tgroup ...")
seus <- split_seurat_by_tgroup(seu,debug2=0)

plist <- lapply(names(seus),function(tgroup) {
	rds_file <- get_out_fpath(args,sprintf("%s_deg_on_treatment.rds",tgroup))
	if (args$reuse==1 & file.exists(rds_file)) {
		deg_dt <- readRDS(rds_file)
	} else {
		seu <- seus[[tgroup]]
		seu$comp <- "ctrl"
		seu$comp[seu$resp=="NR"] <- "expr"
		
		region_col <- "cd4_cd8"
		deg_dt <- diff_from_integration(seu,
																		cluster_col=region_col,
																		min_cell_cnts=3,
																		method=args$dge_method,
																		ncpu=3,
																		debug=0)
		
		saveRDS(deg_dt,file=rds_file)
	}
	
	deg_wo_clusters(deg_dt[cluster_id!="unk",],
									domain1="CD4",
									domain2="CD8",
									ctrl1="R",
									expr2="NR",
									pdf_file=get_out_fpath(args,sprintf("%s_deg_on_treatment.pdf",tgroup)),
									mtitle = sprintf("any[%s]",tgroup),
									apval.co=0.01,
									logFc.co=0.3,
									debug2=0)
})

# ===========
# pseudobulk

pbulk.m = load_pseudobulk(seu)
pmeta <- readRDS(get_meta_rds())

de2in <- list()

de2in$coldata = pmeta[match(colnames(pbulk.m$CD8),sample),]
de2in$coldata = de2in$coldata[tpoint!="D0",] #select only post_inf
de2in$coldata$comp <- "ctrl"
de2in$coldata[outcome=="NR",comp:="expr"]

#resync w/ count matrix
pbulk.m$CD8 = pbulk.m$CD8[,de2in$coldata$sample]

de2in$coldata = as.data.frame(de2in$coldata,row.names = de2in$coldata$sample)

de2in$cts <- as.matrix(pbulk.m$CD8)
design_f_str = "~ tpoint + comp"
deg.dt <- run_deseq2(de2in,args$outd,design_f_str,shrink_lfc=FALSE,alpha=0.1,exptag="CD8_post_inf_resp",debug2=0)

deg.dt = deg.dt[!is.na(padj),]

setnames(deg.dt,'rn','gene')
setnames(deg.dt,'log2FoldChange','logFC')
setnames(deg.dt,'padj','adj.P.Val')

p=EnhancedVolcanoDESeq2(deg.dt[gene %in% get_immreg_genes(extended=1),],
												column2disp="gene",
												AdjustedCutoff=0.05, 
												LabellingCutoff=0.05, 
												LFCCutoff=0.5, 
												main=sprintf("DEG LogFC(%s/%s);%s in pseudo bulk mRNA","NR","R","CD8"),
												debug2=0)

pdf_fn <- file.path(args$outd,"deseq2_deg_by_resp_post_CD8.pdf")
pdf(file=pdf_fn,width=8,height=10)
plot(p)
dev.off()

# --------------
de2in <- list()
pbulk.m = load_pseudobulk(seu)

de2in$coldata = pmeta[match(colnames(pbulk.m$CD8),sample),]
de2in$coldata = de2in$coldata[tpoint=="D0",] #select only post_inf
de2in$coldata$comp <- "ctrl"
de2in$coldata[outcome=="NR",comp:="expr"]

#resync w/ count matrix
pbulk.m$CD8 = pbulk.m$CD8[,de2in$coldata$sample]

de2in$coldata = as.data.frame(de2in$coldata,row.names = de2in$coldata$sample)

de2in$cts <- as.matrix(pbulk.m$CD8)
design_f_str = "~ comp"
deg.dt <- run_deseq2(de2in,args$outd,design_f_str,shrink_lfc=FALSE,alpha=0.1,exptag="CD8_pre_inf_resp",debug2=0)

deg.dt = deg.dt[!is.na(padj),]

setnames(deg.dt,'rn','gene')
setnames(deg.dt,'log2FoldChange','logFC')
setnames(deg.dt,'padj','adj.P.Val')

p=EnhancedVolcanoDESeq2(deg.dt[gene %in% get_immreg_genes(extended=1),],
												column2disp="gene",
												AdjustedCutoff=0.05, 
												LabellingCutoff=0.05, 
												LFCCutoff=0.5, 
												main=sprintf("DEG LogFC(%s/%s);%s in pseudo bulk mRNA","NR","R","CD8"),
												debug2=0)

pdf_fn <- file.path(args$outd,"deseq2_deg_by_resp_pre_CD8.pdf")
pdf(file=pdf_fn,width=8,height=10)
plot(p)
dev.off()

# ========================

message("highlight immreg genes ...")

plist <- lapply(names(seus),function(tgroup) {
	rds_file <- get_out_fpath(args,sprintf("%s_deg_on_treatment.rds",tgroup))
	if (args$reuse==1 & file.exists(rds_file)) {
		deg_dt <- readRDS(rds_file)
	} else {
		seu <- seus[[tgroup]]
		seu$comp <- "ctrl"
		seu$comp[seu$resp=="NR"] <- "expr"
		
		region_col <- "cd4_cd8"
		deg_dt <- diff_from_integration(seu,
																		cluster_col=region_col,
																		min_cell_cnts=3,
																		method=args$dge_method,
																		ncpu=3,
																		debug=0)
		
		saveRDS(deg_dt,file=rds_file)
	}
	
	#fig3.C
	deg_wo_clusters2(deg_dt[cluster_id!="unk",],
									domain1="CD4",
									domain2="CD8",
									genes_incl = get_immreg_genes(),
									ctrl1="R",
									expr2="NR",
									pdf_file=get_out_fpath(args,sprintf("%s_immreg_deg_on_treatment.pdf",tgroup)),
									mtitle = sprintf("immreg[%s]",tgroup),
									apval.co=1e-5,
									logFc.co=0.25,
									min_logfc2disp=0.3,
									topk1 = 15,
									topk2 = 15,
									debug2=0)
	
})
