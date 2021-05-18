# this is from 01_aucell_signature.r but it is for a general usage
source("lib/lib_project.r")

# ==============
FeaturePlotSingle2 <- function(obj, feature, metadata_column, assay2="RNA", red2="umap",minimal=0,maximal=1, sample_cellcnt=30590, ...) {
	
	obj <- obj[, sample(colnames(obj), size=sample_cellcnt, replace=F)]
	
	all_cells<- colnames(obj)
	# groups<- levels(obj@meta.data[, metadata_column])
	groups<- sort(unique(obj@meta.data[[metadata_column]]))
	
	# the minimal and maximal of the value to make the legend scale the same.
	
	ps<- list()
	for (group in groups) {
		subset_indx<- obj@meta.data[, metadata_column] == group
		subset_cells<- all_cells[subset_indx]
		p<- FeaturePlot(obj, features = feature, cells= subset_cells, reduction=red2, raster = FALSE, ...) +
			scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
			ggtitle(group) +
			theme(plot.title = element_text(size = 10, face = "bold"))
		ps[[group]]<- p
	}
	return(ps)
}

# ==============
FeaturePlotSingle3 <- function(obj, feature, assay2="RNA", red2="umap",minimal=0,maximal=1, sample_cellcnt=26213,title2="featureMap", ...) {
	
	obj <- obj[, sample(colnames(obj), size=sample_cellcnt, replace=F)]
	
	# the minimal and maximal of the value to make the legend scale the same.
	p <- FeaturePlot(obj, features = feature, reduction=red2, raster = FALSE, ...) +
		scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
		ggtitle(title2) +
		theme(plot.title = element_text(size = 10, face = "bold"))
	
	return(p)
}

violin_dysfs <- function(sm.melted,title2,cluster_ids,gene_sets,debug2=0) {
	if (debug2==1){browser()}
	my_comps <- list(c("R","NR"))
	
	plst <- list()
	plst[[1]] <- ggplot(sm.melted[(seurat_clusters %in% cluster_ids) & (gene_set %in% gene_sets),],aes(x=resp,y=dysfunction_score,fill=resp)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		stat_compare_means(method="wilcox.test",
											 hide.ns=T,
											 vjust=0.8,
											 label = "p.signif",
											 comparisons = my_comps) +
		facet_grid(rows=vars(gene_set),scales="free_y") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
	
	
	plst[[2]] <-ggplot(sm.melted[(seurat_clusters %in% cluster_ids) & (gene_set %in% gene_sets),],aes(x=resp,y=dysfunction_score,fill=resp)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		stat_compare_means(method="wilcox.test",
											 hide.ns=T,
											 vjust=0.8,
											 label = "p.signif",
											 comparisons = my_comps) +
		facet_grid(cols=vars(seurat_clusters),rows=vars(gene_set),scales="free_y") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
	
	plst
}

violin_dysfs_cellid <- function(sm.melted,title2,gene_sets) {
	my_comps <- list(c("R","NR"))
	ggplot(sm.melted[(gene_set %in% gene_sets),],aes(x=resp,y=dysfunction_score,fill=resp)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		stat_compare_means(method="wilcox.test",
											 hide.ns=T,
											 vjust=0.8,
											 label = "p.signif",
											 comparisons = my_comps) +
		facet_grid(cols=vars(pmid_30726743),rows=vars(gene_set),scales="free_y") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
}

violin_dysfs_tigit_w_pvalue <- function(sm.melted2,title2,debug2=0) {
	if (debug2==1){browser()}
	
	y.pos2 <- sm.melted2[,.(y.position=max(dysfunction_score)),by=c("gene_set")]
	
	astat.test = as.data.table(compare_means(dysfunction_score ~ TIGIT_status,
																					 data = sm.melted2, 
																					 group.by = c("gene_set"),
																					 method="wilcox.test"))
	astat.test= merge(astat.test,
									 y.pos2,
									 by.x=c("gene_set"),
									 by.y=c("gene_set"),
									 all.x=T)
	
	y.pos2 <- sm.melted2[,.(y.position=max(dysfunction_score)),by=c("seurat_clusters","gene_set")]
	
	stat.test = as.data.table(compare_means(dysfunction_score ~ TIGIT_status,
																					data = sm.melted2, 
																					group.by = c("seurat_clusters","gene_set"),
																					method="wilcox.test"))
	
	stat.test= merge(stat.test,
									 y.pos2,
									 by.x=c("seurat_clusters","gene_set"),
									 by.y=c("seurat_clusters","gene_set"),
									 all.x=T)
	
	my_comps=c("OFF","ON")
	
	# gene_sets = c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional")
	
	plst <- list()
	
	plst[[1]]=ggplot(sm.melted2,aes(x=TIGIT_status,y=dysfunction_score,color=TIGIT_status)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		facet_grid(rows=vars(gene_set),scales="free_y") +
		stat_pvalue_manual(astat.test,
											 xmin="group1", 
											 xmax="group2",
											 label = "p.signif") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
	
	
	plst[[2]]=ggplot(sm.melted2,aes(x=TIGIT_status,y=dysfunction_score,color=TIGIT_status)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		facet_grid(cols=vars(seurat_clusters),rows=vars(gene_set),scales="free_y") +
		stat_pvalue_manual(stat.test,
											 xmin="group1", 
											 xmax="group2",
											 label = "p.signif") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
	
	plst
}

violin_dysfs_tigit_cellid_w_pvalue <- function(sm.melted2,title2,debug2=0) {
	
	y.pos2 <- sm.melted2[,.(y.position=max(dysfunction_score)),by=c("pmid_30726743","gene_set")]
	
	stat.test = as.data.table(compare_means(dysfunction_score ~ TIGIT_status,
																					data = sm.melted2, 
																					group.by = c("pmid_30726743","gene_set"),
																					method="wilcox.test"))
	
	stat.test= merge(stat.test,
									 y.pos2,
									 by.x=c("pmid_30726743","gene_set"),
									 by.y=c("pmid_30726743","gene_set"),
									 all.x=T)
	
	my_comps=c("OFF","ON")
	
	# gene_sets = c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional")
	
	ggplot(sm.melted2,aes(x=TIGIT_status,y=dysfunction_score,color=TIGIT_status)) +
		geom_violin(position = position_dodge(width = 0.9)) +
		stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
								 geom="pointrange", color="black") +
		facet_grid(cols=vars(pmid_30726743),rows=vars(gene_set),scales="free_y") +
		stat_pvalue_manual(stat.test,
											 xmin="group1", 
											 xmax="group2",
											 label = "p.signif") +
		theme_bw()+
		theme(axis.text.x = element_blank(),legend.position="bottom") + ggtitle(title2)
	
}

# ==================
fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd=get_wkd(get_projd0(fpath_dt$parentd),fpath_dt$fpath),
									 cluster_col="seurat_clusters",
									 topr = 0.2,
									 ncpu=2,
									 reuse=1,
									 debug=0)

message("loading seus ...")
harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)

# -------------------
#label TIGIT-on/off info to metadata
tigit_expr <- GetAssayData(seu,assay="RNA",slot = "data")['TIGIT',]
seu <- AddMetaData(seu,tigit_expr,col.name = "TIGIT")
indic <- rep("OFF",length(tigit_expr))
indic[tigit_expr>0] <- "ON"
seu <- AddMetaData(seu,indic,col.name = "TIGIT_status")

# 0. select CD8+ cells
seus <- mclapply(split_seurat_by_sample(seu),function(seu) {
	extract_cell_type(seu,ctype = "CD8")
}
,mc.cores=8
)
S <- length(seus)

# -------------
rds_fn <- file.path(args$outd,"cd8_cell_rankings.rds")
if (file.exists(rds_fn)) {
	cells_rankings <- readRDS(rds_fn)
} else {
	cells_rankings <- imap(seus,function(seu,sname) {
		message(sname)
	
		 AUCell_buildRankings(as.matrix(GetAssayData(seu,
		 																						assay="RNA",
		 																						slot="counts")),
		 										 nCores=8,
		 										 plotStats=FALSE,
		 										 verbose=FALSE)
	})
	names(cells_rankings) <- names(seus)
	saveRDS(cells_rankings,file=rds_fn)
}
#rm(seu)
# --------
# collect genes to assess the dysfunction score from multiple studies

topr=0.2
topkj=500

gene_list <- dysfunction_genes_list(topK = topkj)

rds_fn = file.path(args$outd,"sma_dt.rds")
if (args$reuse==1 & file.exists(rds_fn)) {
	sma.dt = readRDS(rds_fn)
} else {
	sma.dt <- rbindlist(imap(cells_rankings,function(cells_ranking,sname) {
		cells_ranking=cells_rankings[[sname]]
		seuj <- seus[[sname]]
		
		sm.dt<-as.data.table(seuj@meta.data,keep.rownames = T)
		for (gsn in names(gene_list)) {
			message(sprintf("%s;%s",sname,gsn))
			cells_AUC <- AUCell_calcAUC(setdiff(gene_list[[gsn]],'TIGIT'),
																	cells_rankings[[sname]],
																	nCores = 8,
																	aucMaxRank = ceiling(topr * nrow(cells_ranking)))
			
			auc_mtrx <- cells_AUC@assays@data$AUC
			
			sm.dt[[gsn]] <- auc_mtrx[1,match(sm.dt$rn,colnames(auc_mtrx))]
		}
		
		sm.dt[,list(rn,long2015.lite,long2015.full,sadefeldman2018.CD8T_exhaust,leun2020.CD8T_dysfunctional,TIGIT,TIGIT_status,pmid_30726743,seurat_clusters,resp,tgroup)]
	}))
	
	saveRDS(sma.dt,file=rds_fn)
}

# ===============
by_tgroups = split_seurat_by_tgroup(extract_cell_type(seu,ctype = "CD8"))

by_tgroups <- lapply(by_tgroups,function(by_tgroup) {
	mdt <- as.data.table(by_tgroup@meta.data,keep.rownames = T)
	
	mdt <- merge(mdt,sma.dt[,list(rn,long2015.lite,long2015.full,sadefeldman2018.CD8T_exhaust,leun2020.CD8T_dysfunctional)],by="rn",all.x=T)
	
	by_tgroup@meta.data <- data.frame(mdt[,!"rn",with=F],row.names = mdt$rn)
	by_tgroup
})

plist<-list()

plist[[1]]<-FeaturePlotSingle3(by_tgroups$pre_inf,feature="TIGIT",red2="tsne",minimal=0,maximal=3.8,title2 = "CD8;pre_inf;TIGIT")

#Fig4.B(left)
plist[[2]]<-FeaturePlotSingle3(by_tgroups$post_inf,feature="TIGIT",red2="tsne",minimal=0,maximal=3.8,title2 = "CD8;post_inf;TIGIT")

#Fig4.C(right)
plist[[3]]<-FeaturePlotSingle3(by_tgroups$pre_inf,feature="long2015.full",red2="tsne",minimal=0,maximal=0.35,title2 = "CD8;pre_inf;long2015.full")

plist[[4]]<-FeaturePlotSingle3(by_tgroups$post_inf,feature="long2015.full",red2="tsne",minimal=0,maximal=0.35,title2 = "CD8;post_inf;long2015.full")

plist[[5]]<-FeaturePlotSingle3(by_tgroups$pre_inf,feature="sadefeldman2018.CD8T_exhaust",red2="tsne",minimal=0,maximal=1,title2 = "CD8;pre_inf;sadefeldman2018.CD8T_exhaust")
plist[[6]]<-FeaturePlotSingle3(by_tgroups$post_inf,feature="sadefeldman2018.CD8T_exhaust",red2="tsne",minimal=0,maximal=1,title2 = "CD8;post_inf;long2015.full")

plist[[7]]<-FeaturePlotSingle3(by_tgroups$pre_inf,feature="leun2020.CD8T_dysfunctional",red2="tsne",minimal=0,maximal=1,title2 = "CD8;pre_inf;leun2020.CD8T_dysfunctional")

plist[[8]]<-FeaturePlotSingle3(by_tgroups$post_inf,feature="leun2020.CD8T_dysfunctional",red2="tsne",minimal=0,maximal=1,title2 = "CD8;post_inf;leun2020.CD8T_dysfunctional")

p=wrap_plots(plist,ncol=4,nrow=2)
pdf_fn <- file.path(args$outd,"cd8_TIGIT_dysfscore_on_tsne.pdf")
pdf(pdf_fn,width = 14, height=6)
plot(p)
dev.off()

sma.dt$rn <- NULL
# ============
sm.melted = melt(sma.dt,
								 id.vars = c("TIGIT_status","pmid_30726743","seurat_clusters","resp","tgroup","TIGIT"),
								 measure.vars = c("long2015.lite","long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional"),variable.name = "gene_set",value.name = "dysfunction_score")

sm.melted[,resp:=factor(resp,level=c("R","NR"))]

plist <- list()
title2 = sprintf("CD8 Tcell dysfunctional score [clusters:0-10]")
plst = violin_dysfs(sm.melted,title2,c(0:10),c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional"))
plist[[1]] <- wrap_plots(plst,ncol=2,widths = c(1,8))

title2 = sprintf("CD8 Tcell dysfunctional score [clusters:11-]")
plst = violin_dysfs(sm.melted,title2,c(11:21),c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional"))
plist[[2]] <- wrap_plots(plst,ncol=2,widths = c(1,8))

p <- wrap_plots(plist,nrow=2)

#Fig3.E
pdf_fn <- file.path(args$outd,"cd8_dysfscore_per_region.pdf")
pdf(pdf_fn,width = 12, height=12)
plot(p)
# ---------------

#SFig3.B
title2 = sprintf("CD8 Tcell dysfunctional score per subtype")
p = violin_dysfs_cellid(sm.melted,title2,c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional"))
plot(p)

# ---------------
# Fig4.C(left)
title2 = sprintf("CD8 Tcell dysfunctional score compared by TIGIT[cluster:0-10]")
plist <- list()

sm.melted2 = sm.melted[(seurat_clusters %in% c(0:10) & gene_set %in% c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional")),]
plst = violin_dysfs_tigit_w_pvalue(sm.melted2,title2)
plist[[1]] <- wrap_plots(plst,ncol=2,widths = c(1,8))

title2 = sprintf("CD8 Tcell dysfunctional score compared by TIGIT[cluster:11-]")
sm.melted2 = sm.melted[(seurat_clusters %in% c(11:20) & gene_set %in% c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional")),]
plst = violin_dysfs_tigit_w_pvalue(sm.melted2,title2)
plist[[2]] <- wrap_plots(plst,ncol=2,widths = c(1,8))

p <- wrap_plots(plist,nrow=2)
plot(p)

# -------------
title2 = sprintf("CD8 Tcell dysfunctional score by TIGIT per subtype")
sm.melted2 <- sm.melted[gene_set %in% c("long2015.full","sadefeldman2018.CD8T_exhaust","leun2020.CD8T_dysfunctional"),]
p = violin_dysfs_tigit_cellid_w_pvalue(sm.melted2,title2)
plot(p)

dev.off()

# =================

#Fig4.B(right)
plist <- list()
plist[[1]] <- annotate_figure(RidgePlot(subset(seu,TIGIT_status=="ON"),
																					assay="RNA",
																					# idents=sname,
																					slot="data",
																					features = "TIGIT", 
																					ncol = 1,
																					group.by = "pmid_30726743"),
																top = text_grob(label = "TIGIT expression level \n in TIGIT[ON] cells", face = 'bold'))


plist[[2]] <- annotate_figure(RidgePlot(subset(seu,PDCD1>0.),
															 assay="RNA",
															 # idents=sname,
															 slot="data",
															 features = "PDCD1", 
															 ncol = 1,
															 group.by = "pmid_30726743"),
										 top = text_grob(label = "PDCD1 expression level \n in PDCD1[ON] cells", face = 'bold'))

p = wrap_plots(plist,ncol=2)
pdf_fn <- file.path(args$outd,"cd8_ridge_per_cellid.pdf")
pdf(file=pdf_fn,width=12,height=8)
plot(p)
dev.off()