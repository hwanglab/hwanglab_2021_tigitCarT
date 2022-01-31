#0. load shared library
source('lib/lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd=get_wkd(get_projd0(fpath_dt$parentd),fpath_dt$fpath),
									 ncpu=2,
									 reuse=0,
									 debug=0)

violin_plots2 <- function(seuj,
													goi,
													my_comps,
													field2cmp,
													pdf_fn,
													comp_label = "p.signif",
													width1=8,
													height1=10,
													title1="expressed",
													title2="dropout or not expressed",
													debug2=0) {
	
	if (debug2==1){browser()}
	
	m <- GetAssayData(seuj,assay = "RNA",slot = "data")
	goi_avail <- goi[goi %in% rownames(m)]
	
	gexpr.dt <- as.data.table(melt(as.matrix(m[goi_avail,]),varnames = c("gene","cbc"),value.name = "expression_level"))
	
	smdt = as.data.table(seuj@meta.data,keep.rownames = T)
	
	gexpr.dt = cbind(gexpr.dt,
									 smdt[match(gexpr.dt$cbc,rn),
									 		 list(orig.ident,tpoint,tgroup,cd4_cd8,resp,seurat_clusters,pmid_30726743)])
	
	# ----------------
	gexpr1.dt<-gexpr.dt[expression_level>0.,]
	
	p =ggplot(gexpr1.dt,
						aes_string(x=field2cmp, y='expression_level', fill=field2cmp)) +
		geom_violin(position = position_dodge(width = 0.5)) +
		stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",size=0.25,fun.args = list(mult = 1)) +
		stat_compare_means(aes_string(group=field2cmp),
											 # method="t.test",
											 method="wilcox.test",
											 hide.ns=T,
											 vjust=0.8,
											 label = comp_label,
											 comparisons = my_comps) +
		
		facet_grid(cols=vars(gene)) +
		theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
		theme_bw() +
		theme(legend.position="bottom") +
		ggtitle(title1)
	
	plist <- list()
	plist[[1]] <- p
	
	# ------------
	gexpr0.dt<-gexpr.dt[expression_level==0.,]
	
	count0.dt <- gexpr0.dt[,.N,by=c("orig.ident","gene",`field2cmp`)]
	denom<-gexpr.dt[,.N,by="orig.ident"]
	count0.dt$denom <- denom[match(count0.dt$orig.ident,orig.ident),N]
	
	count0.dt[,pct0:=100*(N/denom)]
	
	p <- ggplot(count0.dt,
							aes_string(x=field2cmp, y='pct0', fill=field2cmp)) +
		geom_violin(position = position_dodge(width = 0.5)) +
		stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",size=0.25,fun.args = list(mult = 1)) +
		stat_compare_means(aes_string(group=field2cmp),
											 method="t.test",
											 # method="wilcox.test",
											 hide.ns=T,
											 vjust=0.8,
											 label = comp_label,
											 # label = "p.signif",
											 comparisons = my_comps) +
		facet_grid(cols=vars(gene)) +
		theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
		theme_bw() +
		theme(legend.position="bottom") +
		ggtitle(title2)
	
	plist[[2]] <- p
	
	p <- wrap_plots(plist,nrow=2)
	
	pdf(file=pdf_fn,width=width1,height=height1)
	plot(p)
	dev.off()
}

violin_by_clusters_plots2 <- function(seuj,
																			goi,
																			my_comps,
																			field2cmp,
																			field2cmp_order,
																			pdf_fn,
																			width1=8,
																			height1=10,
																			comp_label = "p.signif",
																			title1="expressed",
																			debug2=0) {
	
	if (debug2==1){browser()}
	
	m <- GetAssayData(seuj,assay = "RNA",slot = "data")
	goi_avail <- goi[goi %in% rownames(m)]
	
	gexpr.dt <- as.data.table(melt(as.matrix(m[goi_avail,]),varnames = c("gene","cbc"),value.name = "expression_level"))
	
	smdt = as.data.table(seuj@meta.data,keep.rownames = T)
	
	gexpr.dt = cbind(gexpr.dt,
									 smdt[match(gexpr.dt$cbc,rn),
									 		 list(orig.ident,tpoint,tgroup,cd4_cd8,resp,seurat_clusters,pmid_30726743)])
	
	gexpr.dt[[field2cmp]] <- factor(gexpr.dt[[field2cmp]],level=field2cmp_order)
	
	# ----------------
	
	cluster_group = chunk_by_unitlen(sort(unique(gexpr.dt$seurat_clusters)),unit_len = 5)
	pdf(file=pdf_fn,width=width1,height=height1)
	
	lapply(cluster_group,function(clsgroup){
		p <- ggplot(gexpr.dt[seurat_clusters %in% clsgroup,],
								aes_string(x=field2cmp, y='expression_level', fill=field2cmp)) +
			geom_violin(position = position_dodge(width = 0.5)) +
			stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",size=0.25,fun.args = list(mult = 1)) +
			stat_compare_means(aes_string(group=field2cmp),
												 # method="t.test",
												 method="wilcox.test",
												 hide.ns=T,
												 vjust=0.8,
												 label = comp_label,
												 #label = "p.signif",
												 comparisons = my_comps) +
			
			facet_grid(cols=vars(gene),rows=vars(seurat_clusters)) +
			theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
			theme_bw() +
			theme(legend.position="bottom") +
			ggtitle(title1)
		plot(p)
		NA
	})
	dev.off()
}


violin_by_patient_plots <- function(seuj,
																		goi,
																		my_comps,
																		field2cmp,
																		field2cmp_order,
																		pdf_fn,
																		width1=8,
																		height1=10,
																		title1="expressed",
																		debug2=0) {
	
	if (debug2==1){browser()}
	
	m <- GetAssayData(seuj,assay = "RNA",slot = "data")
	goi_avail <- goi[goi %in% rownames(m)]
	
	gexpr.dt <- as.data.table(melt(as.matrix(m[goi_avail,]),
																 varnames = c("gene","cbc"),
																 value.name = "expression_level"))
	
	smdt = as.data.table(seuj@meta.data,keep.rownames = T)
	
	gexpr.dt = cbind(gexpr.dt,
									 smdt[match(gexpr.dt$cbc,rn),
									 		 list(orig.ident,tpoint,cd4_cd8,resp,patid,tgroup)])
	
	gexpr.dt[[field2cmp]] <- factor(gexpr.dt[[field2cmp]],level=field2cmp_order)
	
	# ----------------
	
	by_resps<-split(gexpr.dt,by="resp")
	
	plist = imap(by_resps,function(clsgroup,rgroup) {
		# clsgroup = by_resps[[1]]
		# rgroup = "R"
		
		#we are interested in pat sample that has both pre and post
		a=acast(clsgroup[,.N,by=c('patid','tgroup')], patid ~ tgroup,value.var = "N")
		a[is.na(a)]=0
		clsgroup = clsgroup[patid %in% rownames(a)[rowSums(a>0)>1],]
		
		p <- ggplot(clsgroup,
								aes_string(x=field2cmp, y='expression_level', fill=field2cmp)) +
			geom_violin(position = position_dodge(width = 0.5)) +
			stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",size=0.25,fun.args = list(mult = 1)) +
			# stat_compare_means(aes_string(group=field2cmp),
			# 									 # method="t.test",
			# 									 method="wilcox.test",
			# 									 hide.ns=T,
			# 									 vjust=0.8,
			# 									 label = "p.signif",
			# 									 comparisons = my_comps) +
			
			facet_grid(cols=vars(patid),rows=vars(gene)) +
			theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
			theme_bw() +
			theme(legend.position="bottom") +
			ggtitle(sprintf("[%s]%s",rgroup,title1))
		p
	})
	
	pdf(file=pdf_fn,width=width1,height=height1)
	p=wrap_plots(plist,ncol=2,widths = c(2,1))
	plot(p)
	dev.off()
}

harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)

# ===========
# load cluster size

csize = get_harmony_seui_smeta()

ctype_pct <- csize$smeta_smry[,.(avg_pct=mean(pct)),by=c("seurat_clusters","cd4_cd8")]
m = acast(ctype_pct,seurat_clusters ~ cd4_cd8,value.var = "avg_pct")
m[is.na(m)]<-0
m_scaled<-get_zscore(m,margin1 = 1)

seurat_clusters <- rownames(m_scaled)
maj_ctype_cluster <- list()
maj_ctype_cluster[['CD4']] <- seurat_clusters[m_scaled[,'CD4']>=m_scaled[,'CD8']]
maj_ctype_cluster[['CD8']] <- seurat_clusters[m_scaled[,'CD4']<m_scaled[,'CD8']]

# ========================
message("split cells into tgroup and violin plot of checkpoint receptor expression at pre/post per each patient...")
my_comps <- list(c("NR","R"))
goi <- c("CTLA4","LAG3","HAVCR2","PDCD1","TIGIT","VSIR","PVRIG")

seus <- split_seurat_by_cellgroup(seu,debug2=0)

plist <- lapply(names(seus),function(ctype) {
	# ctype <- "CD8"
	# tgroup<-"post_inf"
	
	message(sprintf("c:%s",ctype))
	
	seuj <- seus[[ctype]]
	
	#Fig3.F
	violin_by_patient_plots(seuj,
													goi,
													my_comps,
													field2cmp="tgroup",
													field2cmp_order = c("pre_inf","post_inf"),
													pdf_fn=file.path(args$outd,sprintf("%s_checkpoint_inhib_expr_per_patid.pdf",ctype)),
													width1=8,
													height1=10,
													title1=sprintf("[%s] checkpoint inhibitors expression level per patid",ctype),
													debug2=0)
})

# ========================
message("split cells into tgroup ...")
my_comps <- list(c("NR","R"))
goi <- c("CTLA4","LAG3","HAVCR2","PDCD1","TIGIT","VSIR")

seus <- split_seurat_by_tgroup(seu,debug2=0)

plist <- lapply(names(seus),function(tgroup) {
	# tgroup<-"post_inf"
	seu_by_ctypes <- split_seurat_by_cellgroup(seus[[tgroup]],debug2=0)
	lapply(names(seu_by_ctypes),function(ctype){
		# ctype <- "CD8"
		message(sprintf("t:%s;c:%s",tgroup,ctype))
		
		seuj <- seu_by_ctypes[[ctype]]
		
		maj_clusters <- maj_ctype_cluster[[ctype]]
		
		seujm <- subset(seuj,seurat_clusters %in% maj_clusters)
		
		violin_by_clusters_plots2(seujm,
															goi,
															my_comps,
															field2cmp="resp",
															field2cmp_order = c("R","NR"),
															comp_label = "p.signif",
															#comp_label = "p.format",
															pdf_fn=file.path(args$outd,sprintf("%s_%s_checkpoint_inhib_expr_per_resp.pdf",tgroup,ctype)),
															width1=8,
															height1=10,
															title1=sprintf("[%s;%s] checkpoint inhibitors expression level per resp",tgroup,ctype),
															debug2=0)
	})
})

# =================

seu_by_ctypes <- split_seurat_by_cellgroup(seu,debug2=0)

my_comps <- list(c("D0","D14"),c("D0","D30"),c("D14","D30"))

lapply(names(seu_by_ctypes),function(ctype) {
	pdf_fn <- file.path(args$outd,sprintf("%s_checkpoint_inhib_exprlev_per_tpoint.pdf",ctype))
	#Fig2.D
	violin_plots2(seu_by_ctypes[[ctype]],
								goi,
								my_comps,
								field2cmp="tpoint",
								pdf_fn=file.path(args$outd,sprintf("%s_checkpoint_inhib_expr_per_tpoint.pdf",ctype)),
								width1=8,
								comp_label = "p.signif",
								# comp_label = "p.format",
								height1=10,
								title1=sprintf("[%s] checkpoint inhibitors expression level (>0.) per timepoint",ctype),
								title2=sprintf("[%s] cell pct where checkpoint inhibitors never expressed",ctype),
								debug2=0)
})
