#0. load shared library
source('lib/lib_project.r')

# --------------------------
#1. locate PBMC srun1 and 2 QC'ed/cellID seurat objects
# note: collect all cell ID seurat samples from all 3 seqruns;assign a unified orig.ident before concatenating them

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd="out/cart/vcluster_heatmap",
									 cluster_col="seurat_clusters",
									 heatmap_title="CART",
									 ncpu=1,
									 topK=120,
									 req_genes="CD4,CD40LG,TNFRSF4,CD8A,CD8B,TIGIT",
									 pval=0.05,
									 reuse=1,
									 debug=0)

create_dir_if_not_exist(args$outd)

# ==================
harmony_rds=locate_harmony_rds()

args$harmony_rds <-harmony_rds$seu
args$marker_rds <-harmony_rds$deg
vcluster_heatmap(args) #Fig1.B, Fig1.C, Fig2.E, Fig2.F, Fig3.B

# ================
rds_fn <- file.path(args$outd,"CART_harmony_meta.rds")
if (args$reuse==1 & file.exists(rds_fn)) {
	smeta <- readRDS(rds_fn)
} else {
	seui <- readRDS_w_msg(harmony_rds$seu)
	smeta = as.data.table(seui@meta.data)
	rm(seui)
	saveRDS(smeta,file=rds_fn)
}

ret <- cluster_by_pmid_30726743_composition(smeta,debug2=args$debug) #SFig1.C

rds_fn <- file.path(args$outd,"pmid_30726743_pct_by_seurat_clusters.rds")
saveRDS(ret$cellid_comp,file=rds_fn)

pdf_fn = file.path(args$outd,"pmid_30726743_abundance_by_cluster.pdf")
pdf(file=pdf_fn,width=8,height=12)
plot(ret$plt)
dev.off()

# ================
# CD8 cell type comparison by patient response group
pmeta <- readRDS(get_meta_rds())

cellcnt_by_cluster <- smeta[,.N,by=c("orig.ident","pmid_30726743")]
cellcnt_by_sample <- smeta[,.(denom=.N),by=c("orig.ident")]
cellcnt_by_cluster = merge(cellcnt_by_cluster,cellcnt_by_sample,by="orig.ident")
cellcnt_by_cluster[,pct:=100.*N/denom]
cellcnt_by_cluster <- merge(cellcnt_by_cluster,pmeta,by.x="orig.ident",by.y="sample")
cellcnt_by_cluster[,outcome:=factor(outcome,level=c("R","NR"))]

my_comps <- list(c("R","NR"))

#SFig3.A
title2="T cell subtype in fraction and compared by the patient treatment response"
p <- ggplot(cellcnt_by_cluster[!is.na(pmid_30726743),],aes(x=outcome,y=pct,fill=outcome)) +
	geom_boxplot(position=position_dodge2(width=0.75)) +
	# ylim(c(-0.1,y_max))+
	geom_jitter(size=0.5) +
	stat_compare_means(method="wilcox.test",
										 hide.ns=T,
										 label = "p.signif",
										 vjust=0.8,
										 comparisons = my_comps) +
	facet_grid(cols=vars(pmid_30726743)) +
	theme_bw()+
	labs(title = "T cell subtype in fraction and compared by the patient treatment response\nCD4_TE(P.value=0.02733),CD8_EM(P.value=0.0439)",
			 y="Normalized cell fraction (%)") +
	theme(axis.text.x = element_blank(),legend.position="bottom")

pdf_fn = file.path(args$outd,"cell_type_comp_by_resp.pdf")
pdf(file=pdf_fn,width=12,height=5)
plot(p)
dev.off()

# =================
# to precompute p.value stats in the figure above
x=cellcnt_by_cluster[pmid_30726743=="CD8_EM" & outcome=="R",pct]
y=cellcnt_by_cluster[pmid_30726743=="CD8_EM" & outcome=="NR",pct]
ret=rankSumTest_tryCatch(x,y)

x=cellcnt_by_cluster[pmid_30726743=="CD4_TE" & outcome=="R",pct]
y=cellcnt_by_cluster[pmid_30726743=="CD4_TE" & outcome=="NR",pct]
ret=rankSumTest_tryCatch(x,y)


# ==============
if (debug2==1){browser()}
pmeta <- readRDS(get_meta_rds())

cellcnt_by_cluster <- smeta[,.N,by=c("orig.ident","seurat_clusters")]
cellcnt_by_sample <- smeta[,.(denom=.N),by=c("orig.ident")]
cellcnt_by_cluster = merge(cellcnt_by_cluster,cellcnt_by_sample,by="orig.ident")
cellcnt_by_cluster[,pct:=100.*N/denom]

cellcnt_by_cluster <- merge(cellcnt_by_cluster,pmeta,by.x="orig.ident",by.y="sample")

m <- dcast(cellcnt_by_cluster,
					 seurat_clusters + tpoint ~ patient, 
					 fun.aggregate = sum,
					 value.var = "pct")

agg_pct <- as.data.table(melt(m))
colnames(agg_pct) <- c("seurat_clusters","tpoint","patient","sum.pct")
agg_pct[,seurat_clusters:=factor(seurat_clusters)]

agg_denom = agg_pct[,.(denom=sum(sum.pct)),by=c("tpoint")]
agg_pct = merge(agg_pct,agg_denom,by="tpoint")
agg_pct[,norm.pct:=100*(sum.pct/denom)]

#SFig1.A
delim <- list(s1=c(0,11),e2=c(10,19))
plist <- lapply(1:2,function(j) {
	p <- ggplot(agg_pct[seurat_clusters %in% c(delim$s1[j]:delim$e2[j]),], 
							aes(fill=patient, y=norm.pct, x=seurat_clusters)) + 
		geom_bar(position="stack", stat="identity") +
		facet_grid(rows = vars(tpoint)) +
		theme_bw()+
		labs(title = "Patient proportion in each cluster",
				 y="Normalized cell fraction (%)")
	p
})
pdf_fn = file.path(args$outd,"patient_spct_by_cluster.pdf")
pdf(file=pdf_fn,width=10,height=8)
wrap_plots(plist,ncol2=2)
dev.off()

tsv_fn = file.path(args$outd,"patient_spct_by_cluster.tsv")
fwrite(agg_pct,file=tsv_fn,sep="\t")

# ==============
if (debug2==1){browser()}
pmeta <- readRDS(get_meta_rds())

cellcnt_by_cluster <- smeta[,.N,by=c("orig.ident","seurat_clusters","pmid_30726743","Phase")]
cellcnt_by_sample <- smeta[,.(denom=.N),by=c("orig.ident")]
cellcnt_by_cluster = merge(cellcnt_by_cluster,cellcnt_by_sample,by="orig.ident")
cellcnt_by_cluster[,pct:=100.*N/denom]
cellcnt_by_cluster <- merge(cellcnt_by_cluster,pmeta,by.x="orig.ident",by.y="sample")

m <- dcast(cellcnt_by_cluster,
					 seurat_clusters + pmid_30726743 ~ Phase, 
					 fun.aggregate = sum,
					 value.var = "pct")

agg_pct <- as.data.table(melt(m))
colnames(agg_pct) <- c("seurat_clusters","pmid_30726743","cell_cycle","sum.pct")
agg_pct[,seurat_clusters:=factor(seurat_clusters)]

agg_pct[is.na(pmid_30726743),pmid_30726743:="unk"]

agg_pct$cd4_cd8 <- "CD4"
agg_pct$cd4_cd8[grep("CD8",agg_pct$pmid_30726743)] <- "CD8"

agg_denom = agg_pct[,.(denom=sum(sum.pct)),by=c("cd4_cd8")]
agg_pct = merge(agg_pct,agg_denom,by="cd4_cd8")
agg_pct[,norm.pct:=100*(sum.pct/denom)]

#SFig1.B
delim <- list(s1=c(0,11),e2=c(10,19))
plist <- lapply(1:2,function(j) {
	p <- ggplot(agg_pct[seurat_clusters %in% c(delim$s1[j]:delim$e2[j]),], 
							aes(fill=cell_cycle, y=norm.pct, x=seurat_clusters)) + 
		geom_bar(position="stack", stat="identity") +
		facet_grid(rows = vars(cd4_cd8)) +
		theme_bw()+
		labs(title = "Abundance of the cell_cycle across the samples",
				 y="Normalized cell fraction (%)")
	p
})

pdf_fn = file.path(args$outd,"cellcycle_spct_by_cluster.pdf")
pdf(file=pdf_fn,width=10,height=8)
wrap_plots(plist,ncol2=2)
dev.off()

tsv_fn = file.path(args$outd,"cellcycle_spct_by_cluster.tsv")
fwrite(agg_pct,file=tsv_fn,sep="\t")

ifttt_notify(fpath_dt$fbase0)