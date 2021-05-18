#0. load shared library
source('lib/lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd="out/cart/celltype_spiderweb",
									 cluster_col="seurat_clusters,pmid_30726743",
									 dominant_cluster=1,
									 ncpu=2,
									 debug=0)

harmony_rds=locate_harmony_rds()
seui <- readRDS_w_msg(harmony_rds$seu)

# ================================
csf=""
if (args$dominant_cluster==1) {
	csf = "c0_10"
}

# ----------
# read seurat metadata
smeta <- as.data.table(seui@meta.data,keep.rownames = T)
smeta$tgroup <- "post_inf"
smeta[tpoint=="D0",tgroup:="pre_inf"]
smeta$cd4_cd8 <- "CD4"
smeta[grep("CD8",pmid_30726743),cd4_cd8:="CD8"]

# -----------
a2 <- smeta[,.N,by=c("orig.ident","tpoint","pmid_30726743","cd4_cd8")]
tcell <-a2[,.(denom=sum(N)),by=c("orig.ident")]
a2$pct <- 100. * (a2$N / tcell[match(a2$orig.ident,orig.ident),denom])

#Fig2.E
plst <- lapply(c("CD4","CD8"),function(ctype){
	avg_subtype_comp <- acast(a2[cd4_cd8==`ctype` & !is.na(pmid_30726743),], pmid_30726743 ~ tpoint, fun.aggregate = mean, value.var = "pct")
	
	subtype_comp.m <- melt(avg_subtype_comp,varnames = c("pmid_30726743","tpoint"),value.name = "avg_pct")
	
	ggplot(data=subtype_comp.m, aes(x=pmid_30726743, y=avg_pct, color=tpoint,group=tpoint,fill=tpoint)) + 
		geom_point() + 
		geom_line() +
		geom_polygon(size = 1, alpha= 0.2) + 
		facet_grid(cols = vars(tpoint)) +
		ylim(-0.1,max(subtype_comp.m$avg_pct)*1.1) +
		xlab("pmid_30726743") + 
		ylab("avg_pct") + 
		coord_polar() +
		ggtitle(sprintf('CART %s subtype profile per timepoint',ctype)) + 
		theme(text = element_text(size=8))+
		theme_minimal()
})
p = wrap_plots(plst,nrow=2)
pdf_fn <- file.path(args$outd,"pmid_30726743_prof_per_tpoint.pdf")
pdf(file=pdf_fn,width=8,height=6)
plot(p)
dev.off()
