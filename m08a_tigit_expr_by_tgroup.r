#0. load shared library
source('lib/lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 outd="out/cart/tigit_expr_by_tgroup",
									 dge_method="MAST",
									 ncpu=2,
									 reuse=1,
									 debug=0)

message("loading seus ...")

harmony_rds=locate_harmony_rds()
seu <- readRDS(harmony_rds$seu)
seu <- add_meta2(seu)

# ========================
message("split cells into tgroup and violin plot of checkpoint receptor expression at pre/post per each patient...")
my_comps <- list(c("NR","R"))
goi <- c("CTLA4","LAG3","HAVCR2","PDCD1","TIGIT","VSIR")

intra_deg_rds <- file.path(args$outd,"intra_deg_by_tgroup.rds")

by_pats <- split_seurat_by_patient(seu)

if (args$reuse==1 & file.exists(intra_deg_rds)){
	deg_dts <- readRDS(intra_deg_rds)
} else {
	# perform intra DEG between pre vs. post to get logFC/p-value on the genes of interest
	deg_dts=imap(by_pats,function(by_pat,patid){
		deg_dt <- NA
		if (length(unique(by_pat$tgroup))>1){
			by_pat$comp <- "ctrl"
			by_pat$comp[by_pat$tgroup=="post_inf"] <- "expr"
			
			deg_dt <- diff_from_integration(by_pat,
																			cluster_col="cd4_cd8",
																			min_cell_cnts=3,
																			method=args$dge_method,
																			ncpu=2,
																			debug=0)
			
		}
		deg_dt
	})
	deg_dts <- deg_dts[!is.na(deg_dts)]
	saveRDS(deg_dts,file=intra_deg_rds)
}

# ================
intra_deg_dt <- rbindlist(imap(deg_dts,function(deg_dt,patid) {
	deg_dt$patid <- patid
	deg_dt
}))

patids <- unique(intra_deg_dt$patid)

tigit_on_pct.dt <- rbindlist(imap(by_pats[patids],function(by_pat,patid) {
	rbindlist(imap(split_seurat_by_cellgroup(by_pat),function(by_cgroup,cgroup) {
		by_tgroups <- split_seurat_by_tgroup(by_cgroup)
		
		pre_table=as.data.table(table(GetAssayData(by_tgroups$pre_inf)['TIGIT',]>0))
		post_table=as.data.table(table(GetAssayData(by_tgroups$post_inf)['TIGIT',]>0))
		
		pre_table[,pct:=N/pre_table[,sum(N)]]
		post_table[,pct:=N/post_table[,sum(N)]]
		data.table(patid=patid,
							 cgroup=cgroup,
							 pct.1=pre_table[V1==TRUE,pct],
							 pct.2=post_table[V1==TRUE,pct])
		}))
	}))

tigit_on_pct.dt$avg_logFC<-NA
tigit_on_pct.dt$p_val<-1
tigit_on_pct.dt$p_val_adj<-NA
tigit_on_pct.dt$gene <- "TIGIT"
tigit_on_pct.dt$resp<-"R"
setnames(tigit_on_pct.dt,'cgroup','cluster_id')
tigit_on_pct.dt[patid %in% c("P3","P8"),resp:="NR"]

cgroup <- "CD8"

tigit_on_pct.dt2 <- tigit_on_pct.dt[cluster_id=="CD8"]

intra_deg_dt2 <- intra_deg_dt[gene=="TIGIT" & cluster_id=="CD8",]
intra_deg_dt2$resp <- "R"
intra_deg_dt2[patid %in% c("P3","P8"),resp:="NR"]

intra_deg_dt2 = rbind(tigit_on_pct.dt2[!(patid %in% intra_deg_dt2$patid)],
			intra_deg_dt2)

setnames(intra_deg_dt2,'pct.2','pre_inf')
setnames(intra_deg_dt2,'pct.1','post_inf')

intra_deg_dt2[,resp:=factor(resp,level=c("R","NR"))]

on_pct.dt <- melt(intra_deg_dt2,
									id.vars = c("patid"),
									measure.vars = c("pre_inf","post_inf"),
									variable.name = "tgroup",
									value.name = "on_pct")

on_pct.dt$resp <- intra_deg_dt2[match(on_pct.dt$patid,patid),resp]

on_pct.dt2 <- cbind(on_pct.dt,intra_deg_dt2[match(on_pct.dt$patid,patid),list(avg_logFC,p_val_adj)])

on_pct.dt2$annot <- ""
on_pct.dt2[tgroup=="pre_inf",annot:=patid]
on_pct.dt2[tgroup=="post_inf",annot:=sprintf("Log2FC:%3.2g(FDR:%3.2g)",avg_logFC,p_val_adj)]

#fig3.E
p = ggplot(on_pct.dt2,aes(x=tgroup,
										 y=on_pct,
										 group=patid,
										 color=resp,
										 label=annot)) +
	geom_line()+
	geom_point()+
	geom_text_repel(color="black",min.segment.length = 0, seed = 42, box.padding = 0.5)+
	labs(x="timepoint",y="The % of the cells TIGIT expressed",title=sprintf("%s;The %% of cells TIGIT expressed before/after CART infusion",cgroup)) +
	theme_minimal()

pdf_fn <- file.path(args$outd,"TIGIT_status_by_tgroup.pdf")
pdf(file=pdf_fn,width=10,height=8)
plot(p)
dev.off()

tsv_fn <- file.path(args$outd,"TIGIT_status_by_tgroup_CD8.tsv")
fwrite(on_pct.dt2,tsv_fn,sep="\t")
