# author1: hongc2@ccf.org
 
library(data.table)

if (T) {
	library(argparse)
	args_tmp <- commandArgs(trailingOnly = F)
	fname <- normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)]))
	message(sprintf("script_filepath:%s",fname))
	message(sprintf("script_file_directory:%s",dirname(fname)))
	
	prog_ver <- "CAR-T:v0.1"
	parser <- ArgumentParser(description=sprintf('[%s]to assign cell ID given a user cell marker genes; to generate the marker gene matrix, UMAP w/ the best predicted cell ID and marker expression level',prog_ver))
	
	parser$add_argument("-i", "--seu_rd", type="character", required=TRUE,
											dest="seu_rd",
											help="input seurat obj rd file path derived from seurat_solo.r")
	parser$add_argument("-x", "--sample_sheet", type="character", required=TRUE,
											dest="sample_sheet", default="",
											help="sample info excel/csv file where the columns contain sample,condition,Patient,Timepoint,type,comp,and sample_tag")
	
	parser$add_argument("-j", "--sheet_number", type="integer", required=FALSE,
											dest="sheet_number", default=1,
											help="sheet number in the sample sheet")
	parser$add_argument("-m", "--cell2grep", type="character", required=FALSE,
											dest="cell2grep",default="",
											help="cell type to analyze for grep ['']")
	
	parser$add_argument("-g", "--gois", type="character", required=FALSE,
											dest="gois", default ="",
											help="genes of interest for featplot; comma separated")
	
	parser$add_argument("-c", "--pmid_ref", type="character", required=TRUE,
											dest="pmid_ref",
											help="PMID of marker_csv, e.g., pmid_12345")
	
	parser$add_argument("-s", "--singler_ref_rd", type="character", required=TRUE,
											dest="singler_ref_rd",
											help="SingleR reference rd file path that contains ref.sce")
	
	parser$add_argument("-S", "--ref_label", type="character", required=FALSE,
											dest="ref_label",
											help="label.main or [label.fine]?")
	
	parser$add_argument("-r", "--reuse", type="integer", required=FALSE,
											dest="reuse", default=1,
											help="reuse the same prev cell ID if exists? 0:No,[1]:Yes")
	
	parser$add_argument("-o", "--outd", type="character", required=FALSE,
											dest="outd", default="out",
											help="output directory")
	
	parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
											dest="ncpu", default=4,
											help="ncpu? [4]")
	
	parser$add_argument("-d", "--debug", type="integer", required=FALSE,
											dest="debug", default=0,
											help="is debug mode? [0],1")
	
	parser$add_argument("-dd", "--debugd", type="character", required=FALSE,
											dest="debugd", default="",
											help="debug directory")
	
	args <- parser$parse_args()
} else {
	args <- data.table(seu_rd="seu.rds",
										 sample_sheet="doc/seqrun_info/cellranger_samplesheet.xlsx",
										 sheet_number=5,
										 cell2grep="",
										 singler_ref_rd="data/cell_id_markers/pmid_30726743/GSE107011_pbmc.rd",
										 gois="BATF,BTLA,CCR7,CD160,CD226,CD244,CD274,CD96,CREM,CTLA4,CXCR5,EOMES,FASLG,FGL2,FOS,FOSB,FOSL1,FOSL2,FOXO1,FOXP1,FOXP3,GATA3,GZMA,GZMB,GZMC,GZMK,GZMM,HAVCR2,IFNG,IL10,IL21,IL21R,IRF4,JUN,JUNB,JUND,KLRG1,LAG3,LRIG1,NECTIN3,NFAT5,NFATC1,NFATC2,NFATC2IP,NFATC3,NFKBIA,NR4A1,NR4A2,NR4A3,PDCD1,PDCD10,PDCD11,PDCD1LG2,PRDM1,PRF1,PVR,PVRIG,SELL,SIGIRR,SPRY2,TBET,TCF7,TGFB1,TGFB2,TIGIT,TNF,TOX,TOX2,VSIR",
										 ncpu=1,
										 reuse=1,
										 debug=0,
										 pmid_ref="pmid_30726743",
										 ref_label="label.main",
										 outd="out/cart/cell_id_by_singler",
										 debugd="out/cart/cell_id_by_singler/pmid_30726743.tmp")
}

message(str(args))

source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_workflow.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_seurat3.r'))
source('lib_project.r')

library(Seurat)
library(parallel)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(randomcoloR)
library(reshape2)
library(pheatmap)
library(ggplotify)

if (!file.exists(args$outd)) {dir.create(args$outd,showWarnings = FALSE, recursive = TRUE)}

if (args$debugd!="" & !file.exists(args$debugd)) {dir.create(args$debugd,showWarnings = FALSE, recursive = TRUE)}

if (args$sample_sheet=="" | !file.exists(args$sample_sheet)) {
	stop(sprintf("specify [%s]",args$sample_sheet))
} else {
	sdt <- as.data.table(read.xlsx(args$sample_sheet,sheet = args$sheet_number))
	is_all <- all(c("patient","timepoint","sample") %in% colnames(sdt))
	if (!is_all) {
		stop(sprintf("make sure the colummns:patient, timepoint, and sample exist in %s[%d]",args$sample_sheet,args$sheet_number))
	}
}

if (args$debug == 1) {browser()}
message(sprintf("loading seurat rd file [%s] ...",args$seu_rd))
load(args$seu_rd)

message(sprintf("loading singler ref rd file [%s] ...",args$singler_ref_rd))
load(args$singler_ref_rd)
if (args$cell2grep=="") {
	cell_ids_oi <- sort(unique(ref.sce[[args$ref_label]]))
} else {
	j <- grep(args$cell2grep,ref.sce[[args$ref_label]])
	cell_ids_oi <- sort(unique(ref.sce[[args$ref_label]][j]))
	ref.sce <- ref.sce[,j]
}

# --------------
message("preparing palette for cell annotation ...")
message(sprintf("cell_types:%s",paste0(cell_ids_oi,collapse = ",")))

n <- length(cell_ids_oi)
cellid_color <- data.table(ctype=cell_ids_oi,
													 pal=distinctColorPalette(n))
if (args$gois=="") {
	gois <- list(NA)
} else {
	gois <- comma_string_to_list(args$gois)
}

# -------------
# check if cell id is available
snames <- names(seus)

rds_fpath <- file.path(args$outd,"cellid_seus.rds")
if (args$reuse==1 & file.exists(rds_fpath)) {
	seus <- readRDS(rds_fpath)
} else {
	seus <- mclapply(snames,function(sname) {
		seu <- seus[[sname]]
		message(sname)
		seu <- cell_id_by_singler(seu,
															ref.sce,
															pmid=args$pmid,
															label2id=args$ref_label)
		
		message(sprintf("Done[%s]",sname))
		return(seu)
	}
	,mc.cores = args$ncpu
	)
	names(seus) <- snames
	saveRDS(seus,file=rds_fpath)
}

rds_fpath <- file.path(args$outd,"seus.rds")
if (args$reuse==1 & file.exists(rds_fpath)) {
	seus2 <- readRDS(rds_fpath)
} else {
	seus2 <- lapply(seus,function(slist){slist$seu})
	saveRDS(seus2,rds_fpath)
}

rds_fpath <- file.path(args$outd,"cell_ids.rds")
if (args$reuse==1 & file.exists(rds_fpath)){
	cell.ids <- readRDS(rds_fpath)
} else {
	cell.ids <- lapply(seus,function(slist){slist$pred.id})
	saveRDS(cell.ids,rds_fpath)
}

# -----------
cell_membership_dts <- lapply(seus2,function(seu) {
	dt2 <- as.data.table(seu@meta.data[,c(args$pmid_ref,'seurat_clusters','orig.ident')])
	dt2 <- cbind(dt2,sdt[match(dt2$orig.ident,sample),list(patient,timepoint)])
	return(dt2)
})
names(cell_membership_dts) <- names(seus2)

# -------------
wb <- createWorkbook()
message("cell count in seurat cluster by cell_id ...")
cellid_clusters <- lapply(names(cell_membership_dts),function(sname) {
	d <- acast(cell_membership_dts[[sname]],
						 seurat_clusters ~ get(args$pmid_ref), 
						 fun.aggregate = length, 
						 value.var = 'patient')
	message(sprintf("storing seurat_cluster x cell id count [%s]",sname))
	addWorksheet(wb,sheetName = sname)
	writeData(wb,sname,d,rowNames = TRUE)
	return(d)
})
saveWorkbook(wb,file=file.path(args$outd,sprintf('seurat_cluster_by_%s.xlsx',args$pmid_ref)),overwrite=T)
# ------------

wb <- createWorkbook()
message("cell count in timepoint by cell_id ...")
cell_membership_dt <- rbindlist(cell_membership_dts)
by_pats <- split(cell_membership_dt,by="patient")
cellid_tpoints <- lapply(names(by_pats),function(pat) {
	by_pat <- by_pats[[pat]]
	d <- acast(by_pat, timepoint ~ get(args$pmid_ref), fun.aggregate = length, value.var = 'patient')
	message(sprintf("storing seurat_cluster x cell id count [%s]",pat))
	addWorksheet(wb,sheetName = pat)
	writeData(wb,pat,d,rowNames = TRUE)
	return(d)
})
names(cellid_tpoints) <- names(by_pats)
saveWorkbook(wb,file=file.path(args$outd,sprintf('timepoint_by_%s.xlsx',args$pmid_ref)),overwrite=T)

wb <- createWorkbook()
tpoint_cellid <- acast(cell_membership_dt, timepoint ~ get(args$pmid_ref), fun.aggregate = length, value.var = 'timepoint')
addWorksheet(wb,sheetName = "tpoint_cellid")
writeData(wb,"tpoint_cellid",tpoint_cellid,rowNames=TRUE)

pat_by_cellid <- acast(cell_membership_dt, orig.ident ~ get(args$pmid_ref), fun.aggregate = length, value.var = 'orig.ident')
addWorksheet(wb,sheetName = "sname_cellid")
writeData(wb,"sname_cellid",pat_by_cellid,rowNames=TRUE)

norm.pat_by_cellid <- pat_by_cellid/rowSums2(pat_by_cellid)
norm.pat_by_cellid <- norm.pat_by_cellid[,sort(colnames(norm.pat_by_cellid))]
addWorksheet(wb,sheetName = "sname_cellid.norm")
writeData(wb,"sname_cellid.norm",norm.pat_by_cellid,rowNames=TRUE)
saveWorkbook(wb,file=file.path(args$outd,sprintf('all_by_%s.xlsx',args$pmid_ref)),overwrite=T)

# -------------
message("performing cell ID plot ...")
goi_per_samples <- list()
dge_dts <- list()
de_test <- "MAST"

message("change idents by args$pmid_ref ...")
snames <- names(seus2)
seus <- lapply(names(seus2),function(sname){
	seu <- seus2[[sname]]
	DefaultAssay(seu) <- "SCT"
	Idents(seu) <- args$pmid_ref
	return(seu)
})
names(seus) <- snames
rm(seus2)

reuse <- 0
cls_dge_rd <- file.path(args$debugd,"cls_dge.rd")
if (args$reuse==1 & args$debugd!="") {
	if (file.exists(cls_dge_rd)) {reuse <- 1}
}

if (reuse==1) {
	message(sprintf("reuse prev rd [%s] ...",cls_dge_rd))
	load(cls_dge_rd)
} else {
	cluster_marker_dges <- batch_solo_FindAllMarkers(seus,
																									 ident2=args$pmid_ref,
																									 assay2="SCT")
	names(cluster_marker_dges) <- names(seus)
	if (args$debugd!="") {
		save(cluster_marker_dges,file=cls_dge_rd,compress=TRUE)
	}
}

pdf_file <- file.path(args$outd,"dimplot_w_log2fc_heatmap.pdf")
if (args$debug==0) {pdf(pdf_file,width = 14,height = 6)}

for(sname in names(seus)) {
	message(sprintf("generating dimplot + heatmap for markers [%s]...",sname))
	seu <- seus[[sname]]
	cluster_marker_dge <- cluster_marker_dges[[sname]]
	# ------------
	message("setting colors ...")
	message(paste0(cellid_color$ctype,collpase=","))
	message(paste0(unique(seu@meta.data[[args$pmid_ref]]),collapse=","))
	
	cellid_color2 <- cellid_color[ctype %in% unique(seu@meta.data[[args$pmid_ref]]),]
	
	colv <- cellid_color2$pal
	names(colv) <- cellid_color2$ctype
	
	pum <- DimPlot(seu, 
								 reduction="umap",
								 cols = colv,
								 order = names(sort(norm.pat_by_cellid[sname,],decreasing = TRUE)),
								 pt.size=0.5,
								 combine=T)
	
	pum[[1]]$layers[[1]]$aes_params$alpha <- 0.2
	
	pum <- pum + ggtitle(sprintf("%s:%s",sname,args$pmid_ref))
	# ----------
	a1 <- as.data.table(cluster_marker_dge)
	goij <- NA
	phm <- ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
	if (nrow(a1)>0) {
		if (!is.na(gois[[1]])) {a1 <- a1[gene %in% gois,]}
		if (nrow(a1[p_val_adj<0.05])>1) {
			message(sprintf("print umap along with topk DGE heatmap[%s]",sname))
			
			a2 <- acast(a1[p_val_adj<0.05,], gene ~ cluster, value.var='avg_logFC')
			if (dim(a2)[1]==1 | dim(a2)[2]==1) {
				message(sprintf("either 1 feature or 1 cluster is observed and ignore heatmap plotting!"))
			} else {
				a3 <- select_topk_from_each_sample(a2,num_feats = round(60/dim(a2)[2]))
				
				a2 <- a2[rownames(a3),colnames(a3)]
				
				myPal <- get_red_white_blue_for_Heatmap(a2)
				a2[is.na(a2)] <- 0.
				
				cluster_info <- norm.pat_by_cellid[sname,colnames(a2)]
				
				ha1 <- rowAnnotation(frac=anno_barplot(cluster_info),
														 show_annotation_name = TRUE)
				
				phm <- as.ggplot(Heatmap(t(a2),
																 name="logFC",
																 col=myPal,
																 #row_names_side="left",
																 right_annotation = ha1,
																 show_row_names = TRUE,
																 column_title=sprintf("avg_logFC;padj<0.05")))
				
				goij <- rownames(a2)
			}
		} else {
			message(sprintf("No significant DGE observed in [%s]",sname))
		}
	} else {
		message(sprintf("No DGE observed in [%s]",sname))
	}
	p <- pum + phm + plot_layout(ncol = 2,widths = c(1, 2))
	plot(p)
	
	# ------------
	goi_per_samples[[sname]] <- goij
	dge_dts[[sname]] <- a1
	message("Done[dimplot].")
}
if (args$debug==0) {dev.off()}
message("Done.")

# ----------------

rd_file <- file.path(args$outd,"dge.rd")
message("Saving RD file ...")
save(cell_membership_dt,
		 dge_dts,
		 goi_per_samples,
		 file=rd_file,
		 compress=TRUE)