# author1: hongc2@ccf.org

source('lib/lib_project.r')

if (T) {
	library(argparse)
	args_tmp <- commandArgs(trailingOnly = F)
	scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
	
	message(sprintf("script_filepath [%s]",scriptPath))
	
	prog_ver <- "CAR-T:v0.1"
	parser <- ArgumentParser(description=sprintf('[%s]QC control: to obtain a higher quality of seurat objects',prog_ver))
	
	parser$add_argument("-i", "--raw_seu_rd", type="character", required=TRUE,
											dest="raw_seu_rd",
											help="raw_seu_rd file path from seurat_from_10x.r")
	parser$add_argument("-o", "--outd", type="character", required=FALSE,
											dest="outd", default="out",
											help="output directory")
	parser$add_argument("-f", "--sd_filtx", type="double", required=FALSE,
											dest="sd_filtx", default=2.5,
											help="sd_filtx @ mean +/- sd_filtx * sd [2.5]")
	parser$add_argument("-m", "--mt_pct", type="integer", required=FALSE,
											dest="mt_pct", default=-1,
											help="max mitochondrial pct [-1]:mean + 1*sd")
	parser$add_argument("-D", "--doubletr", type="double", required=FALSE,
											dest="doubletr", default=0.,
											help="doublet rate [0.]:auto detect")
	parser$add_argument("-C", "--ncell_for_active_gene", type="integer", required=FALSE,
											dest="ncell_for_active_gene", default=0,
											help="ncell_for_active_gene [0]")
	parser$add_argument("-d", "--debug", type="integer", required=FALSE,
											dest="debug", default=0,
											help="is debug mode? [0],1")
	
	args <- parser$parse_args()
} else {
	args <- data.table(raw_seu_rd="results/raw_seus.rds",
										 outd="results/qc_cells",
										 sd_filtx=2.5,
										 mt_pct=15.,
										 doubletr=0.,
										 ncell_for_active_gene=0,
										 debug=1)
}

#use_condaenv(condaenv="Renv", conda="~/miniconda3/bin/conda")

library(Seurat)
library(openxlsx)
library(ggplot2)
library(parallel)
library(reticulate)

if (!file.exists(args$outd)) {dir.create(args$outd,showWarnings = FALSE, recursive = TRUE)}

# ---------------------
message(sprintf("loading seurat objects from [%s] ...",args$raw_seu_rds))
sc_raws <- read_rd_or_rds(args$raw_seu_rds) #sc_raws
sc_raws <- sc_raws[sort(names(sc_raws))]
# --------------------
message("seurat object metrics before QC ...")
pdf_file <- get_out_fpath(args,"summary_beforeQC.pdf")
message(pdf_file)
print_qc_figures(sc_raws,pdf_file,debug2=0)

dim_report_tsv <- get_out_fpath(args,"cellcounts_dim_QC.tsv")
message(dim_report_tsv)
report_sc3_dim(sc_raws,tsv_fpath=dim_report_tsv)

sample_meta_dt = load_sample_worksheet(sheet_name="patient_meta",debug2=0)

cell_cnt1 <- cell_count_w_meta(sc_raws,sample_meta = sample_meta_dt,qc_comment="N",debug2=0)

# ---------------------
message("gathering single cell read/gene metrics ...")

seus <- mclapply(sc_raws,function(seu) {
	filter_by_reads_features_cnt2(seu,
																sd_filtx=args$sd_filtx,
																mt_pct=args$mt_pct,
																doubletr=args$doubletr,
																debug2=args$debug)
}
,mc.cores = args$ncpu
)
names(seus) <- names(sc_raws)
rm(sc_raws)

# ---------------------
qc_comment <- sprintf("|sd(nC,nF)|<%g,mt<%g%%",args$sd_filtx,args$mt_pct)

cell_cnt2 <- cell_count_w_meta(seus,sample_meta = sample_meta_dt, qc_comment=qc_comment)

cell_cnt <- rbind(cell_cnt1,cell_cnt2)

p <- ggplot(cell_cnt,aes(x=patient,y=cnt,group=qc,fill=qc)) +
	geom_bar(stat="identity",position=position_dodge()) +
	facet_grid(rows=vars(timepoint)) +
	theme_bw()+
	theme(legend.position="bottom")

pdf_file <- get_out_fpath(args,"cell_cnt_after_qc.pdf")
message(pdf_file)
pdf(file=pdf_file,width=6,height = 8)
plot(p)
dev.off()

rd_file <- get_out_fpath(args,"seus.rd")
message(sprintf("saving rd file[%s] ...",rd_file))
save(seus,file=rd_file,compress=TRUE)

message("seurat object metrics after QC ...")

pdf_file <- get_out_fpath(args,"summary_afterQC.pdf")
message(pdf_file)

print_qc_figures(seus,pdf_file,qc_comment=qc_comment,debug2=0)
message(pdf_file)

report_sc3_dim(seus,tsv_fpath=dim_report_tsv,qc_comment=qc_comment,append2=TRUE,debug2=0)

