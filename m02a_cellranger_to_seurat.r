# author1: hongc2@ccf.org

library(data.table)
source('lib/lib_project.r')

if (T) {
  library(argparse)
  args_tmp <- commandArgs(trailingOnly = F)
  scriptPath <- dirname(normalizePath(sub("--file=","",args_tmp[grep("--file",args_tmp)])))
  
  message(sprintf("script_filepath [%s]",scriptPath))
  prog_ver <- "CAR-T:v0.1"
  parser <- ArgumentParser(description=sprintf('[%s]to read cellranger count filtered_feature_bc_matrix and generate seurat3 objects',prog_ver))
  
  parser$add_argument("-i", "--cellranger_countd", type="character", required=TRUE,
                      dest="cellranger_countd",
                      help="cellranger_count d")
  
  parser$add_argument("-x", "--sample_sheet", type="character", required=FALSE,
                      dest="sample_sheet", default="data/sample_sheets/samples_info.xlsx",
                      help="sample info excel/csv file where the columns contain sample,condition,Patient,Timepoint,type,comp,and sample_tag")

  parser$add_argument("-o", "--outd", type="character", required=FALSE,
                      dest="outd", default="out",
                      help="output directory")
  
  parser$add_argument("-s", "--soi", type="character", required=FALSE,
                      dest="soi", default="",
                      help="samples of interest, use a word w/o space that you want to select in sample_tag field.")
  

  parser$add_argument("-n", "--ncpu", type="integer", required=FALSE,
                      dest="ncpu", default=3,
                      help="ncpu? [3]")
  
  parser$add_argument("-d", "--debug", type="integer", required=FALSE,
                      dest="debug", default=0,
                      help="is debug mode? [0],1")
  
  args <- parser$parse_args()
} else {
  fpath_dt <- get_current_script_fpath()
  args <- data.table(cellranger_countd="results/cellranger",
                     pos_quantile=0.6,
                     outd=get_outd(get_proj_out0(fpath_dt$parentd),fpath_dt$fpath),
                     soi="",
                     ncpu=12,
                     debug=0)
}

message(str(args))

library(reticulate)
library(openxlsx)
library(Seurat)

if (!file.exists(args$outd)) {dir.create(args$outd,showWarnings = F,recursive = T)}

sample_dt <- as.data.table(read.xlsx(args$sample_sheet,sheet=1))

if (args$soi!="") {
  if (args$debug==1) {browser()}
  message(paste0(sample_dt$sample,collapse = ","))
  message(sprintf('selecting only %s samples ...',args$soi))
  sample_dt <- sample_dt[grep(args$soi,sample),]
  message("after selection ...")
  message(paste0(sample_dt$sample,collapse = ","))
}

sample_dt <- sample_dt[order(sample),]

reticulate::source_python(file.path(Sys.getenv('R_UTIL'),'lib_scrublet.py'))

sc_raws <- mclapply(1:nrow(sample_dt),function(i) {
  
  if (args$debug==1) {browser()}
  sample <- sample_dt[i,sample]
  message(sample)
  
  bc_matrix_dir <- file.path(args$cellranger_countd,sample,"outs","filtered_feature_bc_matrix")
  
  message(sprintf('start to read 10x count at %s',bc_matrix_dir))
  
  mtrx_gz <- file.path(bc_matrix_dir,'matrix.mtx.gz')
  sc_raw <- NA
  if (file.exists(mtrx_gz)) {
    # stopifnot(file.exists(bc_matrix_dir))
    
    sc.data <- Read10X(data.dir = bc_matrix_dir)
    only_ge <- FALSE
    if (is.null(names(sc.data))) {only_ge <- TRUE}
    if (only_ge) {
      sc_raw <- CreateSeuratObject(counts = sc.data,
                                   project = sample_dt[i,sample],
                                   min.cells = 0,
                                   min.features = 0)
    } else {
      sc_raw <- CreateSeuratObject(counts = sc.data$`Gene Expression`,
                                   project = sample_dt[i,sample],
                                   min.cells = 0,
                                   min.features = 0)
      
      if ('Antibody Capture' %in% names(sc.data)) {
        message("Storing Protein Feature Barcodes ...")
        sc_raw[['Protein']] = CreateAssayObject(counts = sc.data$`Antibody Capture`)
      }
      
      if ('Pathoscope2' %in% names(sc.data)) {
        message("Storing Pathoscope2 ...")
        sc_raw[['Pathoscope2']] = CreateAssayObject(counts = sc.data$`Pathoscope2`)
      }
      sc_raw <- NormalizeData(sc_raw)
      
      assay_names = names(sc.data)
      if ("Custom" %in% assay_names) {
        assay_names[assay_names=="Custom"] <- "HTO"
        names(sc.data) <- assay_names
        sc_raw[['HTO']] <- CreateAssayObject(counts = sc.data$HTO)
        
        sc_raw <- NormalizeData(sc_raw, assay = "HTO", normalization.method = "CLR")
      }
    }
    if (!is.na(sc_raw)) {
      sc_raw[["percent.mt"]] <- PercentageFeatureSet(sc_raw, pattern = "^MT-")
      sc_raw[["doublet_score"]] <- scrublet_annot_batch(bc_matrix_dir,args$debug)
      
    }
    message(sprintf("Done[%s]",bc_matrix_dir))
  } else {
    message(sprintf("the sample[%s] count matrix is not available yet!",sample))
  }
  if (args$debug==1) {browser()}
  return(sc_raw)
}
,mc.cores = args$ncpu
)
names(sc_raws) <- sample_dt$sample
sc_raws <- sc_raws[!is.na(sc_raws)]

rd_file <- file.path(args$outd,"raw_seurats.rd")
message(sprintf("saving RD file [%s] ...",rd_file))
save(sc_raws,file=rd_file,compress = TRUE)

message("reporting cell counts for each sample ...")
if (args$debug == 1) {browser()}
report_sc3_dim(sc_raws,file.path(args$outd,"cellcounts_dim.xlsx"))
