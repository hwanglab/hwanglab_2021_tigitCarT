# combine all QC'ed single sample Seurat object and generate one batch-removal or aligned Seurat object file
#1. load shared library
source('lib/lib_project.r')

# --------------------------
fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 seus_rds="data/jackson_hong_2021CART.rds",
									 ncpu=12,
									 outd=get_wkd(get_projd0(),fpath_dt$fpath),
									 reuse=1,
									 debug=0)

seus <- readRDS(args$seus_rds)

harmony_findmarkers(seus,wkd=args$outd,reuse=args$reuse)
