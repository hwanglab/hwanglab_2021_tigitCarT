## R version used in the manuscript
- R 3.5.1

## Prepapre computational environment and your bash configuration
- Prepare a project directory to work on. For example,
	- `mkdir -p ~/projects`
- Clone the paper reproducible code from this github
	- `cd ~/projects`
	- `git clone https://github.com/hwanglab/hwanglab_2021_tigitCarT`
	
- Add the environment variable to your bash shell script (`~/.bash_profile` or `~/.profile`)
	- `export R_UTIL=~/projects/hwanglab_2021_tigitCarT/lib`
	- relogin to apply the environmen variable to your bash environment

- Setup some required R library modules
	- `Rscript ./m01a_R_module_setup.r`

- Data analysis step
	- Download FASTQ files from EGA and save the files into `hwanglab_2021_tigitCarT/data/fastq`


## Workflow
### Option1: FASTQ to cellranger read count matrix
- Perform cellranger for 1st seqrun. For example,
	- ```cellranger count --id P3D0.c1 --transcriptome=${REFD}/10x_genomics/GRCh38 --fastqs=data/fastq/D1_1,data/fastq/D1_2 --nosecondary --localcores=20 --localmem=196```

- Perform cellranger with feature barcodes for 2nd seqrun. For example,
	- ```cellranger count --id P10D0.c2 --transcriptome=${REFD}/10x_genomics/GRCh38 --libraries=data/sample_sheets/P10_day0_CART.csv --feature-ref=data/sample_sheets/srun2_fr_totalseqB.csv --localcores=12 --localmem=96```

- convert the cellranger read count matrices to Seurat object files: `m02a_cellranger_to_seurat.r`

- perform QC and retain a high quality read count matrix: `m03a_qc_cells.r`

- run SingleR on each sample using the reference database (GSE107011): `m04a_cell_id_by_single.r`

### Option2: Starting from the preprocessed Seurat object files directly
- download the QC'ed Seurat object file (`seus_20210202.rds`) at `https://m-17902d.f0ce50.e229.dn.glob.us/projects/pub_data/2021a.CART_TIGIT/seus_20210202.rds.gpg` and copy the file into `out/cart` directory, i.e (`out/cart/seus_20210202.rds`)


### Data integration to all downstream analyses
- `m05a_harmony_findmarkers.r`
	- input: `out/cart/seus_20210202.rds`
	- process: to integrate each patient CART samples; to perform batch effect removal; to align cells; to generate cluster marker genes in DEG table format
	- output: `out/cart/harmony_findmarkers/{seui.rds,cls_marker.tsv,cls_deg.rds}`

### Heatmap of marker genes and tSNE plots with cell annotation
- `m06a_vcluster_heatmap.r`
	- input: `out/cart/harmony_findmarkers/{seui.rds,cls_deg.rds}`
	- process: select top-k marker genes in each cluster and generate a heatmap
	- output: `out/cart/vcluster_heatmap/*.pdf` Fig1.B, Fig1.C, SFig1.A, SFig1.B, SFig1.C, and SFig3.A(left), Fig2.E(upper), Fig2.F, Fig3.B

### CD4/CD8 T cell profile in spider web
- `m07a_celltype_spiderweb.r`
	- output: Fig2.E (bottom)

### TIGIT expression comparison before and after CART infusion
- `m08a_tigit_expr_by_tgroup.r`
	- output: Fig3.G

### CART manufacture vs. post infusion comparison
- `m09_deg_pre_vs_post_inf.r`
	- output: Fig2.C, Fig2.G, Fig2.B, Fig3.C

### DEG heatmap from resp vs. nonresp patient group
- `m10_deg_by_cluster_inter.r`
	- output: Fig3.D

### Immune checkpoint expression change per patient
- `m11_violin_ridge_gexpr.r`
	- output: Fig2.D, Fig3.F

### TIGIT expression and compared by the treatment outcome
- `m12_dyfscore_TIGIT_on_tsne.r`
	- output: Fig3.E, Fig4.B(left), Fig4.B(right), Fig4.C(left), Fig4.C(right), SFig3.B(right)
	
### DEG between TIGIT-on cells and TIGIT-off cells
- `m13_tigitON_pre_vs_post.r`
	- output: Fig4.A