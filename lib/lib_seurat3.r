# author: hongc2@ccf.org
# input: cellranger output directory
# objective: to compare groups of scRNA
# ref: 
# https://satijalab.org/seurat/v3.0/immune_alignment.html
# https://www.biorxiv.org/content/biorxiv/early/2018/11/02/460147.full.pdf

library(Seurat)
library(data.table)
library(openxlsx)
library(reshape2)
library(dplyr)
library(reticulate)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)
library(parallel)
# library(DoubletFinder)
library(R.utils)
library(ggplot2)
#library(DEsingle)
library(BiocParallel)
library(matrixStats)
library(pheatmap)
library(SingleCellExperiment)
library(SingleR)
library(harmony) #install_github("immunogenomics/harmony")
library(SeuratWrappers) #remotes::install_github('satijalab/seurat-wrappers')
library(monocle3) #devtools::install_github('cole-trapnell-lab/leidenbase');devtools::install_github('cole-trapnell-lab/monocle3');https://cole-trapnell-lab.github.io/monocle3/docs/installation/

# use_condaenv(condaenv="scell",conda="/media/sammy/apps/miniconda3/bin/conda")

#' Convert 10x cellranger count output to seurat3 object
#' more detail
#' 
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'    
#' @param sample_dt A data frame or data.table in the column names (sample,cellranger_count_outd,sgroup). For example,
#'  D1_1, /home/blahblah/1907UNHS-0026_cellranger_count/D1_1, pre_infusion
#' @param min_cell The minimum number of cells to include.
#' @param min_features The min number of features to include.
#' @param ncpu The number of cpus to utilize
#' @param debug 0 (no debug), 1 (debug; test with the first sample only)
#' @return seurat3 object in a list.
#' 
sc10x_to_seurat3 <- function(sample_dt,min_cell=3,min_features=200,ncpu=2,debug=0,mt_pct=10.) {
  
  #make sure that sample_dt has 3 columns
  stopifnot(dim(sample_dt)[2]==3)
  if (debug==1){browser()}
  S <- dim(sample_dt)[1]
  message(sprintf("total number of samples [%d]",S))
  
  if (debug==1) {S <- 1}
  
  seurat_raws <- mclapply(1:S,function(i) {
    if (debug==1){browser()}
    sample <- sample_dt[i,sample]
    crcount_out_fpath <- sample_dt[i,cellranger_count_outd]
    message(sample)
    
    bc_matrix_dir <- sprintf("%s/outs/filtered_feature_bc_matrix",crcount_out_fpath)
    
    message(sprintf('start to read 10x count at [%s]',bc_matrix_dir))
    
    sc.data <- Read10X(data.dir = bc_matrix_dir)
    
    sc_raw <- CreateSeuratObject(counts = sc.data,
                                 project = sample,
                                 min.cells = min_cell,
                                 min.features = min_features)
    
    sc_raw[["percent.mt"]] <- PercentageFeatureSet(sc_raw, pattern = "^MT-")
    sc_raw$sgroup <- sample_dt[i,sgroup]
    if (debug==1){browser()}
    sc_raw <- subset(sc_raw,subset = percent.mt < 10.)
    
    message(sprintf("Done[%s]",sample))
    return(sc_raw)
  },mc.cores = ncpu) #,mc.cores = ncpu
  
  names(seurat_raws) <- sample_dt$sample
  if (debug==1){browser()}
  return(seurat_raws)
}

print_qc_figures <- function(sc_raws,pdf_file,qc_comment="N",debug2=0) {
  
  graphics.off()
  if (debug2==1) {browser()
  } else { pdf(pdf_file) }

  # ---------------------
  message("gathering single cell read/gene metrics ...")
  sc_metrics <- lapply(sc_raws,function(seu) {
    
    sc_metric <- as.data.table(seu@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","doublet_score")])
    return(sc_metric)
  })
  sc_metrics <- rbindlist(sc_metrics)
  qc_tsv <- sprintf("%s_QC.tsv",pdf_file)
  fwrite(sc_metrics[,.(nCount=mean(nCount_RNA),nFeature=mean(nFeature_RNA),mt_pct=mean(percent.mt),double_score=mean(doublet_score)),by="orig.ident"],file=qc_tsv,sep="\t")
  # ---------------------
  message("violin plot (read count) ...")
  p <- ggplot(sc_metrics, aes(x = orig.ident,
                              y = log10(nCount_RNA))) + 
    geom_violin(trim=TRUE)
  
  p <- p + stat_summary(fun.data=mean_sdl,
                        geom="pointrange", color="red") +
    labs(title=sprintf("nCount_RNA[QC:%s]",qc_comment),
         x="samples",
         y="log10(read count)")+
    coord_flip()
  
  plot(p)
  
  # ---------------------
  message("violin plot (gene) ...")
  p <- ggplot(sc_metrics, aes(x = orig.ident,
                              y = log10(nFeature_RNA))) + 
    geom_violin(trim=TRUE)
  
  p <- p + stat_summary(fun.data=mean_sdl,
                        geom="pointrange", color="red") +
    labs(title=sprintf("nFeature_RNA[QC:%s]",qc_comment),
         x="samples",
         y="log10(features)")+
    coord_flip()
  
  plot(p)
  
  # ---------------------
  message("violin plot (mt.pct) ...")
  p <- ggplot(sc_metrics, aes(x = orig.ident,
                              y = percent.mt)) + 
    geom_violin(trim=TRUE) +
    ylim(0,max(sc_metrics$percent.mt))
  
  p <- p + stat_summary(fun.data=mean_sdl,
                        geom="pointrange", color="red") +
    labs(title=sprintf("percent.mt[QC:%s]",qc_comment),
         x="samples",
         y="mitrochondria [%]") +
    coord_flip()
  
  plot(p)
  
  # ---------------------
  message("violin plot (doublet score) ...")
  doublet_avail <- 0
  if (sum(sc_metrics$doublet_score)>0) {
    doublet_avail <- 1
    p <- ggplot(sc_metrics, aes(x = orig.ident,
                                y = doublet_score)) + 
      geom_violin(trim=TRUE)
    
    p <- p + stat_summary(fun.data=mean_sdl,
                          geom="pointrange", color="red") +
      labs(title=sprintf("doublet score by scrublet[QC:%s]",qc_comment),
           x="samples",
           y="doublet score") +
      coord_flip()
    
    plot(p)
  }
  # -------------------------
  message("scatter plots between read count vs. mitochondira counts ...")
  cnt_vs_mito <- list()
  for(sample in names(sc_raws)) {
    p <- FeatureScatter(sc_raws[[sample]],
                        feature1="nCount_RNA",
                        feature2="percent.mt")
    cnt_vs_mito[[sample]] <- p +
      theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
      theme(legend.position='none') +
      ggtitle(sample)
  }
  cnt_vs_mitos <- chunk2(cnt_vs_mito,4)
  if (debug2==1) {browser()}
  ret <- lapply(cnt_vs_mitos,function(cnt_vs_mito){
    p <- wrap_plots(cnt_vs_mito,ncol=round(sqrt(length(cnt_vs_mito))))
    plot(p)
    NA
  })
  
  # -------------------------
  message("scatter plots between read count vs. feature counts ...")
  if (debug2==1) {browser()}
  cnt_vs_feature <- list()
  for(sample in names(sc_raws)) {
    p <- FeatureScatter(sc_raws[[sample]],
                        feature1="nCount_RNA",
                        feature2="nFeature_RNA")
    
    cnt_vs_feature[[sample]] <- p + 
      theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
      theme(legend.position='none') +
      ggtitle(sample)
  }
  
  cnt_vs_features <- chunk2(cnt_vs_feature,4)
  ret <- lapply(cnt_vs_features,function(cnt_vs_feature){
    p <- wrap_plots(cnt_vs_feature,,ncol=round(sqrt(length(cnt_vs_feature))))
    plot(p)
    NA
  })
  
  # -------------------------
  message("scatter plots between read count vs. doublet score ...")
  if (doublet_avail==1) {
    cnt_vs_feature <- list()
    for(sample in names(sc_raws)) {
  
      p <- FeatureScatter(sc_raws[[sample]],
                          feature1="nCount_RNA",
                          feature2="doublet_score")
      
      cnt_vs_feature[[sample]] <- p + 
        theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
        theme(legend.position='none') +
        ggtitle(sample)
      
    }
    
    cnt_vs_features <- chunk2(cnt_vs_feature,4)
    ret <- lapply(cnt_vs_features,function(cnt_vs_feature){
      p <- wrap_plots(cnt_vs_feature,,ncol=round(sqrt(length(cnt_vs_feature))))
      plot(p)
      NA
    })
  }
  if (debug2==1) {browser()
  } else {dev.off()}
}

filter_by_reads_features_cnt2 <- function(seu,
                                         sd_filtx=2.5,
                                         nc_cutoffs=NA,
                                         ncell_for_active_gene=0,
                                         mt_pct=15,
                                         max_mt_pct=25,
                                         min_ncount=500,
                                         min_nfeature=100,
                                         doubletr=-1.,
                                         nfeature_to_ncount_coff=0.1,
                                         debug2=0) {
  
  sc_metric <- as.data.table(seu@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","doublet_score")])
  
  if (debug2==1) {browser()}
  
  if (doubletr == 0.) {
    ncells <- dim(sc_metric)[1]
    doubletr <- get_doublet_rate_in_10x(ncells)
  }
  
  pre_nc_cutoff <- 1
  if (is.na(nc_cutoffs)) {
    pre_nc_cutoff <- 0
    nCount_sdv <-sd(seu$nCount_RNA)
    nCount_meanv <- mean(seu$nCount_RNA)
    
    nFeature_sdv <-sd(seu$nFeature_RNA)
    nFeature_meanv <- mean(seu$nFeature_RNA)
    
    nc_cutoffs <- c(nCount_meanv-sd_filtx*nCount_sdv,nCount_meanv+sd_filtx*nCount_sdv)
    if (nc_cutoffs[[1]]<min_ncount){nc_cutoffs[[1]] <- min_ncount}
  
    nf_cutoffs <- c(nFeature_meanv-sd_filtx*nFeature_sdv,nFeature_meanv+sd_filtx*nFeature_sdv)
    if (nf_cutoffs[[1]]<min_nfeature){nf_cutoffs[[1]] <- min_nfeature}
    
  } else {
    nf_cutoffs <- NA
  }
  
  if (mt_pct < 0) {
    mt_meanv <- mean(seu$percent.mt)
    mt_sdv <- sd(seu$percent.mt)
    mt_pct <- min(c(mt_meanv + mt_sdv,max_mt_pct))
  }
  
  dr_cutoff2 <- 1.
  if (doubletr != -1) {
    dr_cutoff2 <- quantile(seu$doublet_score, c(1.-doubletr))
  }
  
  if (nfeature_to_ncount_coff>0) {
    smeta = seu@meta.data
    nfeature_to_ncount <- smeta$nFeature_RNA/(smeta$nCount_RNA+1)
    smeta$feat_complexity <- "high"
    smeta$feat_complexity[nfeature_to_ncount<nfeature_to_ncount_coff] <- "low"
    seu$feat_complexity <- smeta$feat_complexity
    seu <- subset(seu,feat_complexity=="high")
    seu$feat_complexity <- NULL
    rm(smeta)
  }
  
  cell_ids <- rownames(seu@meta.data)

  if (pre_nc_cutoff==0) {
    celloi <- cell_ids[seu@meta.data$nCount_RNA>=nc_cutoffs[[1]] & 
                         seu@meta.data$nCount_RNA<=nc_cutoffs[[2]] & 
                         seu@meta.data$nFeature_RNA>=nf_cutoffs[[1]] & 
                         seu@meta.data$nFeature_RNA<=nf_cutoffs[[2]] &
                         seu@meta.data$percent.mt<mt_pct &
                         seu@meta.data$doublet_score<dr_cutoff2]
  } else {
    celloi <- cell_ids[seu@meta.data$nCount_RNA>=nc_cutoffs[[1]] & 
                         seu@meta.data$nCount_RNA<=nc_cutoffs[[2]] & 
                         seu@meta.data$percent.mt<mt_pct &
                         seu@meta.data$doublet_score<dr_cutoff2]
  }
  
  nc_genes <- Matrix::rowSums(seu@assays$RNA@counts>0)
  genes <- names(nc_genes)
  active_features <- genes[nc_genes>ncell_for_active_gene]
  
  seu <- SubsetData(seu,cells=celloi,features=active_features)
  
  seu <- NormalizeData(seu)
  
  if (debug2==1) {browser()}
  return(seu)
}

filter_by_reads_features_cnt <- function(seu,
                                         nc_lor=0.02,
                                         nc_hir=0.98,
                                         nf_lor=0.02,
                                         nf_hir=0.98,
                                         ncell_for_active_gene=0,
                                         mt_pct=10.,
                                         doubletr=-1.,
                                         debug2=0) {
  
  sc_metric <- as.data.table(seu@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","doublet_score")])
  if (debug2==1) {browser()}
  # browser()
  
  if (doubletr == 0.) {
    ncells <- dim(sc_metric)[1]
    doubletr <- get_doublet_rate_in_10x(ncells)
  }
  
  nc_cutoffs <- quantile(seu$nCount_RNA, c(nc_lor,nc_hir))
  nf_cutoffs <- quantile(seu$nFeature_RNA, c(nf_lor,nf_hir))
  dr_cutoffs <- quantile(seu$doublet_score, c(1.-doubletr))
  
  cell_ids <- rownames(seu@meta.data)
  
  celloi <- cell_ids[seu@meta.data$nCount_RNA>nc_cutoffs[[1]] & 
                       seu@meta.data$nCount_RNA<nc_cutoffs[[2]] & 
                       seu@meta.data$nFeature_RNA>nf_cutoffs[[1]] & 
                       seu@meta.data$nFeature_RNA<nf_cutoffs[[2]] &
                       seu@meta.data$percent.mt<mt_pct &
                       seu@meta.data$doublet_score<dr_cutoffs[[1]]]
  
  nc_genes <- Matrix::rowSums(seu@assays$RNA@counts>0)
  genes <- names(nc_genes)
  active_features <- genes[nc_genes>ncell_for_active_gene]
  
  seu <- subset(seu,cells=celloi,features=active_features)
  if (debug2==1) {browser()}
  return(seu)
}

get_k_most_distinct_genes_per_cluster <- function(sc3.markers.dt,
                                                  top_k_pos=10,
                                                  max_p_val_adj=0.001) {
  
  sc3.markers.dt2 <- split(sc3.markers.dt,by="cluster")
  
  top_most_vars <- lapply(sc3.markers.dt2,function(sc3m) {
    sc3m <- sc3m[p_val_adj<max_p_val_adj,]
    if (dim(sc3m)[1]<top_k_pos){
      sc3m[order(-avg_logFC)]
    } else {
      sc3m[order(-avg_logFC)][1:top_k_pos,]
    }
  })
  
  top_vars <- rbindlist(top_most_vars)
  
  # 
  # if (top_k_neg>0) {
  # 	top_most_vars <- lapply(sc3.markers.dt2,function(sc3m) {
  # 		sc3m <- sc3m[p_val_adj<max_p_val_adj,]
  # 		if (dim(sc3m)[1]<top_k_pos){
  # 			sc3m[order(avg_logFC)]
  # 		} else {
  # 			sc3m[order(avg_logFC)][1:top_k_pos,]
  # 		}
  # 	})
  # 	top_neg_vars <- rbindlist(top_most_vars)
  # 	top_vars <- distinct(rbindlist(list(top_vars,top_neg_vars)))
  # }
  
  
  return(top_vars)
}

find_clusters <- function(seu,genes=NA){
  seu <- NormalizeData(seu)
  if (is.na(genes)){
    seu <- FindVariableFeatures(seu, selection.method = "vst")
  } else {
    VariableFeatures(seu) <- genes
  }
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:20)
  seu
}

solo_analysis <- function(wkd,sc3,sname,max_cluster=12,vst_feature=2000) {
  
  fig_prefix <- "01_feature_readCnt_mt"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p<-VlnPlot(sc3, pt.size=0.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(p)
  dev.off()
  
  # ============
  
  fig_prefix <- "02_readCnt_vs_mt_feature"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  plot1 <- FeatureScatter(sc3, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sc3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p<-CombinePlots(plots = list(plot1, plot2))
  plot(p)
  dev.off()
  
  fig_prefix <- "03_readCnt_features_hist"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  plot1 <- ggplot(sc3@meta.data,aes(x=nCount_RNA))+geom_histogram(color="black",fill="white")
  plot2 <- ggplot(sc3@meta.data,aes(x=nFeature_RNA))+geom_histogram(color="black",fill="white")
  p <- CombinePlots(plots = list(plot1,plot2))
  plot(p)
  dev.off()
  
  # ..............................................
  # find variable features if necessary. Otherwise, use default
  sc3 <- FindVariableFeatures(sc3, selection.method = "vst", nfeatures = vst_feature)
  
  
  # Scaling the data across cells before dim reduction
  all.genes <- rownames(sc3)
  sc3 <- ScaleData(sc3, features = all.genes)
  
  # Perform PCA
  num_cells <- dim(sc3@assays$RNA@data)[2]
  if (num_cells < 200) {
    npcs <- floor(num_cells/4)
  } else {
    npcs <- 50
  }
  
  sc3 <- RunPCA(sc3, features = VariableFeatures(object = sc3), npcs = npcs)
  
  # ..............................................
  VizDimLoadings(sc3, dims = 1:2, reduction = "pca")
  
  fig_prefix <- "04_pca"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DimPlot(sc3, reduction = "pca")
  plot(p)
  dev.off()
  
  fig_prefix <- "05_pca_heatmap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  DimHeatmap(sc3, dims = 1:6, cells = 500, balanced = TRUE)
  dev.off()
  
  # ..............................................
  
  if (max_cluster==0) {
    # Determine the ‘dimensionality’ of the dataset to cluster cells
    sc3 <- JackStraw(sc3, num.replicate = 100)
    sc3 <- ScoreJackStraw(sc3, dims = 1:20)
    L <- length(sc3@reductions$pca@stdev)
    D <- which(sc3@reductions$pca@stdev < sc3@reductions$pca@stdev[L]*1.2)[1]
  } else {
    # ..............................................
    D <- max_cluster
  }
  
  if (D > npcs) {
    D <- npcs
  }
  
  message(sprintf("the number of max clusters [%d]",D))
  #Cluster the cells
  sc3 <- FindNeighbors(sc3, dims = 1:D) #KNN-graph
  sc3 <- FindClusters(sc3, resolution = 0.5) #Louvain algorithm (granularity: 0.4 ~ 1.2 for 3K cells)
  
  # Non-linear dim reduction for visualization
  # for umap installation
  # reticulate::py_install(packages ='umap-learn')
  
  sc3 <- RunUMAP(sc3, dims = 1:D, umap.method = "umap-learn", metric = "correlation")
  
  fig_prefix <- "06_umap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DimPlot(sc3,reduction="umap")
  plot(p)
  dev.off()
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sc3.markers <- FindAllMarkers(sc3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # sc3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  sc3.markers.dt <- as.data.table(sc3.markers)
  
  top_most_vars <- get_k_most_distinct_genes_per_cluster(sc3.markers.dt,top_k_pos=10)
  
  fig_prefix <- "07_most_vgenes_heatmap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DoHeatmap(sc3, features = top_most_vars$gene) + NoLegend()
  plot(p)
  dev.off()
  
  tab_prefix <- "08_distinct_expr_genes_per_cluster"
  wb <- createWorkbook()
  sheet_name <- sprintf('%s_top%d_pos_expr_genes',sname,10)
  message(sprintf("sheet_name:%s",sheet_name))
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,top_most_vars)
  
  top_most_vars <- get_k_most_distinct_genes_per_cluster(sc3.markers.dt,top_k_pos=1000)
  sheet_name <- sprintf('%s_top%d_pos_expr_genes',sname,1000)
  message(sprintf("sheet_name:%s",sheet_name))
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,top_most_vars)
  
  saveWorkbook(wb,file=file.path(wkd,sprintf('%s_%s.xlsx',sname,tab_prefix)),overwrite=T)
  
  return(sc3)
}

sc3_normalization <- function(sc_raws,norm_method="default",nfeature_rna=300,mt_pct=10,vst_feature=0,norm_scale_factor=10000,ncpu=8) {
  
  sc_norms <- mclapply(names(sc_raws),function(name) {
    
    sc <- sc_raws[[name]]
    if (FALSE) {
      sinfo <- samples[sample == name,]
      sc$sample <- sinfo$sample
      sc$patient <- sinfo$Patient
      sc$timepoint <- sinfo$Timepoint
      sc$type <- sinfo$type
    }
    message(sprintf('computing mitochondria pct [%s] ...',name))
    sc <- PercentageFeatureSet(sc,pattern="^MT-",col.name = "percent.mt")
    
    if (FALSE) {
      fig_tag <- "mitochondria"
      pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',fig_tag,name))
      message(pdf_file)
      pdf(pdf_file)
      hist(sc$percent.mt,xlab = "mitochondria[%]",main = name,breaks = 50)
      dev.off()
    }
    
    message("filtering cells containing high mt ...")
    sc <- subset(sc, cells = which(sc$nFeature_RNA > nfeature_rna & sc$percent.mt < mt_pct))
    
    if (norm_method == "sctf") {
      message(sprintf("norm_method[%s]",norm_method))
      sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)
    } else {
      message(sprintf("norm_method[%s]",norm_method))
      message(sprintf('perform a normalization [%s] ...',name))
      sc <- NormalizeData(sc, scale.factor = norm_scale_factor, verbose = FALSE)
    }
    
    message(sprintf('find variable features [%s] ...',name))
    if (vst_feature>0) {
      sc <- FindVariableFeatures(sc, 
                                 selection.method = "vst", 
                                 nfeatures = vst_feature,
                                 verbose = FALSE)
    }
    message(sprintf('Done[%s]',name))
    return(sc)
  },mc.cores = ncpu) #,mc.cores = ncpu
  
  names(sc_norms) <- names(sc_raws)
  return(sc_norms)
}


runFIA_with_tryCatch <- function(sc.list,sc.features,D,k_filter=200) {
  sc.anchors <- tryCatch(
    {
      sc.anchors <- FindIntegrationAnchors(object.list = sc.list, 
                                           normalization.method = "SCT", 
                                           anchor.features = sc.features,
                                           dims = 1:D,
                                           k.filter = k_filter)
    },
    error=function(cond){
      message(sprintf('k_filter[%d]',k_filter))
      return(NA)
    },
    finally={
      message(sprintf('Done[%d]',k_filter))
    }
  )
  return(sc.anchors)
}

integration_analysis <- function(sc_samples,min.cell=50,num_features=3000,debug2=0) {
  
  if (debug2==1) {browser()}
  
  feat_cnts <- sapply(sc_samples,function(seu){return(dim(seu)[1])})
  mean_feat <- mean(feat_cnts)
  D <- 30
  pca_approx <- TRUE
  if (D > mean_feat) {
    D <- round(mean_feat * 0.7)
    pca_approx <- FALSE
  }
  if (num_features > mean_feat) {num_features <- mean_feat}
  
  if (num_features>0) {
    message(sprintf('Find [%d] variable features ...',num_features))
    sc_samples <- lapply(sc_samples,function(sc_sample){
      FindVariableFeatures(sc_sample, selection.method = "vst", nfeatures = num_features)
    })
  }
  
  message(sprintf('integrate the anchors found to each sample ([%d] anchor features)...',num_features))
  sc.anchors <- FindIntegrationAnchors(object.list = sc_samples, dims = 1:D, anchor.features = num_features)
  rm(sc_samples)
  sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:D)
  rm(sc.anchors)
  # -----------------
  DefaultAssay(sc.integrated) <- "integrated"
  
  message('cell counts per sample')
  summary(sc.integrated@active.ident)
  
  message('data matrix dimension')
  dim(sc.integrated@assays$integrated@data)
  
  # Run the standard workflow for visualization and clustering
  sc.integrated <- ScaleData(sc.integrated, verbose = TRUE)
  sc.integrated <- RunPCA(sc.integrated, npcs = D, verbose = TRUE, approx=pca_approx)
  
  # UMAP and Clustering
  sc.integrated <- RunUMAP(sc.integrated, reduction = "pca", dims = 1:D)
  
  sc.integrated <- FindNeighbors(sc.integrated, reduction = "pca", dims = 1:D)
  sc.integrated <- FindClusters(sc.integrated, resolution = 0.5)
  
  return(sc.integrated)
}

rm_assay <- function(seus,assay2default="RNA",assay2del="SCT") {
  
  message(sprintf("removing assay [%s] ...",assay2del))
  snames<-names(seus)
  seus <- lapply(seus, function(seu){
    avail_assays <- names(seu@assays)
    stopifnot(assay2default %in% avail_assays)
    DefaultAssay(seu) <- assay2default
    if (assay2del %in% avail_assays) {
      seu[[assay2del]] <- NULL
    }
    return(seu)
  })
  names(seus) <- snames
  return(seus)
}

merge_seurats <- function(seus,debug2=0) {
  if (debug2==1){browser()}
  S <- length(seus)
  stopifnot(S>1)
  
  message("appending an index to distinguish cbc ...")
  snames<-names(seus)
  for (s in 1:S) {
    seus[[s]] <- RenameCells(seus[[s]],add.cell.id = s)
  }
  names(seus) <- snames
  
  message(sprintf("merging [%d] seurat objects ...",S))
  seui <- seus[[1]]
  for (s in 2:S) {
    message(s)
    seui <- merge(seui,seus[[s]])
  }
  return(seui)
}

aligned_by_harmony0 <- function(seus,D=30,min.cells=50,debug2=0,reused="") {
  if (debug2==1){browser()}
  rd_file<-file.path(reused,"regress_out_scaled.rd")
  if (file.exists(rd_file)) {
    message(sprintf("reuse [%s] ...",rd_file))
    load(rd_file)
  } else {
    message("running harmony for data integration ...")
    seui <- merge_seurats(seus)
    rm(seus)
    DefaultAssay(seui) <- "RNA"
    seui <- NormalizeData(seui)
    seui <- FindVariableFeatures(seui,selection.method = "vst")
    
    seui <- cellcycle_scoring(seui)
    vars_to_regress <- "orig.ident"
    cmeta <- colnames(seui@meta.data)
    if ('percent.mt' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'percent.mt')}
    if ('S.Score' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'S.Score')}
    if ('G2M.Score' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'G2M.Score')}
    if ('seqrun' %in% cmeta) {
      if (length(unique(seui[['seqrun']]))>1) {
        vars_to_regress <- c(vars_to_regress,'seqrun')
      }
    }
    vars_to_regress <- vars_to_regress[!is.na(vars_to_regress)]
    
    message("--ScaleData ...")
    seui <- ScaleData(seui)
    if (reused!="") {
      save(seui,vars_to_regress,file=rd_file,compress=T)
    }
  }
  message("--RunPCA ...")
  seui <- RunPCA(seui,
                 features = VariableFeatures(object = seui))
  
  message("running harmony ...")
  
  seui <- RunHarmony(seui, group.by.vars = vars_to_regress)
  
  # Dimensional reduction and plotting
  message("--RunUMAP ...")
  seui <- RunUMAP(seui, dims = 1:D, reduction = "harmony")
  
  message("--FindNeighbors ...")
  seui <- FindNeighbors(seui, reduction = "harmony", dims = 1:D)
  
  message("--FindClusters ...")
  seui <- FindClusters(seui, resolution = 0.55)
  
  # message("--DimPlot ...")
  # DimPlot(seui, group.by = c("orig.ident"), ncol = 2)
  
  message("Done[harmony].")
  if (debug2==1){browser()}
  return(seui)
}

#checked 05/14/2021
aligned_by_harmony_w_cc <- function(seus,
                                    D=30,
                                    min.cells=50,
                                    fc.resolution=0.8,
                                    uvars_to_regress="",
                                    do_cluster=1,
                                    harmony=1,
                                    do_tsne=1,
                                    cell_cycle=1,
                                    debug2=0,
                                    reused="") {
  if (debug2==1){browser()}
  rds_file<-file.path(reused,"regress_out_scaled.rds")
  if (file.exists(rds_file)) {
    message(sprintf("reuse [%s] ...",rds_file))
    seui=readRDS(rds_file)
  } else {
    message("deleting SCT assay ...")
    seus <- lapply(seus,function(seu) {
      if ('RNA' %in% names(seu@assays)) {
        DefaultAssay(seu) <- "RNA"
      }
      
      if ('SCT' %in% names(seu@assays)) {
        seu@assays[['SCT']] <- NULL
        seu@reductions <- list()
        seu@graphs <- list()
      }
      seu
    })
    
    message("running harmony for data integration ...")
    seui <- merge_seurats(seus,debug2=debug2)
    rm(seus)
    DefaultAssay(seui) <- "RNA"
    seui <- NormalizeData(seui)
    seui <- FindVariableFeatures(seui,selection.method = "vst")

    message("--ScaleData ...")
    seui <- regress_out_scaledata(seui,seqrun=NA,cell_cycle=cell_cycle)
    if (reused!="") {
      saveRDS(seui,file=rds_file)
    }
  }
  
  if (do_cluster==0 & harmony==1){
    do_cluster=1
  }
  
  if (do_cluster==1) {
    message("--RunPCA ...")
    seui <- RunPCA(seui,
                   features = VariableFeatures(object = seui))
    
    message("running harmony ...")
    vars_to_regress <- "orig.ident"
    cmeta <- colnames(seui@meta.data)
  }
  
  red_method="pca"
  if (harmony==1) {
    if ('seqrun' %in% cmeta) {
      if (length(unique(seui[['seqrun']]))>1) {
        vars_to_regress <- c(vars_to_regress,'seqrun')
      }
    }
    if (uvars_to_regress!="") {
      vars_to_regress <- unique(c(uvars_to_regress,vars_to_regress))
    }
    
    seui <- RunHarmony(seui, group.by.vars = vars_to_regress)
    red_method="harmony"
  }
  
  if (do_cluster==1) {
    # Dimensional reduction and plotting
    message(sprintf("--RunUMAP [%s] ...",red_method))
    seui <- RunUMAP(seui, dims = 1:D, reduction = red_method)
    
    if (do_tsne==1) {
      message("--RunTSNE ...")
      seui <- RunTSNE(seui, dims = 1:D, reduction = red_method)
    }
    
    message("--FindNeighbors ...")
    seui <- FindNeighbors(seui, reduction = red_method, dims = 1:D)
    
    message("--FindClusters ...")
    seui <- FindClusters(seui, resolution = fc.resolution)
    
    # message("--DimPlot ...")
    # DimPlot(seui, group.by = c("orig.ident"), ncol = 2)
  }
  message("Done[harmony].")
  if (debug2==1){browser()}
  return(seui)
}


aligned_by_harmony_wo_cc <- function(seus,D=30,min.cells=50,debug2=0,reused="") {
  if (debug2==1){browser()}
  rd_file<-file.path(reused,"regress_out_scaled.rd")
  if (file.exists(rd_file)) {
    message(sprintf("reuse [%s] ...",rd_file))
    load(rd_file)
  } else {
    message("deleting SCT assay ...")
    seus <- lapply(seus,function(seu) {
      if ('SCT' %in% names(seu@assays)) {seu@assays$SCT <- NULL}
      return(seu)
    })
    message("running harmony for data integration ...")
    seui <- merge_seurats(seus)
    rm(seus)
    DefaultAssay(seui) <- "RNA"
    seui <- NormalizeData(seui)
    seui <- FindVariableFeatures(seui,selection.method = "vst")

    vars_to_regress <- "orig.ident"
    cmeta <- colnames(seui@meta.data)
    
    if ('seqrun' %in% cmeta) {
      if (length(unique(seui[['seqrun']]))>1) {
        vars_to_regress <- c(vars_to_regress,'seqrun')
      }
    }

    message("--ScaleData ...")
    seui <- ScaleData(seui)
    if (reused!="") {
      save(seui,vars_to_regress,file=rd_file,compress=T)
    }
  }
  message("--RunPCA ...")
  seui <- RunPCA(seui,
                 features = VariableFeatures(object = seui))
  
  message("running harmony ...")
  
  seui <- RunHarmony(seui, group.by.vars = vars_to_regress)
  
  # Dimensional reduction and plotting
  message("--RunUMAP ...")
  seui <- RunUMAP(seui, dims = 1:D, reduction = "harmony")
  
  message("--FindNeighbors ...")
  seui <- FindNeighbors(seui, reduction = "harmony", dims = 1:D)
  
  message("--FindClusters ...")
  seui <- FindClusters(seui, resolution = 0.55)
  
  # message("--DimPlot ...")
  # DimPlot(seui, group.by = c("orig.ident"), ncol = 2)
  
  message("Done[harmony].")
  if (debug2==1){browser()}
  return(seui)
}

cellcycle_genes_avail <- function(seu,debug2=0) {
  if (debug2==1){browser()}
  feats.avail <- rownames(GetAssayData(seu))
  ccg <- cc.genes.updated.2019
  ccg$s.genes <- ccg$s.genes[ccg$s.genes %in% feats.avail]
  ccg$g2m.genes <- ccg$g2m.genes[ccg$g2m.genes %in% feats.avail]
  return(ccg)
}

cellcycle_scoring <- function(seu,debug2=0) {
  seu <- tryCatch(
    {
      ccg <- cellcycle_genes_avail(seu,debug2=debug2)
      
      feats.avail <- rownames(GetAssayData(seu))
      
      if (any(ccg$s.genes %in% feats.avail) & any(ccg$g2m.genes %in% feats.avail)) {
        message("assigning cell cycle scores ...")
        seu <- CellCycleScoring(seu,
                                s.features = ccg$s.genes, 
                                g2m.features = ccg$g2m.genes, 
                                set.ident = TRUE)
        return(seu)
      }
    },
    error=function(cond){
      message(sprintf("cellcycle_scoring failed![%s]",cond))
      return(NA)
    },
    warning=function(cond){
      message(sprintf("cellcycle_scoring caused warning![%s]",cond))
      return(seu)
    },
    finally={
      message("Done[cellcycle_score].")
    }
  )
  return(seu)
}

regress_out_scaledata <- function(seu,seqrun="seqrun",cell_cycle=1,do_center=TRUE,split_by=NA,debug2=0) {
  if (debug2==1){browser()}
  
  if (DefaultAssay(seu) != "Protein" & cell_cycle==1) {
    seu <- cellcycle_scoring(seu,debug2=debug2)
  }
  
  cmeta <- colnames(seu@meta.data)
  vars_to_regress <- NA
  
  if ('percent.mt' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'percent.mt')}
  if ('S.Score' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'S.Score')}
  if ('G2M.Score' %in% cmeta) {vars_to_regress <- c(vars_to_regress,'G2M.Score')}
  
  if (!is.na(seqrun)) {
    if (seqrun %in% cmeta) {
      if (dim(unique(seu[[seqrun]]))[1]>1) {
        vars_to_regress <- c(vars_to_regress,seqrun)
      }
    }
  }
  vars_to_regress <- vars_to_regress[!is.na(vars_to_regress)]
  
  message(sprintf("scaling/regressing out [%s]...",paste0(vars_to_regress,collapse = ",")))
  
  if (length(vars_to_regress)>0) {
    
    if (is.na(split_by)) {
      seu <- ScaleData(seu, 
                       do.center=do_center,
                       vars.to.regress = vars_to_regress)
    } else {
      seu <- ScaleData(seu,
                       do.center=do_center,
                       split.by=split_by,
                       vars.to.regress = vars_to_regress)
    }
  } else {
    if (is.na(split_by)) {
      seu <- ScaleData(seu,do.center=do_center)
    } else {
      seu <- ScaleData(seu,do.center=do_center,split.by=split_by)
    }
  }
  return(seu)
}

integration_analysis_sctf <- function(sc.samples,si_features=3000,D=30,min.cells=50,debug2=0,reused="") {
  
  #use_condaenv(condaenv="Renv", conda="~/miniconda3/bin/conda")
  #reticulate::py_install(packages ='umap-learn',method="conda")
  if (debug2==1){browser()}
  
  sc_anchors_rd <- file.path(reused,"sc_anchors.rd")
  if (file.exists(sc_anchors_rd)) {
    load(sc_anchors_rd)
  } else {
    message("finding the number of features ...")
    feat_cnts <- sapply(sc.samples,function(seu){return(dim(seu)[1])})
    if (D > mean(feat_cnts)) {D <- mean(feat_cnts) * 0.5}
    
    cell_cnts <- sapply(sc.samples,function(seu){return(dim(seu@meta.data)[1])})
    sc.samples <- sc.samples[cell_cnts>min.cells]
    if (length(sc.samples)<2) {return(NA)}
    message("selectIntegrationFeatures ...")
    sc.features <- SelectIntegrationFeatures(object.list = sc.samples,
                                             nfeatures = si_features)

    feats.avail <- unique(Reduce(intersect,lapply(sc.samples,function(sc3){rownames(sc3)})))
    
    ccg <- cc.genes.updated.2019
    ccg$s.genes <- ccg$s.genes[ccg$s.genes %in% feats.avail]
    ccg$g2m.genes <- ccg$g2m.genes[ccg$g2m.genes %in% feats.avail]
    
    sc.features <- unique(c(sc.features,ccg$s.genes,ccg$g2m.genes))
    message("PrepSCTIntegration ...")
    sc.list <- PrepSCTIntegration(object.list = sc.samples,
                                  anchor.features = sc.features)
    
    rm(sc.samples)
    
    k_filterj <- 200
    sc.anchors <- NA
    j <- 0
    while(is.na(sc.anchors) & (j<8)) {
      message(sprintf("FindIntegrationAnchors[%d] ...",j))
      sc.anchors <- runFIA_with_tryCatch(sc.list,sc.features,D,k_filter=(k_filterj-25*j))
      j <- j + 1
    }
    rm(sc.list)
    # rm(sc.features)
    
    stopifnot(!is.na(sc.anchors))
    
    if (reused!="") {
      if (!file.exists(reused)) {dir.create(reused,showWarnings = F,recursive = T)}
      save(sc.anchors,file=sc_anchors_rd,compress=T)
    }
  }
  
  message(sprintf("--IntegrateData ..."))
  sc.integrated <- IntegrateData(anchorset = sc.anchors,
                                 normalization.method = "SCT",
                                 dims = 1:D)
  rm(sc.anchors)
  
  sc.integrated <- regress_out_scaledata(sc.integrated)

  sc.integrated <- RunPCA(sc.integrated, 
                          features = VariableFeatures(object = sc.integrated))

  message(sprintf("--RunUMAP ..."))
  sc.integrated <- RunUMAP(object = sc.integrated,
                           umap.method="umap-learn",
                           metric="correlation",
                           dims = 1:D)
  
  message(sprintf("--FindNeighbors ..."))
  sc.integrated <- FindNeighbors(sc.integrated, 
                                 dims = 1:D, 
                                 verbose = TRUE)
  
  message(sprintf("--FindClusters ..."))
  sc.integrated <- FindClusters(sc.integrated, 
                                verbose = TRUE)
  
  if (debug2==1){browser()}
  
  return(sc.integrated)
}

print_integ_result <- function(sci,
                               pdf_prefix,
                               split_by,
                               group_by="seurat_clusters",
                               ncol = 3,
                               pdfwinch=0,
                               pdfhinch=0,
                               debug=0) {
  
  pdf_file <- sprintf("%s_dimPlot.pdf",pdf_prefix)
  orig_sample_cnt <- length(unique(sci@meta.data$orig.ident))
  
  if (debug==1) {browser()} 
  # else {pdf(pdf_file,height=5.5,width=11)}
  else {
    if (pdfwinch>0){
      pdf(pdf_file,width = pdfwinch,height = pdfhinch)
    } else {
      pdf(pdf_file)
    }
  }
  
  p1 <- DimPlot(sci, 
                reduction = "umap", 
                group.by = split_by,
                label = FALSE,
                combine=TRUE)
  
  plot(p1,asp=1.)

  split_by_ids <- unique(sci[[split_by]])
  
  N <- dim(split_by_ids)[1]
  M <- ceiling(N/4)
  for (m in 1:M) {
    st1 <- (m-1)*4 + 1
    ed2 <- m*4
    if (ed2>N) {ed2 <- N}
    
    seu <- subset(sci,cells = Cells(sci)[sci@meta.data[,split_by] %in% split_by_ids[st1:ed2,]])

    p2 <- DimPlot(seu,
                  reduction = "umap", 
                  split.by = split_by, 
                  group.by = "seurat_clusters",
                  ncol=2,
                  label = TRUE,
                  combine=TRUE)
    plot(p2)
  }

  if (group_by != "seurat_clusters") {
    if (group_by %in% colnames(sci@meta.data)) {
      for (m in 1:M) {
        st1 <- (m-1)*4 + 1
        ed2 <- m*4
        if (ed2>N) {ed2 <- N}
        
        seu <- subset(sci,cells = Cells(sci)[sci@meta.data[,split_by] %in% split_by_ids[st1:ed2,]])
        
        p3 <- DimPlot(seu,
                      reduction = "umap", 
                      split.by = split_by, 
                      group.by = group_by,
                      ncol=2,
                      label = FALSE,
                      combine=TRUE)
        plot(p3)
      }
      
    } else {
      message(sprintf("group_by[%s] does not exist in sci@meta.data",group_by))
    }
  }

  if (debug==0) {
    dev.off()
  }
}

print_integ_featurePlot <- function(sc.integrated,pdf_prefix,group_by,goi,debug2=0) {
  
  library(randomcoloR)
  if (debug2==1){browser()}
  DefaultAssay(sc.integrated) <- "SCT"

  features <- rownames(sc.integrated@assays$SCT@data)
  goi_eff <- sort(goi[goi %in% features])
  
  L <- length(goi_eff)
  if (L == 0) {
    return(FALSE)
  }
  
  S <- dim(unique(sc.integrated[[group_by]]))[1]
  
  if (L > 6) {
    goi_effs <- split(goi_eff, ceiling(seq_along(goi_eff)/(L/ceiling(L/6))))
  } else {
    goi_effs <- list(goi_eff)
  }
  
  for (i in 1:length(goi_effs)) {
    if (debug2==0){
      pdf(sprintf("%s_featurePlot_%d.pdf",pdf_prefix,i),width=5.5,height=11)
    }
    
    if (debug2==1){browser()}
    my_pal <- distinctColorPalette(S)
    
    p <- FeaturePlot(sc.integrated, 
                     features = goi_effs[[i]], 
                     split.by = group_by)
    plot(p)
    
    plots <- VlnPlot(sc.integrated, 
                     features = goi_effs[[i]],
                     split.by = group_by,
                     pt.size = 0, 
                     combine = FALSE)
        
    p <- CombinePlots(plots = plots, ncol = 1)
    plot(p)
    
    
    p <- RidgePlot(sc.integrated, 
                   features = goi_effs[[i]], 
                   ncol = 1,
                   group.by = group_by)
    plot(p)
    
    p <- DoHeatmap(sc.integrated, 
                   features = goi_effs[[i]], 
                   size = 3,
                   group.by = group_by,
                   assay ="SCT",
                   slot = "scale.data")
    plot(p)
    if (debug2==1){browser()}
    if (debug2==0){dev.off()}
  }
  return(TRUE)
}

#' https://divingintogeneticsandgenomics.rbind.io/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
#' p_list<- FeaturePlotSingle(pbmc, feature= "MS4A1", metadata_column = "samples", pt.size = 0.05, order =TRUE)
#' 
FeaturePlotSingle<- function(obj, feature, metadata_column, assay2="RNA", red2="umap", ...) {
  
  all_cells<- colnames(obj)
  # groups<- levels(obj@meta.data[, metadata_column])
  groups<- unique(obj@meta.data[[metadata_column]])
  
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[[assay2]]@data[feature, ])
  maximal<- max(obj[[assay2]]@data[feature, ])
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


markedFeaturePlot_solo <- function(goi,sc3,pdf_prefix,plotName) {

  DefaultAssay(sc3) <- "RNA"
  features <- rownames(sc3@assays$RNA@data)
  goi_eff <- goi[goi %in% features]
  
  pdf_file <- sprintf("%s.pdf",pdf_prefix)
  if (plotName == "FeaturePlot") {
    plot_list <- FeaturePlot(sc3,features=goi_eff,pt.size=0.2,coord.fixed=TRUE,combine = FALSE)
  } else if (plotName == "DoHeatmap") {
    plot_list <- DoHeatmap(sc3, features=goi_eff,combine = FALSE)
  }
  
  multi.page <- ggarrange(plotlist=plot_list,
                          ncol = 2,
                          nrow = 2)
  message(pdf_file)
  ggexport(multi.page, filename = pdf_file)
  message("Done.")
}

depreciated_markedFeaturePlot_solo <- function(goi,sc3,pdf_prefix,plotName) {
  source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
  
  # browser()
  DefaultAssay(sc3) <- "RNA"
  features <- rownames(sc3@assays$RNA@data)
  goi_eff <- goi[goi %in% features]
  
  L <- length(goi_eff)
  if (L == 0) {
    return(FALSE)
  }
  goi_eff <- sort(goi_eff)
  
  if (plotName=="FeaturePlot") {
    M <- 4
  } else if (plotName == "DoHeatmap") {
    M <- 100
  }
  
  if (L > M) {
    goi_effs <- chunk2(goi_eff,M)
  } else {
    goi_effs <- list(goi_eff)
  }
  
  pdf_file <- sprintf("%s.pdf",pdf_prefix)
  message(pdf_file)
  pdf(pdf_file)
  
  for (i in 1:length(goi_effs)) {
    if (length(goi_effs[[i]]) > 0) {
      if (plotName == "FeaturePlot") {
        p <- FeaturePlot(sc3,features=goi_effs[[i]],ncol=2,coord.fixed=TRUE)
      } else if (plotName == "DoHeatmap") {
        p <- DoHeatmap(sc3, features = goi_effs[[i]]) + NoLegend()
      }
      plot(p)
    }
  }
  dev.off()
}

diff_analy_prep <- function(sci) {
  sci$cluster_sgroup <- paste0(Idents(sci),'_',sci$sgroup)
  sci$cluster <- Idents(sci)
  Idents(sci) <- "cluster_sgroup"
  return(sci)
}


diff_analy_anchors <- function(sci,comps,gsample='group_sample',min_cell_cnts=3) {
  
  #comps: a pair of sgroups to compare with
  #comps <- data.table(ctrl=c('pre-infusion','day14','pre-infusion'),
  #	expr=c('day14','day30','day30'))
  # browser()
  DefaultAssay(sci) <- "RNA"
  
  cluster_labels <- sort(unique(sci@meta.data$seurat_clusters))
  
  sci2 <- diff_analy_prep(sci)
  diffMarkers <- list()
  for (clab in cluster_labels) {
    
    message(sprintf("cluster_label:%s",clab))
    for (r in 1:dim(comps)[1]) {
      message(sprintf("sgroup:%s vs. %s",comps$ctrl[r],comps$expr[r]))
      ctrl_lab <- sprintf("%s_%s",clab,comps$ctrl[r])
      expr_lab <- sprintf("%s_%s",clab,comps$expr[r])
      
      if ((length(which(sci2$cluster_sgroup==ctrl_lab)) >= min_cell_cnts) & (length(which(sci2$cluster_sgroup==expr_lab)) >= min_cell_cnts)) {
        message('enough samples for diff analysis.')
        comp_label <- sprintf("%s_%s_%d",gsample,clab,r)
        fm_result <- FindMarkers(sci2,
                                 ident.1=expr_lab,
                                 ident.2=ctrl_lab,
                                 verbose=TRUE)
        diffMarkers[[comp_label]] <- as.data.table(fm_result)
        diffMarkers[[comp_label]]$gsample <- gsample
        diffMarkers[[comp_label]]$cluster <- clab
        diffMarkers[[comp_label]]$gene <- rownames(fm_result)
        diffMarkers[[comp_label]]$sgroups <- sprintf("%s_%s",comps$ctrl[r],comps$expr[r])
      }
    }
  }
  
  return(rbindlist(diffMarkers))
  
}


tmp_split_multilanes_crcount <- function(gene_names,seurat3_raws,total_sidx) {
  
  rcnts <- lapply(names(seurat3_raws),function(sname) {
    
    message(sname)
    srd1 <- seurat3_raws[[sname]]
    features <- rownames(srd1@assays$RNA@counts)
    cell_bc <- colnames(srd1@assays$RNA@counts)
    
    bc_start_idx <- sapply(1:total_sidx,function(i){which(str_detect(cell_bc,as.character(i)))[1]})
    bc_start_idx <- c(bc_start_idx,length(cell_bc)+1)
    rcnts2 <- list()
    for (i in 1:(length(bc_start_idx)-1)) {
      message(i)
      from_bc_idx <- bc_start_idx[i]
      to_bc_idx <- bc_start_idx[i+1]-1
      
      rcnts2[[i]] <- as.data.table(Matrix::rowSums(srd1@assays$RNA@counts[,(from_bc_idx:to_bc_idx)]))
    }
    
    # browser()
    rcnt <- do.call(cbind,rcnts2)
    colnames(rcnt) <- sprintf("%d",1:total_sidx)
    rcnt2 <- rcnt[match(gene_names,features),]
    rcnt2$gene <- gene_names
    return(rcnt2)
  })
  
  names(rcnts) <- names(seurat3_raws)
  
  return(rcnts)
}


report_sc3_dim <- function(seus,tsv_fpath=NA,qc_comment="N",append2=FALSE,debug2=0) {
  if (debug2==1){browser()}
  rc_dims <- sapply(names(seus),function(sname) {
    return(dim(seus[[sname]]@assays$RNA@counts))
  }, simplify = TRUE)
  
  rc_dims <- as.data.table(t(rc_dims),keep.rownames=TRUE)
  colnames(rc_dims) <- c("sample","features","cells")
  
  if (!is.na(tsv_fpath)) {
    rc_dims$qc_comment <- qc_comment
    fwrite(rc_dims,file=tsv_fpath,append=append2,sep="\t")
  }
  if (debug2==1){browser()}
  return(rc_dims)
}

cell_typing <- function(cell_marker,dge_per_cluster) {

  marker_genes <- rownames(cell_marker)
  ctypes <- colnames(cell_marker)
  
  dge_ctypes <- list()
  
  for (j in 1:dim(cell_marker)[2]) {
    # j <- 2
    ctype <- ctypes[j]
    mg_indicator <- cell_marker[,j]>0
    
    mprof <- get_mprof(marker_genes[mg_indicator])
    
    sample_tags <- names(dge_per_cluster)
    dgeis <- list()
    for (sample_tag in sample_tags) {
      message(sprintf("j:%d/sample_tag:%s",j,sample_tag))
      
      dgei <- as.data.table(dge_per_cluster[[sample_tag]])
      dgei <- dgei[order(avg_logFC,-p_val)]
      dgei$sample <- sample_tag
      dgei$ctype <- ctype
      
      dgei$rscore <- 1:dim(dgei)[1]/dim(dgei)[1]
      dgei_up <- dgei[gene %in% mprof$up & avg_logFC>=0.,]
      
      if (length(mprof$down)>0) {
        dgei$rscore <- (1-1:dim(dgei)[1]/dim(dgei)[1])
        dgei_down <- dgei[gene %in% mprof$down & avg_logFC<=0.,]
        dgei <- rbind(dgei_up,dgei_down)
      } else {
        dgei <- dgei_up
      }
      dgeis[[sample_tag]] <- dgei
    }
    dge_ctypes[[j]] <- rbindlist(dgeis)
  }
  dge_dt <- rbindlist(dge_ctypes)
  
  dge_by_samples <- split(dge_dt,by=c("sample"))
  
  sample_cluster_ctype_profs <- list()
  
  for (sname in names(dge_by_samples)) {
    dge_by_clusters <- split(dge_by_samples[[sname]],by='cluster')
    dge_by_clusters <- dge_by_clusters[!isEmpty(dge_by_clusters)]
    
    for (cid in names(dge_by_clusters)) {
      
      name2 <- sprintf("%s_%s",sname,cid)
      message(name2)
      
      dge_cluster <- dge_by_clusters[[cid]]
      
      dge_df <- dcast(dge_cluster,ctype ~ gene, value.var = 'rscore')
      
      ctypes <- dge_df$ctype
      dge_df <- dge_df[,2:dim(dge_df)[2]]
      
      dge_df[is.na(dge_df)] <- 0.
      
      topK <- rowSums(as.data.table(dge_df) > 0)
      topK[topK>1] <- 2
      topK[topK==1] <- 1
      
      dge_df <- cbind(dge_df,topK)
      
      rownames(dge_df) <- ctypes
      
      D <- dim(dge_df)[2]
      
      topKscore <- apply(dge_df,1,function(cvals) {
        rscores <- cvals[1:(D-1)]
        topK <- cvals[D]
        meanScore <- mean(sort(rscores,decreasing = T)[1:topK])
        return(meanScore)
      })
      
      ctype2 <- names(topKscore)
      cscore <- as.data.table(topKscore)
      cscore$ctype <- ctype2
      cscore$cluster <- cid
      cscore$sample <- sname
      
      cscore <- cscore[order(-topKscore)]
      sample_cluster_ctype_profs[[name2]] <- cscore
    }
  }
  
  sample_cluster_ctype_profs <- rbindlist(sample_cluster_ctype_profs)
  
  #dge_by_sample_cids <- dge_by_sample_cids[!isEmpty(dge_by_sample_cids)]
  screp <- as.data.table(dcast(sample_cluster_ctype_profs,sample + cluster ~ ctype, fun.aggreate = max, value.var = "topKscore"))
  
  screps <- split(screp,by="sample")
  
  for (sname in names(screps)) {
    screpj <- screps[[sname]]
    M <- dim(screpj)[2]
    cluster_id <- screpj$cluster
    screpj <- screpj[,3:M]
    ctypes <- colnames(screpj)
    
    rownames(screpj) <- paste0('cluster_',cluster_id)
    colnames(screpj) <- tcell_to_fullname(ctypes)
    pdf_file <- file.path(wkd,sprintf('%s_tcellSubtype_x_cluster.pdf',sname))
    message(pdf_file)
    screpj[is.na(screpj)] <- 0.
    pheatmap(screpj,main = sprintf("sctf_%s[T-cell subtypes vs. cluster]",sname),file=pdf_file)
    
  }
  
  return(screps)
}


integs_to_dgeTable <- function(sc.integs,comps,wkd,min_cell_cnts=3) {
  
  # comps <- data.table(ctrl=c('pre-infusion','day14','pre-infusion'),
  # 										expr=c('day14','day30','day30'))
  
  message("identify differentially expressed genes across conditions ...")
  
  group_by <- "sgroup"
  
  group_samples <- names(sc.integs)
  
  diffMarkers <- list()
  
  for (gsample in group_samples) {
    
    message(sprintf("gsample:%s",gsample))
    
    sci <- sc.integs[[gsample]]
    if (!is.na(sci)) {
      DefaultAssay(sci) <- "RNA"
      sci <- NormalizeData(sci, verbose = FALSE)
      
      cluster_labels <- sort(unique(sci@meta.data$seurat_clusters))
      
      sci2 <- seurat3_diff_analy_prep(sci)
      
      for (clab in cluster_labels) {
        message(sprintf("cluster_label:%s",clab))
        for (r in 1:dim(comps)[1]) {
          message(sprintf("sgroup:%s vs. %s",comps$ctrl[r],comps$expr[r]))
          ctrl_lab <- sprintf("%s_%s",clab,comps$ctrl[r])
          expr_lab <- sprintf("%s_%s",clab,comps$expr[r])
          
          if ((length(which(sci2$cluster_sgroup==ctrl_lab)) >= min_cell_cnts) & (length(which(sci2$cluster_sgroup==expr_lab)) >= min_cell_cnts)) {
            message('enough samples for diff analysis.')
            comp_label <- sprintf("%s_%s_%d",gsample,clab,r)
            fm_result <- FindMarkers(sci2,
                                     ident.1=expr_lab,
                                     ident.2=ctrl_lab,
                                     verbose=TRUE)
            diffMarkers[[comp_label]] <- as.data.table(fm_result)
            diffMarkers[[comp_label]]$gsample <- gsample
            diffMarkers[[comp_label]]$cluster <- clab
            diffMarkers[[comp_label]]$gene <- rownames(fm_result)
            diffMarkers[[comp_label]]$sgroups <- sprintf("%s_%s",comps$ctrl[r],comps$expr[r])
          }
        }
      }
    }
  }
  
  diffMarker_dt <- rbindlist(diffMarkers)
  warnings()
  save(diffMarker_dt,file=file.path(wkd,'seurat3_diff.rd'),compress=TRUE)
  
  wb <- createWorkbook("seurat3_diff_analysis")
  dts <- split(diffMarker_dt,by=c("gsample","sgroups"))
  
  for (i in 1:length(dts)) {
    sheet_name <- names(dts)[i]
    addWorksheet(wb,sheetName = sheet_name)
    writeData(wb,sheet_name,dts[[i]])
  }
  saveWorkbook(wb,file=file.path(wkd,"seurat3_diff.xlsx"),overwrite=T)
}


export_sc3_to_text <- function(seu,wkd,exp_tag="sample",slot="norm",debug=0) {
  
  if (debug==1){browser()}
  
  if (!file.exists(wkd)) {dir.create(wkd,showWarnings = FALSE, recursive=TRUE)}
  
  if (slot == "norm") {
    sc3_dt <- as.data.table(seu@assays$SCT@data)
    genes <- rownames(seu@assays$SCT@data)
  } else {
    sc3_dt <- as.data.table(seu@assays$RNA@counts)
    genes <- rownames(seu@assays$RNA@counts)
  }
  
  cluster_avail <- 0
  if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    cluster_avail <- 1
  }
  
  cb_mat <- colnames(sc3_dt)
  
  meta_info <- data.table(cb=rownames(seu@meta.data),
                          orig.ident=seu@meta.data$orig.ident)
  
  if (cluster_avail==1){
    meta_info$cluster_id <- seu@meta.data$seurat_clusters
  }
  sample_cb <- paste0(meta_info[match(cb_mat,meta_info$cb),orig.ident],'.',meta_info$cb)
  
  colnames(sc3_dt) <- sample_cb
  
  sc3_dt$gene <- genes
  
  N <- dim(sc3_dt)[2]
  sc3_dt <- cbind(sc3_dt[,N,with=F],sc3_dt[,1:(N-1),with=F])
  
  tsv_fpath <- file.path(wkd,sprintf("%s_%s_data_matrix.tsv",exp_tag,slot))
  fwrite(sc3_dt,file=tsv_fpath,sep="\t",row.names=FALSE,col.names = TRUE)
  tsv_gz <- file.path(wkd,sprintf("%s_%s_matrix.tsv.gz",exp_tag,slot))
  if (file.exists(tsv_gz)) {unlink(tsv_gz)}
  gzip(tsv_fpath,destname=tsv_gz)
  
  if (cluster_avail==1) {
    tsv_fpath <- file.path(wkd,sprintf("%s_%s_cluster_info.tsv",exp_tag,slot))
    meta_info$sample_cb <- paste0(meta_info$orig.ident,'.',meta_info$cb)
    fwrite(meta_info,file=tsv_fpath,sep="\t",row.names=FALSE,col.names = TRUE)
    tsv_gz <- file.path(wkd,sprintf("%s_%s_cluster_info.tsv.gz",exp_tag,slot))
    if (file.exists(tsv_gz)) {unlink(tsv_gz)}
    gzip(tsv_fpath,destname=tsv_gz)
  }
}

export_sc3_to_loom <- function(seu,wkd,exp_tag,slot="raw") {
  
  # library(devtools)
  # devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  library(hdf5r)
  library(loomR)
  
  if (slot=="norm") {
    seu_mtrx <- seu@assays$RNA@data
  } else {
    seu_mtrx <- seu@assays$RNA@counts
  }
  
  loom_fpath <- file.path(wkd,sprintf("%s_data_matrix.loom",exp_tag))
  
  genes<- list(rownames(seu_mtrx))
  names(genes) <- genes
  
  cids<- list(colnames(seu_mtrx))
  names(cids) <- cids
  
  loomR::create(filename=loom_fpath,
                data=seu_mtrx,
                gene.attrs=genes,
                cell.attrs=cids,
                do.transpose=TRUE,
                overwrite=TRUE,
                verbose=FALSE)
  
  return(loom_fpath)
}

export_to_tsv <- function(seu,out_pref,assay2="RNA",slot2="data") {

  mtrx <- as.matrix(GetAssayData(object = seu,
                                 assay=assay2,
                                 slot = slot2))
  
  mtrx <- as.data.table(mtrx,keep.rownames=TRUE)
  setnames(mtrx,'rn','gene')
  fwrite(mtrx,
         file=sprintf("%s_mtrx.tsv.gz",out_pref),
         sep="\t",
         col.names=TRUE,
         quote=FALSE)
  
  smeta <- as.data.table(seu@meta.data,keep.rownames=TRUE)
  setnames(smeta,'rn','cell')
  fwrite(smeta,
         file=sprintf("%s_meta.tsv.gz",out_pref),
         sep="\t",
         col.names=TRUE,
         quote=FALSE)
}

seurat3_rename_to_sample_w_cluster <- function(sci) {
  sci$cluster_sgroup <- paste0(Idents(sci),'_',sci$sgroup)
  sci$cluster <- Idents(sci)
  Idents(sci) <- "cluster_sgroup"
  return(sci)
}

diff_from_integration <- function(sci,
                                   cluster_col="seurat_clusters",
                                  assay2="RNA",
                                   min_cell_cnts=3,
                                  method="MAST",
                                   ncpu=3,
                                   debug=0) {

  diffMarker_dt <- NA
  
  if (debug==1){browser()}
  stopifnot(!is.null(method) & !is.na(method))
  
  if (!is.na(sci)) {

    DefaultAssay(sci) <- assay2
    
    cluster_labels <- sort(unique(sci@meta.data[[cluster_col]]))
    sci$cluster_sample <- paste0(sci@meta.data[[cluster_col]],'_',sci@meta.data$comp)
    
    Idents(sci) <- "cluster_sample"
    
    diffMarker2 <- mclapply(cluster_labels,function(clab) {
      comp_label <- sprintf("cluster_label [%s]",clab)
      
      ctrl_lab <- sprintf("%s_ctrl",clab)
      expr_lab <- sprintf("%s_expr",clab)
      
      if (debug==1){browser()}
      
      if ((length(which(sci$cluster_sample==ctrl_lab)) >= min_cell_cnts) & 
          (length(which(sci$cluster_sample==expr_lab)) >= min_cell_cnts)) {

        # message(sprintf("assay[%s];test.method[%s]",assay2,method))
        
        fm_result <- FindMarkers(sci,
                                 test.use=method,
                                 assay=assay2,
                                 ident.1=expr_lab,
                                 ident.2=ctrl_lab,
                                 only.pos=FALSE,
                                 verbose=FALSE) #NOTE that ident.1 is of our interest
        
        diffMarker <- as.data.table(fm_result)
        diffMarker$cluster_id <- clab
        diffMarker$gene <- rownames(fm_result)
      } else {
        message(sprintf('not enough samples for diff analysis [%s].',comp_label))
        diffMarker <- NA
      }
      return(diffMarker)
    }
    ,mc.cores = ncpu
    )
    
    if (debug==1){browser()}
    
    diffMarker2 <- diffMarker2[!is.na(diffMarker2)]
    
    diffMarker_dt <- rbindlist(diffMarker2)
  }
  
  if (debug==1){browser()}
  
  warnings()
  return(diffMarker_dt)
}

#' Perform DEsingle on Seurat Integrated Data
#' NOTE: this is under development. DGE postprocess (table/visualization) will be affected.
#' 
#'    Changjin Hong, hongc2@ccf.org

#'    ref: DEsingle for detecting three types of differential expression in single-cell RNA-seq data
#'    
#' @param sci seurat integrated object
#' @param cluster_col meta.data column to use for clustering/group
#' @param min_cell_cnts the min number of cells to avoid DGE comp error due to insufficient samples
#' @param ncpu The number of cpus to utilize
#' @param debug 0 (no debug), 1 (debug; test with the first sample only)
#' @return DGE report dt.table
#' 
#' 
DEsingle_from_integration <- function(sci,
                                  cluster_col="seurat_clusters",
                                  min_cell_cnts=3,
                                  ncpu=3,
                                  debug=0) {
  
  diffMarker_dt <- NA
  
  if (debug==1){browser()}
  
  if (!is.na(sci)) {
    param <- MulticoreParam(workers = ncpu, progressbar = TRUE)
    register(param)
    
    assay2 <- "SCT"
    DefaultAssay(sci) <- assay2
    cluster_labels <- sort(unlist(unique(sci[[cluster_col]]),use.names=F))
    
    if (is.na(cluster_col)) {
      stopifnot(length(unique(sci$comp)),2)
      sci$cluster_sample <- sci$comp
    } else {
      sci$cluster_sample <- paste0(unlist(sci[[cluster_col]],use.names=F),'_',sci$comp)
    }
    
    Idents(sci) <- "cluster_sample"
    
    diffMarker2 <- lapply(cluster_labels,function(clab) {
      comp_label <- sprintf("cluster_label [%s]",clab)
      
      ctrl_lab <- sprintf("%s_ctrl",clab)
      expr_lab <- sprintf("%s_expr",clab)
      
      if ((length(which(sci$cluster_sample==ctrl_lab)) >= min_cell_cnts) & 
          (length(which(sci$cluster_sample==expr_lab)) >= min_cell_cnts)) {
        
        message(sprintf('enough samples for diff analysis [%s].',comp_label))
        
        message(sprintf("assay[%s]",assay2))
        cbc.ctrl <- Cells(sci)[sci$cluster_sample==ctrl_lab]
        cbc.expr <- Cells(sci)[sci$cluster_sample==expr_lab]
    
        group2 <- factor(c(rep(ctrl_lab,length(cbc.ctrl)),
                           rep(expr_lab,length(cbc.expr))))
        
        if (debug==1){browser()}
        
        sci.mtx <- as.matrix(GetAssayData(sci,slot = "counts"))[,c(cbc.ctrl,cbc.expr)]

        # Detecting the DE genes in parallelization with ncpu cores
        fm_result <- DEsingle(counts = sci.mtx, group = group2, parallel = TRUE, BPPARAM = param)
        rm(sci.mtx)
        fm_result.classified <- DEtype(results = fm_result, threshold = 0.05)
        
        diffMarker <- as.data.table(fm_result.classified)
        diffMarker$cluster_id <- clab
        diffMarker$gene <- rownames(fm_result.classified)
      } else {
        diffMarker <- NA
      }
      return(diffMarker)
    }
    )
    
    if (debug==1){browser()}
    
    diffMarker2 <- diffMarker2[!is.na(diffMarker2)]
    
    diffMarker_dt <- rbindlist(diffMarker2)
  }
  
  if (debug==1){browser()}
  
  warnings()
  return(diffMarker_dt)
}

ctrl_vs_expr_for_comp <- function(snames,comp_sheet,sheet_idx=1,ctrl_by_or=NA,expr_by_or=NA,debug2=0) {
  if (debug2==1) {browser()}
  
  if (!is.na(comp_sheet) & file.exists(comp_sheet)) {
    if (endsWith(comp_sheet,"xlsx")|endsWith(comp_sheet,"xls")) {
      message(sprintf("reading xlsx[%s],sheet_idx[%d]",comp_sheet,sheet_idx))
      comp2_dt <- as.data.table(read.xlsx(comp_sheet,sheet=sheet_idx))
    } else {
      comp2_dt <- fread(comp_sheet)
    }
    
    message(sprintf("column headers avail at %s: %s",comp_sheet,paste0(colnames(comp2_dt),collapse = ',')))
    if (any(snames %in% comp2_dt$sample)) {
      comp2_dt <- comp2_dt[(sample %in% snames),]
      sfield <- "sample"
    } else if (any(snames %in% comp2_dt$sample_tag)) {
      comp2_dt <- comp2_dt[(sample_tag %in% snames),]
      sfield <- "sample_tag"
    } else {
      stop("verify comp_sheet if either sample or sample_tag does match with seurat object ...")
    }
    
    stopifnot(nrow(comp2_dt[comp=="ctrl",])>0 & nrow(comp2_dt[comp=="expr",])>0)
    
    comp_dt <- data.table(ctrl=comp2_dt[comp=="ctrl",paste0(get(sfield),collapse=",")],
                          expr=comp2_dt[comp=="expr",paste0(get(sfield),collapse=",")])
  } else {
    comp_dt <- as.data.table(t(combn(snames,2,simplify = TRUE)))
    colnames(comp_dt) <- c("ctrl","expr")
  }
  if (debug2==1) {browser()}
  return(comp_dt)
}

batch_dge_from_seu_integs_objs <- function(sc.integs,
                                           comp_sheet="",
                                           cluster_col="seurat_clusters",
                                           dge_method="MAST",
                                           ncpu=3,
                                           sheet_idx=1,
                                           debug=0) {
  if (debug==1){browser()}
  
  diffMarker_dts <- list()
  sci_names <- names(sc.integs)
  N <- length(sci_names)
  for (gsi in 1:N) {
    gsample <- sci_names[[gsi]]
    # gsample <- names(sc.integs)[[2]]
    message(gsample)
    sci <- sc.integs[[gsample]]
    snames <- sort(unique(sci@meta.data$orig.ident))
    if (length(snames)>1) {
      comp_dt <- ctrl_vs_expr_for_comp(snames,comp_sheet,sheet_idx=sheet_idx,debug2 = debug)
      
      diffMarker_dt2 <- list()
      R <- dim(comp_dt)[1]
      for (r in 1:R) {
        # r <- 1
        job_notice <- sprintf("%s[%d/%d],comp[%d/%d]",gsample,gsi,N,r,R)
        message(job_notice)
        samples_to_comp <- unlist(apply(comp_dt[r,],2,comma_string_to_list),use.names = F)
        
        seu <- subset(sci,cells = Cells(sci)[sci@meta.data$orig.ident %in% samples_to_comp])
        seu@meta.data$comp <- "comp"
        seu@meta.data$comp[which(seu@meta.data$orig.ident %in% comma_string_to_list(comp_dt[r,'ctrl']))] <- "ctrl"
        seu@meta.data$comp[which(seu@meta.data$orig.ident %in% comma_string_to_list(comp_dt[r,'expr']))] <- "expr"
        
        if (dge_method=="desingle") {
          diffMarker_dt <- DEsingle_from_integration(seu,
                                                     cluster_col=cluster_col,
                                                     min_cell_cnts=3,
                                                     ncpu=ncpu,
                                                     debug=debug)
        } else {
          diffMarker_dt <- diff_from_integration(seu,
                                                 cluster_col=cluster_col,
                                                 min_cell_cnts=3,
                                                 method=dge_method,
                                                 ncpu=ncpu,
                                                 debug=debug)
        }
        if (!is.na(diffMarker_dt)) {
          diffMarker_dt$ctrl <- comp_dt[r,'ctrl']
          diffMarker_dt$expr <- comp_dt[r,'expr']
          diffMarker_dt$comp <- sprintf("comp_%d",r)
          diffMarker_dt$gsample <- gsample
          diffMarker_dt$group1 <- diffMarker_dt$gsample
          diffMarker_dt$group2 <- NA
          if (!isEmpty(grep("\\.",diffMarker_dt$gsample))) {
            diffMarker_dt$group1 <- tstrsplit(diffMarker_dt$gsample,"\\.")[[1]]
            diffMarker_dt$group2 <- tstrsplit(diffMarker_dt$gsample,"\\.")[[2]]
          }
          
          cmpname <- diffMarker_dt[1,paste0(ctrl,"_",expr)]
          diffMarker_dt2[[cmpname]] <- diffMarker_dt
        }
        message(sprintf("Done[%s]",job_notice))
      }
      diffMarker_dts[[gsample]] <- diffMarker_dt2
    }
  }
  return(diffMarker_dts)
}


findConservedMarkers_with_tryCatch <- function(sci,expr_ident,group_by,method="MAST",debug2=0,exp_tag="fcm") {
  ge_conserved_df <- tryCatch(
    {
      DefaultAssay(sci) <- "RNA"
      if (debug2==1){browser()}
      
      if (method=="desingle") {
        cbc.ctrl <- names(Idents(sci))[Idents(sci)!=expr_ident]
        cbc.expr <- names(Idents(sci))[Idents(sci)==expr_ident]
        
        group2 <- factor(c(rep(sprintf("no_%s",expr_ident),length(cbc.ctrl)),
                           rep(expr_ident,length(cbc.expr))))
        
        sci.mtx <- as.matrix(GetAssayData(sci,slot = "counts"))[,c(cbc.ctrl,cbc.expr)]
        
        # Detecting the DE genes in parallelization with ncpu cores
        fm_result <- DEsingle(counts = sci.mtx, group = group2)
        rm(sci.mtx)
        ge_conserved_df <- DEtype(results = fm_result, threshold = 0.05)
      } else {
        ge_conserved_df <-FindConservedMarkers(sci,
                                   test.use=method,
                                   ident.1 = expr_ident,
                                   grouping.var = group_by,
                                   verbose = FALSE)
      }
    },
    error=function(cond) {
      # browser() #debug
      message(sprintf('FindConservedMarkers(fcm=%s,expr_ident=%s,group_by=%s) is failed!',exp_tag,expr_ident,group_by))
      return(NA)
    },
    finally={
      # message('Done.')
    }
  )
  return(ge_conserved_df)
}

conserved_from_integration <- function(sc.integs,cluster_col="seurat_clusters",group_by="orig.ident",dge_method="MAST",ncpu=4,debug2=0) {
  
  # cons_rd_file <- sprintf("%s.rd",out_pref)
  
  conservedMarkers <- list()
  # for (scomp in c("CART.P12")) {
  for (scomp in names(sc.integs)) { #debug

    message(sprintf("scomp[%s]",scomp))
    if (debug2==1) {browser()}
    
    sci <- sc.integs[[scomp]]
    DefaultAssay(sci) <- "RNA"
    Idents(sci) <- cluster_col
    
    if (!is.na(sci)) {
      ident_cluster_df <- table(sci$orig.ident,unlist(sci[[cluster_col]]))
      if (dge_method=="desingle") {
        comm_clusters <- colnames(ident_cluster_df)[colSums(ident_cluster_df)>=10]
      } else {
        comm_clusters <- colnames(ident_cluster_df)[which(colSums(ident_cluster_df>2)>=2)]
      }
      
      cluster_markers <- mclapply(comm_clusters, function(identj) {
        message(identj)
        ret <-findConservedMarkers_with_tryCatch(sci,identj,group_by,method=dge_method,debug2=debug2,exp_tag=scomp)
        return(ret)}
        ,mc.cores=ncpu
        )
      
      names(cluster_markers) <- comm_clusters
      cluster_markers <- cluster_markers[!is.na(cluster_markers)]
      
      conMarkers <- lapply(names(cluster_markers),function(identj) {
        message(identj)
        cmarker <- cluster_markers[[identj]]
        if (dim(cmarker)[1]==0) {
          cmarker <- NA
        } else {
          cmarker$cluster_id<-identj
          cmarker$gene<-rownames(cmarker)
          cmarker <- as.data.table(cmarker)
        }
        return(cmarker)
      })
      names(conMarkers) <- names(cluster_markers)
      conMarkers <- conMarkers[!is.na(conMarkers)]
      
      if (debug2==1) {browser()}
      
      consMarker_dt <- rbindlist2(conMarkers)
      
      conservedMarkers[[scomp]] <- list(shared=consMarker_dt)
      
    }
  }
  if (debug2==1) {browser()}
  # save(conservedMarkers,file=cons_rd_file,compress=TRUE)
  return(conservedMarkers)
}
#' To get DGE tables from seurat3 SCT integration where DGE is defined from one cluster vs. all the other clusters in the same sample
#' 
#' To obtain DGE tables from seurat3 SCT integration
#' more detail
#' 
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'    
#' @param sc.integs a list of Seurat3 SCT integration objects

#' @return list of data.table of DGE.
#' 
dge_within_sample_from_integration <- function(sc.integs) {
  
  dgei <- list()
  for (cmp_name in names(sc.integs)) {
    # browser()
    sci <- sc.integs[[cmp_name]]
    
    DefaultAssay(sci) <- "RNA"
    
    sci <- NormalizeData(sci, verbose = FALSE)
    
    if (!is.na(sci)) {
      Idents(object = sci) <- "orig.ident"
      
      dgei[[cmp_name]] <- mclapply(unique(Idents(sci)),function(sample) {
        
        message(sprintf('DGE analy from sc3 integ [%s,%s]',cmp_name,sample))
        
        sci_j <- subset(x = sci, idents = sample)
        Idents(sci_j) <- sci_j$seurat_clusters
        
        dge_j <- FindAllMarkers(object = sci_j,
                                only.pos = FALSE)
        return(dge_j)
      },mc.cores = length(unique(Idents(sci))))
      
      names(dgei[[cmp_name]]) <- unique(Idents(sci))
    }
  }
  
  return(dgei)
}

# --------------------
#' To annotate sc3 object integ clusters by cellassign
#' 
#' To obtain ca_fits and ca_fits_matrix_info, refer to cellassign_pbmc.ipynb
#' 
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'    
#' @param sc.integs a list of Seurat3 SCT integration objects
#' @param ca_fits a list of cellassign() output variables
#' @param ca_fits_matrix_info: a list of matrix row(genes)/col(cell_ids) information used in ca_fits
#' 
#' @return list of Seurat3 integration objects where celltype is added into $meta.data
#' 
#' 
annotate_integ_cluster_by_cellassign <- function(sc.integs,annot_cell,ncpu=3,debug2=0) {
  #ca_fits,ca_fits_matrix_info
  #
  message("annotating a cell type predicted by cellassign ...")
  if (debug2==1) {browser()}
  sc.integs.annot <- lapply(sc.integs, function(sc.integ) {
    if (debug2==1) {browser()}
    sm <- sc.integ@meta.data
    sm$uuid <- paste0(sm$orig.ident,':',tstrsplit(rownames(sm),"_")[[1]])
    
    anns <- list()
    for (sname in unique(sm$orig.ident)) {
      anns[[sname]] <- data.table(uuid=paste0(sname,':',annot_cell$mtx_meta[[sname]]$cellid),
                                  cell_type=annot_cell$fits[[sname]]$cell_type)
    }
    ann <- rbindlist(anns)
    
    sm$cellassign <- ann[match(sm$uuid,ann$uuid),cell_type]
    sm$uuid <- NULL
    sc.integ@meta.data <- sm
    return(sc.integ)
  }) #,mc.cores = ncpu
  
  return(sc.integs.annot)
}


# --------------------
#' To annotate sc3 object integ clusters by cellassign
#' 
#' To obtain ca_fits and ca_fits_matrix_info, refer to cellassign_pbmc.ipynb
#' 
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'    NOTE: using table(), the implementation can be simpler!
#'    
#' @param sc.integ a Seurat3 SCT integration object
#' @param group_by2 column variable of sc.integ meta.data
#' 
#' @return a data.frame of count matrix in the format of (group_by2,sample)
#' 
cell_counts_by <- function(sc.integ,
                           out_prefix,
                           cell_id="seurat_clusters",
                           plot_title="cell_count",
                           debug=0) {

  # meta <- data.table(sample = sc.integ$orig.ident,
  #                    cid = sc.integ$seurat_clusters)
  if (debug==1) {browser()}
  meta <- data.table(sample = sc.integ$orig.ident,
                     ucell_id = sc.integ@meta.data[[cell_id]],
                     cid = unlist(sc.integ[["seurat_clusters"]]))
  
  meta_by_sample <- split(meta,by="sample")
  samples <- names(meta_by_sample)
  
  nn_sizes <- list()
  nn_size.pcts <- list()
  cellid_sizes <- list()
  cellid_size.pcts <- list()
  
  wb <- createWorkbook(plot_title)
  
  pdf_file <- sprintf('%s.pdf',out_prefix)
  message(sprintf("generating pdf_file[%s]",pdf_file))
  pdf(pdf_file)
  
  for (sname in samples) {
    a.ldf <- meta_by_sample[[sname]][,.N,by=c("cid","ucell_id")]
    a.wdf <- dcast(a.ldf,cid ~ ucell_id, value.var = 'N')
    a.mat <- as.matrix(a.wdf[,2:dim(a.wdf)[2]])
    rownames(a.mat) <- a.wdf[,'cid']
    a.mat[is.na(a.mat)] <- 0
    
    sheet_name <- sname
    addWorksheet(wb,sheetName = sheet_name)
    writeData(wb,sheet_name,a.mat,colNames=TRUE,rowNames=TRUE)
    
    nn_size <- data.table(sample=sname,
                          cid=rownames(a.mat),
                          val=rowSums2(a.mat))
    
    nn_size.pct <- data.table(sample=sname,
                              cid=rownames(a.mat),
                              val=100.*nn_size$val/sum(nn_size$val))
    
    cellid_size <- data.table(sample=sname,
                              cid=colnames(a.mat),
                              val=colSums2(a.mat))
    
    cellid_size.pct <- data.table(sample=sname,
                                  cid=colnames(a.mat),
                                  val=100.*cellid_size$val/sum(cellid_size$val))
    
    nn_sizes[[sname]] <- nn_size
    nn_size.pcts[[sname]] <- nn_size.pct
    cellid_sizes[[sname]] <- cellid_size
    cellid_size.pcts[[sname]] <- cellid_size.pct
    
    if (debug==1) {browser()}
    r1c2 <- 1
    a.norm_mat <- apply(a.mat,r1c2,function(vec) {
      S <- sum(vec)
      if (S>0) {
        norm_vec <- vec/S
      } else {
        norm_vec <- vec
      }
      return(norm_vec)
    })
    if (r1c2==1) {
      a.norm_mat <- t(a.norm_mat)
    }
    pheatmap(a.norm_mat,angle_col=45,main=sname,cluster_rows=FALSE,cluster_cols=FALSE)
  }
  if (debug==1) {browser()}
  
  p2 <- list()
  nn_sizes <- rbindlist(nn_sizes)
  p2[['nn_sizes']] <- ggplot(data=nn_sizes, aes(x=sample, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: total cells",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("sample") +
    ylab("cell counts")
  
  nn_size.wdf <- dcast(nn_sizes,cid ~ sample, value.var = 'val')
  sheet_name <- "nn_sizes"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,nn_size.wdf)
  
  # ----------
  nn_size.pcts <- rbindlist(nn_size.pcts)
  
  p2[['nn_size.pcts']] <- ggplot(data=nn_size.pcts, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell counts[%%]",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("sample") +
    ylab("cell counts [%]")
  
  nn_size.pct.wdf <- dcast(nn_size.pcts,cid ~ sample, value.var = 'val')
  sheet_name <- "nn_size.pcts"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,nn_size.pct.wdf)
  
  p <- CombinePlots(p2,ncol=1)
  plot(p)
  # -----------
  p2 <- list()
  cellid_sizes <- rbindlist(cellid_sizes)
  
  p2[['cellid_sizes']] <- ggplot(data=cellid_sizes, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell ids",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("cell counts") +
    ylab("cell ids") # + coord_flip()
  
  cellid_size.wdf <- dcast(cellid_sizes,cid ~ sample, value.var = 'val')
  sheet_name <- "cellid_sizes"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,cellid_size.wdf)
  # ------------
  cellid_size.pcts <- rbindlist(cellid_size.pcts)
  
  p2[['cellid_size.pcts']] <- ggplot(data=cellid_size.pcts, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell ids[%%]",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("cell counts") +
    ylab("cell ids [%]")# + coord_flip()
  
  cellid_size.pct.wdf <- dcast(cellid_size.pcts,cid ~ sample, value.var = 'val')
  sheet_name <- "cellid_size.pcts"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,cellid_size.pct.wdf)
  
  xlsx_file <- sprintf("%s.xlsx",out_prefix)
  message(sprintf("creating an excel file [%s]",xlsx_file))
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  
  p <- CombinePlots(p2,ncol = 1)
  plot(p)
  dev.off()
  if (debug==1) {browser()}
  message(sprintf("Done [%s].",plot_title))
}

scrublet_annot_batch <- function(bc_matrix_dir,debug2=0) {
  if (debug2==1){browser()}
  mtx_file <- file.path(bc_matrix_dir,'matrix.mtx')
  if (file.exists(sprintf("%s.gz",mtx_file))) {
    cmd <- sprintf("gunzip -fc %s.gz > %s",mtx_file,mtx_file)
    system(cmd)
  } else {
    stop(sprintf("check if %s.gz exists",mtx_file))
  }
  
  fea_file <- file.path(bc_matrix_dir,'features.tsv')
  if (file.exists(sprintf("%s.gz",fea_file))) {
    cmd <- sprintf("gunzip -fc %s.gz > %s",fea_file,fea_file)
    system(cmd)
  }
  
  scrublet_annot <- annotate_doublets(mtx_fpath=normalizePath(mtx_file),
                                      feature_fpath=normalizePath(fea_file))
  
  if (file.exists(mtx_file)) {
    unlink(mtx_file)
  }
  
  if (file.exists(fea_file)) {
    unlink(fea_file)
  }
  
  return(as.numeric(scrublet_annot[[1]]))
}


get_doublet_rate_in_10x <- function(ncells) {
  # ref: https://www.biotech.wisc.edu/services/gec/services/rochelightcyclerservices
  
  # ~0.8% per 1,000 cells; ~1.6% per 2,000 cells; ~2.3% per 3,000 cells; ~3.1% per 4,000 cells; ~3.9% per 5,000 cells; etc.
  
  library(splines)
  x <- c(1000,  2000, 3000, 4000, 5000, 8000)
  y <- c(0.008,0.016,0.023,0.031,0.039,0.070)
  fit2 <- lm( y~ns(x, 3) )
  # plot(x,y, xlim=c(1000,8001), ylim=c(0,0.1))
  # xx <- seq(1000,8001, length.out=250)
  # lines(xx, predict(fit2, data.frame(x=xx)), col='orange')
  return(predict(fit2,data.frame(x=ncells)))
}


rbindlist2 <- function(dt_list) {
  
  ucname <- sort(unique(unlist(lapply(dt_list,colnames))))
  head_col<-c("gene","minimump_p_val","max_pval","cluster_id")
  ucname <- c(head_col,setdiff(ucname,head_col))
  
  dts <- lapply(dt_list,function(dt){
    delta_col <- setdiff(ucname,colnames(dt))
    if (!isEmpty(delta_col)) {
      dt[,c(delta_col):=NA]
    }
    setcolorder(dt,ucname)
    return(dt)
  })
  return(rbindlist(dts))
}

solo_FindAllMarkers <- function(seu,ident2=NA,assay2="SCT",de_test="MAST",only_pos=FALSE,rd_fpath=NA) {
  reuse <- 0
  if (!is.na(rd_fpath)) {
    if (file.exists(rd_fpath)) {reuse <- 1}
  }
  
  if (reuse==1) {
    message(sprintf("use prev result [%s]",rd_fpath))
    load(rd_fpath)
  } else {
    if (assay2 %in% names(seu@assays)) {
      my_assay <- assay2
    } else {
      my_assay <- "SCT"
    }
    
    if (!is.na(ident2)) {
      Idents(seu) <- ident2
    }
    seu.markers.full <- FindAllMarkers(object = seu,
                                       assay = my_assay,
                                       test.use = de_test,
                                       only.pos = FALSE)
    if (reuse==1) {
      save(seu.markers.full,file=rd_fpath,compress=T)
    }
  
  }
  return(seu.markers.full)
}

#checked 05/14/2021
batch_solo_FindAllMarkers <- function(seus,ident2=NA, assay2="SCT",only_pos=FALSE,ncpu=3,debug2=0) {
  if (debug2==1){browser()}
  snames <- names(seus)
  
  cluster_marker_dges <- mclapply(snames, function(sname) {
    message(sprintf("running FindAllMarkers[%s];Idents[%s]",sname,ident2))
    
    assays_avail <- names(seus[[sname]]@assays)
    if (!(assay2 %in% assays_avail)) {
      if (("SCT" %in% assays_avail)) {
        assay2 <- "SCT"
      } else {
        if ("RNA" %in% assays_avail) {
          assay2 <- "RNA"
        } else {
          stop(sprintf("no %s,SCT,RNA is available from the seurat object",assay2))
        }
      }
    }
   
    ret <- solo_FindAllMarkers(seus[[sname]],ident2=ident2,assay2=assay2,only_pos=only_pos)
    message(sprintf("Done[%s]",sname))
    return(ret)
  }
  ,mc.cores = ncpu
  )
  names(cluster_marker_dges) <- snames
  if (debug2==1){browser()}
  return(cluster_marker_dges)
}


cell_id_by_singler <- function(seu,ref.sce,pmid,label2id="label.fine",debug2=0) {
  if (debug2==1){browser()}
  
  query.sce <- as.SingleCellExperiment(seu)

  shared_feat <- intersect(rownames(ref.sce),rownames(query.sce))
  
  message(sprintf("ref feat# [%d]",length(rownames(ref.sce))))
  message(sprintf("test feat# [%d]",length(rownames(query.sce))))
  message(sprintf("shared feat# [%d]",length(shared_feat)))
  message("cell ID annotation in progress ...")
  pred.id <- SingleR(test=query.sce,
                     ref=ref.sce,
                     method="single",
                     labels=ref.sce[[label2id]])
  message("Done.")
  seu[[pmid]] <- pred.id$pruned.labels
  if (debug2==1){browser()}
  return(list(seu=seu,pred.id=pred.id))
}


single_seurat_to_sce <- function(seu,assay="SCT") {
  
  DefaultAssay(seu) <- assay
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(seu,slot = "counts"),
                                            logcounts=GetAssayData(seu,slot = "data")),
                              colData = as.data.frame(seu@meta.data))
  
  reducedDims(sce) <- SimpleList(PCA=seu@reductions$pca@cell.embeddings[,1:2],
                                 UMAP=seu@reductions$umap@cell.embeddings)
  return(sce)
}

seurat3_to_monocle3 <- function(seurat,clid="seurat_clusters",debug2=0) {
  if (debug2==1){browser()}
  # part one, gene annotations
  
  # gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), 
  #                                  row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
  gene_annotation <- as.data.frame(rownames(seurat),
                                   row.names = rownames(seurat))
  colnames(gene_annotation) <- "gene_short_name"
  
  # part two, cell information
  
  cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]],
                                 row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
  
  colnames(cell_metadata) <- "barcode"
  
  cell_metadata <- cbind(cell_metadata,seurat@meta.data)
  
  # part three, counts sparse matrix
  New_matrix <- seurat@assays[["RNA"]]@counts
  #New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
  New_matrix <- New_matrix[rownames(seurat), ]
  expression_matrix <- New_matrix
  
  
  ### Construct the basic cds object
  cds_from_seurat <- new_cell_data_set(expression_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)
  
  
  ### Construct and assign the made up partition
  
  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  
  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
  
  ### Assign the cluster info
  
  list_cluster <- seurat@meta.data[[clid]]
  names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
  
  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  
  
  ### Could be a space-holder, but essentially fills out louvain parameters
  
  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  
  
  ### Assign UMAP coordinate
  cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings
  
  ### Assign feature loading for downstream module analysis
  
  cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings
  
  
  ### Learn graph, this step usually takes a significant period of time for larger samples
  
  print("Learning graph, which can take a while depends on the sample")
  cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)
  return(cds_from_seurat)
  
  
}
get_cluster_by_expr <- function(seu2,assay2,slot2="data",debug2=0) {
  if(debug2==1) {browser()}
  avgExpr <- AverageExpression(seu2,assay=assay2,slot = slot2)
  avgExpr[[assay2]]
}


evaluate_metrics_by_cellranger_recommend<-function() {
  lower_bnd=data.table(`Valid Barcodes`=75,
                       `Q30 Bases in RNA Read`=65,
                       `Estimated Number of Cells`=500,
                       `Fraction Reads in Cells`=70,
                       `Mean Reads per Cell`=20000,
                       `Reads Mapped Confidently to Transcriptome`=30,
                       `Reads Mapped Antisense to Gene`=0)
  
  upper_bnd=data.table(`Valid Barcodes`=100,
                       `Q30 Bases in RNA Read`=100,
                       `Estimated Number of Cells`=10000,
                       `Fraction Reads in Cells`=100,
                       `Mean Reads per Cell`=Inf,
                       `Reads Mapped Confidently to Transcriptome`=100,
                       `Reads Mapped Antisense to Gene`=10)
  
  a = rbind(lower_bnd,upper_bnd)
  a$type = c("lower","upper")
  a
}

cell_count_w_meta <- function(sc_raws,sample_meta=NA,qc_comment="N",debug2=0) {
  if (debug2==1) {browser()}
  cell_cnts <- as.data.table(sapply(sc_raws,function(seu) {dim(seu)[2]}),keep.rownames = T)
  colnames(cell_cnts) <- c("sample","cnt")
  
  if (!is.na(sample_meta)) {
    cell_cnts<-merge(cell_cnts,sample_meta,by="sample")
  }
  cell_cnts$qc <- qc_comment
  cell_cnts
}


scina_with_tryCatch <- function(seu,markers,allow_unknown,outd) {
  pred_cell_ids <- tryCatch(
    {
      pred_cell_ids <- SCINA(GetAssayData(seu,slot = "counts"),
                             markers,
                             max_iter = 100,
                             convergence_n = 10,
                             convergence_rate = 0.999,
                             sensitivity_cutoff = 0.9,
                             rm_overlap=FALSE,
                             allow_unknown=allow_unknown,
                             log_file=file.path(outd,sprintf('SCINA_%s.log',sname)))
    },
    error=function(cond){
      message('scina failed in handling input!')
      return(NULL)
    },
    finally={
      message('scina processes input!')
    }
  )
  return(pred_cell_ids)
}
