#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(Rsubread)
library(Rsamtools)
library(GenomicRanges)
library(parallel)
library(jsonlite)
library(reshape2)
# library(circlize)
library(Matrix)
library(uwot)
library(pheatmap)

# library(gamlss)
source(file.path(Sys.getenv('R_UTIL'),'lib_workflow.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))

add_chr_contig_header <- function(bam){
  out_bam <- sprintf('%s_chr.bam',bam)
  
  cmd <- sprintf("samtools view -H %s | sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - %s > %s",bam,bam,out_bam)
  if (!file.exists(out_bam)){
    system(cmd)
  }

  cmd <- sprintf("samtools index %s",out_bam)
  bai <- sprintf("%s.bai",out_bam)
  if (!file.exists(bai)){
    system(cmd)
  }
  
  return(out_bam)
}

get_gunzip_fn <- function(file) {
	S = str_length(file)
	L = length(grep('\\.gz$',file,value=TRUE))
	if (L>0){
		file2 = substr(file,1,S-3)
		system(sprintf('gunzip -fc %s > %s',file,file2))
	} else {
		file2 = file
	}
	return(file2)
}

run_dreme_motif_search <- function(fasta_gz,control_fasta_gz,outdir) {
	message(sprintf('run dreme motif analysis on [%s;%s]',fasta_gz,control_fasta_gz))
	if (!dir.exists(outdir)){
		dir.create(outdir, showWarnings = FALSE)
	}

	fasta = get_gunzip_fn(fasta_gz)
	control_fasta = get_gunzip_fn(control_fasta_gz)

	cmd = sprintf("dreme -oc %s -p %s -n %s -png -dna -norc",outdir,fasta,control_fasta)
	message(cmd)
	system(cmd)
	message('Done.')
	if (fasta_gz!=fasta){unlink(fasta)}
	if (control_fasta_gz!=control_fasta){unlink(control_fasta)}
}

macs2_peak <- function(bam,sample_name,out_dir,peak_type='sharp',out_prefix='narrow',in_bam=NA,macs2_path=NA) {
  
  
  
  if (is.na(macs2_path)){
    macs2 <- "macs2"
  } else {
    macs2 <- macs2_path
  }
  
  cmd <- sprintf("%s callpeak -t %s -f BAM -g hs -n %s -B --outdir %s",macs2,bam,sample_name,out_dir)
  
  #cmd <- sprintf("%s callpeak -t %s -f BAM -g hs -n %s -B --keep-dup all -m 3 100 --outdir %s",macs2,bam,sample_name,out_dir)
  
  if (!is.na(in_bam)){cmd <- sprintf("%s -c %s",cmd,in_bam)}
  
  if (peak_type=='broad'){
    cmd <- sprintf("%s --broad",cmd)
  }
  
  message(cmd)
  system(cmd) #debug
  
  if (out_prefix=='narrow'){
    peak_tsv <- file.path(out_dir,sprintf('%s_peaks.narrowPeak',sample_name))
  } else {
    peak_tsv <- file.path(out_dir,sprintf('%s_peaks.xls',sample_name))
  }
  
  return(peak_tsv)
}

macs2_merge <- function(macs2_peaks) {
  
  dts <- list()
  if (startsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  dts[[1]] <- fread(sprintf("%s -v ^# %s",grep_cmd,macs2_peaks[1]))
  
  if (startsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  dts[[2]] <- fread(sprintf("%s -v ^# %s",grep_cmd,macs2_peaks[2]))
  
  peaks <- lapply(dts,function(x){
    GRanges(seqnames = x$chr,
            ranges = IRanges(start = x$start,
                             end = x$end,
                             names = x$name))
  })
  merged <- union(peaks[[1]],peaks[[2]])
  merged_dt <- data.table(GeneID=paste0('macs2_',1:length(merged)),
                          Chr=as.character(seqnames(merged)),
                          Start=start(merged),
                          End=end(merged),
                          Strand='*')
  return(merged_dt)
}

macs2_narrowpeak_concat <- function(macs2_peaks) {
  
  dts <- list()
  
  tmp_out1 <- sprintf('%s.tmp',macs2_peaks[1])
  
  if (endsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n > %s",grep_cmd,macs2_peaks[1],tmp_out1)

  message(cmd)
  system(cmd)
  
  tmp_out2 <- sprintf('%s.tmp',macs2_peaks[2])
  if (endsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n > %s",grep_cmd,macs2_peaks[2],tmp_out2)
  
  message(cmd)
  system(cmd)

  cmd <- sprintf('cat %s %s | sort -k1,1V -k2,2n -k3,3n',tmp_out1,tmp_out2)
  
  merged <- fread(cmd,col.names=c('Chr','Start','End'))
  merged$GeneID <- paste0('macs2_',1:dim(merged)[1])
  merged$Strand <- '*'

  unlink(tmp_out1)
  unlink(tmp_out2)
  
  return(merged)
}

macs2_narrowpeak_merge <- function(macs2_peaks) {
  
  dts <- list()
  if (endsWith(macs2_peaks[1], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n",grep_cmd,macs2_peaks[1])
  message(cmd)
  dts[[1]] <- fread(cmd)
  
  if (endsWith(macs2_peaks[2], ".gz")) {grep_cmd <- 'zgrep'
  } else {grep_cmd <- 'grep'}
  cmd <- sprintf("%s -v ^# %s | cut -f1,2,3 | sort -k1,1V -k2,2n -k3,3n",grep_cmd,macs2_peaks[2])
  message(cmd)
  dts[[2]] <- fread(cmd)
  
  peaks <- lapply(dts,function(x){
    GRanges(seqnames = x$V1,
            ranges = IRanges(start = x$V2,
                             end = x$V3))
  })
  merged <- union(peaks[[1]],peaks[[2]])
  merged_dt <- data.table(GeneID=paste0('macs2_',1:length(merged)),
                          Chr=as.character(seqnames(merged)),
                          Start=start(merged),
                          End=end(merged),
                          Strand='*')
  return(merged_dt)
}

get_total_reads_from_bam <- function(bam){
  cmd <- sprintf("samtools flagstat %s | head -n1 | cut -d\' \' -f1",bam)
  message(cmd)
  out <- system(cmd,intern = TRUE)
  num_reads <- as.integer(out)
  return(num_reads)
}

samtools_flagstat <- function(bam,samtools_bin='samtools',reuse=T){
  
  flag_stat_file <- paste0(bam,'.flagstat')
  cmd <- sprintf("samtools flagstat %s > %s",bam,flag_stat_file)
  message(cmd)
  
  if (!reuse | !file.exists(flag_stat_file)) {
    out <- system(cmd,intern = TRUE)
  }
  
  con  <- file(flag_stat_file, open = "r")
  
  stat.item <- list()
  stat.val <- list()
  
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    val <- as.numeric(str_extract(oneLine, pattern = "[0-9]+(?=\\s\\+)"))
    
    if (grepl('total', oneLine)) {
      stat.item <- c(stat.item,'total')
      stat.val <- c(stat.val,val)
    } else if (grepl('secondary',oneLine)) {
      stat.item <- c(stat.item,'secondary')
      stat.val <- c(stat.val,val)
    } else if (grepl('supplementary',oneLine)) {
      stat.item <- c(stat.item,'supplementary')
      stat.val <- c(stat.val,val)
    } else if (grepl('duplicates',oneLine)) {
      stat.item <- c(stat.item,'duplicates')
      stat.val <- c(stat.val,val)
    } else if (grepl('with itself and mate mapped',oneLine)) {
      stat.item <- c(stat.item,'both_mapped')
      stat.val <- c(stat.val,val)
    } else if (grepl('paired in sequencing',oneLine)) {
      stat.item <- c(stat.item,'paired_reads')
      stat.val <- c(stat.val,val)
    } else if (grepl('read1',oneLine)) {
      stat.item <- c(stat.item,'read1')
      stat.val <- c(stat.val,val)
    } else if (grepl('read2',oneLine)) {
      stat.item <- c(stat.item,'read2')
      stat.val <- c(stat.val,val)
    } else if (grepl('properly paired',oneLine)) {
      stat.item <- c(stat.item,'properly_paired')
      stat.val <- c(stat.val,val)
    } else if (grepl('singletons',oneLine)) {
      stat.item <- c(stat.item,'singletons')
      stat.val <- c(stat.val,val)
    } else if (grepl('with mate mapped to a different chr \\(mapQ>=5\\)',oneLine)) {
      stat.item <- c(stat.item,'mate_in_diffchr_mapq5')
      stat.val <- c(stat.val,val)
    } else if (grepl('with mate mapped to a different chr',oneLine)) {
      stat.item <- c(stat.item,'mate_in_diffchr')
      stat.val <- c(stat.val,val)
    }  else if (grepl('mapped',oneLine)) {
      stat.item <- c(stat.item,'mapped')
      stat.val <- c(stat.val,val)
    } else {
      message(sprintf('WARNING[%s]',oneLine))
    }
  } 
  
  close(con)
  
  stat_dt <- data.table(item=stat.item,val=stat.val)
  if (file.exists(flag_stat_file)) {
    if (!reuse) {
      unlink(flag_stat_file)
    }
  }
  return(stat_dt)
}


macs2_bdgdiff <- function(samp,
                          oprefix,
                          out_dir,
                          loglr_cutoff=1,
                          mode='TF',
                          macs2_bin='macs2'){
  # samp <- data.table(bam=c(bam1,bam2),
  # bdg=c(bdg1,bdg2),
  # bdgi=c(bdgi1,bdgi2),
  # nread=c(nread1,nread2),
  # nreadi=c(nreadi1,nreadi2),
  # group = c("ctrl","experiment"),
  # name=strsplit(args$name,",")[[1]])

  # ref: https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands
  
  if (mode == 'TF'){
    minlen <- 120
    maxgap <- 60
  } else if (mode =='histone'){
    minlen <- 120
    maxgap <- 60
  }
  
  cmd <- sprintf('%s bdgdiff --t1 %s --t2 %s --c1 %s --c2 %s --d1 %d --d2 %d -l %d -g %d --o-prefix %s --outdir %s',
                 macs2_bin,samp$bdg[1],samp$bdg[2],
                 samp$bdgi[1],samp$bdgi[2],
                 as.integer(round(samp$nread[1]/1e6)),
                 as.integer(round(samp$nread[2]/1e6)),
                 minlen,maxgap,oprefix,out_dir)
  message(cmd)
  system(cmd)
  
}

encode_narrowPeak <- function(encode_narrow_peak_file){
  #contig_name, start, end, peak_name, -10log10(qval), ?, fold-change, -log10(pval), -log10(qval), relative_summit_position
  if (file.exists(encode_narrow_peak_file)) {
    encode_npeak <- fread(encode_narrow_peak_file,header = F,col.names = c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
    #browser()
    encode_npeak <- encode_npeak[order(chr,start,end),]
    
    return(encode_npeak)
  } else {
    stop(sprintf('the input file[%s] does not exist!',encode_narrow_peak_file))
  }
}

prep_deseq2 <- function(bam_dt,region2eval,countMultiMappingReads=F,allowMultiOverlap=F,fraction=F,ncore=8) {
  
  R <- dim(region2eval)[1]
  C <- dim(bam_dt)[1]
  
  message(sprintf('R=%d',R)) #debug
  message(sprintf('C=%d',C))
  
  fc <- mclapply(bam_dt$bam,
                 function(bam){
                   Rsubread::featureCounts(bam,
                                 annot.ext = region2eval,
				 countMultiMappingReads=countMultiMappingReads,
				 allowMultiOverlap=allowMultiOverlap,
				 fraction=fraction)
                   },
                 mc.cores=ncore)
  
  cts <- matrix(NA, nrow = R, ncol = C)
  
  rownames(cts) <- paste0(fc[[1]]$annotation$Chr,':',
                          fc[[1]]$annotation$Start,'-',
                          fc[[1]]$annotation$End)
  colnames(cts) <- bam_dt$name
  
  for (i in 1:C){
    cts[,i] <- fc[[i]]$counts
  }

  cts <- round(cts)
  
  coldata <- DataFrame(group=bam_dt$group)
  rownames(coldata) <- bam_dt$name
  coldata$group <- factor(coldata$group)
  return(list(cts=cts,coldata=coldata))
}


picard_dedup <- function(bam,bam_out_dir,min_mapq=0,picard_jar_path='/media/sammy/projects/apps/picard/latest/picard.jar', reuse=TRUE){
  
  
  bam_fname <- basename(bam)
  mdup_bam <- file.path(bam_out_dir,sprintf('%s.mdup.bam',bam_fname))
  out_metrics <- file.path(bam_out_dir,sprintf('%s.mdup.metrics',bam_fname))
  
  if (is.na(picard_jar_path)){
    cmd_prefix <- "java -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT"
  } else {
    cmd_prefix <- sprintf("java -jar %s MarkDuplicates VALIDATION_STRINGENCY=SILENT",picard_jar_path)
  }
  if(file.exists(mdup_bam)&reuse) {
    message('reuse the previous result ...')
  } else {
    cmd <- sprintf('%s I=%s O=%s M=%s',cmd_prefix,bam,mdup_bam,out_metrics)
    message(cmd)
    system(cmd)
  }
  
  if (min_mapq>0) {
    mdup_mapq_bam <- file.path(bam_out_dir,sprintf('%s.mapq%d.bam',bam_fname,min_mapq))
    ncpu <- 4
    #cmd <- sprintf('samtools view -@ %d -b -h -F0x404 %s > %s',ncpu,mdup_bam,dedup_bam)
    if(file.exists(mdup_mapq_bam)&reuse){
      message('reuse the previous result ...')
    } else {
      cmd <- sprintf('samtools view -@ %d -b -h -q %d %s > %s',ncpu,min_mapq,mdup_bam,mdup_mapq_bam)
      message(cmd)
      system(cmd)
    }
    out_bam <- mdup_mapq_bam
  } else {
    out_bam <- mdup_bam
  }
  
  build_bai(out_bam,reuse)
  return(out_bam)
}

picard_insertSize <- function(wks,outd,debug2=0) {
	
	if (debug2==1) {browser()}
	if (!file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	
	java <- get_prog_info_yaml('java')
	java_cmd <- sprintf("%s %s",java$bin_path,java$default_opt)
	
	picard <- get_prog_info_yaml("picard")
	
	insertSizeFiles <- sapply(1:nrow(wks),function(r) {
		bam <- wks$bam[r]
		sample <- wks$sample[r]
		cmd <- java_cmd
		cmd <- sprintf("%s %s",cmd,picard$bin_path)
		cmd <- sprintf("%s %s",cmd,"CollectInsertSizeMetrics")
		cmd <- sprintf("%s I=%s",cmd,bam)
		output1 <- file.path(outd,sprintf("%s_picard_insertsize.txt",sample))
		cmd <- sprintf("%s O=%s",cmd,output1)
		pdf1 <- file.path(outd,sprintf("%s_picard_insertsize.pdf",sample))
		cmd <- sprintf("%s H=%s",cmd,pdf1)
		cmd <- sprintf("%s M=0.5",cmd)
		message(cmd)
		if (debug2==0) {system(cmd)}
		return(output1)
	})
	if (debug2==1) {browser()}
	return(insertSizeFiles)
}

build_bai <- function(out_bam,reuse=TRUE) {
  
  cmd <- sprintf('samtools index %s',out_bam)
  bai <- sprintf('%s.bai',out_bam)
  if (file.exists(bai) & reuse) {
    message('reuse the previous result ...')
  } else {
    system(cmd)
  }
  return(bai)
}

log2fc <- function(ctrl_cnt,expl_cnt){
  if (ctrl_cnt>0){
    if (expl_cnt>0){
      fcL2 <- log2(expl_cnt/ctrl_cnt)
    } else {
      fcL2 <- -Inf
    }
  } else {
    fcL2 <- Inf
  }
  return(fcL2)
}

bam_to_bdg <- function(bam,out_bdg,bin_size=5,bamCoverage_bin='bamCoverage'){
  
  cmd <- sprintf('%s --bam %s --binSize %d --outFileFormat bedgraph --outFileName %s',bamCoverage_bin,bam,bin_size,out_bdg)
  
  message(cmd)
  
  if (!file.exists(out_bdg)) {
    system(cmd)
  }
}

bedtools_bamtobed_qc <- function(in_bam,tag_bed,bedtools_bin="bedtools"){
  
  cmd <- sprintf("samtools view -hb -f 0x404 %s | %s bamtobed -i stdin > %s",
                 in_bam,
                 bedtools_bin,
                 tag_bed)
  message(cmd)
  system(cmd)
}

bedtools_bamtobed <- function(in_bam,out_bed,bedtools_bin="bedtools",optStr=NA){
  cmd <- sprinf("%s bamtobed",bedtools_bin)
  if (is.na(optStr)){
    cmd <- sprintf("%s %s",cmd,optStr)
  }
  cmd <- sprintf("%s -i %s > %s",cmd,in_bam,out_bed)
  message(cmd)
  system(cmd)
}


parse_encode_summary_json <- function(json_fn,row2print){
  
  js1 <- fromJSON(json_fn)

  head <- strsplit(js1$qc_files$head[row2print],'\t')[[1]]

  content <- strsplit(js1$qc_files$contents[row2print],'\t')[[1]]
  
  
  return (data.table(head=head,content=content))
}

diffReps <- function(bed_co,bed_tr,outfile,
                     genome="hg19",
                     test_model='nb',
                     window=200,
                     frag=150,
                     diffReps_bin='diffReps.pl'){
  cmd <- diffReps_bin
  
  cmd <- sprintf("%s -tr %s -co %s -gn hg19 -re %s -me %s --frag %d --nproc 8 --window %d",diffReps_bin,bed_tr,bed_co[1],outfile,test_model,frag,window)
  
  message(cmd)
  system(message)
}


macs2_read_narrow_peak <- function(narrow_peak_fn){
  macs_narrow_peak_cols <- c('chr','start','end','name','ucsc_score','strand','fc','pval_nlog10','qval_nlog10','summit')
  return(fread(narrow_peak_fn,col.names=macs_narrow_peak_cols))
}

sort_matrix_by_rowsum <- function(m2) {
  L <- dim(m2)[2]
  m2.row.sum <- cbind(m2,rowSums(m2))
  o2 <- rev(order(m2.row.sum[,L+1]))
  m2.row.sum <- m2.row.sum[o2,]
  return(m2.row.sum[,1:L])
}

getMeanDepthFromBam <- function(bam,get,minMapQuality=30,debug=F)
{
  if(debug){browser()}
  message(bam)
  param <- ApplyPileupsParam(which=get, 
                             what=c("seq", "qual"),
                             yieldBy="position",
                             yieldAll=TRUE,
                             minMapQuality=minMapQuality,
                             maxDepth=1000)
  
  fls <- PileupFiles(bam, param=param)
  calcInfo <- function(x)
  {
    if(debug){browser()}
    info <- apply(x[["seq"]], 2, function(y) {
      y <- y[c("A", "C", "G", "T"), , drop=FALSE]
      cvg <- colSums(y)
    })
    info
  }
  res <- applyPileups(fls, calcInfo, param=param)
  genecov <- t(do.call(cbind,res))
  return(list(mean_dp=mean(genecov),cov_vec=genecov))
}

# --------------------------------------------------------------------
grange_to_bed_file <- function(gr, file, name=NULL, offset=0)
{
  if(file.exists(file)){file.remove(file)}
  
  df <- data.frame(chr=as.character(seqnames(gr)),
                   start=start(gr)-1+offset,
                   end=end(gr)-offset,
                   anot=sprintf("%d_%s_%s",gr$test_id,gr$type,gr$gene_name))
  if(!is.null(name))
  {
    df$name <- values(gr)[,name]
  }
  #df <- df[df$chr %in% c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"),]
  write.table(df, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
}

parse_apa_bed_file <- function(mosdepth_region_fpath) {
  
  dt <- fread(cmd=sprintf("zcat %s",mosdepth_region_fpath),
              col.names=c('chr','start','end','annot','mean_depth'))
  
  tmp <- tstrsplit(dt$annot,'_',fixed=TRUE)
  dt$test_id <- as.numeric(tmp[[1]])
  dt$type <- tmp[[2]]
  dt$gene_name <- tmp[[3]]
  dt$annot <- NULL
  return(dt)
}

compPauRatioFromBam_mosdepth <- function(bam,apa_test_gr,minMapQuality=30, min_depth=0.,ncpu=4,mosdepth_bin="mosdepth", comp_ratio=T,debug=F)
{
  #apa_test_reg can be multiple regions
  # test_id chrom    start      end   type gene_name
  # 1:       1  chr4  8160406  8160455 shared    ABLIM2
  # 2:       1  chr4  8147865  8148132 unique    ABLIM2
  # 3:       2 chr16 89371613 89371752 shared   ANKRD11
  # 4:       2 chr16 89365283 89367370 unique   ANKRD11
  # 5:       3  chr2  9347148  9347349 shared     ASAP2
  # 6:       3  chr2  9347350  9349629 unique     ASAP2
  # 7:       4 chr12 58022496 58022686 shared  B4GALNT1
  # 8:       4 chr12 58022210 58022495 unique  B4GALNT1
  # 9:       5 chr12 58022496 58022686 shared  B4GALNT1
  # 10:       5 chr12 58019917 58020744 unique  B4GALNT1
  # 11:       6  chr8 41525785 41526082 shared      ANK1
  # 12:       6  chr8 41529872 41530430 shared      ANK1
  # 13:       6  chr8 41510744 41513272 unique      ANK1
  
  if(debug){browser()}
  message(bam)
  
  apa_test_bed_fpath <- paste0(bam,'.goi.bed')
  
  out_prefix <- paste0(bam,'.md')
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  if (!file.exists(region_bed_fpath)) {

    if (file.exists(apa_test_bed_fpath)) {
      unlink(apa_test_bed_fpath)
    }
    
    grange_to_bed_file(apa_test_gr,apa_test_bed_fpath)
  
    if (!file.exists(apa_test_bed_fpath)){
      stop('%s is not available and probably goi_gr to bed conversion is failed',apa_test_bed_fpath)
    }
    
    cmd <- sprintf('%s -t %d -b %s -n -Q %d %s %s',mosdepth_bin,ncpu,apa_test_bed_fpath,minMapQuality,out_prefix,bam)
    message(cmd)
    
    system(cmd)
  }
  
  bnd <- parse_apa_bed_file(region_bed_fpath)
  bnds <- split(bnd,by=c('test_id'),sorted=FALSE,drop=FALSE)
  
  eps2 <- 1e-7
  rm(bnd)
  
  pau_ratios <- sapply(1:length(bnds), function(i){
    reg <- bnds[[i]]
    
    width <- reg[type=='shared',end] - reg[type=='shared',start]+ 1
    shared_mdp0 <- sum(reg[type=='shared',mean_depth] * width) / sum(width)
    
    if (!is.finite(shared_mdp0)) {
      shared_mdp0 <- 0
    }
    shared_mdp <- shared_mdp0 + eps2 #shared
    
    width <- reg[type=='unique',end] - reg[type=='unique',start]+ 1
    unique_mdp <- sum(reg[type=='unique',mean_depth] * width) / sum(width)
    
    if (!is.finite(unique_mdp)) {
      unique_mdp <- 0
    }
    
    if (comp_ratio) {
      if (shared_mdp > 0) {
        paur <- c(unique_mdp/shared_mdp,shared_mdp,unique_mdp)
      } else {
        paur <- c(NA,shared_mdp,unique_mdp)
      }
    } else {
      paur <- c(-1.,shared_mdp,unique_mdp)
    }
    
    return(paur)
  })
  colnames(pau_ratios)<-names(bnds)
  rownames(pau_ratios) <-c('ratio','common','unique')
  message(sprintf("Done mean depth calculation on [%s]",bam))

  unlink(apa_test_bed_fpath)

  return(pau_ratios)
}

compPauRatioFromBam <- function(bam,apa_test_gr,minMapQuality=30, min_depth=0.,comp_ratio=T,debug=F)
{
  #apa_test_reg can be multiple regions
  # test_id chrom    start      end   type gene_name
  # 1:       1  chr4  8160406  8160455 shared    ABLIM2
  # 2:       1  chr4  8147865  8148132 unique    ABLIM2
  # 3:       2 chr16 89371613 89371752 shared   ANKRD11
  # 4:       2 chr16 89365283 89367370 unique   ANKRD11
  # 5:       3  chr2  9347148  9347349 shared     ASAP2
  # 6:       3  chr2  9347350  9349629 unique     ASAP2
  # 7:       4 chr12 58022496 58022686 shared  B4GALNT1
  # 8:       4 chr12 58022210 58022495 unique  B4GALNT1
  # 9:       5 chr12 58022496 58022686 shared  B4GALNT1
  # 10:       5 chr12 58019917 58020744 unique  B4GALNT1
  # 11:       6  chr8 41525785 41526082 shared      ANK1
  # 12:       6  chr8 41529872 41530430 shared      ANK1
  # 13:       6  chr8 41510744 41513272 unique      ANK1
  
  if(debug){browser()}
  message(bam)
  param <- ApplyPileupsParam(which=apa_test_gr,
                             what=c("seq", "qual"),
                             yieldBy="position",
                             yieldAll=TRUE,
                             minMapQuality=minMapQuality,
                             maxDepth=5000)
  
  fls <- PileupFiles(bam, param=param)
  
  calcInfo <- function(x) {
    info <- apply(x[["seq"]], 2, function(y) {
      y <- y[c("A", "C", "G", "T"), , drop=FALSE] #TODO: what about indels?
      cvg <- colSums(y)
    })
    info
  }
  res <- applyPileups(fls, calcInfo, param=param)
  
  if(debug){browser()}

  cvg_end <- cumsum(width(apa_test_gr))
  L <- length(apa_test_gr)

  apa_test_gr<-as.data.table(apa_test_gr)
  apa_test_gr$cvg_start <- c(1,cvg_end[1:(L-1)]+1)
  apa_test_gr$cvg_end <- cvg_end
  
  bnd <- apa_test_gr[,.(start=min(cvg_start),end=max(cvg_end)),by=c('test_id','type')]
  
  bnds <- split(bnd,by=c('test_id'),sorted=FALSE,drop=FALSE)
  
  genecov <- t(do.call(cbind,res))
  
  eps2 <- 1e-7
  rm(bnd)
  rm(res)
  
  pau_ratios <- sapply(1:length(bnds), function(i){
    reg <- bnds[[i]]
    shared_mdp0 <- mean(genecov[reg[type=='shared',start]:reg[type=='shared',end]])
    if (!is.finite(shared_mdp0)) {
      shared_mdp0 <- 0
    }
    shared_mdp <- shared_mdp0 + eps2 #shared
    
    unique_mdp <- mean(genecov[reg[type=='unique',start]:reg[type=='unique',end]]) #unique
    if (!is.finite(unique_mdp)) {
      unique_mdp <- 0
    }
    
    if (comp_ratio) {
      if (unique_mdp > shared_mdp) {
        if (unique_mdp >= min_depth){
          paur <- c(1.,shared_mdp,unique_mdp)
        } else {
          paur <- c(NA,shared_mdp,unique_mdp)
        }
      } else {
        if (shared_mdp >= min_depth){
          paur <- c(unique_mdp/shared_mdp,shared_mdp,unique_mdp)
        } else {
          paur <- c(NA,shared_mdp,unique_mdp)
        }
      }
    } else {
      paur <- c(-1.,shared_mdp,unique_mdp)
    }
    
    return(paur)
  })
  colnames(pau_ratios)<-names(bnds)
  rownames(pau_ratios) <-c('ratio','common','unique')
  message(sprintf("Done meand depth calculation on [%s]",bam))
  return(pau_ratios)
}

gamlss_general_with_tryCatch <- function(dat2) {
  fitspm <- tryCatch(
    {
      fitspm <- gamlss(x ~ y, data=dat2, family=NO)
    },
    error=function(cond){
      message('gamlss failed in handling input dat2!')
      return(NULL)
    },
    finally={
      message('gamlss processes input dat2!')
    }
  )
  return(fitspm)
}

gamlss_with_tryCatch <- function(dat2,mode=0) {
  fitspm <- tryCatch(
    {
      if (mode==1){
        fitspm <- gamlss(pau ~ 1, data=dat2)
      } else {
        fitspm <- gamlss(nbeta ~ pau, data=dat2, family=NO)
      }
    },
    error=function(cond){
      message(sprintf('gamlss failed in handling input dat2![%s]',cond))
      return(NA)
    },
    finally={
      message('gamlss processes input dat2!')
    }
  )
  return(fitspm)
}

regression_fit <- function(dat2,method='spline_monotonic',debug=False) {
  
  if (debug){browser()}
  
  if (method=='spline_monotonic'){
    fitspm_up <- gamlss(pau ~ pbm(nbeta,mono='up'), data=dat2)
    fitspm_dn <- gamlss(pau ~ pbm(nbeta,mono='down'), data=dat2)
    
    if (fitspm_up$sbc < fitspm_dn$sbc){
      fitspm <- fitspm_up
    } else {
      fitspm <- fitspm_dn
    }
  } else { #use a simple linear regression
    fitspm <- gamlss_with_tryCatch(dat2)
  }
  
  if (data.class(fitspm)!="gamlss") {
    regfit_info <- list(rsq=0.,pval=1.)
  } else {
    if (F){
      plot(pau ~ nbeta, data=dat2, col='lightblue')
      lines(fitted(fitspm)[order(dat2$nbeta)] ~ dat2$nbeta[order(dat2$nbeta)])
      summary(fitspm)
    }
  
    fitspmi <- list()
    fitspmi$R2 <- Rsq(fitspm)
    fitspmi$rss <- sum(residuals(fitspm, what='mu',type='simple')^2)
    
    fit0 <- gamlss_with_tryCatch(dat2,mode=1)
    if (data.class(fit0)!="gamlss") {
      regfit_info <- list(rsq=0.,pval=1.)
    } else {
      if (F){
        plot(pau ~ nbeta, data=dat2, col='lightblue')
        lines(fitted(fit0)~dat2$nbeta)
        summary(fit0)
      }
  
      fit0i <- list()
      fit0i$R2 <- Rsq(fit0)
      fit0i$rss <- sum(residuals(fit0, what='mu',type='simple')^2)
      
      n <- dim(dat2)[1]
      
      if (method=='spline_monotonic') {
        d <- round(edf(fitspm)[[1]][1])
        if (d==1) {d <- 2}
      } else {
        d <- 2
      }
      
      Fvalue <- (fit0i$rss-fitspmi$rss)/(d-1)/fitspmi$rss*(n-d)
      regfit_info <- list(rsq=fitspmi$R2,pval=pf(Fvalue,2,7,lower.tail=FALSE)[1])
    }
  }
  return(regfit_info)
}

regression_fit_general <- function(dat2,method='spline_monotonic',debug=False) {
  
  if (debug){browser()}
  
  if (method=='spline_monotonic'){
    fitspm_up <- gamlss(y ~ pbm(x,mono='up'), data=dat2)
    fitspm_dn <- gamlss(y ~ pbm(x,mono='down'), data=dat2)
    
    if (fitspm_up$sbc < fitspm_dn$sbc){
      fitspm <- fitspm_up
    } else {
      fitspm <- fitspm_dn
    }
  } else { #use a simple linear regression
    fitspm <- gamlss_general_with_tryCatch(dat2)
  }
  
  if (is.null(fitspm)) {
    regfit_info <- list(rsq=0.,pval=1.)
  } else {
    if (F){
      plot(y ~ x, data=dat2, col='lightblue')
      lines(fitted(fitspm)[order(dat2$x)] ~ dat2$x[order(dat2$x)])
      summary(fitspm)
    }
    
    fitspmi <- list()
    fitspmi$R2 <- Rsq(fitspm)
    fitspmi$rss <- sum(residuals(fitspm, what='mu',type='simple')^2)
    
    fit0 <- gamlss(y ~ 1, data=dat2)
    if (F){
      plot(y ~ x, data=dat2, col='lightblue')
      lines(fitted(fit0)~dat2$x)
      summary(fit0)
    }
    
    fit0i <- list()
    fit0i$R2 <- Rsq(fit0)
    fit0i$rss <- sum(residuals(fit0, what='mu',type='simple')^2)
    
    n <- dim(dat2)[1]
    
    if (method=='spline_monotonic') {
      d <- round(edf(fitspm)[[1]][1])
      if (d==1) {d <- 2}
    } else {
      d <- 2
    }
    
    Fvalue <- (fit0i$rss-fitspmi$rss)/(d-1)/fitspmi$rss*(n-d)
    regfit_info <- list(rsq=fitspmi$R2,pval=pf(Fvalue,2,7,lower.tail=FALSE)[1])
  }
  return(regfit_info)
}

pathoscope2b <- function(args,wks,action="script",debug2=1) {
	if (debug2==1){browser()}
	message('pathocope2 ...')
	# if (args$debug) {browser()}
	
	cmethod <- get_prog_info_yaml('pathoscope2')
	if (args$ps2_tgt_cfg != "") {
	  stopifnot(file.exists(args$ps2_tgt_cfg))
	  cmethod$config_path_opt <- args$ps2_tgt_cfg
	}
	message(cmethod)
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	message(sprintf("outputd:%s",args$outd))
	outd0 <- args$outd
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		if (debug2==1){browser()}
		if (isEmpty(wks$cmd[i])) {
			cmd <- ""
		} else {
			cmd <- sprintf("%s\n",wks$cmd[i])
		}
		
		reqcmd_str <- ""
		if (!is.null(cmethod$reqcmd)) {
			reqcmd_str <- comma_string_to_list(cmethod$reqcmd,sep=";")
			cmd <- sprintf("%s\n%s;python",cmd,reqcmd_str)
		} else {
			cmd <- sprintf("%s\npython",cmd)
		}
		
		cmd <- sprintf("%s %s",cmd,cmethod$bin_path)
		if (!is.null(cmethod$config_path_opt)) {cmd <- sprintf("%s -c %s",cmd,cmethod$config_path_opt)}
		if ("default_opt" %in% names(cmethod)) {
			if (!is.null(cmethod$default_opt)) {
				cmd <- sprintf("%s %s",cmd,cmethod$default_opt)
			}
		}
		
		cmd <- sprintf("%s -1 %s",cmd,wks$read1[i])
		cmd <- sprintf("%s -2 %s",cmd,wks$read2[i])
		
		outd <- file.path(normalizePath(outd0),wks$sample[i])
		message(sprintf("[%d]outputd:%s",i,outd))
		
		message(sprintf("creating [%s]...",outd))
		if (!file.exists(outd)) {
			dir.create(outd)
		}
		cmd <- sprintf("%s -o %s",cmd,outd)
		
		# cmd <- sprintf("%s --get_consensus_contig",cmd)
		
		if (args$runopt_pa2!="") {
			cmd <- sprintf("%s %s",cmd,args$runopt_pa2)
		}
		output_file <- file.path(outd,"ps2rank_n0_ti1_s100_top0.txt")
		
		if (file.exists(output_file)) {
			cmd <- sprintf("## reuse [%s]",output_file)
		}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		
		# browser()
		return(list(sample=wks$sample[i],output=output_file,cmd=cmd))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														output=sapply(1:M,function(m) {ret[[m]]$output}))
	
	if (debug2==1){browser()}
	
	return(cmethod)
}

pathoscope2 <- function(ps2_info,reads,outdir,readL=0,max_num_ti_report=0,debug=F) {
  if (debug) {browser()}
  
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  
  if (!file.exists(ps2_info$conf_path)) {
    stop(sprintf('check PS2 conf file [%s] exists!',ps2_info$conf_fpath))
  }
  
  cmd <- sprintf("python %s -c %s -r Illumina -a bt2 -t %d -o %s -p %d --adjreflen -b very-sensitive-local",ps2_info$prog_path,ps2_info$conf_path,max_num_ti_report,outdir,ps2_info$ncpu)
  
  if (readL > 0) {
    cmd <- sprintf("%s -L %d",cmd,readL)
  }
  
  if (length(reads)>1) {
    if (file.exists(reads[[1]]) & file.exists(reads[[2]])) {
      cmd <- sprintf('%s -1 %s -2 %s',cmd,reads[[1]],reads[[2]])
    } else {
      if (!debug) {
	      stop('check if both reads [%s] and [%s] exist!',reads[[1]],reads[[2]])
      }
    }
  } else {
    if (file.exists(reads[[1]])) {
      cmd <- sprintf('%s -u %s',cmd,reads[[1]])
    } else {
      if (!debug) {
      	stop('check if read [%s] exists!',reads[[1]])
      }
    }
  }
  message(cmd)
  if (!debug) {
  	system(cmd)
  }
}

pathoscope2_tiLen_db <- function() {
	cmethod <- get_prog_info_yaml('pathoscope2')
	cmethod$resource$header_path <- run_system(cmd=sprintf("ls %s",cmethod$resource$header_path))
	tiLen_rd <- sprintf("%s.rd",cmethod$resource$header_path)
	
	if (file.exists(tiLen_rd)) {
		load(tiLen_rd)
	} else {
		cmd_str <- sprintf("zcat %s | cut -d '|' -f1,3 | uniq",cmethod$resource$header_path)
		tiLen_dt <- fread(cmd=cmd_str,col.names = c('ti','tiLen'))
		tiLen_dt <- tiLen_dt[,lapply(.SD, function(x) {as.numeric(tstrsplit(x,":")[[2]])})]
		save(tiLen_dt,file=tiLen_rd,compress=T)
	}
	return(tiLen_dt)
}

pathoscope2report_to_dt <- function(wks) {
	
	#wks data.table(sample=sample_name_wo_space,
	#								output=ps2_report_file_path)
	
	rank.dts <- lapply(1:nrow(wks),function(r) {
		message(r)
		fread_between_two_lines(wks$output[r],from_word = '>pathogen_frequency',to_word = '>contigs',debug2=0)
	})
	rank.dt <- rbindlist(rank.dts[!is.na(rank.dts)])
	
	contig.dts <- lapply(1:nrow(wks),function(r) {
		message(r)
		a=fread_between_two_lines(wks$output[r],from_word = '>contigs',to_word = '>total_number_of_input_reads',debug2=0)
		ret<-NA
		if (!is.na(a)) {if (ncol(a)==4) {
			a$sample <- wks$sample[r]
			ret <- a
		}}
		ret
	})
	contig.dt <- rbindlist(contig.dts[!is.na(contig.dts)])
	contig.dt[,idx:=1:nrow(contig.dt)]
	
	list(rank=rank.dt,contig=contig.dt)
}

#this function will be obsolete and replaced by pathoscope2report_to_dt
pathoscope2_get_rankings <- function(ps2_report_file,python_path="python",sample=NA) {
  # browser()
  ps2_dir <- file.path(Sys.getenv("PPLINE"),"pathoscope")
  pybin2litereport <- file.path(ps2_dir,'pathoreport','postprocess_ps2_report_lite.py')
  stopifnot(file.exists(pybin2litereport))
  
  ps2_lite_tmp_fn <- paste0(ps2_report_file,'.tmp')
  cmd <- sprintf('%s %s -i %s -o %s',
  							 python_path,
                 pybin2litereport, 
                 ps2_report_file, 
                 ps2_lite_tmp_fn)
  
  if (!is.na(sample)) {
    cmd <- sprintf('%s -s %s',cmd,sample)
  }
  
  message(cmd)
  system(cmd)
  ps2_ranking_dt <- fread(ps2_lite_tmp_fn)
  unlink(ps2_lite_tmp_fn)
  return(ps2_ranking_dt)
}

pathoqc <- function(pathoqc_bin,raw_reads,read_base,outdir,ncpu=4,debug=F) {

	qc_read1 <- file.path(outdir,paste0(read_base,'_R1.fastq'))
	
	if (!file.exists(qc_read1)){
		
		if (length(raw_reads)>1) {
			qc_reads <- file.path(outdir,paste0(read_base,c('_R1.fastq','_R2.fastq')))
			cmd <- sprintf("python %s -1 %s -2 %s -t 33 -m 25 -e 50 -g 1 -a Y -a2 Y -d 0 -q 1 -p %d -o %s",
				pathoqc_bin,
				raw_reads[[1]],
				raw_reads[[2]],
				ncpu,
				outdir)
			message(cmd)
			if (!debug){
			  system(cmd)
			}
		} else {
				qc_reads <- file.path(outdir,paste0(read_base,'_R1.fastq'))
				cmd <- sprintf("python %s -1 %s -t 33 -m 25 -e 50 -g 1 -a Y -d 0 -q 1 -p %d -o %s",
				pathoqc_bin,
				raw_reads[[1]],
				ncpu,
				outdir)
				
				message(cmd)
				
				if (!debug) {
				  system(cmd)
				}
				
		}
	} else {
		if (length(raw_reads)>1) {
			qc_reads <- file.path(outdir,paste0(read_base,c('_R1.fastq','_R2.fastq')))
		} else {
			qc_reads <- file.path(outdir,paste0(read_base,'_R1.fastq'))
		}
	}
	return(qc_reads)
}

samtools_extract_by_bed <- function(bam,bed,out_bam,reuse=F,samtools_bin='samtools'){
  if (reuse & file.exists(out_bam)) {
    message(sprintf('reuse prev out_bam [%s] ...',out_bam))
  } else {
    cmd <- sprintf("%s view -O BAM -q 1 -F 4 -@ 4 -L %s -o %s %s",samtools_bin,bed,out_bam,bam)
    message(cmd)
    system(cmd)
    bai <- build_bai(out_bam)
  }
}


samtools_extract_by_contigs <- function(bam,contigs,out_bam,reuse=F,samtools_bin='samtools') {
  
  if (reuse & file.exists(out_bam)) {
    message(sprintf('reuse prev out_bam [%s] ...',out_bam))
  } else {
    cmd <- sprintf("%s view -h -f2 -@2 %s %s | %s sort -n -o %s -O BAM -@2 -",
                   samtools_bin,
                   bam,
                   paste0(contigs,collapse=" "),
                   samtools_bin,
                   out_bam)
    jmessage(cmd)
    system(cmd)
  }
}


bedtools_genomecov <- function(bam,out_file,reuse=F,bedtools_bin='bedtools') {
  if (reuse & file.exists(out_file)) {
    message(sprintf('reuse prev out_file [%s] ...',out_file))
  } else {
    cmd <- sprintf("%s genomecov -bga -split -ibam %s | gzip -fc > %s",bedtools_bin,bam,out_file)
    message(cmd)
    system(cmd)
  }
}

mosdepth_md_in_bed <- function(bam,roi_bed,mosdepth_bin,ncpu=1,minMapQuality=1) {
  
  out_prefix <- paste0(bam,'.md')
  
  if (!file.exists(roi_bed)){
    stop('%s is not available and probably goi_gr to bed conversion is failed',roi_bed)
  }
  
  cmd <- sprintf('%s -t %d -b %s -n -Q %d %s %s',mosdepth_bin,ncpu,roi_bed,minMapQuality,out_prefix,bam)
  message(cmd)
  system(cmd)
  
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  if (!file.exists(region_bed_fpath)) {
    stop('%s is not available and probably mosdepth run on [%s] failed',region_bed_fpath,bam)
  }
  
  bnd <- parse_apa_bed_file(region_bed_fpath)
  
  if (file.exists(roi_bed)) {
    unlink(roi_bed)
  }
  
}

mosdepth_md_by_window <- function(bam,out_prefix=NA,mosdepth_bin="mosdepth",window=30,ncpu=1,minMapQuality=1) {
	
	if (is.na(out_prefix)) {
		out_prefix <- paste0(bam,'.md')
	}

  cmd <- sprintf('%s -t %d -b %d -n -Q %d %s %s',mosdepth_bin,ncpu,window,minMapQuality,out_prefix,bam)
  message(cmd)
  system(cmd)
  
  region_bed_fpath <- sprintf("%s.regions.bed.gz",out_prefix)
  
  if (!file.exists(region_bed_fpath)) {
    stop('%s is not available and probably mosdepth run on [%s] failed',region_bed_fpath,bam)
  }
	return(region_bed_fpath)
}


get_mount_dir<-function(){
  mount_dir <- data.table(input='',
                          output='')
  
  hostname <- Sys.getenv('HOSTNAME')
  message(sprintf('hostname:%s',hostname))
  return(mount_dir)
}

get_limit_from_quantile <- function(x,target_r) {
  stopifnot(target_r>=0. & target_r<=1.)
  quantile(x,target_r)
}

build_minimap_ref_index <- function(ref,k=15,w=10) {
	stopifnot(file.exists(ref))
	ref_idx <- tag_file(ref,tag=sprintf('k%dw%d',k,w),new_ext='mmi')
	if (!file.exists(ref_idx)) {
		cmd <- sprintf("minimap2 -k%d -w%d -d %s %s",k,w,ref_idx,ref)
		run_system(cmd)
	}
	return(ref_idx)
}

run_minimap <- function(ref,r1,outd=NA,r2=NA,runopt=NA) {
  
  inpath <- decompose_fname(r1)
  if (is.na(outd)) {
    outd <- file.path(inpath$dir,'minimap2')
    dir.create(outd,recursive = TRUE, showWarnings = FALSE)
  }
  
  out_bam <- file.path(outd,sprintf('%s.bam',inpath$fbase0))
  
  cmd <- sprintf("minimap2 -aL -m200 %s %s | samtools sort -@ 4 -O bam -o %s -", ref, r1, out_bam)

  run_system(cmd)
  
  build_bai(out_bam,reuse=FALSE)
  return(out_bam)
}


catcmd_samtools_sort <- function(cmd,bam_prefix,tmpdir="") {
	
	output <- sprintf("%s_so.bam",bam_prefix)
	
	cmd <- sprintf("%s\nsamtools sort -o %s",cmd,output)
	
	if (tmpdir!="") {
		cmd <- sprintf("%s -T %s",cmd,tmpdir)
	}
	cmd <- sprintf("%s %s.bam\n\n",cmd,bam_prefix)

	cmd <- sprintf("%s\nsamtools index %s\n\n",cmd,output)

	return(list(cmd=cmd,output=output))
}


is_per_bam <- function(bam_fpath,samtools_bin="samtools") {
  # browser()
  cmd <- sprintf("%s view -h %s | head -n20000 | %s view -c -f 1 -",samtools_bin,bam_fpath,samtools_bin)
  jmessage(cmd)
  out <- system(cmd,intern = TRUE)
  per_count <- as.integer(out)
  if (per_count>0) {
    is_per=TRUE
  } else {
    is_per=FALSE
  }
  return(is_per)
}

is_valid_bam <- function(bam_fpath,samtools_bin="samtools") {
  
  cmd <- sprintf("%s quickcheck %s",samtools_bin,bam_fpath)
  jmessage(cmd)
  out <- system(cmd,intern = FALSE)
  valid_bam <- TRUE
  if (out != 0) {
    valid_bam <- FALSE
  }
  return(valid_bam)
}


bam_to_fastq <- function(per,bam,fqs,bamToFastq_bin="bedtools") {
  # browser()
  if (per) {
    cmd <- sprintf("%s bamtofastq -i %s -fq %s -fq2 %s",bamToFastq_bin,bam,fqs$R1,fqs$R2)
  } else {
    cmd <- sprintf("%s bamtofastq -i %s -fq %s",bamToFastq_bin,bam,fqs$R1)
  }
  jmessage(cmd)
  system(cmd)
  
}

unmapped_reads_from_bam <- function(bam_fpath,unmapped_dir,bamToFastq_bin="bedtools",per=TRUE) {
  
  unmapped_bam <- file.path(unmapped_dir,'unmapped.bam')
  
  cmd <- sprintf("samtools view -h -b -F0x102 -@2 %s | samtools sort -n -o %s -O BAM -@2",bam_fpath,unmapped_bam)
  if (!file.exists(unmapped_bam)) {
    jmessage(cmd)
    system(cmd)
  }
  
  read_base <- file.path(unmapped_dir,'unmapped_reads')
  
  fqs <- data.table(R1=sprintf('%s_R1.fastq',read_base),
                    R2=NA)
  if (per) {
    fqs$R2 <- sprintf('%s_R2.fastq',read_base)
  }
  
  bam_to_fastq(per,unmapped_bam,fqs,bamToFastq_bin = bamToFastq_bin)
  
  unlink(unmapped_bam)
  return(fqs)
}


get_hla_contigs_from_bam <- function(bam_fpath,
                                     fastq_dir,
                                     per=TRUE) {
  # browser()
  tmp_fn <- file.path(fastq_dir,"chr6_contigs.tmp")
  cmd <- sprintf("samtools view -H %s | grep -P \'(SN:chr6)|(SN:HLA)|(SN:6)\' > %s",bam_fpath,tmp_fn)
  jmessage(cmd)
  system(cmd)
  if (file.exists(tmp_fn)) {
    bam_part_head <- fread(tmp_fn,header=FALSE)
    hla_contigs <- tstrsplit(bam_part_head$V2,":")[[2]]
    unlink(tmp_fn)
  } else {
    hla_contigs <- NA
  }
  return(hla_contigs)
}


potential_hla_reads_from_bam <- function(bam_fpath,
                                         fastq_dir,
                                         bamToFastq_bin="bedtools",
                                         per=TRUE) {
  
  hla_contigs <- get_hla_contigs_from_bam(bam_fpath,fastq_dir,per=per)
  hla_fqs <- NA
  
  if (!is.na(hla_contigs)) {
    out_bam <- file.path(fastq_dir,'potential_hla.bam')
    samtools_extract_by_contigs(bam_fpath,hla_contigs,out_bam) 
    
    hla_fqs <- data.table(R1=file.path(fastq_dir,'hla_R1.fastq'),
                          R2=NA)
    if (per) {hla_fqs$R2 <- file.path(fastq_dir,'hla_R2.fastq')}
    
    bam_to_fastq(per,out_bam,hla_fqs,bamToFastq_bin=bamToFastq_bin)
    if (file.exists(out_bam)) {
      unlink(out_bam)
    }
  }
  return(hla_fqs)
}


samtools_consensus <- function(bamFn,refFn,refConsFq_out,samtools_dir=NA) {
	
	message('collecting consensus genomic segment sequences ...')
	
	cmd <- sprintf('bcftools mpileup -f %s %s | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq | gzip -fc > %s', refFn, bamFn, refConsFq_out)
	
	message(cmd)
	
	system(cmd)
	
}

numeric_wide_df_to_mat <- function(wide_df) {
	rowname2 <- wide_df[,1]
	M <- as.matrix(wide_df[,-1])
	rownames(M) <- rowname2
	M[is.na(M)] <- 0.
	return(M)
}

fastp <- function(args,wks,action="script",debug2=0) {
	if (debug2==1){browser()}
	message('fastp ...')
	cmethod <- get_prog_info_yaml('fastp')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		cmd <- cmethod$bin_path
		if (!is.null(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
		
		r1fpath <- wks$read1[i]
		r2fpath <- wks$read2[i]
		
		cmd <- sprintf("%s -i %s",cmd,r1fpath)
		cmd <- sprintf("%s -I %s",cmd,r2fpath)
		
		if (("outd" %in% colnames(wks))) {
			outd <- normalizePath(wks$outd[i])
		} else {
			stopifnot(args$outd!="")
			outd <- normalizePath(args$outd)
		}
		
		if (!file.exists(outd)) {
			message(sprintf("creating [%s]...",outd))
			dir.create(outd,recursive = T)
		}
		
		qcreadd <- file.path(outd,wks$sample[i])
		if (!file.exists(qcreadd)) {
			dir.create(qcreadd,recursive = T)
		}
		
		qr1 <- file.path(normalizePath(qcreadd),"qc_R1.fastq.gz")
		qr2 <- file.path(normalizePath(qcreadd),"qc_R2.fastq.gz")
		
		cmd <- sprintf("%s -o %s",cmd,qr1)
		cmd <- sprintf("%s -O %s",cmd,qr2)
		
		if (args$trim_polya=="Y"){cmd <- sprintf("%s -g",cmd)} #trim poly
		if (args$min_len>0) {
			cmd <- sprintf("%s -l %d",cmd,args$min_len) #min length required
		}
		
		if (args$min_complexity>0) {
			cmd <- sprintf("%s --low_complexity_filter",cmd)
			cmd <- sprintf("%s -Y %d",cmd,args$min_complexity) #low complexity filtering
		}
		
		cmd <- sprintf("%s -w %d",cmd,cmethod$ncpu) #multithread
		cmd <- sprintf("%s --json %s",cmd,file.path(normalizePath(qcreadd),"fastp.json"))
		cmd <- sprintf("%s --html %s",cmd,file.path(normalizePath(qcreadd),"fastp.html"))
		
		if (args$runopt_fastp!="") {
			cmd <- sprintf("%s %s",cmd,args$runopt_fastp)
		}
		
		if (file.exists(qr1) & args$reuse==1) {cmd <- sprintf("## reuse [%s]",qr1)}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		# browser()
		return(list(sample=wks$sample[i],read1=qr1,read2=qr2,cmd=cmd,outputd=qcreadd))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														outputd=sapply(1:M,function(m) {ret[[m]]$outputd}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

rnaseq_filto_mouse <- function(args,wks,action="script",debug2=0) {
	if (debug2==1){browser()}
	message('fastp ...')
	cmethod <- get_prog_info_yaml('hisat2')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	cmethod2 <- get_prog_info_yaml('rna_wo_mouse')
	hisat2_ss_site <- cmethod2$resource$hisat2_ss_site
	hisat2_refd <- cmethod2$resource$hisat2_refd
	sam_to_fastq_perl <- sprintf("perl %s",cmethod2$bin_path)
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		
		cmd <- cmethod$bin_path
		if (!is.null(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
		
		r1fpath <- wks$read1[i]
		r2fpath <- wks$read2[i]
		
		if (("outd" %in% colnames(wks))) {
			outd <- normalizePath(wks$outd[i])
		} else {
			stopifnot(args$outd!="")
			outd <- file.path(normalizePath(args$outd),wks$sample[i])
		}
		
		aligned <- file.path(outd,'hisat')
		if (!file.exists(aligned)) {
			message(sprintf("creating [%s]...",aligned))
			dir.create(aligned,recursive = T)
		}
		
		opref <- file.path(aligned,sprintf("%s_hisat",wks$sample[i]))
		# ------------
		cmd <- sprintf("%s -q --known-splicesite-infile %s -x %s --summary-file %s_summary.tsv --met-file %s_met.tsv -p %d -1 %s -2 %s -S %s.sam",
									 cmd,
									 hisat2_ss_site,
									 hisat2_refd,
									 opref,
									 opref,
									 args$ncpu,
									 r1fpath,
									 r2fpath,
									 opref)
		
		cmd <- sprintf("%s\nsamtools view -H %s.sam > %s_minus_reads_uniquely_mapped_to_mouse.sam",cmd,opref,opref)
		
		# cmd <- sprintf("%s\nsamtools view -H %s.sam > %s_reads_uniquely_mapped_to_mouse.sam",cmd,opref,opref)
		
		cmd <- sprintf("%s\nawk '$2 < 256' %s.sam | awk '$3 ~ /human.*/' >> %s_minus_reads_uniquely_mapped_to_mouse.sam",cmd,opref,opref)
		
		cmd <- sprintf("%s\nawk '$2 < 256' %s.sam | awk '$3 ~ /mouse.*/ && $21 !~ /^NH:i:1$/' >> %s_minus_reads_uniquely_mapped_to_mouse.sam",cmd,opref,opref)
		
		# cmd <- sprintf("%s\nawk '$2 < 256' %s.sam | awk '$3 ~ /mouse.*/ && $21 ~ /^NH:i:1$/' >> %s_reads_uniquely_mapped_to_mouse.sam",cmd,opref,opref)
		
		# cmd <- sprintf("%s\nsamtools view -bS %s_minus_reads_uniquely_mapped_to_mouse.sam > %s_minus_reads_uniquely_mapped_to_mouse.bam",cmd,opref,opref)
		# cmd <- sprintf("%s\nsamtools view -bS %s_reads_uniquely_mapped_to_mouse.sam > %s_reads_uniquely_mapped_to_mouse.bam",cmd,opref,opref)
		
		opref2 <- file.path(outd,wks$sample[i])
		cmd <- sprintf("%s\n%s %s_minus_reads_uniquely_mapped_to_mouse.sam %s",cmd,sam_to_fastq_perl,opref,opref2)
		
		cmd <- sprintf("%s\nrm -rf %s\n\n",cmd,aligned)
		
		human_R1_fastqz = sprintf("%s_R1.fastq.gz",opref2)
		human_R2_fastqz = sprintf("%s_R2.fastq.gz",opref2)
		if (file.exists(human_R1_fastqz) & file.exists(human_R2_fastqz) & args$reuse==1) {cmd <- sprintf("## reuse [%s]",human_R1_fastqz)}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		# browser()
		return(list(sample=wks$sample[i],read1=human_R1_fastqz,read2=human_R2_fastqz,cmd=cmd,outputd=opref2))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														outputd=sapply(1:M,function(m) {ret[[m]]$outputd}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

# rnabam2salmon <- function(args,wks,action,debug2=0) {
# 	if (debug2==1){browser()}
# 
# 	cmethod <- get_prog_info_yaml("salmon")
# 	
# 	cmethod <- set_user_option(cmethod,args)
# 	cmethod$action <- action
# 	
# 	#check_consistency(cmethod,wks)
# 	if ('out1' %in% colnames(wks)) {
# 		wks$in1 <- wks$out1
# 		wks$out1 <- NULL
# 	}
# 	
# 	if (args$outd=="") {
# 		outd <- "rnabam2salmon"
# 	} else {
# 		outd <- args$outd
# 	}
# 	if (!file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
# 	
# 	M <- dim(wks)[1]
# 	
# 	ret <- lapply(1:M,function(i) {
# 		# browser()
# 		sample <- wks$sample[i]
# 		bam <- wks$in1[i]
# 		outdi <- file.path(outd,sample)
# 		
# 		if (!file.exists(outdi)) {dir.create(outdi,showWarnings = F,recursive = T)}
# 		
# 		cmd <- sprintf("%s -a %s -o %s",cmethod$bin_path, bam, outdi)
# 		
# 		if (!isEmpty(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
# 		if (!isEmpty(cmethod$config_path_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$config_path_opt)}
# 
# 		# cmd <- sprintf("%s %s",cmd,args$runopt)
# 		
# 		message(cmd)
# 		
# 		# ./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant
# 		out1 <- file.path(outdi,"quant.sf")
# 		return(list(out1=out1,cmd=cmd))
# 	})
# 	
# 	wks$cmd <- as.data.table(sapply(1:M,function(m) {ret[[m]]$cmd}))
# 	wks$out1 <- as.data.table(sapply(1:M,function(m) {ret[[m]]$out1}))
# 	cmethod$wks <- wks
# 	if (debug2==1){browser()}
# 	return(cmethod)
# 	
# }

#handling multiple bam files
rnabam2featcnt <- function(args,wks,debug2=0) {
	if (debug2==1){browser()}
	
	if ('out1' %in% colnames(wks)) {
		wks$in1 <- wks$out1
		wks$out1 <- NULL
	}

	method_tag <-  "rnabam2featcnt"
	if (args$outd=="") {
		outd <- method_tag
	} else {
		outd <- args$outd
	}
	if (!file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	
	M <- dim(wks)[1]
	
	# load(file.path(outd,"rsubread_fc.rd")) #fc
	
	if (debug2==1){browser()}
	
	fc <- Rsubread::featureCounts(wks$in1,
													annot.ext = args$encode_gtf,
													isGTFAnnotationFile = TRUE,
													GTF.featureType = "exon",
													GTF.attrType="gene_id",
													isPairedEnd = TRUE,
													strandSpecific = 0, 
													# countChimericFragments=FALSE,
													fraction=TRUE,
													nthreads=args$ncpu)
	
	cmethod <- list(exptag=method_tag,
									wks=fc)
	
	#save(cmethods,file=file.path(outd,"rsubread_fc.rd"),compress=T)
	
	# wks$cmd <- as.data.table(sapply(1:M,function(m) {ret[[m]]$cmd}))
	# wks$out1 <- as.data.table(sapply(1:M,function(m) {ret[[m]]$out1}))
	# cmethod$wks <- wks
	if (debug2==1){browser()}
	return(cmethod)
	
}

#handling single bam
read_count <- function(bam,region2eval,
											 countMultiMappingReads=F,
											 allowMultiOverlap=F,
											 fraction=F,
											 ignoreDup=F) {
	
	R <- dim(region2eval)[1]
	
	message(sprintf('R=%d',R)) #debug
	
	fc <- Rsubread::featureCounts(bam,
																annot.ext = region2eval,
																countMultiMappingReads=countMultiMappingReads,
																allowMultiOverlap=allowMultiOverlap,
																fraction=fraction,
																ignoreDup=ignoreDup)
	
	cts <- matrix(NA, nrow = R, ncol = 1)
	
	rownames(cts) <- paste0(fc$annotation$Chr,':',
													fc$annotation$Start,'-',
													fc$annotation$End)
	
	colnames(cts) <- 'read_count'
	cts[,1] <- round(fc$counts)
	
	return(cts)
}

salmon2 <- function(args,wks,action="script",debug2=0) {
	if (debug2==1){browser()}
	message('salmon ...')
	cmethod <- get_prog_info_yaml('salmon')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	M <- dim(wks)[1]
	
	stopifnot(args$outd!="")
	outd0 <- normalizePath(args$outd)
	if (!file.exists(outd0)) {
		message(sprintf("creating [%s]...",outd0))
		dir.create(outd0,recursive = T)
	}
	
	ret <- lapply(1:M,function(i) {
		if (debug2==1){browser()}
		cmd <- sprintf("%s quant",cmethod$bin_path)
		if (!is.null(cmethod$config_path_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$config_path_opt)}
		if (!is.null(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
		
		if (args$inputd!="") {
			r1fpath <- file.path(args$inputd,wks$read1[i])
			r2fpath <- file.path(args$inputd,wks$read2[i])
		} else {
			r1fpath <- wks$read1[i]
			r2fpath <- wks$read2[i]
		}
		cmd <- sprintf("%s -1 %s",cmd,r1fpath)
		cmd <- sprintf("%s -2 %s",cmd,r2fpath)
		
		outd <- file.path(outd0,wks$sample[i])
		if (!file.exists(outd)) {
			message(sprintf("creating [%s]...",outd))
			dir.create(outd,recursive = T)
		}
		cmd <- sprintf("%s -o %s",cmd,outd)
		
		out1 <- file.path(outd,'quant.sf')
		
		if (args$runopt!="") {
			cmd <- sprintf("%s %s",cmd,args$runopt)
		}
		
		if (action == "run") {
			run_system(cmd)
		}
		# browser()
		return(list(sample=wks$sample[i],read1=r1fpath,read2=r2fpath,cmd=cmd,out1=out1))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														out1=sapply(1:M,function(m) {ret[[m]]$out1}))
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

cutadapt2 <- function(cmethod) {
	
	reads <- cmethod$input
	
	ret <- lapply(1:dim(reads)[1],function(i){
		
		cmd <- sprintf("%s\n\n",reads$cmd[i])
		cmd <- sprintf("%s%s",cmd,cmethod$bin_path)
		
		if (!is.null(cmethod$option1) & !is_empty(cmethod$option1)) {
			cmd <- sprintf("%s %s",cmd,cmethod$option1)
		}
		
		outd <- file.path(reads$dir[i],'cutadapt')
		dir.create(outd, recursive=TRUE, showWarnings = FALSE)
		
		out1 <- file.path(outd,basename(reads$r1[i]))
		cmd <- sprintf("%s -o %s",cmd,out1)
		if (cmethod$is_paired==1){
			out2 <- file.path(outd,basename(reads$r2[i]))
			cmd <- sprintf("%s -p %s",cmd,out2)
		} else {
			out2 <- NULL
		}
		
		cmd <- sprintf("%s %s",cmd,reads$r1[i])
		
		if (cmethod$is_paired==1){
			cmd <- sprintf("%s %s",cmd,reads$r2[i])
		}
		
		ret2 <- gen_cmd(cmethod,cmd,out1)
		cmd <- ret2$cmd
		status <- ret2$status
		
		return(list(out1,out2,status,cmd))
	})
	
	ret <- rbindlist(ret)
	reads$out1 <- ret$V1
	if (cmethod$is_paired==1) {
		reads$out2 <- ret$V2
	}
	reads$status <- ret$V3
	reads$cmd <- ret$V4
	cmethod$input <- reads
	
	return(cmethod)
}


norm_by_geometric_mean <- function(Cgs) {
	
	message('function of geometric norm ...')
	#Cgs <- sc.integs$`CAR-T.before_infusion`@assays$RNA@counts
	
	S <- dim(Cgs)[2]
	
	# Compute the geometric mean over the line
	
	message('compute geometric mean ...')
	gm.mean  <-  apply(Cgs, 1, function(x) {
		exp(sum(log(x+1.))/S)
	})
	
	Cgs <- sweep(Cgs, 1, gm.mean, FUN="/")
	
	message('computing scale factors ...')
	scale_factors <- apply(Cgs, 2, function(my_cell){median(my_cell[my_cell>0])})
	
	message('normalizing ...')
	Ngs <- sweep(Cgs, 2, scale_factors, FUN="/")
	
	return(Ngs)
	
}

tf_idf <- function(Ngs) {
	
	message('function of aggreated marker scores ...')
	
	totalReadCntPerCell <- Matrix::colSums(Ngs)
	totalReadCntPerGene <- Matrix::rowSums(Ngs)
	S <- dim(Ngs)[2]
	Tgs <- Ngs/totalReadCntPerCell * log(1 + S/totalReadCntPerGene)
	
	return(Tgs)
}

gatk_rnabam2vcf <- function(args,wks,action,debug2=0) {
	#NOTE: we use GATK4 here
	if (debug2==1){browser()}
	java <- get_prog_info_yaml('java')
	if (args$tmpdir!="") {
		tmpdir <- str_extract_all(java$default_opt,
															"(?<=tmpdir\\=)\\S+(?=\\s+)",
															simplify = TRUE)
		
		java$default_opt <- str_replace(java$default_opt,
																		sprintf("tmpdir=%s",tmpdir),
																		sprintf("tmpdir=%s",args$tmpdir))
	}
	
	java_cmd <- sprintf("%s %s",java$bin_path,java$default_opt)
	
	cmethod <- get_prog_info_yaml("gatk4")
	
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	#check_consistency(cmethod,wks)
	if ('out1' %in% colnames(wks)) {
		wks$in1 <- wks$out1
		wks$out1 <- NULL
	}
	
	if (args$outd=="") {
		outd <- "rnabam2vcf"
	} else {
		outd <- args$outd
	}
	if (!file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	
	star <- get_prog_info_yaml("star")
	ref_fpath <- star$resource$ref_path
	# check_consistency(cmethod,wks)
	M <- dim(wks)[1]
	ret <- lapply(1:M,function(i) {
		# browser()
		sample <- wks$sample[i]
		#'SplitNCigarReads------------------'
		if (wks$cmd[i]!="") {
			jcmd <- sprintf("%s\n\n%s %s",wks$cmd[i],java_cmd,cmethod$bin_path)
		} else {
			jcmd <- sprintf("%s %s",java_cmd,cmethod$bin_path)
		}
		
		cmd <- sprintf("%s SplitNCigarReads",jcmd)
		cmd <- sprintf("%s -R %s",cmd,ref_fpath)
		cmd <- sprintf("%s -I %s",cmd,wks$in1[i])
		
		if (args$runopt!="") {
			cmd <- sprintf("%s -L %s",cmd,args$runopt)
		}
		
		outprefix <- file.path(outd,sample)
		cmd <- sprintf("%s -O %s.split.bam\n\n",cmd,outprefix)
		
		#'HaploytypeCaller------------------'
		cmd <- sprintf("%s %s HaplotypeCaller",cmd,jcmd)
		cmd <- sprintf("%s -R %s",cmd,ref_fpath)
		cmd <- sprintf("%s -I %s.split.bam",cmd,outprefix)
		cmd <- sprintf("%s --dont-use-soft-clipped-bases -stand-call-conf 20.0",cmd)
		cmd <- sprintf("%s -O %s.split.vcf\n\n",cmd,outprefix)
		
		#'VariantFiltration------------------'
		cmd <- sprintf("%s %s VariantFiltration",cmd,jcmd)
		cmd <- sprintf("%s -R %s",cmd,ref_fpath)
		cmd <- sprintf("%s -V %s.split.vcf",cmd,outprefix)
		cmd <- sprintf("%s --filter-name DP5 --filter-expression 'DP < 5'",cmd)
		
		out1 <- sprintf("%s.split_dp5.vcf",outprefix)
		
		cmd <- sprintf("%s -O %s\n\n",cmd,out1)

		if (action == "run") {
			run_system(cmd)
		}
		if (debug2==1){browser()}
		return(list(out1=out1,cmd=cmd))
	})
	
	wks$cmd <- as.data.table(sapply(1:M,function(m) {ret[[m]]$cmd}))
	wks$out1 <- as.data.table(sapply(1:M,function(m) {ret[[m]]$out1}))
	cmethod$wks <- wks
	if (debug2==1){browser()}
	return(cmethod)
}

merge_svcfs_to_mvcf <- function(args,vcf_dt,unified_sname_in_svcf) {
	
	stopifnot(all(c('sname','vcf_fpath') %in% colnames(vcf_dt)))
	stopifnot(all(c('ncpu','ref_fa','reuse','outd') %in% colnames(args)))
	
	vcf_dt$usname_vcf_fpath <- unlist(mclapply(1:nrow(vcf_dt),function(r2) {
		out_vcf <- get_out_fpath(args,sprintf("%s.vcf",vcf_dt$sname[r2]),subd = 'varscan2_w_uname')
		vcf_assign_unique_sample_name(vcf_dt$vcf_fpath[r2],
																	unified_sname_in_svcf,
																	vcf_dt$sname[r2],
																	out_vcf,
																	reuse=args$reuse)
		out_vcf
	}
	,mc.cores = args$ncpu
	))
	
	merged_vcf <- get_out_fpath(args,'varscan2_normal.vcf')

	in_vcfs_str <- sprintf("--variant %s",paste0(vcf_dt$usname_vcf_fpath,collapse=" --variant "))
	
	# message(in_vcfs_str)
	
	if (args$reuse==0 | !file.exists(merged_vcf)) {
		m.gatk <- get_prog_info_yaml('gatk')
		m.java <- get_prog_info_yaml('java')
		
		gatk_cmd <- sprintf("%s -jar %s -T CombineVariants %s -R %s -nt %d -o %s",
												m.java$bin_path,m.gatk$bin_path,in_vcfs_str,args$ref_fa,args$ncpu,merged_vcf)
		# message(gatk_cmd)
		run_system(gatk_cmd)
	}
	merged_vcf
}

varscan2_vcall <- function(args,wks,action="script",debug2=0) {
	if (debug2==1){browser()}
	message('varscan2 ...')
	cmethod <- get_prog_info_yaml('varscan2')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	M <- dim(wks)[1]
	
	java <- get_prog_info_yaml('java')
	if (args$tmpdir!="") {
		tmpdir <- str_extract_all(java$default_opt,
															"(?<=tmpdir\\=)\\S+(?=\\s+)",
															simplify = TRUE)
		
		java$default_opt <- str_replace(java$default_opt,
																		sprintf("tmpdir=%s",tmpdir),
																		sprintf("tmpdir=%s",args$tmpdir))
	}
	
	java_cmd <- sprintf("%s %s",java$bin_path,java$default_opt)
	
	stopifnot(args$outd!="")
	outd <- normalizePath(args$outd)
	if (!file.exists(outd)) {
		message(sprintf("creating [%s]...",outd))
		dir.create(outd,recursive = T)
	}
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		
		out_vcf <- file.path(outd,sprintf("%s_varscan2.vcf",wks$sample[i]))
		
		cmd <- "samtools mpileup -d 1000 -q 15 -Q 15 -A"
		cmd <- sprintf("%s -f %s",cmd,args$ref)
		
		cmd <- sprintf("%s %s",cmd,wks$bam[i])
		
		cmd <- sprintf("%s | %s %s",cmd,java_cmd,cmethod$bin_path)
		cmd <- sprintf("%s mpileup2cns",cmd)
		cmd <- sprintf("%s - --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.02 --output-vcf 1",cmd)

		if (!is.null(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
		
		cmd <- sprintf("%s > %s",cmd,out_vcf)
		
		if (args$reuse==1 & file.exists(out_vcf)) {cmd <- sprintf("## reuse [%s]",out_vcf)}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		# browser()
		
		return(list(sample=wks$sample[i],input1=wks$bam[i],cmd=cmd,output1=out_vcf))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														input1=sapply(1:M,function(m) {ret[[m]]$input1}),
														output1=sapply(1:M,function(m) {ret[[m]]$output1}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

vcf_assign_unique_sample_name <- function(in_vcf,old_sname,new_sname,out_vcf,reuse=0) {
	
	stopifnot(file.exists(in_vcf))
	
	sed_cmd<-sprintf("sed 's/%s/%s/' %s > %s",old_sname,new_sname,in_vcf,out_vcf)
	if (reuse==0|!file.exists(out_vcf)) {
		run_system(sed_cmd)
	}
}
	
shrink_mtrx <- function(a,num_feats=60,debug2=0) {
	if (debug2==1) {browser()}
	if (dim(a)[1]>num_feats) {
		N2 <- num_feats + 5
		if (N2>dim(a)[1]) {N2 <- dim(a)[1]}
		feat_max <- apply(a,1,function(r){max(abs(r[!is.na(r)]))})
		ret <- get_limit_from_quantile(feat_max,(1 - N2/dim(a)[1]))
		a <- a[feat_max>=ret,]
	}
	return(a)
}

select_topk_from_each_sample <- function(L10.mtrx,num_feats=20,debug2=0,incl_feats='TIGIT',excl_feat_pat=NA) {
	if (debug2==1){browser()}
	
	if (!is.na(excl_feat_pat)) {
		L10.mtrx <- L10.mtrx[grep(excl_feat_pat,rownames(L10.mtrx),invert = T),]
	}
	feats <- rownames(L10.mtrx)
	
	N <- length(feats)
	if (N<num_feats){num_feats<-N}
	if (num_feats == 1) {
		return(L10.mtrx)
	} else {
		topk_feats <- apply(L10.mtrx,2,function(col2) {
			dt2 <-data.table(gene=feats,L10=col2)
			return(dt2[order(-L10),gene][1:num_feats])
		})
		g1s <- unique(sort(array(topk_feats)))
		if (!is.na(incl_feats)) {
			incl_feats <- incl_feats[incl_feats %in% feats]
			g1s <- unique(c(g1s,incl_feats))
		}
		L10.mtrx.shrinked <- as.matrix(L10.mtrx[g1s,])
		colnames(L10.mtrx.shrinked) <- colnames(L10.mtrx)
		return(L10.mtrx.shrinked)
	}
	if (debug2==1){browser()}
}

get_red_white_blue_for_heatmap <- function(a,paletteLength=50,debug2=0) {
	
	if (debug2==1) {browser()}
	
	def_min_neg <- -1e-3
	def_max_pos <- 1e-3
	
	myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	
	myBreaks <- c(seq(min(min(a[!is.na(a)]),def_min_neg), 
										0, 
										length.out=ceiling(paletteLength/2) + 1),
								seq(max(a[!is.na(a)])/paletteLength, 
										max(max(a[!is.na(a)]),def_max_pos), 
										length.out=floor(paletteLength/2)))
	
	dupindc <- duplicated(myBreaks)
	myBreaks <- myBreaks[!dupindc]

	return(list(colors=myColor,breaks=myBreaks))
}


get_red_white_blue_for_Heatmap <- function(a,midv=0,debug2=0) {
	
	if (debug2==1) {browser()}

	if (midv==0) {
		myMin <- min(min(a[!is.na(a)]),-1e-3)
		myMax <- max(max(a[!is.na(a)]),1e-3)
		
	} else {
		myMin <- min(min(a[!is.na(a)]))
		myMax <- max(max(a[!is.na(a)]))
	}
	return(circlize::colorRamp2(c(myMin,midv,myMax),c("blue", "white", "red")))
}


heatmap_reduced_features_logfc <- function(logfc_mtrx,
																					 max_feats=100,
																					 fig_label="heatmap",
																					 row_names_max_width=4,
																					 row_names_gp=12,
																					 cell_name="score",
																					 row_km = 1, 
																					 column_km = 1,
																					 debug2=0) {
	
	if (debug2==1){browser()}
	
	rownames(logfc_mtrx) <- sapply(rownames(logfc_mtrx),function(long_str) {strtrim(long_str, 75)})

	a2 <- select_topk_from_each_sample(abs(logfc_mtrx),num_feats = round(max_feats/dim(logfc_mtrx)[2]))
	
	a <- as.matrix(logfc_mtrx[rownames(a2),colnames(a2)])
	colnames(a) <- colnames(a2)
	
	p <- NA
	if (dim(a)[2]>0) {
		a[is.na(a)] <- 0.
		
		col_fun = get_red_white_blue_for_Heatmap(a)
		
		if (row_km>=dim(a)[1]) {
			row_km <- round(dim(a)[1]*0.4)
			if (row_km<0) {row_km=1}
		}
		
		if (column_km>=dim(a)[2]) {
			column_km <- round(dim(a)[2]*0.4)
			if (column_km<0) {column_km=1}
		}
		
		p <- as.ggplot(Heatmap(a,
													 col=col_fun,
													 column_title=fig_label,
													 name = cell_name,
													 row_km = row_km, 
													 column_km = column_km,
													 border = TRUE,
													 row_names_max_width = unit(row_names_max_width, "cm"),
													 row_names_gp = gpar(fontsize = row_names_gp)))
		
	}
	p
}

rankSumTest_tryCatch <- function(x,y,test="wilcox") {
	test_res <- tryCatch(
		{
			if (test=="wilcox") {
				wilcox.test(x,y)
			}
			else if (test=="ttest") {
				t.test(x,y)
			} else if (test=="kruskal") {
				kruskal.test(x,y)
			}
		},
		error=function(cond) {
			# message('wil.test failed')
			return(NA)
		},
		finally={
			# message('Done.')
		}
	)
	return(test_res)
}

wilcox_tryCatch <- function(x,y) {
	wilcox_res <- tryCatch(
		{
			wilcox.test(x,y)
		},
		error=function(cond){
			# message('wil.test failed')
			return(NA)
		},
		finally={
			# message('Done.')
		}
	)
	return(wilcox_res)
}


ttest_tryCatch <- function(y1,y2) {
	ttest_res <- tryCatch(
		{
			ttest_res <- t.test(y1,y2)
		},
		error=function(cond){
			# message('t.test failed')
			return(NA)
		},
		finally={
			# message('Done.')
		}
	)
	return(ttest_res)
}

read_gtf_in_dt <- function(ensembl_gtf_fpath) {
	
	rd_fpath <- sprintf("%s.gene.rd",ensembl_gtf_fpath)
	
	message(rd_fpath)
	
	if (file.exists(rd_fpath)) {
		message(sprintf("reuse the stored result[%s]",rd_fpath))
		load(rd_fpath)
	} else {
		if (endsWith(ensembl_gtf_fpath,'.gz')) {
			grep_cmd <- "zgrep"
		} else {
			grep_cmd <- "grep"
		}
		cmd <- sprintf("%s -Pv \'^##\' %s | awk \'$3 == \"gene\" {print $0}\'",grep_cmd,ensembl_gtf_fpath)
		message(cmd)
		gtf_dt <- fread(cmd=cmd)
		gtf_dt2 <- gtf_dt[,list(V1,V2,V3,V4,V5,V7)]
		colnames(gtf_dt2) <- c("chr","source","mtype","start","end","strand")
		
		gtf_dt2$ensembl_with_ver <- str_extract(gtf_dt$V9,"(?<=gene_id \")\\S+(?=\";)")
		gtf_dt2$ensembl <- tstrsplit(gtf_dt2$ensembl_with_ver,"\\.")[[1]]
		gtf_dt2$type <- str_extract(gtf_dt$V9,"(?<=gene_type \")\\S+(?=\";)")
		gtf_dt2$symbol <- str_extract(gtf_dt$V9,"(?<=gene_name \")\\S+(?=\";)")
		save(gtf_dt2,file=rd_fpath,compress = TRUE)
	}
	return(gtf_dt2)
}


mds_plot <- function(feat_by_sample,sample_dt,group.color,group.shape=NA,title2="mds_plot",xrange=NA,yrange=NA,debug2=0) {
	#feat_by_sample: assume that it is scaled across the features
	if (debug2==1){browser()}
	stopifnot(group.color %in% colnames(sample_dt))
	if (!is.na(group.shape)) {
		stopifnot(group.shape %in% colnames(sample_dt))
	}
	message("computing distance ...")
	d <- dist(t(feat_by_sample))
	
	message("performing MDS ...")
	fit <- cmdscale(d, eig = TRUE, k = 2)
	x <- fit$points[,1]
	y <- fit$points[,2]
	mds_pca <- data.frame(C1=x,C2=y)
	rownames(mds_pca) <- colnames(feat_by_sample)
	mds_pca <- cbind(mds_pca,sample_dt)
	mds_pca[[group.color]] <- as.factor(mds_pca[[group.color]])
	if (is.na(group.shape)) {
		p <- ggplot(data=mds_pca,aes_string(x="C1",y="C2",color=group.color))
	} else {
		mds_pca[[group.shape]] <- as.factor(mds_pca[[group.shape]])
		p <- ggplot(data=mds_pca,aes_string(x="C1",y="C2",color=group.color,shape=group.shape))
	}
	
	# p <- p + scale_shape_manual(values=c(20,6,3,4))
	p <- p + geom_point(alpha=0.4,size=2)
	p <- p + ggtitle(sprintf("%s;%s;%s",title2,group.color,group.shape))
	p <- p + theme_bw()
	
	if (!is.na(xrange)) {p <- p + xlim(xrange)}
	if (!is.na(yrange)) {p <- p + ylim(yrange)}
	
	if (debug2==1){browser()}
	return(list(df=mds_pca,fig_ptr=p))
}

umap_plot <- function(sample_by_feat,sample_dt,group.color,group.shape=NA,title2="umap_plot",n_neighbors=12,debug2=0) {
	
	if (debug2==1){browser()}
	
	#feat_by_sample: assume that it is scaled across the features
	M <- dim(sample_by_feat)[1]
	if (M < n_neighbors) {n_neighbors <- round(M * 0.3)}
	
	message(sprintf("target n_neighbors[%d]",n_neighbors))
	umap_out <- umap(sample_by_feat, n_neighbors=n_neighbors, init="spca",scale=TRUE)
	umap_out <- data.frame(umap_out)
	colnames(umap_out) <- c('C1','C2')
	rownames(umap_out) <- rownames(sample_by_feat)
	
	umap_out <- cbind(umap_out,sample_dt)
	umap_out[[group.color]] <- as.factor(umap_out[[group.color]])
	if (is.na(group.shape)) {
		p <- ggplot(data=umap_out,aes_string(x="C1",y="C2",color=group.color))
	} else {
		umap_out[[group.shape]] <- as.factor(umap_out[[group.shape]])
		p <- ggplot(data=umap_out,aes_string(x="C1",y="C2",color=group.color,shape=group.shape))
	}
	#scale_shape_manual(values=c(3, 16, 17)) : to change point shapes 

	p <- p + geom_point(alpha=0.4,size=2)
	p <- p + ggtitle(sprintf("%s;%s;%s",title2,group.color,group.shape))
	if (debug2==1){browser()}
	return(list(df=umap_out,fig_ptr=p))
}


hlahd <- function(args,wks,action="script",debug2=1) {
	if (debug2==1){browser()}
	progname="hlahd"
	message(sprintf('prep %s cmd...',progname))
	cmethod <- get_prog_info_yaml(progname)
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	outd <- normalizePath(args$outd)
	if (!file.exists(outd)) {
		message(sprintf("creating [%s]...",outd))
		dir.create(outd,recursive = T)
	}
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		
		if (args$inputd!="") {
			r1fpath <- file.path(args$inputd,wks$read1[i])
			r2fpath <- file.path(args$inputd,wks$read2[i])
		} else {
			r1fpath <- wks$read1[i]
			r2fpath <- wks$read2[i]
		}
		
		cmd <-""
		if (F) {
		to_delete <- list()
		
		if (endsWith(r1fpath,'.gz')) {
			r1_unzipped <- gsub("\\.gz$","",r1fpath)
			cmd <- sprintf("gunzip -fc %s > %s",r1fpath,r1_unzipped)
			r1fpath <- r1_unzipped
			to_delete[[1]] <- r1_unzipped
			
			r2_unzipped <- gsub("\\.gz$","",r2fpath)
			cmd <- sprintf("%s\ngunzip -fc %s > %s",cmd,r2fpath,r2_unzipped)
			r2fpath <- r2_unzipped
			to_delete[[2]] <- r2_unzipped
		}
		}
		cmd <- sprintf("%s\n%s",cmd,cmethod$bin_path)
		
		if (!is.null(cmethod$default_opt)) {cmd <- sprintf("%s %s",cmd,cmethod$default_opt)}
		
		cmd <- sprintf("%s -t %d",cmd,cmethod$ncpu) #multithread
		cmd <- sprintf("%s -m %d",cmd,args$min_rlen) #read shorter than this will be ignored
		cmd <- sprintf("%s -f %s",cmd,cmethod$resource$freqdb) #Use information of allele frequencies from hlahd db
		if (args$runopt!="") {
			cmd <- sprintf("%s %s",cmd,args$runopt)
		}
		
		cmd <- sprintf("%s %s",cmd,r1fpath)
		cmd <- sprintf("%s %s",cmd,r2fpath)
		
		cmd <- sprintf("%s %s",cmd,cmethod$resource$split)
		cmd <- sprintf("%s %s",cmd,cmethod$resource$dict)
		cmd <- sprintf("%s %s",cmd,wks$sample[i])
		cmd <- sprintf("%s %s",cmd,outd)
		
		if (F) {
		if (length(to_delete)>0) {
			for (j in 1:length(to_delete)) {
				cmd <- sprintf("%s\nrm -rf %s",cmd,to_delete[[j]])
			}
		}
		}
		
		result_fn <- file.path(outd,wks$sample[i],"result",sprintf("%s_final.result.txt",wks$sample[i]))
		if (debug2==1){browser()}
		if (file.exists(result_fn)) {cmd <- sprintf("## reuse [%s]",result_fn)}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		# browser()
		return(list(sample=wks$sample[i],read1=r1fpath,read2=r2fpath,cmd=cmd,output=result_fn))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														output=sapply(1:M,function(m) {ret[[m]]$output}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}

get_zscore <- function(feat_by_sample_mtx,margin1=1,debug2=0) {
	# same as "ret = t(scale(feat_by_sample_mtx,center=TRUE,scale=TRUE))" but it can hanlde NA/NAN
	if (debug2==1) {browser()}
	zscore2 <- apply(feat_by_sample_mtx,margin1,function(c2) {
		if (debug2==1) {browser()}
		(c2-mean(c2, na.rm=TRUE))/sd(c2, na.rm=TRUE)
	})
	if (margin1==1) {
		zscore2 = t(zscore2)
	}
	zscore2[is.nan(zscore2)] <- 0.
	zscore2
}

set_NA_to_mean <- function(sample_by_feat,margin1=1,debug=0) {
	sample_by_feat <- apply(sample_by_feat,margin1,function(r2) {
		if (debug==1){browser()}
		nz_mean <- mean(r2[!is.na(r2)])
		r2[is.na(r2)]<-nz_mean
		r2
	})
	if (margin1==1) {
		ret = t(sample_by_feat)
	} else {
		ret = sample_by_feat
	}
	ret
}


heatmap_feat_by_cluster_avgexprscaled <- function(normAvgExpr, pdf_fpath, width2=8, height2=6, title2 = "heatmap avgExpr scaled", debug2=0) {
	if (debug2==1){browser()}
	
	normAvgExpr_scaled <- t(apply(normAvgExpr,1,scale))
	colnames(normAvgExpr_scaled) <- colnames(normAvgExpr)
	
	clrArr <- get_red_white_blue_for_heatmap(normAvgExpr_scaled)
	
	# message(pdf_fpath)
	# pdf(file=pdf_fpath,width=width2,height=height2)
	
	pheatmap(normAvgExpr_scaled,color = clrArr$colors,breaks = clrArr$breaks,main = title2)
	
}

build_blat_index <- function(ref_fa_fpath,faToTwoBit_bin) {
	ref_2bit <- tag_file(ref_fa_fpath,new_ext = '2bit')
	if (!file.exists(ref_2bit)) {
		cmd <- sprintf("%s %s %s",faToTwoBit_bin,ref_fa_fpath,ref_2bit)
		run_system(cmd)
	}
	ref_2bit
}

pblat <- function(query_fa,ref2bit,pblat_bin,out_fpath,query_type="dna",ncpu=4) {
	cmd <- sprintf("%s %s %s -q=%s -out=psl -threads=%d -noHead %s",pblat_bin,ref2bit,query_fa,query_type,ncpu,out_fpath)
	run_system(cmd)
	
	out_fpath
}

purecn_mod_seqnames <- function(res,is_chr=TRUE) {
	centromeres=res$input$centromeres
	centromeres=renameSeqlevels(centromeres, sub("chrM","chrMT",seqlevels(centromeres)))
	res$input$centromeres=renameSeqlevels(centromeres, sub("chr", "", seqlevels(centromeres)))
	res
}

wrap_plot_grid1 <- function(plst,nrow1=3,ncol1=3,title2="my plots",pdf_fn=NA,width=8,height=6) {
	plist2b<-chunk_by_unitlen(plst,unit_len = nrow1*ncol1)
	
	if (!is.na(pdf_fpath)) {pdf(pdf_fpath,width=width,height=height)}
	
	lapply(plist2b,function(plist2bj){
		p<-wrap_plots(plist2bj,nrow=nrow1,ncol=ncol1)
		p<- p & theme(legend.position = "bottom") #shared legend in patchwork
		p <- p + plot_layout(guides = "collect")
		p <- p + plot_annotation(title2)
		plot(p)
		NA
	})
	if (!is.na(pdf_fpath)) {dev.off()}
}

ncbi_blast <- function(args,contig_wks,action="script",debug2=0) {
	if (debug2==1){browser()}
	
	progname="blast"
	message(sprintf('prep %s cmd...',progname))
	cmethod <- get_prog_info_yaml(progname)
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	outd <- create_dir_if_not_exist(normalizePath(args$outd))
	
	by_samples <- split(contig_wks,by="sample")
	
	queryd <-create_dir_if_not_exist(file.path(outd,"query_fasta"))
	mapped <-create_dir_if_not_exist(file.path(outd,"sam"))
	
	ret <- lapply(names(by_samples),function(sname) {
		
		bys<-by_samples[[sname]]
		
		faObjs <- lapply(1:nrow(bys),function(r) {
			myseq <- s2c(string(bys$contig[r]))
			header1=sprintf("%d_%s",bys$idx[r],sname)
			as.SeqFastadna(myseq, name = header1, Annot=header1)
		})
		names(faObjs) <- sprintf("%d_%s",bys$idx,sname)
		
		query_fa <- file.path(queryd,sprintf("%s.fasta",sname))
		
		write.fasta(faObjs,names=names(faObjs),file.out=query_fa,as.string = FALSE)
		
		sam_fn <- file.path(mapped,sprintf("%s.sam",sname))
		
		cmd_str <- sprintf("%s -db %s -query %s -num_threads %d -out %s",cmethod$bin_path,args$ntdb,query_fa,args$ncpu,sam_fn)
		
		if (cmethod$run_opt!=""){
			cmd_str <- sprintf("%s %s",cmd_str,cmethod$run_opt)
		}
		
		if ((!"default_opt" %in% colnames(args))) {args$default_opt = ""}
		if (args$default_opt!="") {
			cmd_str <- sprintf("%s %s",cmd_str,args$default_opt)
		} else if (cmethod$default_opt!="") {
			cmd_str <- sprintf("%s %s",cmd_str,cmethod$default_opt)
		}
		
		cmd_str <- sprintf("%s\ngzip -f %s",cmd_str,sam_fn)
		message(cmd_str)
		
		return(list(sample=sname,
								input1=query_fa,
								output1=sprintf("%s.gz",sam_fn),
								cmd=cmd_str))
	})
	
	M = length(ret)
	message(sprintf("total number of cmds generated[%d]",M))
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														input1=sapply(1:M,function(m) {ret[[m]]$input1}),
														output=sapply(1:M,function(m) {ret[[m]]$output1}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}


convert_genelist_to_gmtfile <- function(member_genes,gmt_fpath,shared_url="https://www") {
	fpw<-file(gmt_fpath,open = 'w')
	imap(member_genes,function(genes,gs_name){
		line2 = sprintf("%s\t%s\t%s",gs_name,shared_url,paste0(genes,collapse='\t'))
		writeLines(line2,fpw)
	})
	close(fpw)
}

gsmtx_to_gctfile <- function(gsmtx,gct_fpath) {
	# gsmtx <- matus$`PanCK+`$`Non-Tumor`$rc$adjcnt
	dim_str <- paste0(dim(gsmtx),collapse="\t")
	genes <- rownames(gsmtx)
	gct.dt <- data.table(NAME=genes,Description=genes)
	
	gexpr.dt<-as.data.table(gsmtx)
	gct.dt = cbind(gct.dt,gexpr.dt)
	
	fpw <- file(gct_fpath,open = 'w')
	writeLines(c("#1.2",dim_str),fpw)
	close(fpw)
	message(gct_fpath)
	fwrite(gct.dt,file=gct_fpath,sep="\t",append=TRUE,col.names = TRUE)
}