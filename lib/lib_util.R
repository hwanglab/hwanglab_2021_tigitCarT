#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(tools)

jmessage <- function(script=NA,func='main',msg='msg1',action='info') {
	now = Sys.time()
	if (action == 'error') {
		message(sprintf("[%s:%s]%s:%s|%s ...",now,action,script,func,msg))
		stop(msg2display)
	} else {
		message(sprintf("[%s:%s]%s:%s|%s ...",now,action,script,func,msg))
	}
}

print_table <- function(u_title,u_fontSize,ms3a) {
	
	title = textGrob(u_title,gp=gpar(fontsize=12))
	padding = unit(5,'mm')
	table = gtable_add_rows(tableGrob(ms3a,rows=NULL,theme=ttheme_default(base_size=u_fontSize)),
													heights = grobHeight(title) + padding,
													pos=0)
	
	table = gtable_add_grob(table,
													title,
													1,1,1,ncol(table),
													clip = 'off')
	grid.newpage()
	grid.draw(table)
}

file_filter_rows <- function(fn,colidx1,key) {
	fn_head = sprintf('%s.head',fn)
	system(sprintf('head -n1 %s > %s',fn,fn_head))
	fn_tmp = sprintf('%s.tmp',fn)
	system(sprintf('awk \'$%d == \"%s\" {print $0}\' %s > %s',colidx1,key,fn,fn_tmp))
	system(sprintf('cat %s >> %s',fn_tmp,fn_head))
	unlink(fn_tmp)
	fn_head
}

get_max_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  maxVals <- rep(-1e7,M)
  for (i in 1:M){
    maxVals[i] <- max(my_list[[i]],na.rm = T)
  }
  return(max(maxVals))
}


get_min_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  minVals <- rep(1e7,M)
  for (i in 1:M){
    minVals[i] <- min(my_list[[i]],na.rm = T)
  }
  return(min(minVals))
}

get_mount_dir <- function(){
  hostname <- Sys.getenv('HOSTNAME')
  if (hostname == 'lri-107577'){
    mount_prefix <- '/home/hongc2/orange_passport/apa_atingLab2019'
  } else {
    mount_prefix <- file.path(Sys.getenv('HOME'),'projects','apa_atingLab2019')
  }
  return(mount_prefix)
}

comma1k <- function(my_integers){
  
  my_integer_dc <- data.class(my_integers[1])
  
  if (my_integer_dc=='numeric'){
    char_integers_with_comma <- format(my_integers, big.mark=",", scientific=FALSE)
  } else if (my_integer_dc == 'character') {
    char_integers_with_comma <- format(as.numeric(my_integers), big.mark=",", scientific=FALSE)
  } else {
    stop(sprintf('cannot support the data type [%s]',my_integer_dc))
  }
  return(char_integers_with_comma)
}

wo_bgnd_ggplot <- function() {
	wo_bgnd <- theme_bw() + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(),
											 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
	return(wo_bgnd)
}

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
# 
# b37_to_hg19_bam <- function(bam){
# 	chr_bam <- sprintf('%s.chr.bam',bam)
# 	cmd <- sprintf("samtools view -H %s | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | grep -vP '(SN:GL)|(SN:hs)' | samtools reheader - %s > %s",bam,bam,chr_bam)
# 	system(cmd)
# 	return(chr_bam)
# }

decompose_fname <- function(filepath2) {
  
  fbase <- basename(filepath2)
  if (endsWith(fbase,'.gz')) {
    fbase <- sub(".gz","",fbase)
  }
  
  #finally strip off the original file extension
  fbase0 <- tools::file_path_sans_ext(fbase)
  
  dirsegs <- tstrsplit(dirname(filepath2),.Platform$file.sep)
  D <- length(dirsegs)
  
  fname <- data.table(dir=dirname(filepath2),
  										parentd=dirsegs[[D]],
                      fbase0=file_path_sans_ext(fbase),
                      ext=file_ext(fbase))
  
  return(fname)
}

get_file_base <- function(filepath2) {
  # message(filepath2)
  fbase <- basename(filepath2)
  if (endsWith(fbase,'.gz')) {
    fbase <- sub(".gz","",fbase)
  }
  
  #finally strip off the original file extension
  fbase0 <- tools::file_path_sans_ext(fbase)
  return(fbase0)
}

tag_file <- function(filepath2,tag=NA,new_ext=NA){
  
  message(sprintf('in:%s',filepath2))
  fname <- decompose_fname(filepath2)
  if (is.na(tag) & is.na(new_ext)) {
  	stop('assign a value to either tag or new_ext!')
  }
  
  if (is.na(new_ext)) {
    tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s_%s.%s',fname$fbase0,tag,fname$ext)))
  } else {
  	if (is.na(tag)) {
  		tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s.%s',fname$fbase0,new_ext)))
  	} else {
  		tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s_%s.%s',fname$fbase0,tag,new_ext)))
  	}
  }
  
  message(sprintf('in:%s',tagged_fname))
  return(tagged_fname)
}

page_to_print <- function(range_str) {
  
  range2 <- strsplit(range_str,',')[[1]]
  range2a <- strsplit(range2,'-')
  range0 <- list()
  for (i in length(range2a)) {
    range3 <- range2a[[i]]
    range4 <- as.numeric(unlist(range3))
    if (length(range4)==2) {
      range0[[i]] <- range4[1]:range4[2]
    } else {
      range0[[i]] <- range4
    }
  }
  page_to_print <- unique(unlist(range0))
  return(page_to_print)
}

annotate_reads <- function(read_dt,is_paired) {
  read_dt$sample_label1 <- sapply(read_dt$r1,function(read_fpath) {
    fname <- decompose_fname(read_fpath)
    return(substr(fname$fbase0, 1, nchar(fname$fbase0)-3))
  })
  if (is_paired==1) {
    read_dt$sample_label2 <- sapply(read_dt$r2,function(read_fpath) {
      fname <- decompose_fname(read_fpath)
      return(substr(fname$fbase0, 1, nchar(fname$fbase0)-3))
    })
  }
  
  read_dt$dir1 <- sapply(read_dt$r1,function(read_fpath) {
    fname <- decompose_fname(read_fpath)
    return(fname$dir)
  })
  
  if (is_paired==1) {
    read_dt$dir2 <- sapply(read_dt$r2,function(read_fpath) {
      fname <- decompose_fname(read_fpath)
      return(fname$dir)
    })
  }
  
  if (is_paired==1) {
    is_valid <- all(read_dt[,sample_label1==sample_label2]) & all(read_dt[,dir1==dir2])
    if (is_valid) {
      read_dt$dir <- read_dt$dir1
      read_dt$sample_label <- read_dt$sample_label1
    } else {
      read_dt$dir <- NA
      read_dt$sample_label <- NA
    }
    read_dt[,sample_label1:=NULL]
    read_dt[,sample_label2:=NULL]
    read_dt[,dir1:=NULL]
    read_dt[,dir2:=NULL]
  } else {
    is_valid <- TRUE
    read_dt$dir <- read_dt$dir1
    read_dt$sample_label <- read_dt$sample_label1
  }
  
  return(list(is_valid=is_valid,read_dt=read_dt))
}

get_reads_from_prefix <- function(dir_pref,is_paired=1,suffix='fastq.gz') {
  message(sprintf('%s/*%s',dir_pref,suffix))
  ret <- sort(Sys.glob(sprintf('%s*/*%s',dir_pref,suffix)))
  reads <- list()
  if (is_paired==1) {
    reads[[1]] <- ret[endsWith(ret,'R1.fastq.gz')]
    reads[[2]] <- ret[endsWith(ret,'R2.fastq.gz')]
    
    read_dt <- data.table(r1=reads[[1]],r2=reads[[2]])
    ret <- annotate_reads(read_dt,is_paired)
    stopifnot(ret$is_valid)
    read_dt <- ret$read_dt
  } else {
    read_dt <- data.table(r1=ret)
    ret <- annotate_reads(read_dt,is_paired)
    read_dt <- ret$read_dt
  }
  return(read_dt)
}


run_system <- function(cmd,syslog_fn=NA,intern=TRUE,verbose=0,debug2=0) {
  if (debug2==1){browser()}
	if (verbose>0) {message(cmd)}
	
  sys_out <- system(cmd,intern=intern)
  
  if (intern & 'Error' %in% sys_out) {
    stop('error occurs while runing %s\nRefer %s',cmd,sys_out)
  }
  if (!is.na(syslog_fn)) {
    fwrite(sys_out,file=syslog_fn)
  }
  return(sys_out)
}


comma_string_to_list <- function(comma_string,sort_flag="",sep=",",debug2=0) {
	if (debug2==1) {browser()}
	items <- unique(str_trim(tstrsplit(comma_string,sep)))
	if (sort_flag=="ascending") {
		items <- sort(items,decreasing = FALSE)
	} else if (sort_flag=="descending") {
		items <- sort(items,decreasing = TRUE)
	}
	if (debug2==1) {browser()}
	return(items)
}

try2 <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(structure(msg, class = "try-error"))
  })
}

concat_txt_files <- function(text_files,out_fn) {
  cmd <- sprintf("cat %s | gzip -fc > %s",paste0(text_files,sep=" ",collapse=""),out_fn)
  system(cmd)
}

copy_to_local <- function(remote_fpath,local_cached) {
  
	if (!file.exists(local_cached)) {dir.create(local_cached,showWarnings = F,recursive = T)}
  local_cpy <- file.path(local_cached,basename(remote_fpath))
  cmd <- sprintf("rsync -ravP %s %s/",remote_fpath,local_cached)
  system(cmd)
  message(sprintf("copied [%s] to [%s]",remote_fpath,local_cpy))
	
  return(local_cpy)
}

chunk2 <- function(my_list,num_batches=3) {
	split(my_list, rep_len(1:num_batches, length(my_list)))
}

chunk_by_unitlen <- function(list2split,unit_len=10) {
	A <- length(list2split)
	if (A <= unit_len) {
		mbatches <- list(list2split)
	} else {
		num_batches <- ceiling(A/unit_len)
		b <- lapply(1:num_batches, function(j) {
			rep(j,unit_len)
		})
		f <- unlist(b)
		f <- f[1:A]
		mbatches <- split(list2split, f)
	}
	message(sprintf("total [%d] batches w/ [%d] elements!",length(mbatches),unit_len))
	mbatches
}


dt_to_xlsx <- function(input_dt,col1,out_fbase,sort_field="padj",debug2=0) {
	if(debug2==1){browser()}
	library(openxlsx)
	
	dts1 <- split(input_dt[order(get(col1)),],by=col1)
	
	wb <- createWorkbook(col1)
	for (name1 in names(dts1)) {
		sheet_name <- name1
		addWorksheet(wb,sheetName = sheet_name)
		freezePane(wb, sheet=sheet_name,firstRow = TRUE, firstCol = TRUE)
		writeData(wb,sheet_name,dts1[[name1]][order(get(sort_field))],withFilter=TRUE)
	}
	saveWorkbook(wb, file=sprintf("%s.xlsx",out_fbase), overwrite = TRUE)
	if(debug2==1){browser()}
}

matdf_to_scina_marker_fmt <- function(matdf_tsv,out_csv) {
	dt <- fread(matdf_tsv)
	N <- dim(dt)[2]
	genes <- dt$marker
	dt$marker <- NULL
	M <- max(colSums(dt))
	
	markers2 <- apply(dt,2,function(mycol) {
		genes2 <- genes[mycol>0]
		M2 <- M - length(genes2)
		return(c(genes2,rep(NA,M2)))
	})
	
	fwrite(as.data.table(markers2),file=out_csv,sep=',')
}


fread_between_two_lines <- function(txt_file,from_word,to_word=NA,debug2=0) {
	if (debug2==1){browser()}
	dt2 <- NA
	if (file.exists(txt_file)) {
		#check if 'from_word' exists
		cmd<-sprintf("grep '%s' %s",from_word,txt_file)
		exit_status <- run_system(cmd,intern = FALSE)
		if (exit_status==0) {
			cmd <- sprintf("sed -n '/%s/=' %s",from_word,txt_file)
			line_str <- run_system(cmd)
			message(line_str)
			line_from <- as.numeric(line_str)
			
			if (!is.na(to_word)) {
				cmd<-sprintf("grep '%s' %s",to_word,txt_file)
				exit_status <- run_system(cmd,intern = FALSE)
				if (exit_status==1) {
					to_word <- NA
				}
			}
			
			if (is.na(to_word)) {
				dt <- fread(file=txt_file,skip=line_from-1)
			} else {
				cmd <- sprintf("sed -n '/%s/=' %s",to_word,txt_file)
				
				line_str <- run_system(cmd)
				message(line_str)
				line_to <- as.numeric(line_str) - 1
				# browser()
				dt <- fread(file=txt_file,skip=line_from-1,nrows=(line_to-line_from))
			}
			dt2<-as.data.table(dt)
		}
	}
	return(dt2)
}

rbindlist2 <- function(res_stat) {
	if (all(is.na(res_stat)|is.null(res_stat))) {
		res_stats <- NA
	} else {
		res_stats <- rbindlist(res_stat[!is.null(res_stat) & !is.na(res_stat)])
	}
	return(res_stats)
}



cut_nice <- function(x, lower = 0, upper, by,
										 sep = "-", above.char = "+",debug2=0) {
	
	if(debug2==1){browser()}
	labs <- c(paste(sprintf("%3.4f",seq(lower, upper - by, by = by)),
									sprintf("%3.4f",seq(lower + by, upper, by = by)),
									sep = sep),
						paste(sprintf("%3.3f",upper), above.char, sep = ""))
	
	cut(x, breaks = c(seq(lower, upper, by = by), Inf),
			right = TRUE, labels = labs)
}

get_current_script_fpath <- function() {
	
	if (rstudioapi::isAvailable()) {
		script_fpath <- rstudioapi::getSourceEditorContext()$path
	} else {
		args = commandArgs()
		script_fpath = args[substr(args,1,7) == '--file=']
		script_fpath <- substr(script_fpath, 8, nchar(script_fpath))
	}
	
	script_fpath <- normalizePath(script_fpath)
	
	message(script_fpath)
	
	fpath_loc <- data.table(fpath=script_fpath)
	fpath_loc <- cbind(fpath_loc,decompose_fname(script_fpath))
	
	return(fpath_loc)
}

create_outdirectory <- function(args,script_loc) {
	outd<-file.path(args$baseoutd,script_loc$fbase0)
	if (!file.exists(outd)){dir.create(outd,showWarnings = F,recursive = T)}
	message("=====================")
	message(sprintf("output directory [%s]",outd))
	message("=====================")
	return(outd)
}

check_fpath <- function(fpath) {
	stopifnot(file.exists(fpath))
}

create_dir_if_not_exist <- function(newd,active=T) {
	if (!file.exists(newd)){
		if (active==T) {
			message(sprintf("creating %s ...",newd))
			dir.create(newd,showWarnings = F,recursive = T)
		}
	}
	message(newd)
	newd
}


get_outd <- function(outd0,script_fpath,create=1,debug2=0) {
	if (debug2==1){browser()}
	decfn <- decompose_fname(script_fpath)

	outd<-file.path(outd0,decfn$fbase0)
	
	if (create==1 & !file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	message(sprintf("output directory[%s]",outd))
	outd
}

get_wkd <- function(projd,script_fpath,create=1) {
	
	decfn <- decompose_fname(script_fpath)
	
	wkd<-file.path(projd,decfn$fbase0)
	
	if (create==1 & !file.exists(wkd)) {dir.create(wkd,showWarnings = F,recursive = T)}
	message(sprintf("output directory[%s]",wkd))
	wkd
}

get_pipeline_step <- function(script_fpath) {
	decfn <- decompose_fname(script_fpath)
	step_name <- tstrsplit(decfn$fbase0,"_")[[1]]
	message(sprintf("step_name[%s]",step_name))
	step_name
}

get_out_fpath <- function(args,fname,subd=NA,verbose=0) {
	outd <- args$outd
	if (!is.na(subd)) {
		outd <- file.path(args$outd,subd)
	}
	
	if (!file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	
	out_fpath <- file.path(outd,sprintf("%s_%s",args$step,fname))
	if (verbose>0) {
		message(sprintf("prep/writing %s ...",out_fpath))
	}
	message(out_fpath)
	out_fpath
}

get_proj_out0 <- function(proj_name) {
	file.path(Sys.getenv('HOME'),'projects',proj_name)
}

read_rd_or_rds <- function(rds_fpath) {
	message(sprintf("loading %s",rds_fpath))
	if (endsWith(rds_fpath,'.rd')) {
		ret <- get(load(rds_fpath))
	} else if (endsWith(rds_fpath,'.rds')) {
		ret <- readRDS(rds_fpath)
	} else {
		stop(sprintf("cannot support the file format %s",rds_fpath))
	}
	ret
}

check_prog_in_path <- function(bin_fpath) {
	cmd_str=sprintf("which %s",bin_fpath)
	ret <- run_system(cmd_str,intern = FALSE)
	bin_in_path = TRUE
	if (ret==1){
		bin_in_path=FALSE
	}
	bin_in_path
}
