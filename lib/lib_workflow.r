#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 2/1/2019

library(data.table)
library(stringr)
library(yaml)
library(configr)
library(httr)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))

get_project_out_dir <- function(unique_project_name) {
	hostname <- Sys.getenv('HOSTNAME')
	message(sprintf("hostname:%s",hostname))
	
	HOMED <- Sys.getenv('HOME')
	
	if (hostname == "cc-dclrilog52") {
		proj_outd<-file.path(HOMED,"lustre","projects","out",unique_project_name)
	} else {
		proj_outd<-file.path(HOMED,"projects","out",unique_project_name)
	}
	
	if (!file.exists(proj_outd)) {dir.create(proj_outd,recursive = T)}
	message(sprintf("project base output directory[%s]",proj_outd))
	return(proj_outd)
}

get_mount_dir <- function() {
  hostname <- Sys.getenv('HOSTNAME')
  message(sprintf("hostname:%s",hostname))
  username <- Sys.getenv('USER')
  message(sprintf("username:%s",username))
  
  mount_dir <- "<hidden>"
}

get_prog_config_filename <- function() {
  hsuffix <- Sys.getenv('HOSTNAME')

  if (is_empty(hsuffix)) {
    hsuffix <-'lri-gen-p-04'
  }
  prog_config_fname <- sprintf('rworkflow_%s.ini',hsuffix)
  return(prog_config_fname)
}

get_prog_info <- function(method,is_paired=1,reuse=0,debug=0,action='script') {
  if (debug==1) {browser()}

  #prog_config_fname <- get_prog_config_filename()
	rwkf.ini <- Sys.getenv('PROG_CONF')
	
	message(sprintf('prog_config[%s]',rwkf.ini))
	stopifnot(file.exists(rwkf.ini))

  prog <- list()
  prog$method <- method
  prog$exptag <- method
  prog$user_email <- ''
  prog$bin_path <- eval.config(file = rwkf.ini, config = method, value = "bin_path")
  prog$option1 <- eval.config(file = rwkf.ini, config = method, value = "default_opt")
  prog$option2 <- eval.config(file = rwkf.ini, config = method, value = "config_path_opt")
  prog$reqcmd <- eval.config(file = rwkf.ini, config = method, value = "reqcmd")
  prog$preamble <- eval.config(file = rwkf.ini, config = method, value = "preamble")
  prog$reuse <- reuse
  prog$debug <- debug
  prog$action <- action
  prog$is_paired <- is_paired
  return(prog)
}

set_user_option <- function(cmethod,args){
	
	opts <- names(args)
	
	message(opts)
	
	stopifnot('exptag' %in% opts)
	stopifnot('user_email' %in% opts)
	stopifnot('outd' %in% opts)

	if (!('reuse' %in% opts)) {args$reuse=0}
	if (!('is_paired' %in% opts)) {args$is_paired=1}
	if (!('ncpu' %in% opts)) {args$ncpu=2}
	if (!('hours' %in% opts)) {args$hours=24}
	if (!('reuse' %in% opts)) {args$reuse=0}
	if (!('memG' %in% opts)) {args$memG=48}
	if (!('run_opt' %in% opts)) {args$run_opt=""}
	
	cmethod$exptag = args$exptag
	cmethod$user_email = args$user_email
	cmethod$outd = args$outd
	cmethod$reuse=  args$reuse
	cmethod$is_paired = args$is_paired
	cmethod$ncpu = args$ncpu
	cmethod$hours = args$hours
	cmethod$reuse = args$reuse
	cmethod$memG = args$memG
	cmethod$run_opt = args$run_opt
	
	return(cmethod)
}

locate_prog_config <- function(config.ini=NA) {
  
  if (is.na(config.ini)) {
  	#mount_prefix <- get_mount_dir()
    if (Sys.getenv("PROG_CONF")=="") {
  	  conf_file <- file.path(Sys.getenv("PPLINE"),'configs','local',sprintf('hwanglab_prog_conf_%s.yml',Sys.getenv("HOSTNAME")))
    } else {
      conf_file <- Sys.getenv("PROG_CONF")
    }
    if (length(conf_file)==0) {
      stop('check if [%s] exists and is valid',conf_file)
    }
  } else {
    conf_file <- config.ini
  }
  
  if (!file.exists(conf_file)) {
    stop(sprintf('check if [%s] exists !',conf_file))
  }
  return(conf_file)
}

load_prog_config <- function(prog_conf_file=NA) {
	
	prog_conf_file <- locate_prog_config()
	message(sprintf('prog_config[%s]',prog_conf_file))
	prog_yaml <- read_yaml(prog_conf_file)

	return(prog_yaml)
}

prep_update_worksheet <- function(in_wksheet_file,out_wksheet_file) {
  
  jmessage('prepare updated worksheet template ...')
  cmd <- sprintf("grep -vP '^#' %s > %s",in_wksheet_file,out_wksheet_file)
  message(cmd)
  system(cmd)
}

replace_env_vars <- function(str1) {
  # browser()
  str2 <- gsub(sprintf("\\$\\{%s\\}","APP"),Sys.getenv('APP'),str1)
  str2 <- gsub(sprintf("\\$\\{%s\\}","LABSHARED"),Sys.getenv('LABSHARED'),str2)
  str2 <- gsub(sprintf("\\$\\{%s\\}","TMP"),Sys.getenv('TMP'),str2)
  
  return(str2)
}

get_prog_info_yaml <- function(method,debug2=0) {
	if (debug2==1) {browser()}
	prog_yaml <- load_prog_config()
	stopifnot(method %in% names(prog_yaml))
	cmethod <- prog_yaml[[method]]
	cmethod$required_input_columns <- comma_string_to_list(replace_env_vars(cmethod$required_input_columns))
	cmethod$optional_input_columns <- comma_string_to_list(replace_env_vars(cmethod$optional_input_columns))
	
	if (!is.null(cmethod[['bin_path']])) {
		if (!startsWith(cmethod$bin_path,'/')) {
			subd_componets <- tstrsplit(replace_env_vars(cmethod$bin_path),'/')
			if (length(subd_componets)>1) {
				cmethod$bin_path <- replace_env_vars(cmethod$bin_path)
			}
		}
	}
	
	if (!is.null(cmethod[['resource']])) {
		for (entry in names(cmethod$resource)) {
			if (!is.null(cmethod$resource[[entry]])){
				if (endsWith(entry,'_path') & !startsWith(cmethod$resource[[entry]],'/')) {
					subd_componets <- tstrsplit(replace_env_vars(cmethod$resource[[entry]]),'/')
					if (length(subd_componets)>1) {
						cmethod$resource[[entry]] <- replace_env_vars(cmethod$resource[[entry]])
					}
				}
			}
		}
	}
	return(cmethod)
}

load_method_in_wksheet <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -P '^#' %s | cut -f2,3,4,5",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=F, col.names = c('step', 'item','value','condition'))
  return(sdt)
}


load_method_in_wksheet_old <- function(worksheet_file,stepi) {
  # browser()
  cmd_str <- sprintf("grep -P '^#' %s | cut -f2,3,4,5 | grep -w '^%d'",worksheet_file,stepi)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=F, col.names = c('step', 'item','value','condition'))
  return(sdt)
}

load_samples_in_wksheet <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -vP '^#' %s",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=T)
  return(sdt)
}

check_consistency <- function(cmethod,worksheet) {
	stopifnot(all(cmethod$required_input_columns %in% colnames(worksheet)))
}

check_consistency_in_wksheet <- function(mdt,sdt) {
  
  # browser()
  
  mdt_req <- mdt[condition=='required' & is.na(value),]
  
  if (dim(mdt_req)[1]>0) {
    stop('empty value in a required field in the worksheet file [method section]')
  }
  
  sample_cols <- colnames(sdt)
  mdt_sreq <- mdt[condition=='s.required',]
  
  for (i in 1:dim(mdt_sreq)[1]) {
    
    sreq_row <- mdt_sreq[i,]
    vals <- strsplit(sreq_row$value,';')[[1]]
    if (length(vals) > 0) {
      
      for (j in 1:length(vals)) {
        if (vals[j] %in% sample_cols) {
          message(sprintf('column[%s] verified',vals[j]))
        } else {
          stop(sprintf('check if column[%s] exist in the worksheet file',
                       vals[j]))
        }
      }
    } else {
      stop(sprintf('empty value at the row[%s] in the worksheet file:method',
                   sreq_row$item))
    }
  }
  
}

bean_leslie_project_ps2 <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -P '^#' %s | cut -f2,3,4",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=F, col.names = c('item','value','condition'))
  return(sdt)
}

bean_leslie_project_samples <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -vP '^#' %s",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=T)
  return(sdt)
}

bean_leslie_project_check_consistency <- function(mdt,sdt) {
  
  sample_cols <- colnames(sdt)
  mdt_sreq <- mdt[condition=='s.required',]
  
  for (i in 1:dim(mdt_sreq)[1]) {
    sreq_row <- mdt_sreq[i,]
    vals <- strsplit(sreq_row$value,';')[[1]]
    if (length(vals) > 0) {
      
      for (j in 1:length(vals)) {
        if (vals[j] %in% sample_cols) {
          message(sprintf('column[%s] verified',vals[j]))
        } else {
          stop(sprintf('check if column[%s] exist in the worksheet file',
                       vals[j]))
        }
      }
    } else {
      stop(sprintf('empty value at the row[%s] in the worksheet file:method',
                   sreq_row$item))
    }
  }
  
}

launch_pipeline <- function(runinst,sample,inputs,output_dir,debug=0,reuse=0) {
  # browser()
  if (debug==0) {
    message(output_dir)
    if (!dir.exists(output_dir)){
      dir.create(output_dir,recursive=TRUE)
    }
  } else {
    # browser()
    x<-1
  }
  
  if (endsWith(runinst$prog_path,'.py')){
    cmd <- sprintf("python %s",runinst$prog_path)
  } else if (endsWith(runinst$prog_path,'.r')) {
    cmd <- sprintf("Rscript %s",runinst$prog_path)
  } else {
    cmd <- runinst$prog_path
  }
  
  outputs <- list()
  
  if (!is.null(runinst$common_opt)) {
    cmd <- sprintf("%s %s",cmd,runinst$common_opt)
  }
  
  message(sprintf('generating a commandline[%s]',runinst$rstep$method))
  
  if (runinst$rstep$method %in% c('pathoqc','pathoscope2','salmon')) {
    I <- length(inputs)
    cmd <- sprintf("%s -1 %s",cmd,inputs[[1]])
    if (I==2) {
      cmd <- sprintf("%s -2 %s",cmd,inputs[[2]])
    }
    cmd <- sprintf("%s -o %s",cmd,output_dir)
    
  } else if (runinst$rstep$method == 'fastqc') {
    cmd <- sprintf("%s %s",cmd,paste(inputs,collapse = ' '))
    cmd <- sprintf("%s -o %s",cmd,output_dir)
  } else if (runinst$rstep$method == 'cutadapt') {
    
    cmd <- sprintf("%s -o %s",cmd,file.path(output_dir,basename(inputs[[1]])))
    if (length(inputs)==2) {
      cmd <- sprintf("%s -p %s",cmd,file.path(output_dir,basename(inputs[[2]])))
    }
    cmd <- sprintf("%s %s",cmd,paste(inputs,collapse = ' '))
    
  } else if (runinst$rstep$method == 'cellranger') {
    run_dir <- dirname(inputs[[1]])
    cmd <- sprintf("%s --id=%s",cmd,sample)
    cmd <- sprintf("%s --run=%s",cmd,run_dir)
    cmd <- sprintf("%s --samplesheet=%s",cmd,inputs[[1]])
    cmd <- sprintf("%s --output-dir %s",cmd,output_dir)
    
  } else {
    stop(sprintf('method [%s] is not supported yet',runinst$rstep$method))
  }
  
  if (!is.null(runinst$add_opt)) {
    cmd <- sprintf("%s %s",cmd,runinst$add_opt)
  }
  
  outfpat <- runinst$rstep$mdt[item=='output_file_suffix',value]
  fpats <- strsplit(outfpat,';')[[1]]
  
  for (j in 1:length(fpats)) {
    outfpat_n <- fpats[j]
    if (startsWith(outfpat_n,'*')) { #check if this is for a suffix type
      outputs[[j]] <- file.path(output_dir,outfpat_n)
      message(outputs[[j]])
    } else {
      outputs[[j]] <- outfpat
    }
  }
  
  log_file <- file.path(output_dir,sprintf("%s.log",sample))
  if (F) {
    cmd <- sprintf("%s 2>&1 | tee %s",cmd,log_file)
  }
  
  message(cmd)
  if (debug==0) {
    out_fpath <- file.path(output_dir,outputs[[1]])
    message(sprintf('expected output file[%s]',out_fpath))
    if (reuse==1 & file.exists(out_fpath)) {
      message('reuse prev result ...')
    } else {
      #message(sprintf('************check <%s>**************',out_fpath))
      system(cmd) #debug
    }
  }
  return(outputs)
}

gen_cmd <- function(cmethod,cmd,out1) {
  if (cmethod$action=="run") {
    if (cmethod$reuse==1) {
      if (file.exists(out1)) {
        message(sprintf('The expected output file [%s] already exists and skip running',out1))
      } else {
        message(cmd)
        system(cmd)
      }
      cmd <- ""
      status <- 'successful'
    }
  } else {
    status <- cmethod$action
    if (file.exists(out1)) {
      cmd <- ""
      status <- 'successful'
    }
  }
  return(list(cmd=cmd,status=status))
}


run_module <- function(cmethod,args) {
  # browser()
  if (cmethod$method == "cutadapt") {
    cmethod <- cutadapt2(cmethod)
  } else if (cmethod$method == "salmon") {
    cmethod <- salmon2(cmethod)
  } else {
    message('not impleted yet...')
  }
  return(cmethod)
  
}

prep_pipeline_shortreads <- function(method,args,action) {
  
	if (args$debug==1){browser()}
  cmethod <- get_prog_info(method,
                           args$is_paired,
                           reuse=args$reuse,
                           debug=args$debug,
                           action=action)

  cmethod <- set_user_runoption(args$runopt,cmethod)
  
  cmethod$input <- get_reads_from_prefix(args$in_pat,
                                         args$is_paired,
                                         args$suffix)
  
  cmethod$input$cmd <- ""
  if (args$debug==1){browser()}
  return(cmethod)
}

salmon_quant <- function(args,wks,action,debug2=args$debug) { #prev_cmethod,method,action
	if (debug2==1){browser()}
	
	#$SALMON quant -i $SALMON_IDX -l IU -p 20 -1 $SALMON_IN/$mMDSC_1_p1 -2 $SALMON_IN/$mMDSC_1_p2 --numBootstraps 100 -o $SALMON_OUT/mMDSC_1
	#
	message('salmon_quant ...')
	cmethod <- get_prog_info_yaml('salmon')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	refidx <- cmethod$config_path_opt
	if (args$refidx!="") {refidx <- args$refidx}
	
	runopt <- cmethod$default_opt
	if (args$runopt!="") {runopt <- args$runopt}
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		cmd <- cmethod$bin_path
		cmd <- sprintf("%s quant",cmd)
		
		cmd <- sprintf("%s -i %s",cmd,refidx)
		
		if (runopt!="") {
			cmd <- sprintf("%s %s",cmd,runopt)
		}
		
		cmd <- sprintf("%s -p %d",cmd,cmethod$ncpu)
		
		r1fpath <- wks$read1[i]
		r2fpath <- wks$read2[i]
		
		cmd <- sprintf("%s -1 %s",cmd,r1fpath)
		cmd <- sprintf("%s -2 %s",cmd,r2fpath)
		
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
		
		quant_outd <- file.path(outd,'salmon_quant',wks$sample[i])
		if (!file.exists(quant_outd)) {
			dir.create(quant_outd,recursive = T)
		}
		
		cmd <- sprintf("%s -o %s",cmd,quant_outd)
		
		out1 <- file.path(quant_outd,"quant.sf")
		
		if (cmethod$reuse==1 & file.exists(out1)) {
			cmd <- sprintf("## reuse [%s]",out1)
		}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
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

parse_run_option <- function(arg_runopt) {
  runopts <- strsplit(arg_runopt,';')[[1]]
  runopt_list <- lapply(runopts,function(runopt){
    items <- strsplit(runopt,':')[[1]]
    method <- items[[1]]
    user_opt <- trimws(items[[2]])
    return(list(method,user_opt))
    })
  
  runopt_dt <- rbindlist(runopt_list)
  colnames(runopt_dt) <- c('method','option1')
  return(runopt_dt)
}

set_user_runoption <- function(user_runopt_str,cmethod) {
	if (user_runopt_str!="") {
	  user_opt_dt <- parse_run_option(user_runopt_str)
	  option1 <- user_opt_dt[method==cmethod$method,option1]
	  if (!is_empty(option1)) {
	    cmethod$option1 <- option1
	  }
	}
  return(cmethod)
}

sbatch_generator <- function(cmethod,nscripts,scriptd,target_hostname="",debug2=0) {
  
  if (debug2>0) {browser()}
	
	W <- dim(cmethod$wks)[1]
	if (nscripts > W){
		nscripts <- W
	}
  cmethod$wks$sbatch_id <- (1:W)%%nscripts
  inputs <- split(cmethod$wks,by=c('sbatch_id'))

  for (j in 1:nscripts) {
  	if (debug2>1) {browser()}
    batch_script_tag <- sprintf("%d_%s",j,cmethod$exptag)
    script_fpath <- file.path(scriptd,sprintf('%s.sh',batch_script_tag))
    
    logdir <- file.path(scriptd,'logs')
    dir.create(logdir,showWarnings = FALSE, recursive = TRUE)
    
    message(sprintf('generating a script file [%s] ...',script_fpath))
    fp <- file(script_fpath,'w')
    
    cmd_str <- sprintf("#!/bin/bash -l")
    if (cmethod$user_email!="") {
      cmd_str <- sprintf("%s
#SBATCH -J %s
#SBATCH -o %s/%s.%%N.%%j.out
#SBATCH -e %s/%s.%%N.%%j.err
#SBATCH -p defq
#SBATCH -c %d
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%s
#SBATCH --mem=%dG",
                         cmd_str,
                         batch_script_tag,
      									 logdir,batch_script_tag,
      									 logdir,batch_script_tag,
      									 cmethod$ncpu,
                         cmethod$user_email,
      									 cmethod$memG)
      if (cmethod$hours>0) {
      	cmd_str <- sprintf("%s\n#SBATCH --time=%d:00:00\n",cmd_str,cmethod$hours)
      }
    }
    
    inputj <- inputs[[j]]
    
    if (!is.null(cmethod$preamble)) {
    	cmd_str <- sprintf("%s\n%s\n",cmd_str,cmethod$preamble)
    }
    
    for (i in 1:dim(inputj)[1]) {
      
      if (target_hostname=='cc-dclrilog52') {
        #change mount directory
        cmd_i <- gsub(Sys.getenv("LABSHARED"),
             '/mnt/isilon/w_QHS/hwangt-share',
             inputj$cmd[i])
        } else {
        cmd_i <-inputj$cmd[i]
        }
      cmd_str <- sprintf("%s\n%s\n",cmd_str,cmd_i)
    }
    writeLines(cmd_str,con=fp)
    close(fp)
    message('Done.')
  }
  if (debug2>0) {browser()}
}


sbatch_generator_osc <- function(cmethod,nscripts,scriptd,jobq="<hidden>",debug2=0) {
  
  if (debug2>0) {browser()}
  
  W <- dim(cmethod$wks)[1]
  if (nscripts > W){
    nscripts <- W
  }
  
  cmethod$wks$sbatch_id <- (1:W)%%nscripts
  inputs <- split(cmethod$wks,by=c('sbatch_id'))
  
  for (j in 1:nscripts) {
    if (debug2>1) {browser()}
    batch_script_tag <- sprintf("%d_%s",j,cmethod$exptag)
    script_fpath <- file.path(scriptd,sprintf('%s.sh',batch_script_tag))
    
    logdir <- file.path(scriptd,'logs')
    if (!file.exists(logdir)) {
    	dir.create(logdir,showWarnings = FALSE, recursive = TRUE)
    }
    
    message(sprintf('generating a script file [%s] ...',script_fpath))
    fp <- file(script_fpath,'w')
    
    cmd_str <- sprintf("#!/bin/bash -l")
    if (cmethod$user_email!="") {
      cmd_str <- sprintf("%s
#SBATCH -J %s
#SBATCH -o %s/%s.%%N.%%j.out
#SBATCH -e %s/%s.%%N.%%j.err
#SBATCH --account=%s
#SBATCH --ntasks=%d
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%s
#SBATCH --mem=%dG",
                         cmd_str,
                         batch_script_tag,
                         logdir,batch_script_tag,
                         logdir,batch_script_tag,
                         jobq,cmethod$ncpu,
                         cmethod$user_email,
                         cmethod$memG)
      
      if (cmethod$hours>0) {
        cmd_str <- sprintf("%s\n#SBATCH --time=%d:00:00\n",cmd_str,cmethod$hours)
      }
    }
    
    if ('bash_pre' %in% names(cmethod)) {
    	pre_cmds = comma_string_to_list(cmethod$bash_pre,sep=";")
    	for (pre_cmd in pre_cmds) {
    		cmd_str <- sprintf("%s\n%s\n",cmd_str,pre_cmd)
    	}
    }
    
    inputj <- inputs[[j]]
    
    for (i in 1:dim(inputj)[1]) {
      cmd_i <-inputj$cmd[i]
      cmd_str <- sprintf("%s\n%s\n",cmd_str,cmd_i)
    }
    writeLines(cmd_str,con=fp)
    close(fp)
    message('Done.')
  }
  if (debug2>0) {browser()}
}


get_job_msg <- function(jobname,state="F") {
	uuid <- format(Sys.time(), "%Y%m%d.%H%M%S")
	job_msg <- sprintf("[%s-%s]%s",uuid,jobname,state)
}

ifttt_notify <- function(jobname="hwang_lab",state="F") {
	key = "<hidden>"
	job_msg <- get_job_msg(jobname,state)
	
	# tmp_logd = file.path(Sys.getenv('TMPDIR'),'logd')
	# sys_logd = file.path(Sys.getenv('PPLINE'),'logd')
	url2post <- sprintf("<hidden>",key,job_msg)
	POST(url2post,body=NULL,verbose())
}

assign_batch_run <- function(tn4_dt,num_scripts=8,batch_run_prefix="br") {
	tn4_dt$sidx <- 1:nrow(tn4_dt)
	run_batches <- chunk_by_unitlen(1:nrow(tn4_dt),unit_len = ceiling(nrow(tn4_dt)/num_scripts))
	names(run_batches) <- sprintf("%s.%s",batch_run_prefix,names(run_batches))
	run_schedule = rbindlist(lapply(names(run_batches),function(brun){
		data.table(brun,sidx=run_batches[[brun]])
	}))
	tn4_dt$brun <- run_schedule[match(tn4_dt$sidx,run_schedule$sidx),brun]
	tn4_dt$sidx <- NULL
	tn4_dt
}
