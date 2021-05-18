source('lib/lib_project.r')

#--------------
args=data.table(outd="./data",
								pat_clin="./data/CASE_2417_patient_info.xlsx",
								ega_template="./data/ega_sample_template.csv",
								seqrun1907_info_sheet="~/hwangt-share/Datasets/TaeHyun_Hwang_1907UNHS-0026/1907UNHS-0026/doc/sample_info.xlsx",
								egacryptor="~/hwangt-share/apps/ega/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar")

# --------------
sample_sheets <- list()
fastqs <- list()

# ---------------
#get EGA sample template
ega_sample.dt <- fread(args$ega_template)
ega_cols <- colnames(ega_sample.dt)

#[1] "title"         "alias"         "description"   "subjectId"     "bioSampleId"   "caseOrControl"
# [7] "gender"        "organismPart"  "cellLine"      "region"        "phenotype"

# ------------
# load pat clin
patclin = as.data.table(read.xlsx(args$pat_clin))
patclin[,patid:=sprintf("P%d",Trial.Number)]
patclin[Gender=="F",Gender:="female"]
patclin[Gender=="M",Gender:="male"]

###################################
srun=1
sample_sheets[["seqrun1907"]]=as.data.table(read.xlsx(args$seqrun1907_info_sheet))
fastqs[["seqrun1907"]] = Sys.glob('~/hwangt-share/Datasets/TaeHyun_Hwang_1907UNHS-0026/1907UNHS-0026/*/*R1.fastq.gz')
fq1=fastqs[["seqrun1907"]]

fqdt2=cbind(fq1,decompose_fname(fq1))
fqdt2[,lane_info:=tstrsplit(fbase0,"_")[[1]]]
# fqdt2[,index_info:=tstrsplit(fbase0,"_")[[2]]]

swks=sample_sheets$seqrun1907
swks[Timepoint=="pre-infusion",tpoint:="D0"]
swks[Timepoint=="day14",tpoint:="D14"]
swks[Timepoint=="day30",tpoint:="D30"]

swks[type=="CAR-T",stype:="c"]
swks[type=="pbmc",stype:="p"]

swks[,alias:=sprintf("P%s.%s.%s%d",Patient,tpoint,stype,srun)]
swks[,subjectId:=sprintf("P%d",Patient)]

swks$gender = patclin[match(swks$Patient,Trial.Number),Gender]
swks$phenotype = patclin[match(swks$Patient,Trial.Number),Disease]

ega_samplej.dt = merge(fqdt2[,list(lane_info,fq1,parentd)],swks,by.x="lane_info",by.y="sample")
ega_samplej.dt[,c('title','description','bioSampleId','caseOrControl','cellLine','region'):=""]
ega_samplej.dt$organismPart <- "blood"

ega_cols2 <- c(ega_cols,"fq1","parentd")

ega_sample_to_upload.dt = ega_samplej.dt[,..ega_cols2]
ega_sample_to_upload.dt[,title:=alias]
ega_sample_to_upload.dt[,alias:=parentd]
ega_sample_to_upload.dt$parentd=NULL
# ega_sample_to_upload.dt[,alias:=sprintf("%s.i%s",alias,index_info)]

swks_tsv <- file.path(args$outd,sprintf("samples%d_w_flink.csv",srun))
fwrite(ega_sample_to_upload.dt,file=swks_tsv,sep=",")

swks_rds <- file.path(args$outd,sprintf("samples%d_w_flink.rds",srun))
saveRDS(ega_sample_to_upload.dt,file=swks_rds)

swks_tsv2 <- file.path(args$outd,sprintf("samples%d.csv",srun))
ega_sample_to_upload.dt$fq1<-NULL
fwrite(ega_sample_to_upload.dt,file=swks_tsv2,sep=",")

#----------------
# link files
# Sample alias	 "First Fastq File"	 "First Checksum"	 "First Unencrypted checksum"	 "Second Fastq File"	 "Second Checksum"	 "Second Unencrypted checksum"

swk.dt=readRDS(file=swks_rds)
swk.dt[,fq2:=gsub("R1","R2",fq1)]

lapply(1:nrow(swk.dt),function(r) {
	
	java_cmd <- sprintf("java -jar %s -i %s",args$egacryptor, swk.dt$fq1[r])
	message(java_cmd)
	run_system(java_cmd)
	
	browser()
	java_cmd <- sprintf("java -jar %s -i %s",args$egacryptor, swk.dt$fq2[r])
	message(java_cmd)
	run_system(java_cmd)
	browser()
	
	tools::md5sum(swk.dt$fq1[r])
	browser()
	
	tools::md5sum(swk.dt$fq2[r])
	browser()
	
})

# ------------
# to generate md5 for unencrypted FASTQ files

