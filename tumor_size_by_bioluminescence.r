#0. load shared library
source('lib_project.r')

fpath_dt <- get_current_script_fpath()

args <- data.table(step=get_pipeline_step(fpath_dt$fpath),
									 input_xlsx="out/cart/m01_3seqruns/cart/09a_bioluminescence/stab5_bioluminescence.xlsx",
									 outd="out/cart/m01_3seqruns/cart/09a_bioluminescence",
									 #test.method2="wilcox.test",
									 test.method2="t.test",
									 ncpu=2,
									 reuse=0,
									 debug=0)

create_dir_if_not_exist(args$outd)

biolum.dt = as.data.table(read.xlsx(args$input_xlsx))

biolum.dt = biolum.dt[sample!="mouse1" | !(timepoint %in% c("Week 3","Week 5")),]

biolum.dt[sgroup=="CART & TIGIT",sgroup:="CART & TIGIT-ab"]

my_comps = list(ctrl_IgG=c("No Treatment","CART & IgG"),
								IgG_TIGIT=c("CART & IgG","CART & TIGIT-ab"),
								ctrl_TIGIT=c("No Treatment","CART & TIGIT-ab"))

biolum.dt[,biolum_log10 := log10(bioluminescence)]

biolum.dt[,timepoint2 := timepoint]

biolum.dt[timepoint %in% c("Week 1", "Week 2", "Week 2.5"),timepoint2:="Week 1-2.5"]
biolum.dt[timepoint %in% c("Week 3","Week 4","Week 5"),timepoint2:="Week 3-5"]

a=biolum.dt[,.N,by=c("timepoint2","sgroup")]
a[order(timepoint2,sgroup),]

biolum.dt[,sgroup:=factor(sgroup,levels = c("No Treatment","CART & IgG","CART & TIGIT-ab"))]

pdf_fn = file.path(args$outd,sprintf("sfig5H_TMB_by_bioluminescence_wk_%s_by_treatment.pdf",args$test.method2))
pdf(file=pdf_fn,width=9,height=6)

p=ggplot(biolum.dt,
			 aes(x=sgroup,y=biolum_log10,color=sgroup)) +
	# geom_boxplot()+
	geom_violin()+
	geom_quasirandom(alpha = 1, dodge.width = 0.2, size=1 , shape = 1) +
	stat_summary(fun = mean, na.rm = TRUE, 
							 geom = "point", color = "black", 
							 size = 3, shape = "diamond") +
	stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
							 geom="errorbar", color="black", width=0.2) +
	stat_compare_means(method=args$test.method2,
										 hide.ns=T,
										 vjust=0.8,
										 label = "p.signif",
										 # label = "p.format",
										 comparisons = my_comps) +
	facet_grid(cols=vars(timepoint),scales="free_y", rows = vars(sample)) +
	theme_bw()+
	theme(axis.text.x = element_blank(),legend.position="right") +
	labs(title=args$test.method2,x="timepoint",y="log10(bioluminescence)")
plot(p)

p = ggplot(biolum.dt,
			 aes(x=sgroup,y=biolum_log10,color=sgroup)) +
	# geom_boxplot()+
	geom_violin()+
	geom_quasirandom(alpha = 1, dodge.width = 0.2, size=1 , shape = 1) +
	stat_summary(fun = mean, na.rm = TRUE, 
							 geom = "point", color = "black", 
							 size = 3, shape = "diamond") +
	stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
							 geom="errorbar", color="black", width=0.2) +
	stat_compare_means(method=args$test.method2,
										 hide.ns=T,
										 vjust=0.8,
										 # label = "p.signif",
										 label = "p.format",
										 comparisons = my_comps) +
	facet_grid(cols=vars(timepoint),scales="free_y", rows = vars(sample)) +
	theme_bw()+
	theme(axis.text.x = element_blank(),legend.position="right") +
	labs(title=args$test.method2,x="timepoint",y="log10(bioluminescence)")
plot(p)
dev.off()

##########################

pdf_fn = file.path(args$outd,sprintf("sfig5H_TMB_by_bioluminescence_group_%s_by_treatment.pdf",args$test.method2))
pdf(file=pdf_fn,width=6.5,height=7.5)

p=ggplot(biolum.dt,
									 aes(x=sgroup,y=biolum_log10,color=sgroup)) +
	# geom_boxplot()+
	geom_violin()+
	geom_quasirandom(alpha = 1, dodge.width = 0.2, size=1 , shape = 1) +
	stat_summary(fun = mean, na.rm = TRUE, 
							 geom = "point", color = "black", 
							 size = 3, shape = "diamond") +
	stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
							 geom="errorbar", color="black", width=0.2) +
	stat_compare_means(method=args$test.method2,
										 hide.ns=T,
										 vjust=0.8,
										 label = "p.signif",
										 # label = "p.format",
										 comparisons = my_comps) +
	facet_grid(cols=vars(timepoint2),scales="free_y", rows = vars(sample)) +
	theme_bw()+
	theme(axis.text.x = element_blank(),legend.position="right") +
	labs(title=args$test.method2,x="timepoint in group",y="log10(bioluminescence)")
plot(p)

p=ggplot(biolum.dt,
			 aes(x=sgroup,y=biolum_log10,color=sgroup)) +
	# geom_boxplot()+
	geom_violin()+
	geom_quasirandom(alpha = 1, dodge.width = 0.2, size=1 , shape = 1) +
	stat_summary(fun = mean, na.rm = TRUE, 
							 geom = "point", color = "black", 
							 size = 3, shape = "diamond") +
	stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
							 geom="errorbar", color="black", width=0.2) +
	stat_compare_means(method=args$test.method2,
										 hide.ns=T,
										 vjust=0.8,
										 # label = "p.signif",
										 label = "p.format",
										 comparisons = my_comps) +
	facet_grid(cols=vars(timepoint2),scales="free_y", rows = vars(sample)) +
	theme_bw()+
	theme(axis.text.x = element_blank(),legend.position="right") +
	labs(title=args$test.method2,x="timepoint in group",y="log10(bioluminescence)")
plot(p)
dev.off()
