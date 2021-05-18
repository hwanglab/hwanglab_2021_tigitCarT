if(!require(BiocManager)){
	install.packages('BiocManager')
	library('BiocManager')
}

req_pkg=read.table('./req_pkg.tsv')
for (pkg in req_pkg$V1) {
	browser()
	if( !is.element(pkg, .packages(all.available = TRUE)) ) {
		message(sprintf("installing %s ...",pkg))
		BiocManager::install(pkgs=pkg)
		message("Done.")
	}
	message(sprintf("loading %s ...",pkg))
	library(pkg,character.only = TRUE)
}

#######################
#install pkgs directly from github
github_pkg <- c("immunogenomics/harmony","satijalab/seurat-wrappers")
names(github_pkg) <- c("harmony","SeuratWrappers")
library(remotes)

for (pkg in names(github_pkg)) {
	if( !is.element(pkg, .packages(all.available = TRUE)) ) {
		remotes::install_github(github_pkg[[pkg]])
	}
	message(sprintf("loading %s ...",pkg))
	library(pkg,character.only = TRUE)
}