#!/bin/bash -l
set -ex

projd0=${HOME}/projects

if [[ ! -d ${projd0} ]]; then
	echo "creating ${projd0} ..."
	mkdir -p ${projd0}

cd ${projd0}

# git clone https://github.com/hwanglab/hwanglab_2021_tigitCarT ${projd0}

projd=${projd0}/hwanglab_2021_tigitCarT
echo "project directory:${projd}"

cd ${projd}

echo "Rscript ${projd}/m01_R_module_setup.r"
# Rscript ${projd}/m01_R_module_setup.r

# the step m05 takes too long to finish and the output files are provided directly into the github!
# echo "Rscript ${projd}/m05a_harmony_findmarkers.r"
# Rscript ${projd}/m05a_harmony_findmarkers.r

echo "Rscript ${projd}/m06a_vcluster_heatmap.r"
# Rscript ${projd}/m06a_vcluster_heatmap.r

echo "Rscript ${projd}/m07a_celltype_spiderweb.r"
# Rscript ${projd}/m07a_celltype_spiderweb.r

echo "Rscript ${projd}/m08a_tigit_expr_by_tgroup.r"
# Rscript ${projd}/m08a_tigit_expr_by_tgroup.r
# 
echo "Rscript ${projd}/m09_deg_pre_vs_post_inf.r"
# Rscript ${projd}/m09_deg_pre_vs_post_inf.r
# 
echo "Rscript ${projd}/m10_deg_by_cluster_inter.r"
# Rscript ${projd}/m10_deg_by_cluster_inter.r
# 
echo "Rscript ${projd}/m11_violin_ridge_gexpr.r"
# Rscript ${projd}/m11_violin_ridge_gexpr.r
# 
echo "Rscript ${projd}/m12_dyfscore_TIGIT_on_tsne.r"
# Rscript ${projd}/m12_dyfscore_TIGIT_on_tsne.r
# 
echo "Rscript ${projd}/m13_tigitON_pre_vs_post.r"
# Rscript ${projd}/m13_tigitON_pre_vs_post.r