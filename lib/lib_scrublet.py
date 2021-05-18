# %matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

def annotate_doublets(mtx_fpath,feature_fpath,expected_doublet_rate2=0.06):
	if False:
		plt.rcParams['font.family'] = 'sans-serif'
		plt.rcParams['font.sans-serif'] = 'Arial'
		plt.rc('font', size=14)
		plt.rcParams['pdf.fonttype'] = 42
	
	counts_matrix = scipy.io.mmread(mtx_fpath).T.tocsc()
	genes = np.array(scr.load_genes(feature_fpath, delimiter='\t', column=1))
	
	print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
	print('Number of genes in gene list: {}'.format(len(genes)))
	
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate2)

	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
	                                                          min_cells=3,
	                                                          min_gene_variability_pctl=85,
	                                                          n_prin_comps=30)
	
	if False:
		scrub.plot_histogram()
		
		print('Running UMAP...')
		scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
		print('Done.')
		
		scrub.plot_embedding('UMAP', order_points=True)
	
	return([doublet_scores,predicted_doublets])
