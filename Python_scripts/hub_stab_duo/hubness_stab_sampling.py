from tqdm import tqdm
from skhubness import Hubness
import matplotlib.pyplot as plt
import pandas as pd
import anndata
import scipy
import numpy as np
import scanpy as sc
import time
import warnings
import os
import rpy2.robjects as ro
import anndata2ri

warnings.filterwarnings("ignore")
anndata2ri.activate()
readRDS = ro.r['readRDS']
path_rds = "/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/DimRedPaper/DuoClustering2018/sce_full/sce_full_"
path_res = "/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_duo/"
get_res_path = lambda path_res: f'{path_res}hub_stab_sampling_norm{do_norm}_scale{norm_scale}_{metric}_{clustering_algo}'

#fixed params
seed = 0
n_iter = 10
weighted = True
n_neighbors = 10
norm_scale = True
do_norm = 'seurat'
bootstrap_size = 0.9
metric = 'cosine'
clustering_algo = 'leiden'
fnames = ['Koh', 'KohTCC', 'Kumar', 'KumarTCC', #'SimKumar4easy', 'SimKumar4hard', 'SimKumar8hard',
          'Trapnell', 'TrapnellTCC']
colors = ['red', 'gold', 'chartreuse', 'deepskyblue']

#vary euclidean
n_comps = (25, 50, 100, 500)

def load_data(fname):
    adata = readRDS(path_rds+fname+'.rds')
    adata = anndata.AnnData(X=adata.X)
    return adata

def resampling(adata, bootstrap_size=bootstrap_size, n_iter=n_iter):
    adata_sampled = dict()
    cell_iter = np.zeros((n_iter, adata.n_obs))
    for iter in tqdm(range(n_iter)):
        cell_bootstrap = np.random.uniform(0, 1, size=adata.n_obs)
        cell_bootstrap[cell_bootstrap <= bootstrap_size] = 2
        cell_bootstrap = cell_bootstrap == 2
        cell_iter[iter, :] = cell_bootstrap
        adata_sampled[iter] = anndata.AnnData(X=adata.X[cell_bootstrap],
                                              obsm={'X_pca': adata.obsm['X_pca'][cell_bootstrap]})
    return adata_sampled, cell_iter

def hub_retrieval(X, metric, k=n_neighbors):
    n_obs = X.shape[0]
    hub = Hubness(k=k, metric=metric, random_state=seed, return_value='hubs')
    hub.fit(X)
    hub_id = np.repeat('norm', n_obs)
    hub_id[hub.score()] = 'hubs'
    return hub_id

def overlap_hubs(hubs_ref, hubs_sampled, cell_iter):
    hubs_ref_match = hubs_ref[cell_iter == 1]
    size_overlap = [100*np.sum([(hubs_ref_match == 'hubs')[loop] and
                                       (hubs_sampled == 'hubs')[loop] for loop in range(len(hubs_ref_match))])/
                            np.sum(hubs_ref == 'hubs')][0]
    return size_overlap


# run functions
size_overlap_tot = np.zeros(shape=(len(fnames), len(n_comps)))
for dim in n_comps:
    index = 0
    size_overlap = np.zeros(len(fnames))
    print(metric, dim, clustering_algo)
    #if np.array([get_res_path(path_res).split('/')[-1] in elt for elt in os.listdir(path_res)]).any(): #setting already computed
    #        continue
    for fname in fnames:
        adata = load_data(fname)
        if dim < adata.shape[0]:
            start=time.time()
            ### preprocess###
            sc.pp.log1p(adata)
            sc.pp.normalize_total(adata, target_sum=1e4)
            exprsn = np.array(adata.X.mean(axis=0)).reshape(-1)
            #adata2 = adata.copy()
            #adata2._inplace_subset_var(np.argsort(exprsn)[::-1][:10000])
            adata._inplace_subset_var(np.argsort(exprsn)[::-1][:10000])
            sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, dim)))
            if scipy.sparse.issparse(adata.X):
                adata.X = adata.X.toarray()
            X = adata.obsm['X_pca']
            print('\t\t\tPreprocessing done:', round((time.time()-start)/60, 2), 'mn')
            start = time.time()
            ### get hubs ref ###
            if dim < adata.shape[1]:
                # Make the lists of ref hubs
                hubs_ref = hub_retrieval(X, metric=metric)
                print('\t\t\tRef hubs done:', round((time.time()-start)/60, 2), 'mn')
                start = time.time()
                ### get sampled data###
                # Make the lists of sampled hubs
                adata_sampled, cell_iter = resampling(adata)
                X_sampled = dict()
                for key in adata_sampled.keys():
                    if scipy.sparse.issparse(adata_sampled[key].X):
                        adata_sampled[key].X = adata_sampled[key].X.toarray()
                    #sc.pp.log1p(adata_sampled[key])
                    #sc.pp.normalize_total(adata_sampled[key], target_sum=1e4)
                    #exprsn = np.array(adata_sampled[key].X.mean(axis=0)).reshape(-1)
                    #keep = np.argsort(exprsn)[::-1][:10000]
                    #sc.tl.pca(adata_sampled[key], n_comps=min(adata_sampled[key].X.shape[1]-1, min(len(adata_sampled[key].X)-1, dim)))
                    X_sampled[key] = adata_sampled[key].obsm['X_pca']
                hubs_s = dict()
                for iter in adata_sampled.keys():
                    hubs_s[iter] = hub_retrieval(X_sampled[iter], metric=metric)
                print('\t\t\tSampled hubs done:', round((time.time()-start)/60, 2), 'mn')
                start = time.time()
                ### get overlap###
                size_overlap[index] = np.mean([overlap_hubs(hubs_ref, hubs_s[i], cell_iter[i]) for i in range(n_iter)])
                print('\t\t\tOverlap calculated:', round((time.time()-start)/60, 2), 'mn')
        else:
            size_overlap[np.argwhere(np.array(fnames)==fname)] = None
        index += 1
    size_overlap_tot[:, np.argwhere(np.array(n_comps)==dim)[0][0]] = size_overlap

boxplot = plt.boxplot(size_overlap_tot, labels=['25', '50', '100', '500'], patch_artist=True)
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)
for dim in range(len(n_comps)):
    plt.scatter(np.repeat(dim+1, size_overlap_tot.shape[0]), size_overlap_tot[:, dim], c=colors[::-1][dim])
plt.savefig(get_res_path(path_res)+'2.png')
#plt.show()

# get the number of hubs per dataset
