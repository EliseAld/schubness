from itertools import product
import os
import anndata2ri
import anndata
import scipy
from skhubness.neighbors import kneighbors_graph
import numpy as np
import scanpy as sc
import time
import warnings
import sklearn
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

anndata2ri.activate()

def recipe_duo(adata, do_log, renorm):
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)
    exprsn = np.array(adata.X.mean(axis=0)).reshape(-1)
    keep = np.argsort(exprsn)[::-1][:5000]
    adata._inplace_subset_var(keep)
    if renorm:
        sc.pp.normalize_total(adata,target_sum=1e4)

def recipe_seurat(adata, do_log, norm_scale):
    #same as scanpy clustering tutorial except initial cells/genes prefiltering /!\ sc.pp.recipe_seurat filters non log data ?
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata._inplace_subset_var(adata.var.highly_variable)
    if norm_scale:
        sc.pp.scale(adata, max_value=10)

def embedding_tsne(adata, do_norm, norm_scale, do_log,do_pca,
                metric,
                seed, n_comps, perplexities):
    hubness_methods = {'nothing': (None, None),
                       'mp_normal': ('mp', {'method': 'normal'}),
                       'ls': ('ls', None),
                       'ls_nicdm': ('ls', {'method': 'nicdm'}),
                       'dsl': ('dsl', None)}
    ### preprocess and PCA step ###
    start0 = time.time()
    if type(do_norm) is str:
        adata.X = scipy.sparse.csr_matrix(adata.X)
        if do_norm == 'seurat':
            recipe_seurat(adata, do_log, norm_scale)
            print(f'\t\tseurat norm retained {adata.X.shape[1]} genes')
        elif do_norm == 'duo':
            recipe_duo(adata, do_log, renorm=norm_scale)
            print(f'\t\tduo norm retained {adata.X.shape[1]} genes')
        else:
            raise ValueError("do_norm not in duo, seurat")
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    if do_log and not(type(do_norm) is str):
        print('\t\tlog_transformed data')
        sc.pp.log1p(adata)
    if do_pca:
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, n_comps)))
        X = adata.obsm['X_pca']
    else:
        print('pca not done!')
        X = adata.X
    n_neighbors = int(np.sqrt(X.shape[0]))
    print('\t\t\tPreprocessing and PCA pipeline:', round((time.time()-start0)/60, 2), 'mn')
    ### Hub reduction ###
    start = time.time()
    all_adata = dict()
    for kernel in ['gauss', 'umap']:
        all_adata[kernel] = adata.copy()
        try:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors + 1, metric=metric, use_rep='X', method=kernel)
        except:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors + 1, metric=metric, use_rep='X', method=kernel,
                            knn=False)
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsp['connectivities'] = kneighbors_graph(X,
                                                                         n_neighbors=n_neighbors,
                                                                         hubness=hubness,
                                                                         hubness_params=hubness_params,
                                                                         metric=metric,
                                                                         mode="connectivity")
        all_adata[method_name].obsp['distances'] = kneighbors_graph(X,
                                                                    n_neighbors=n_neighbors,
                                                                    hubness=hubness,
                                                                    hubness_params=hubness_params,
                                                                    metric=metric,
                                                                    mode="distance")
        all_adata[method_name].uns['neighbors'] = {'connectivities_key': 'connectivities',
                                                   'distances_key': 'distances',
                                                   'params': {'n_neighbors': n_neighbors,
                                                              'method': 'umap',
                                                              'metric': metric}}
    print('\t\t\tkNN graph pipeline:', round((time.time() - start) / 60, 2), 'mn')
    ### tSNE embedding ###
    start = time.time()
    tsne = [sklearn.manifold.TSNE(n_components=2, metric='precomputed', random_state=seed, perplexity=k) for k in perplexities]
    tsne_viz = dict()
    for idx0, method_name in enumerate(['nothing', 'umap', 'gauss', 'mp_normal', 'ls', 'ls_nicdm', 'dsl']):
        tsne_viz[method_name] = dict()
        for idx1, perplex in enumerate(perplexities):
            tsne_viz[method_name][perplex] = tsne[idx1].fit_transform(all_adata[method_name].obsp['distances'].toarray())
    print('\t\t\ttSNE embedding pipeline:', round((time.time() - start) / 60, 2), 'mn')
    ### plot ###
    plt.clf()
    fig = plt.figure(figsize=(15, 25 ))
    k = 0
    for j, method_name in enumerate(tsne_viz.keys()):
        for i, perplex in enumerate(perplexities):
            k += 1
            ax = fig.add_subplot(7, 4, k)
            ax.set_title('perplexity='+str(perplex)+', card='+str(adata.n_obs))
            for group in np.unique(adata.uns['Order']):
                mask = adata.uns['Order'] == group
                ax.scatter(tsne_viz[method_name][perplex][mask, 0],
                            tsne_viz[method_name][perplex][mask, 1], c=colors2[int(group)])
    plt.savefig(get_res_path(fname)+'.png', bbox_inches='tight', pad_inches=0)
    #plt.show()
    print('\t\t\tFull pipeline:', round((time.time()-start0)/60, 2), 'mn')


##### USAGE ######

#fixed params
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_embedding/tsne_test_perplex/"
#func to save results with proper name
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'
do_log = True
do_pca = True
weighted = True
norm_scale = True
clustering_algo = 'leiden'
do_norm = 'seurat'
n_comps = 50
seed = 0
fnames = ['GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9', 'GSE60783_10', 'GSE52529_15', 'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33']
colors2 = ["lightcoral", "red", "gold", "yellow", "yellowgreen", "lightgreen", "forestgreen",
           "aquamarine", "mediumturquoise", "deepskyblue", "royalblue", "blue", "blueviolet",
           "violet", "hotpink"]


#vary params
metrics = ('cosine', 'euclidean')
perplexities = (40, 60, 80, 100)
#params_list = list(product(metric, n_comps))

path_h5ad = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/raw/raw_"
for metric in metrics:
    print(metric, n_comps, do_norm, clustering_algo)
    for fname in fnames:
        if np.array([get_res_path(fname).split('/')[-1] in elt for elt in os.listdir(resdir)]).any():
            continue
        print('\t\t'+fname)
        adata = anndata.read_h5ad(path_h5ad+fname+'.h5ad')

        embedding_tsne(
            adata,
            do_norm=do_norm,
            norm_scale=norm_scale,
            do_log=do_log,
            do_pca=do_pca,
            metric=metric,
            seed=seed,
            n_comps=n_comps,
            perplexities=perplexities)


