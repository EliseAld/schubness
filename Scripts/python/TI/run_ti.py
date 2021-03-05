from itertools import product
import os
from natsort import natsorted
import pandas as pd
import anndata2ri
import anndata
import scipy
import louvain
import leidenalg
from skhubness.neighbors import kneighbors_graph
import numpy as np
import scanpy as sc
import time
import warnings
from tqdm import tqdm
import csv
from benchmark_common_functions import getNclusters, generate_clustering_inputs, recipe_duo, recipe_seurat, load_h5ad, load_rds
warnings.filterwarnings("ignore")

anndata2ri.activate()


def ti_analysis(adata, true_labels, do_norm, norm_scale, do_log, do_pca,
                n_clusters, metric, weighted,  # weighted adjmat for louvain/leiden clustering ?
                seed, n_comps, clustering_algo, n_iter, bootstrap_size):
    hubness_methods = {'nothing': (None, None),
                       'mp_normal': ('mp', {'method': 'normal'}),
                       'ls': ('ls', None),
                       'ls_nicdm': ('ls', {'method': 'nicdm'}),
                       'dsl': ('dsl', None)}
    start = time.time()
    ### preprocess, prepare clustering input ###
    if type(do_norm) is str:
        adata.X = scipy.sparse.csr_matrix(adata.X)
        if do_norm == 'seurat':
            recipe_seurat(adata, do_log, norm_scale)
            # print(f'\t\tseurat norm retained {adata.X.shape[1]} genes')
        elif do_norm == 'duo':
            recipe_duo(adata, do_log, renorm=norm_scale)
            # print(f'\t\tduo norm retained {adata.X.shape[1]} genes')
        else:
            raise ValueError("do_norm not in 'duo', seurat'")
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    if do_log and not(type(do_norm) is str):
        # print('\t\tlog_transformed data')
        sc.pp.log1p(adata)
    if do_pca:
        use_rep = 'X_pca'
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, n_comps)))
        X = adata.obsm['X_pca']
    else:
        # print('pca not done!')
        use_rep = 'X'
        X = adata.X
    n_neighbors = int(np.sqrt(X.shape[0]))
    print('\t\t\tPreprocessing done:', round((time.time()-start)/60, 2), 'mn')
    start = time.time()
    ### clustering and PAGA step ###
    all_adata = dict()
    for kernel in ['umap', 'gauss']:
        all_adata[kernel] = adata.copy()
        try:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric, use_rep=use_rep, method=kernel)
        except:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric, use_rep=use_rep, method=kernel, knn=False)
        G, weights = generate_clustering_inputs(X=X,
                                                metric=metric,
                                                n_neighbors=n_neighbors,
                                                weighted=weighted,
                                                seed=seed,
                                                hubness=None,
                                                hubness_params=None)
        resol, weighted = getNclusters(all_adata[kernel], G, n_clusters=n_clusters, seed=seed,
                                       clustering_algo=clustering_algo, flavor='scanpy', weights=weights)
        if clustering_algo == "leiden":
            sc.tl.leiden(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[kernel], groups="leiden")
        elif clustering_algo == "louvain":
            sc.tl.louvain(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[kernel], groups="louvain")
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
                                                              'method': method_name,
                                                              'metric': metric}}
        G, weights = generate_clustering_inputs(X=X,
                                                metric=metric,
                                                n_neighbors=n_neighbors,
                                                weighted=weighted,
                                                seed=seed,
                                                hubness=hubness,
                                                hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=n_clusters, seed=seed,
                                       clustering_algo=clustering_algo, flavor='base', weights=weights)
        if clustering_algo == "louvain":
            clus = np.array(louvain.find_partition(graph=G,
                                                   partition_type=louvain.RBConfigurationVertexPartition,
                                                   weights=weights,
                                                   resolution_parameter=resol, seed=seed).membership)
            all_adata[method_name].obs['louvain'] = pd.Categorical(values=clus.astype('U'),
                                                                   categories=natsorted(map(str, np.unique(clus))),)
            sc.tl.paga(all_adata[method_name], groups="louvain", neighbors_key='neighbors')
        elif clustering_algo == "leiden":
            clus = np.array(leidenalg.find_partition(graph=G,
                                                     partition_type=leidenalg.RBConfigurationVertexPartition,
                                                     weights=weights,
                                                     resolution_parameter=resol,
                                                     seed=seed).membership)
            all_adata[method_name].obs['leiden'] = pd.Categorical(values=clus.astype('U'),
                                                                  categories=natsorted(map(str, np.unique(clus))),)
            sc.tl.paga(all_adata[method_name], groups="leiden")
    print('\t\t\tHubness and PAGA full pipeline:', round((time.time()-start)/60, 2), 'mn')
    start = time.time()
    ### PAGA stab ###
    all_iter = dict()
    cell_iter = dict()
    feat_iter = dict()
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_iter[method_name] = dict()
        cell_iter[method_name] = np.zeros((n_iter, adata.n_obs))
        feat_iter[method_name] = np.zeros((n_iter, adata.n_vars))
        for iter in tqdm(range(n_iter)):
            feat_bootstrap = np.random.uniform(0, 1, size=adata.n_vars)
            feat_bootstrap[feat_bootstrap <= bootstrap_size] = 0
            feat_bootstrap[feat_bootstrap > bootstrap_size] = 1
            feat_bootstrap = feat_bootstrap == 0
            cell_bootstrap = np.random.uniform(0, 1, size=adata.n_obs)
            cell_bootstrap[cell_bootstrap <= bootstrap_size] = 0
            cell_bootstrap[cell_bootstrap > bootstrap_size] = 1
            cell_bootstrap = cell_bootstrap == 0
            cell_iter[method_name][iter, :] = cell_bootstrap
            feat_iter[method_name][iter, :] = feat_bootstrap
            uns = {'Order': true_labels[cell_bootstrap]}
            adata_sampled = anndata.AnnData(adata.X[cell_bootstrap][:, feat_bootstrap],
                                            uns=uns)
            n_clusters2 = len(np.unique(adata_sampled.uns['Order']))
            if do_pca:
                sc.tl.pca(adata_sampled, n_comps=min(adata_sampled.X.shape[1]-1, min(len(adata_sampled.X)-1, n_comps)))
                X2 = adata_sampled.obsm['X_pca']
            else:
                X2 = adata_sampled.X
            adata_sampled.obsp["connectivities"] = kneighbors_graph(X2,
                                                                    n_neighbors=n_neighbors,
                                                                    hubness=hubness,
                                                                    hubness_params=hubness_params,
                                                                    metric=metric,
                                                                    mode="connectivity")
            adata_sampled.obsp["distances"] = kneighbors_graph(X2,
                                                               n_neighbors=n_neighbors,
                                                               hubness=hubness,
                                                               hubness_params=hubness_params,
                                                               metric=metric,
                                                               mode="distance")
            adata_sampled.uns['neighbors'] = {'connectivities_key': 'connectivities',
                                              'distances_key': 'distances',
                                              'params': {'n_neighbors': n_neighbors,
                                                         'method': method_name,
                                                         'metric': metric}}
            G2, weights2 = generate_clustering_inputs(X=X2,
                                                      metric=metric,
                                                      n_neighbors=n_neighbors,
                                                      weighted=weighted,
                                                      seed=seed,
                                                      hubness=hubness,
                                                      hubness_params=hubness_params)
            resol2, weighted2 = getNclusters(adata_sampled, G2, n_clusters=n_clusters2, seed=seed,
                                             clustering_algo=clustering_algo, flavor='base', weights=weights2)
            if clustering_algo == "leiden":
                clus = np.array(leidenalg.find_partition(graph=G2,
                                                         partition_type=leidenalg.RBConfigurationVertexPartition,
                                                         weights=weights2,
                                                         resolution_parameter=resol2,
                                                         seed=seed).membership)
                adata_sampled.obs['leiden'] = pd.Categorical(values=clus.astype('U'),
                                                             categories=natsorted(map(str, np.unique(clus))),)
                sc.tl.paga(adata_sampled, groups="leiden")
            elif clustering_algo == "louvain":
                clus = np.array(louvain.find_partition(graph=G2,
                                                       partition_type=louvain.RBConfigurationVertexPartition,
                                                       weights=weights2,
                                                       resolution_parameter=resol2, seed=seed).membership)
                adata_sampled.obs['louvain'] = pd.Categorical(values=clus.astype('U'),
                                                              categories=natsorted(map(str, np.unique(clus))), )
                sc.tl.paga(adata_sampled, groups="louvain")
            all_iter[method_name]['iter'+str(iter)] = adata_sampled.uns["paga"]["connectivities_tree"]
    print('\t\t\tPAGA stability pipeline:', round((time.time()-start)/60, 2), 'mn')
    for method_name in all_adata.keys():
        if method_name == "nothing":
            all_adata[method_name] = anndata.AnnData(X=all_adata[method_name].X,
                                                     uns={'Order': all_adata[method_name].uns['Order'],
                                                          'paga': all_adata[method_name].uns['paga']},
                                                     obs=all_adata[method_name].obs)
        else:
            all_adata[method_name] = anndata.AnnData(X=all_adata[method_name].X[:, :2],
                                                     uns={'Order': all_adata[method_name].uns['Order'],
                                                          'paga': all_adata[method_name].uns['paga']},
                                                     obs=all_adata[method_name].obs)
        all_adata[method_name].write_h5ad(filename=get_res_path(fname)+'_'+method_name+".h5ad")
        if method_name not in ["umap", "gauss"]:
            w = csv.writer(open(get_res_path(fname)+'_'+method_name+"_stab.csv", "w"))
            for key, val in all_iter[method_name].items():
                w.writerow([key, val])
            np.savetxt(get_res_path(fname)+'_'+method_name+"_stab_cell.csv", cell_iter[method_name], delimiter=',', fmt='%d')
            np.savetxt(get_res_path(fname)+'_'+method_name+"_stab_feat.csv", feat_iter[method_name], delimiter=',', fmt='%d')


##### USAGE ######
resdir = "/path/to/results/"
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'

# fixed params
do_log = True  # already done with do_norm
do_pca = True
weighted = True
norm_scale = True
seed = 0
n_iter = 10
bootstrap_size = 0.95
fnames = ['dataset1', 'dataset2', ...]
# fnames list of datasets


# varying params
metric = ('cosine', 'euclidean')
n_comps = (25, 50, 100, 500)
do_norm = ('duo', 'seurat')
clustering_algo = ('louvain', 'leiden')
params_list = list(product(metric, n_comps, do_norm, clustering_algo))

path_data = "/path/to/data/"
for metric, n_comps, do_norm, clustering_algo in params_list:
    print(metric, n_comps, do_norm, clustering_algo)
    for fname in fnames:
        if np.array([get_res_path(fname).split('/')[-1] in elt for elt in os.listdir(resdir)]).any():  # dataset already computed for these settings
            continue
        print('\t\t'+fname)
        # if data is RDS
        adata = load_rds(path_data, fname)
        # if data is h5ad
        adata = load_h5ad(path_data, fname)

        ti_analysis(
            adata,
            true_labels=adata.uns['Order'],
            do_norm=do_norm,
            norm_scale=norm_scale,
            do_log=do_log,
            do_pca=do_pca,
            n_clusters=len(np.unique(adata.uns['Order'])),
            metric=metric,
            weighted=weighted, 
            seed=seed,
            n_comps=n_comps,
            clustering_algo=clustering_algo,
            n_iter=n_iter,
            bootstrap_size=bootstrap_size
        )
