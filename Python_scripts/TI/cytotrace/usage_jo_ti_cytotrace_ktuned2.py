from itertools import product
import os
import anndata
import scipy
import louvain
import leidenalg
from skhubness.neighbors import kneighbors_graph
import numpy as np
import skhubness
import scanpy as sc
import gc
import time
from tqdm import tqdm
import csv
import warnings
warnings.filterwarnings("ignore")

def getNclusters(adata, G, n_clusters, seed, clustering_algo, flavor, weights, range_min=0, range_max=3, max_steps=20):
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    weighted = weights is None
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        if clustering_algo == 'louvain':
            if flavor == 'scanpy':
                sc.tl.louvain(adata, resolution=this_resolution, random_state=seed, use_weights=weighted)
                clus = np.array(adata.obs['louvain']).astype(int)
                this_clusters = adata.obs['louvain'].nunique()
            elif flavor == 'base':
                part_louvain = louvain.find_partition(graph = G,
                                   partition_type = louvain.RBConfigurationVertexPartition,
                                   weights = weights,
                                   resolution_parameter=this_resolution, seed=seed)
                clus = np.array(part_louvain.membership)
                this_clusters = len(np.unique(clus))
        elif clustering_algo == 'leiden':
            if flavor == 'scanpy':
                sc.tl.leiden(adata, resolution=this_resolution, random_state=seed, use_weights=weighted)
                clus = np.array(adata.obs['leiden']).astype(int)
                this_clusters = adata.obs['leiden'].nunique()
            elif flavor == 'base':
                part_leiden = leidenalg.find_partition(graph = G,
                                         partition_type=leidenalg.RBConfigurationVertexPartition,
                                         weights=weights,
                                         resolution_parameter=this_resolution,
                                         seed=seed)
                clus = np.array(part_leiden.membership)
                this_clusters = len(np.unique(clus))
        else:
            raise ValueError("incorrect cluster_func, choose 'leiden' or 'louvain'")
        if this_clusters > n_clusters:
            this_max = this_resolution
        elif this_clusters < n_clusters:
            this_min = this_resolution
        else:
            return(this_resolution, weighted)
        this_step += 1
    return this_resolution, weighted

def generate_clustering_inputs(X, metric, n_neighbors, weighted, seed, hubness, hubness_params):
    hub = skhubness.Hubness(k=n_neighbors, metric=metric, hubness=hubness, hubness_params=hubness_params,
                            random_state=seed, store_k_occurrence=True, return_value='all').fit(X)
    del hub.X_train_; gc.collect()
    knn = hub.nn_index_
    if weighted:
        adjmat = knn.kneighbors_graph(mode='distance')
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        try:
            weights = np.array(G.es["weight"]).astype(np.float64)
        except:
            weights = None
        if weights is not None:
            if not np.isfinite(np.sum(weights)): #weights failed
                adjmat = knn.kneighbors_graph(mode='connectivity')
                G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
                weights = None
    else:
        adjmat = knn.kneighbors_graph(mode='connectivity')
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        weights = None
    del hub;gc.collect()
    return G, weights

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

def recipe_zheng17(adata, do_log, norm_scale, n_top_genes=1000):
    #same as sc.pp.recipe_zheng17 except initial cells/genes prefiltering
    sc.pp.normalize_total(adata, key_added='n_counts_all')
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=min(adata.X.shape[1], n_top_genes), log=False)
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    sc.pp.normalize_total(adata)  # renormalize after filtering
    if do_log:
        sc.pp.log1p(adata)  # log transform: X = log(X + 1)
    if norm_scale:
        sc.pp.scale(adata)

def ti_analysis(adata,true_labels,do_norm,norm_scale, do_log,do_pca,
                n_clusters,metric,weighted,  #weighted adjmat for louvain/leiden clustering ?
                seed,n_comps,clustering_algo,n_iter,bootstrap_size):
    hubness_methods = {'nothing': (None, None)#,
                       #'mp_normal': ('mp', {'method': 'normal'}),
                       #'ls': ('ls', None),
                       #'ls_nicdm': ('ls', {'method': 'nicdm'}),
                       #'dsl': ('dsl', None)
                       }
    start=time.time()
    ### preprocess, prepare clustering input ###
    if type(do_norm) is str:
        adata.X = scipy.sparse.csr_matrix(adata.X)
        if do_norm == 'seurat':
            recipe_seurat(adata, do_log, norm_scale)
            print(f'\t\tseurat norm retained {adata.X.shape[1]} genes')
        elif do_norm == 'zheng17':
            recipe_zheng17(adata, do_log, norm_scale, n_top_genes=5000)
            print(f'\t\tzheng norm retained {adata.X.shape[1]} genes')
        elif do_norm == 'duo':
            recipe_duo(adata, do_log, renorm=norm_scale)
            print(f'\t\tduo norm retained {adata.X.shape[1]} genes')
        else:
            raise ValueError("do_norm not in 'duo', seurat', 'zheng17'")
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    if do_log and not(type(do_norm) is str):
        print('\t\tlog_transformed data')
        sc.pp.log1p(adata)
    if do_pca:
        use_rep = 'X_pca'
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, n_comps)))
        X = adata.obsm['X_pca']
    else:
        print('pca not done!')
        use_rep = 'X'
        X = adata.X
    n_neighbors = int(np.sqrt(X.shape[0]))
    print('\t\t\tPreprocessing done:',round((time.time()-start)/60, 2), 'mn')
    start=time.time()
    ### clustering and PAGA step ###
    all_adata = dict()
    for kernel in ['umap','gauss']:
        all_adata[kernel] = adata.copy()
        try:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric,use_rep='X', method=kernel)
        except:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric,use_rep='X', method=kernel, knn=False)
        G, weights = generate_clustering_inputs(X=X,
                                                metric=metric,
                                                n_neighbors=n_neighbors,
                                                weighted=weighted,
                                                seed=seed,
                                                hubness=None,
                                                hubness_params=None)
        resol, weighted = getNclusters(all_adata[kernel], G, n_clusters=n_clusters, seed=seed, clustering_algo=clustering_algo,
                                        flavor='scanpy', weights=weights)
        if clustering_algo == "leiden":
            sc.tl.leiden(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[kernel], groups="leiden")
        elif clustering_algo == "louvain":
            sc.tl.louvain(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[kernel], groups="louvain")
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsm[method_name] = kneighbors_graph(X,
                                                            n_neighbors=n_neighbors,
                                                            hubness=hubness,
                                                            hubness_params=hubness_params,
                                                            metric=metric)
        sc.pp.neighbors(all_adata[method_name], n_neighbors=n_neighbors+1, use_rep=method_name)
        G, weights = generate_clustering_inputs(X=X,
                                                 metric=metric,
                                                 n_neighbors=n_neighbors,
                                                 weighted=weighted,
                                                 seed=seed,
                                                 hubness=hubness,
                                                 hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=n_clusters, seed=seed, clustering_algo=clustering_algo,
                                       flavor='base', weights=weights)
        if clustering_algo == "leiden":
            sc.tl.leiden(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[method_name], groups="leiden")
        elif clustering_algo == "louvain":
            sc.tl.louvain(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=seed)
            sc.tl.paga(all_adata[method_name], groups="louvain")
    print('\t\t\tHubness and PAGA full pipeline:',round((time.time()-start)/60,2),'mn')
    start=time.time()
    ### PAGA stab ###
    #all_iter = dict()
    #cell_iter = dict()
    #feat_iter = dict()
    #for method_name, (hubness, hubness_params) in hubness_methods.items():
    #    all_iter[method_name] = dict()
    #    cell_iter[method_name] = np.zeros((n_iter, adata.n_obs))
    #    feat_iter[method_name] = np.zeros((n_iter, adata.n_vars))
    #    for iter in tqdm(range(n_iter)):
    #        feat_bootstrap = np.random.uniform(0, 1, size=all_adata[method_name].n_vars)
    #        feat_bootstrap[feat_bootstrap <= bootstrap_size] = 0
    #        feat_bootstrap[feat_bootstrap > bootstrap_size]=1
    #        feat_bootstrap = feat_bootstrap == 0
    #        cell_bootstrap = np.random.uniform(0, 1, size=all_adata[method_name].n_obs)
    #        cell_bootstrap[cell_bootstrap <= bootstrap_size] = 0
    #        cell_bootstrap[cell_bootstrap > bootstrap_size] = 1
    #        cell_bootstrap = cell_bootstrap == 0
    #        cell_iter[method_name][iter, :] = cell_bootstrap
    #        feat_iter[method_name][iter, :] = feat_bootstrap
    #        uns = {'Order': true_labels[cell_bootstrap]}
    #        adata_sampled = anndata.AnnData(adata.X[cell_bootstrap][:, feat_bootstrap],
    #                                        uns=uns)
    #        n_clusters2 = len(np.unique(adata_sampled.uns['Order']))
    #        if do_pca:
    #            sc.tl.pca(adata_sampled, n_comps=min(adata_sampled.X.shape[1]-1, min(len(adata_sampled.X)-1, n_comps)))
    #            X2=adata_sampled.obsm['X_pca']
    #        else:
    #            X2=adata_sampled.X
    #        adata_sampled.obsm[method_name] = kneighbors_graph(X2,
    #                                                        n_neighbors=n_neighbors,
    #                                                        hubness=hubness,
    #                                                        hubness_params=hubness_params,
    #                                                           metric=metric)
    #        sc.pp.neighbors(adata_sampled, n_neighbors=n_neighbors+1, use_rep=method_name)
    #        G2, weights2 = generate_clustering_inputs(X=X2,
    #                                             metric=metric,
    #                                             n_neighbors=n_neighbors,
    #                                             weighted=weighted,
    #                                             seed=seed,
    #                                             hubness=hubness,
    #                                             hubness_params=hubness_params)
    #        resol2, weighted2 = getNclusters(adata_sampled, G2, n_clusters=n_clusters2, seed=seed, clustering_algo=clustering_algo,
    #                                   flavor='base', weights=weights2)
    #        if clustering_algo == "leiden":
    #            sc.tl.leiden(adata_sampled, resolution=resol2, use_weights=weighted2, random_state=seed)
    #            sc.tl.paga(adata_sampled, groups="leiden")
    #        elif clustering_algo == "louvain":
    #            sc.tl.louvain(adata_sampled, resolution=resol2, use_weights=weighted2, random_state=seed)
    #            sc.tl.paga(adata_sampled, groups="louvain")
    #        all_iter[method_name]['iter'+str(iter)] = adata_sampled.uns["paga"]["connectivities_tree"]
    #print('\t\t\tPAGA stability pipeline:',round((time.time()-start)/60,2),'mn')
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
        #w = csv.writer(open(get_res_path(fname)+'_'+method_name+"_stab.csv", "w"))
        #for key, val in all_iter[method_name].items():
        #    w.writerow([key, val])
        #np.savetxt(get_res_path(fname)+'_'+method_name+"_stab_cell.csv", cell_iter[method_name], delimiter=',', fmt='%d')
        #np.savetxt(get_res_path(fname)+'_'+method_name+"_stab_feat.csv", feat_iter[method_name], delimiter=',', fmt='%d')

##### USAGE ######
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned2/"
#func to save results with proper name
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'

#fixed params
do_log = True #already done with do_norm
do_pca = True
weighted=True
norm_scale = True
n_neighbors = 10
seed = 0
n_iter = 10
bootstrap_size = 0.95
fnames = ['GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9', 'GSE60783_10', 'GSE67602_11', 'GSE70245_12', 'GSE52529_15', 'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33']

#vary params
metric = ('cosine', 'euclidean')
n_comps = (25, 50, 100, 500)
#do_norm = ('duo', 'seurat', 'zheng17', None)
do_norm = ('duo', 'seurat')
clustering_algo = ('louvain', 'leiden')
params_list = list(product(metric, n_comps, do_norm, clustering_algo))

path_h5ad = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/raw/raw_"
for metric, n_comps, do_norm, clustering_algo in params_list:
    print(metric, n_comps, do_norm, clustering_algo)
    for fname in fnames:
        if np.array([get_res_path(fname).split('/')[-1] in elt for elt in os.listdir(resdir)]).any(): #dataset already computed for these settings
            continue
        print('\t\t'+fname)
        adata = anndata.read_h5ad(path_h5ad+fname+'.h5ad')
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
