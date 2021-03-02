from itertools import product
import os
import rpy2.robjects as ro
from natsort import natsorted
import pandas as pd
import anndata2ri
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
import warnings
import csv
from sklearn import metrics
import math
warnings.filterwarnings("ignore")

anndata2ri.activate()

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
    del hub; gc.collect()
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

def Hbeta(D, beta):
    P = np.exp(-D * beta)
    sumP = np.sum(P)
    if sumP == 0:
        H = 0
        P = D * 0
    else:
        H = np.log(sumP) + beta * np.sum(np.matmul(D,P)) /sumP
        P = P/sumP
    return H, P

def get_sigma(X, perplexity = 15, tol = 1e-5, verbose = True):
    D = X
    n = X.shape[0]
    P = np.zeros((n, n))
    beta = np.ones((n))
    logU = np.log(perplexity)
    for i in range(n):
        if verbose:
            print(i)
        betamin = -math.inf
        betamax = math.inf
        Di = D[i, :][[idx for idx, e in enumerate(range(n)) if e != i]]
        H, thisP = Hbeta(Di, beta[i])
        Hdiff = H - logU
        tries = 0
        while np.abs(Hdiff) > tol and tries < 50:
            if Hdiff > 0:
                betamin = beta[i]
                if betamax == math.inf:
                    beta[i] = beta[i] * 2
                else:
                    beta[i] = (beta[i] + betamax) /2
            else:
                betamax = beta[i]
                if betamin == -math.inf:
                    beta[i] = beta[i]/ 2
                else:
                    beta[i] = (beta[i] + betamin) / 2
            H, thisP = Hbeta(Di, beta[i])
            Hdiff = H - logU
            tries = tries + 1
        P[i, :] = np.insert(thisP, i, 0)
    sigma = np.sqrt(1/beta)
    return sigma, P

def ti_analysis(adata,true_labels,do_norm,norm_scale, do_log,do_pca,
                n_clusters,metric,weighted,  #weighted adjmat for louvain/leiden clustering ?
                seed,n_comps,clustering_algo,n_iter,bootstrap_size):
    hubness_methods = {'nothing': (None, None),
                       'mp_normal': ('mp', {'method': 'normal'}),
                       'ls': ('ls', None),
                       'ls_nicdm': ('ls', {'method': 'nicdm'}),
                       'dsl': ('dsl', None)}
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
    print('\t\t\tPreprocessing done:', round((time.time()-start)/60, 2), 'mn')
    start=time.time()
    ### clustering and PAGA step ###
    all_adata = dict()
    for kernel in ['umap', 'gauss']:
        all_adata[kernel] = adata.copy()
        try:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric, use_rep='X', method=kernel)
        except:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors+1, metric=metric, use_rep='X', method=kernel, knn=False)
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
        sc.pl.paga(all_adata[kernel], show=False, random_state=seed)
        sc.tl.umap(all_adata[kernel], init_pos="paga", random_state=seed)
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsp['connectivities'] = kneighbors_graph(X,
                                                            n_neighbors = n_neighbors,
                                                            hubness = hubness,
                                                            hubness_params = hubness_params,
                                                            metric = metric,
                                                            mode="connectivity")
        all_adata[method_name].obsp['distances'] = kneighbors_graph(X,
                                                            n_neighbors = n_neighbors,
                                                            hubness = hubness,
                                                            hubness_params = hubness_params,
                                                            metric = metric,
                                                            mode="distance")
        all_adata[method_name].uns['neighbors'] = {'connectivities_key': 'connectivities',
                                                   'distances_key': 'distances',
                                                   'params': {'n_neighbors': n_neighbors,
                                                              'method': 'umap',
                                                              'metric': metric}}
        G, weights = generate_clustering_inputs(X=X,
                                                 metric=metric,
                                                 n_neighbors=n_neighbors,
                                                 weighted=weighted,
                                                 seed=seed,
                                                 hubness=hubness,
                                                 hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=n_clusters, seed=seed, clustering_algo=clustering_algo,
                                       flavor='base', weights=weights)
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
        sc.pl.paga(all_adata[method_name], show=False, random_state=seed)
        sc.tl.umap(all_adata[method_name], init_pos="paga", random_state=seed)
    print('\t\t\tHubness and PAGA full pipeline:', round((time.time()-start)/60, 2), 'mn')
    start = time.time()
    ### Evaluate goodness of fit ###
    original = adata.X
    dist_original = metrics.pairwise_distances(original)
    sigma, P = get_sigma(dist_original, perplexity=30, verbose=False)
    cost = []
    for method_name in all_adata.keys():
        dist_paga = metrics.pairwise_distances(all_adata[method_name].obsm['X_umap'])
        H, Q = Hbeta(dist_paga, 2)
        for i in range(Q.shape[0]):
            Q[i, i] = 0
        tot = 0
        for i in range(Q.shape[0]):
            for j in range(Q.shape[0]):
                if P[i, j] > 0 and Q[i, j] != 0:
                    tot += P[i, j]*np.log(P[i, j]/Q[i, j])
        cost.append(tot)
    print('\t\t\tCompute match between original and PAGA:', round((time.time() - start) / 60, 2), 'mn')
    np.savetxt(get_res_path(fname)+"_cost.csv", cost, delimiter=',', fmt='%d')

##### USAGE ######
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned5/"
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned5/"
#func to save results with proper name
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'

#fixed params
do_log = True #already done with do_norm
do_pca = True
weighted = True
norm_scale = True
seed = 0
n_iter = 10
bootstrap_size = 0.95
fnames = ['gold_hematopoiesis-gates_olsson', 'gold_stimulated-dendritic-cells-LPS_shalek',
          'gold_developing-dendritic-cells_schlitzer', 'gold_germline-human-male_guo',
          'gold_aging-hsc-young_kowalczyk']
fnames = ['gold_germline-human-female-weeks_li',
          'gold_germline-human-female_guo', 'gold_mESC-differentiation_hayashi',
          'gold_germline-human-male-weeks_li', 'gold_pancreatic-beta-cell-maturation_zhang', 'gold_human-embryos_petropoulos',
          'gold_myoblast-differentiation_trapnell', 'gold_pancreatic-alpha-cell-maturation_zhang',
          'gold_aging-hsc-old_kowalczyk']
fnames = ['GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9', 'GSE60783_10', 'GSE67602_11']
fnames = ['GSE70245_12', 'GSE52529_15', 'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33']


#vary params
metric = ('cosine', 'euclidean')
n_comps = (25, 50, 100, 500)
do_norm = ('duo', 'seurat')
clustering_algo = ('louvain', 'leiden')
params_list = list(product(metric, n_comps, do_norm, clustering_algo))

readRDS = ro.r['readRDS']
path_rds = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/rds/"
path_h5ad = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/raw/raw_"
for metric, n_comps, do_norm, clustering_algo in params_list:
    print(metric, n_comps, do_norm, clustering_algo)
    for fname in fnames:
        if np.array([get_res_path(fname).split('/')[-1] in elt for elt in os.listdir(resdir)]).any(): #dataset already computed for these settings
            continue
        print('\t\t'+fname)
        #raw = readRDS(path_rds+fname+'.rds')
        #adata_dict = {k:v for k, v in raw.items()}
        #adata = anndata.AnnData(X=adata_dict['counts'])
        #adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
        #del raw, adata_dict
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


