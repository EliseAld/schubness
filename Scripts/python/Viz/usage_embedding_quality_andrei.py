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
from sklearn import metrics
import sklearn
from sklearn.neighbors import NearestNeighbors
from umap import UMAP
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

def natural_pca(data, metric):
    pairwise_dist = sklearn.metrics.pairwise_distances(data, metric=metric)
    nPC = pairwise_dist.shape[0]
    natural_PC = np.zeros(shape=(2, nPC), dtype=int)
    for i in range(nPC):
        if i==0:
            natural_PC[:, i] = [np.where(pairwise_dist == np.amax(pairwise_dist))[0].tolist()[0], np.where(pairwise_dist == np.amax(pairwise_dist))[1].tolist()[0]]
            pairwise_dist[natural_PC[0, i], natural_PC[1, i]] = 0
            pairwise_dist[natural_PC[1, i], natural_PC[0, i]] = 0
        else:
            dist_to_consider = pairwise_dist[[k.tolist()[0] for k in natural_PC[:, :i]]]
            candidates_i = [np.where(pairwise_dist == np.amax(dist_to_consider))[0].tolist()[0], np.where(pairwise_dist == np.amax(dist_to_consider))[1].tolist()[0]]
            natural_PC[0, i] = [k for k in candidates_i if k not in [t.tolist()[0] for t in natural_PC[:, :i]]][0]
            dist_to_consider0 = [max(pairwise_dist[natural_PC[0, i], natural_PC[0, k]],
                                        pairwise_dist[natural_PC[0, k], natural_PC[0, i]]) for k in range(i)]
            dist_to_consider1 = [max(pairwise_dist[natural_PC[0, i], natural_PC[1, k]],
                                        pairwise_dist[natural_PC[1, k], natural_PC[0, i]]) for k in range(i)]
            dist_to_consider = np.append(dist_to_consider0, dist_to_consider1)
            idx = np.where(dist_to_consider == np.amin(dist_to_consider))[0].tolist()[0]
            if idx < i:
                natural_PC[1, i] = natural_PC[0, idx]
            else:
                natural_PC[1, i] = natural_PC[1, idx-i]
            for k in range(i):
                pairwise_dist[natural_PC[0, k], natural_PC[0, i]] = 0
                pairwise_dist[natural_PC[1, k], natural_PC[0, i]] = 0
                pairwise_dist[natural_PC[0, i], natural_PC[0, k]] = 0
                pairwise_dist[natural_PC[0, i], natural_PC[1, k]] = 0
                pairwise_dist[natural_PC[0, k], natural_PC[1, i]] = 0
                pairwise_dist[natural_PC[1, k], natural_PC[1, i]] = 0
                pairwise_dist[natural_PC[1, i], natural_PC[0, k]] = 0
                pairwise_dist[natural_PC[1, i], natural_PC[1, k]] = 0
            pairwise_dist[natural_PC[0, i], natural_PC[1, i]] = 0
            pairwise_dist[natural_PC[1, i], natural_PC[0, i]] = 0
    return natural_PC

def QDM(input, embedding, metric):
    useful_points = natural_pca(input, metric)
    dis_i = sklearn.metrics.pairwise_distances(input, metric=metric)
    dis_o = sklearn.metrics.pairwise_distances(embedding, metric=metric)
    ravel_i = []
    ravel_o = []
    for i in range(useful_points.shape[1]):
        ravel_i.append(dis_i[useful_points[0, i], useful_points[1, i]])
        ravel_o.append(dis_o[useful_points[0, i], useful_points[1, i]])
    corr = np.corrcoef(ravel_i, ravel_o)
    return corr[0, 1]

def QNP(input, embedding, metric, k):
    NN = NearestNeighbors(n_neighbors=k, metric=metric)
    neigh_i = NN.fit(input).kneighbors(return_distance=False)
    neigh_o = NN.fit(embedding).kneighbors(return_distance=False)
    qnp = np.sum([len(np.intersect1d(neigh_i[i, :], neigh_o[i, :])) for i in range(neigh_i.shape[0])])/k/input.shape[0]
    return qnp

def QGC(input, embedding, metric, label, k):
    NN = NearestNeighbors(n_neighbors=k, metric=metric)
    qgc = np.zeros((2, len(np.unique(label))), dtype=int)
    neigh_i = NN.fit(input).kneighbors(return_distance=False)
    neigh_o = NN.fit(embedding).kneighbors(return_distance=False)
    for idx, clust in enumerate(np.unique(label)):
        n_clust = np.sum(label == clust)
        c_i = []
        c_o= []
        for i in range(input.shape[0]):
            if label[i] == clust:
                c_i.append(np.sum(label[neigh_i[i, :]] == clust))
                c_o.append(np.sum(label[neigh_o[i, :]] == clust))
        qgc[0, idx] = np.sum(c_i)/k/n_clust
        qgc[1, idx] = np.sum(c_o)/k/n_clust
    return qgc

def embedding_analysis(adata, do_norm, norm_scale, do_log,do_pca,
                n_clusters, metric, weighted,
                seed, n_comps, clustering_algo):
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
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, 500)))
        original1 = adata.obsm['X_pca']
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, n_comps)))
        X = adata.obsm['X_pca']
    else:
        print('pca not done!')
        X = adata.X
    n_neighbors = int(np.sqrt(X.shape[0]))
    print('\t\t\tPreprocessing and PCA pipeline:', round((time.time()-start0)/60, 2), 'mn')
    ### Hub reduction and clustering ###
    start = time.time()
    all_adata = dict()
    for kernel in ['gauss', 'umap']:
        all_adata[kernel] = adata.copy()
        try:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors + 1, metric=metric, use_rep='X', method=kernel)
        except:
            sc.pp.neighbors(all_adata[kernel], n_neighbors=n_neighbors + 1, metric=metric, use_rep='X', method=kernel,
                            knn=False)
        G, weights = generate_clustering_inputs(X=X,
                                                metric=metric,
                                                n_neighbors=n_neighbors,
                                                weighted=weighted,
                                                seed=seed,
                                                hubness=None,
                                                hubness_params=None)
        resol, weighted = getNclusters(all_adata[kernel], G, n_clusters=n_clusters, seed=seed,
                                       clustering_algo=clustering_algo,
                                       flavor='scanpy', weights=weights)
        if clustering_algo == "leiden":
            sc.tl.leiden(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
        elif clustering_algo == "louvain":
            sc.tl.louvain(all_adata[kernel], resolution=resol, use_weights=weighted, random_state=seed)
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
        G, weights = generate_clustering_inputs(X=X,
                                                metric=metric,
                                                n_neighbors=n_neighbors,
                                                weighted=weighted,
                                                seed=seed,
                                                hubness=hubness,
                                                hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=n_clusters, seed=seed,
                                       clustering_algo=clustering_algo,
                                       flavor='base', weights=weights)
        if clustering_algo == "louvain":
            clus = np.array(louvain.find_partition(graph=G,
                                                   partition_type=louvain.RBConfigurationVertexPartition,
                                                   weights=weights,
                                                   resolution_parameter=resol, seed=seed).membership)
            all_adata[method_name].obs['louvain'] = pd.Categorical(values=clus.astype('U'),
                                                                   categories=natsorted(map(str, np.unique(clus))), )
        elif clustering_algo == "leiden":
            clus = np.array(leidenalg.find_partition(graph=G,
                                                     partition_type=leidenalg.RBConfigurationVertexPartition,
                                                     weights=weights,
                                                     resolution_parameter=resol,
                                                     seed=seed).membership)
            all_adata[method_name].obs['leiden'] = pd.Categorical(values=clus.astype('U'),
                                                                  categories=natsorted(map(str, np.unique(clus))), )
    original0 = adata.X
    original2 = adata.obsm['X_pca'][:, :2]
    print('\t\t\tkNN graph & clustering pipeline:', round((time.time() - start) / 60, 2), 'mn')
    ### tSNE embedding ###
    start = time.time()
    tsne = sklearn.manifold.TSNE(n_components=2, metric='precomputed', random_state=seed, perplexity=60.0)
    q_tsne = np.empty((6, len(all_adata.keys())))
    for idx, method_name in enumerate(all_adata.keys()):
        all_adata[method_name].obsm['X_tsne'] = tsne.fit_transform(all_adata[method_name].obsp['distances'].toarray())
        q_tsne[0, idx] = QDM(original0, all_adata[method_name].obsm['X_tsne'], metric)
        q_tsne[1, idx] = QDM(original1, all_adata[method_name].obsm['X_tsne'], metric)
        q_tsne[2, idx] = QDM(original2, all_adata[method_name].obsm['X_tsne'], metric)
        q_tsne[3, idx] = QNP(original0, all_adata[method_name].obsm['X_tsne'], metric, n_neighbors)
        q_tsne[4, idx] = QNP(original1, all_adata[method_name].obsm['X_tsne'], metric, n_neighbors)
        q_tsne[5, idx] = QNP(original2, all_adata[method_name].obsm['X_tsne'], metric, n_neighbors)
    print('\t\t\ttSNE embedding pipeline:', round((time.time() - start) / 60, 2), 'mn')
    ### UMAP embedding ###
    start = time.time()
    umap = UMAP(n_components=2, metric='precomputed', random_state=seed)
    q_umap = np.empty((6, len(all_adata.keys())))
    for idx, method_name in enumerate(all_adata.keys()):
        all_adata[method_name].obsm['X_umap_'] = umap.fit_transform(all_adata[method_name].obsp['distances'].toarray())
        q_umap[0, idx] = QDM(original0, all_adata[method_name].obsm['X_umap_'], metric)
        q_umap[1, idx] = QDM(original1, all_adata[method_name].obsm['X_umap_'], metric)
        q_umap[2, idx] = QDM(original2, all_adata[method_name].obsm['X_umap_'], metric)
        q_umap[3, idx] = QNP(original0, all_adata[method_name].obsm['X_umap_'], metric, n_neighbors)
        q_umap[4, idx] = QNP(original1, all_adata[method_name].obsm['X_umap_'], metric, n_neighbors)
        q_umap[5, idx] = QNP(original2, all_adata[method_name].obsm['X_umap_'], metric, n_neighbors)
    print('\t\t\tUMAP embedding pipeline:', round((time.time() - start) / 60, 2), 'mn')
    ### PAGA embedding ###
    start = time.time()
    q_paga_umap = np.empty((6, len(all_adata.keys())))
    for idx, method_name in enumerate(all_adata.keys()):
        sc.tl.paga(all_adata[method_name], groups=clustering_algo)
        sc.pl.paga(all_adata[method_name], show=False, random_state=seed, plot=False)
        sc.tl.umap(all_adata[method_name], init_pos="paga", random_state=seed)
        q_paga_umap[0, idx] = QDM(original0, all_adata[method_name].obsm['X_umap'], metric)
        q_paga_umap[1, idx] = QDM(original1, all_adata[method_name].obsm['X_umap'], metric)
        q_paga_umap[2, idx] = QDM(original2, all_adata[method_name].obsm['X_umap'], metric)
        q_paga_umap[3, idx] = QNP(original0, all_adata[method_name].obsm['X_umap'], metric, n_neighbors)
        q_paga_umap[4, idx] = QNP(original1, all_adata[method_name].obsm['X_umap'], metric, n_neighbors)
        q_paga_umap[5, idx] = QNP(original2, all_adata[method_name].obsm['X_umap'], metric, n_neighbors)
    print('\t\t\tPAGA+UMAP embedding pipeline:', round((time.time()-start)/60, 2), 'mn')
    ### Save ###
    np.savetxt(get_res_path(fname)+"_tsne_q.csv", q_tsne, delimiter=',')
    np.savetxt(get_res_path(fname)+"_umap_q.csv", q_umap, delimiter=',')
    np.savetxt(get_res_path(fname)+"_paga_q.csv", q_paga_umap, delimiter=',')
    print('\t\t\tFull pipeline:', round((time.time()-start0)/60, 2), 'mn')

##### USAGE ######

#fixed params
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/h5_jo_ktuned7/"
#resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/cytotrace/h5_jo_ktuned7/"
#func to save results with proper name
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'
do_log = True #already done with do_norm
do_pca = True
weighted = True
norm_scale = True
clustering_algo = 'leiden'
do_norm = 'seurat'
seed = 0
fnames = ['gold_hematopoiesis-gates_olsson', 'gold_stimulated-dendritic-cells-LPS_shalek',
          'gold_developing-dendritic-cells_schlitzer', 'gold_germline-human-male_guo',
          'gold_aging-hsc-young_kowalczyk', 'gold_germline-human-female-weeks_li',
          'gold_germline-human-female_guo', 'gold_mESC-differentiation_hayashi',
          'gold_germline-human-male-weeks_li', 'gold_pancreatic-beta-cell-maturation_zhang', 'gold_human-embryos_petropoulos',
          'gold_myoblast-differentiation_trapnell', 'gold_pancreatic-alpha-cell-maturation_zhang',
          'gold_aging-hsc-old_kowalczyk']
#fnames = ['GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9', 'GSE60783_10', 'GSE67602_11', 'GSE70245_12', 'GSE52529_15', 'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33']


#vary params
metric = ('cosine', 'euclidean')
#metric = ('euclidean', 'cosine')
n_comps = (50, 500)
#n_comps = (500, 50)
params_list = list(product(metric, n_comps))

readRDS = ro.r['readRDS']
path_rds = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/rds/"
path_h5ad = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/raw/raw_"
for metric, n_comps in params_list:
    print(metric, n_comps, do_norm, clustering_algo)
    for fname in fnames:
        if np.array([get_res_path(fname).split('/')[-1]+"_tsne_t.csv" in elt for elt in os.listdir(resdir)]).any(): #dataset already computed for these settings
            continue
        print('\t\t'+fname)
        raw = readRDS(path_rds+fname+'.rds')
        adata_dict = {k:v for k, v in raw.items()}
        adata = anndata.AnnData(X=adata_dict['counts'])
        adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
        del raw, adata_dict
        #adata = anndata.read_h5ad(path_h5ad+fname+'.h5ad')

        embedding_analysis(
            adata,
            do_norm=do_norm,
            norm_scale=norm_scale,
            do_log=do_log,
            do_pca=do_pca,
            n_clusters=len(np.unique(adata.uns['Order'])),
            metric=metric,
            weighted=weighted, 
            seed=seed,
            n_comps=n_comps,
            clustering_algo=clustering_algo)


