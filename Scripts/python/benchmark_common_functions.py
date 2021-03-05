import louvain
import leidenalg
import scanpy as sc
import gc
import numpy as np
import skhubness
import anndata
import rpy2.robjects as ro
import anndata2ri
import warnings
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
                part_louvain = louvain.find_partition(graph=G,
                                                      partition_type=louvain.RBConfigurationVertexPartition,
                                                      weights=weights,
                                                      resolution_parameter=this_resolution, seed=seed)
                clus = np.array(part_louvain.membership)
                this_clusters = len(np.unique(clus))
        elif clustering_algo == 'leiden':
            if flavor == 'scanpy':
                sc.tl.leiden(adata, resolution=this_resolution, random_state=seed, use_weights=weighted)
                clus = np.array(adata.obs['leiden']).astype(int)
                this_clusters = adata.obs['leiden'].nunique()
            elif flavor == 'base':
                part_leiden = leidenalg.find_partition(graph=G,
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
            return this_resolution, weighted
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
            if not np.isfinite(np.sum(weights)):  # weights failed
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
    # same as scanpy clustering tutorial except initial cells/genes prefiltering /!\ sc.pp.recipe_seurat filters non log data ?
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata._inplace_subset_var(adata.var.highly_variable)
    if norm_scale:
        sc.pp.scale(adata, max_value=10)


def load_rds(path, fname):
    readRDS = ro.r['readRDS']
    raw = readRDS(path + fname + '.rds')
    adata_dict = {k: v for k, v in raw.items()}
    adata = anndata.AnnData(X=adata_dict['counts'])
    adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
    del raw, adata_dict
    return adata


def load_h5ad(path, fname):
    return anndata.read_h5ad(path + fname + '.h5ad')