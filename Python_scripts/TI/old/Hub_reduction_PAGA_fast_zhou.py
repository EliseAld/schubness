from skhubness.neighbors import kneighbors_graph
import scanpy as sc
import anndata
import os
import skhubness
import gc
import numpy as np
import rpy2.robjects as ro
import anndata2ri
anndata2ri.activate()

def getNclusters(adata, G, n_clusters, seed, cluster_func, flavor, weights,range_min=0,range_max=3,max_steps=20):
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    weighted = weights is None
    while this_step < max_steps:
#         print('step ' + str(this_step))
        this_resolution = this_min + ((this_max-this_min)/2)
        if cluster_func == 'louvain':
            if flavor == 'scanpy':
                sc.tl.louvain(adata,resolution=this_resolution,random_state=seed,use_weights=weighted)
                clus = np.array(adata.obs['louvain']).astype(int)
                this_clusters = adata.obs['louvain'].nunique()
            elif flavor == 'base':
                louvain.set_rng_seed(seed)
                part_louvain = louvain.find_partition(graph = G, 
                                   partition_type = louvain.RBConfigurationVertexPartition, 
                                   weights = weights,
                                   resolution_parameter=this_resolution)
                clus = np.array(part_louvain.membership)
                this_clusters = len(np.unique(clus))
        elif cluster_func == 'leiden':
            if flavor == 'scanpy':
                sc.tl.leiden(adata,resolution=this_resolution,random_state=seed,use_weights=weighted)
                clus = np.array(adata.obs['leiden']).astype(int)
                this_clusters = adata.obs['leiden'].nunique()
            elif flavor == 'base':
                part_leiden = leidenalg.find_partition(graph = G, 
                                         partition_type = leidenalg.RBConfigurationVertexPartition, 
                                         weights = weights,
                                         resolution_parameter=this_resolution,
                                         seed=seed)
                clus = np.array(part_leiden.membership)
                this_clusters = len(np.unique(clus))
        else:
            raise ValueError("incorrect cluster_func, choose 'leiden' or 'louvain'")
#         print('got ' + str(this_clusters) + ' at resolution ' + str(this_resolution))
        if this_clusters > n_clusters:
            this_max = this_resolution
        elif this_clusters < n_clusters:
            this_min = this_resolution
        else:
            return(this_resolution, weighted)
        this_step += 1
    #print('Cannot find the number of clusters')
    #print('Clustering solution from last iteration is used:' + str(this_clusters) + ' at resolution ' + str(this_resolution))
    return(this_resolution, weighted)

def generate_clustering_inputs(X, metric, n_neighbors, weighted, seed, hubness, hubness_params):
    hub = skhubness.Hubness(k = n_neighbors, metric = metric, hubness=hubness, hubness_params=hubness_params,random_state=seed,store_k_occurrence=True,return_value='all').fit(X)
    #scores = hub.score()
    del hub.X_train_;gc.collect()
    knn = hub.nn_index_
    if weighted:
        adjmat = knn.kneighbors_graph(mode='distance')
        #affinity_matrix = 0.5 * (adjmat + adjmat.T)
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        try:
            weights = np.array(G.es["weight"]).astype(np.float64)
        except:
            weights = None
        if weights is not None:
            if not np.isfinite(np.sum(weights)): #weights failed
                adjmat = knn.kneighbors_graph(mode='connectivity')
                #affinity_matrix = 0.5 * (adjmat + adjmat.T)
                G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
                weights = None
    else:
        adjmat = knn.kneighbors_graph(mode='connectivity')
        #affinity_matrix = 0.5 * (adjmat + adjmat.T)
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        weights = None
    del hub;gc.collect()
    return(G, weights)

def recipe_seurat(adata,do_log=True,norm_scale=True):
    #same as scanpy clustering tutorial except initial cells/genes prefiltering /!\ sc.pp.recipe_seurat filters non log data ?
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata._inplace_subset_var(adata.var.highly_variable)
    if norm_scale:
        sc.pp.scale(adata, max_value=10)

def load(path):
    raw = readRDS(path)
    adata_dict = {k:v for k, v in raw.items()}
    adata = anndata.AnnData(X=adata_dict['counts'])
    adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
    return(adata)

def preprocess(adata, n_dim):
    recipe_seurat(adata) # ! remove datasets with less cells than genes! # change to seurat preprocess (Jo)
    if n_dim < adata.shape[1]:
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_dim)
        print("We have "+str(adata.obsm['X_pca'].shape[0])+" cells over "+str(adata.obsm['X_pca'].shape[1])+" PCs")
        return(adata)
    else:
        return(adata)

def hub_paga(adata,path_out):
    hubness_methods = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'}),
                   'ls':('ls',None),
                   'ls_nicdm':('ls',{'method': 'nicdm'}),
                   'dsl':('dsl',None)
                   }
    all_adata = dict()
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsm[method_name] = kneighbors_graph(all_adata.get(method_name).obsm['X_pca'],
                                                            n_neighbors=10,
                                                            hubness=hubness,
                                                            hubness_params=hubness_params)
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        sc.pp.neighbors(all_adata[method_name], n_neighbors=10, use_rep=method_name)
        G, weights = generate_clustering_inputs(X=adata.obsm['X_pca'],
                                     metric='euclidean',
                                     n_neighbors=10,
                                     weighted=True,
                                     seed=0,
                                     hubness=hubness,
                                     hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=len(np.unique(adata.uns['Order'])), seed=0, cluster_func='leiden', flavor='scanpy', weights=weights)
        sc.tl.leiden(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=1)
        sc.tl.paga(all_adata[method_name], groups="leiden")
    for method_name in all_adata.keys():
        del(all_adata[method_name].obsm[method_name])
    for method_name in all_adata.keys():
        all_adata[method_name].write_h5ad(filename=path_out+method_name+"2.h5ad") # v2 since iter=1 and fitted param for clustering

def all_fast(dataset_folder, n_dim):
    path = main_path+dataset_folder
    print('reading '+dataset_folder)
    adata = load(path)
    if n_dim < adata.shape[0]:
        adata = preprocess(adata, n_dim)
        if n_dim < adata.shape[1]:
            path_out = path.replace('csv/', 'h5/')+'_'+str(n_dim)+'dims_'
            hub_paga(adata, path_out)
        else:
            print('remove '+dataset_folder)
    else:
        print('remove '+dataset_folder)

if __name__ == '__main__':
    main_path = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/"
    dataset_set = [f.name for f in os.scandir(main_path+'/rds/') if f.name.endswith(".rds")]
    methods_set = ['nothing','mp_normal','ls','ls_nicdm','dsl']
    n_dim = [25, 50, 100, 500]
    readRDS = ro.r['readRDS']
    for dim in n_dim:
        for set in dataset_set:
            all_fast(set, dim)
