from skhubness.neighbors import kneighbors_graph
import scanpy as sc
import anndata
import os
import csv
from tqdm import tqdm
import skhubness
import gc
import numpy as np

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
    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' + str(this_clusters) + ' at resolution ' + str(this_resolution))
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

def load_preprocess(path, n_dim):
    adata = anndata.read_h5ad(path)
    sc.pp.recipe_zheng17(adata)
    if n_dim <= adata.shape[0]:
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_dim)
        print("We have "+str(adata.obsm['X_pca'].shape[0])+" cells over "+str(adata.obsm['X_pca'].shape[1])+" PCs")
    else:
        print('less cells than dim')
    return(adata)

def prune_dict(dictionary, mask):
    result = {}
    for key in dictionary:
        if isinstance(dictionary[key], dict):
            newmask = [maskname[mask.find(".") + 1:] for maskname in mask if maskname.startswith(key + ".")]
            result[k] = filterstage1(dictionary[key], newmask)
        elif key in mask:
            result[key] = dictionary[key]
    return result

def hub_paga(adata, methods_set,path_out, n_dim, n_iter=10):
    hubness_choice = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'}),
                   'ls':('ls',None),
                   'ls_nicdm':('ls',{'method': 'nicdm'}),
                   'dsl':('dsl',None)
                   }
    hubness_methods = prune_dict(hubness_choice,methods_set)
    if n_dim <= adata.shape[0]:
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
            resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=len(np.unique(adata.uns['Order'])), seed=1, cluster_func='leiden', flavor='scanpy', weights=weights)
        #all_paga = dict()
            sc.tl.leiden(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=1)
            sc.tl.paga(all_adata[method_name], groups="leiden")
            #all_paga[method_name] = dict()
            #for i in tqdm(range(n_iter)):
            #    sc.tl.paga(all_adata[method_name], groups="leiden")
            #    all_paga[method_name]["iter"+str(i)] = all_adata[method_name].uns["paga"]["connectivities_tree"]
        for method_name in all_adata.keys():
            del(all_adata[method_name].obsm[method_name])
        for method_name in all_adata.keys():
            all_adata[method_name].write_h5ad(filename=path_out+method_name+"2.h5ad")
            #w = csv.writer(open(path_out+method_name+".csv", "w"))
            #for key, val in all_paga[method_name].items():
            #    w.writerow([key, val])

def all_fast(dataset_folder,n_dim,methods_set):
    path = main_path+dataset_folder+'/data.h5ad'
    adata = load_preprocess(path, n_dim)
    path_out = main_path+dataset_folder+'/'+str(n_dim)+'_dims_'
    hub_paga(adata,methods_set,path_out, n_dim)

if __name__ == '__main__':
    main_path = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/"
    print('Available data sets:')
    print([f.name for f in os.scandir(main_path) if f.is_dir()])
    dataset_choice = ['GSE36552','GSE45719','GSE52529','GSE52583','GSE52583b','GSE59114','GSE60783','GSE64447',
                      'GSE67123','GSE67602','GSE69761','GSE70245','GSE74767','GSE74767b','GSE75330','GSE75330b',
                      'GSE75748','GSE76408','GSE85066','GSE86146','GSE87375','GSE87375b','GSE90047','GSE90860',
                      'GSE92332','GSE92332b','GSE93421','GSE94641','GSE95753','GSE95753b','GSE97391','GSE97391b',
                      'GSE98451','GSE98664','GSE99933','GSE102066','GSE103633','GSE106587','GSE107122','GSE107910',
                      'GSE109774','GSE109774b']
    choice1 = input("Do you mant to test all datasets? y/n\n")
    if choice1=='y':
        dataset_set = dataset_choice
    else:
        dataset_set = input("Which datasets? \n")
    methods_choice = ['nothing','mp_normal','ls','ls_nicdm','dsl']
    choice2 = input("Do you mant to test all methods? y/n\n")
    if choice2=='y':
        methods_set = methods_choice
    else:
        methods_set = input("Which methods? \n")
    n_dim = int(input("How many PCs do you mant? \n"))
    for set in dataset_choice:
        all_fast(set,n_dim,methods_set)
