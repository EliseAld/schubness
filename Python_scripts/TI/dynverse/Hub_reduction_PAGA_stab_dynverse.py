from skhubness.neighbors import kneighbors_graph
import scanpy as sc
import anndata
import os
import skhubness
import gc
import numpy as np
import pandas as pd
from tqdm import tqdm
import csv

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
    adata = anndata.read_csv(path+'_data.csv')
    order = pd.read_csv(path+'_order.csv')
    adata.uns['Order'] = order.iloc[:,1].values
    return(adata)
def preprocess(adata, n_dim):
    recipe_seurat(adata) # ! remove datasets with less cells than genes! # change to seurat preprocess (Jo)
    if n_dim < adata.shape[1]:
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_dim)
        print("We have "+str(adata.obsm['X_pca'].shape[0])+" cells over "+str(adata.obsm['X_pca'].shape[1])+" PCs")
        return(adata)
    else:
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

def hub_paga_stab(adata, methods_set, path_out, n_iter=10, bootstrap=0.95):
    hubness_choice = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'}),
                   'ls':('ls',None),
                   'ls_nicdm':('ls',{'method': 'nicdm'}),
                   'dsl':('dsl',None)
                   }
    hubness_methods = prune_dict(hubness_choice,methods_set)
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
        sc.tl.leiden(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=0)
    all_iter = dict()
    cell_iter = dict()
    feat_iter = dict()
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_iter[method_name]=dict()
        cell_iter[method_name]=np.zeros((n_iter,adata.n_obs))
        feat_iter[method_name]=np.zeros((n_iter,adata.n_vars))
        for iter in tqdm(range(n_iter)):
            feat_bootstrap=np.random.uniform(0, 1, size=all_adata[method_name].n_vars)
            feat_bootstrap[feat_bootstrap<=bootstrap]=0
            feat_bootstrap[feat_bootstrap>bootstrap]=1
            feat_bootstrap = feat_bootstrap==0
            cell_bootstrap=np.random.uniform(0, 1, size=all_adata[method_name].n_obs)
            cell_bootstrap[cell_bootstrap<=bootstrap]=0
            cell_bootstrap[cell_bootstrap>bootstrap]=1
            cell_bootstrap = cell_bootstrap==0
            cell_iter[method_name][iter, :]=cell_bootstrap
            feat_iter[method_name][iter, :]=feat_bootstrap
            leiden = all_adata[method_name].obs['leiden'][cell_bootstrap]
            uns = {'Order':all_adata[method_name].uns['Order'][cell_bootstrap]}
            obsm = {'X_pca':adata.obsm['X_pca'][cell_bootstrap]}
            adata_sampled = anndata.AnnData(adata.X[cell_bootstrap][:, feat_bootstrap],
                                            obs=pd.DataFrame(leiden),
                                            uns=uns,
                                            obsm=obsm)
            adata_sampled.obsm[method_name] = kneighbors_graph(adata_sampled.obsm['X_pca'],
                                                            n_neighbors=10,
                                                            hubness=hubness,
                                                            hubness_params=hubness_params)
            sc.pp.neighbors(adata_sampled, n_neighbors=10, use_rep=method_name)
            sc.tl.paga(adata_sampled, groups="leiden")
            all_iter[method_name]['iter'+str(iter)]=adata_sampled.uns["paga"]["connectivities_tree"]
    for method_name in all_adata.keys():
        w = csv.writer(open(path_out+method_name+"_stab.csv", "w"))
        for key, val in all_iter[method_name].items():
            w.writerow([key, val])
        np.savetxt(path_out+method_name+"_stab_cell.csv", cell_iter[method_name], delimiter=',', fmt='%d')
        np.savetxt(path_out+method_name+"_stab_feat.csv", feat_iter[method_name], delimiter=',', fmt='%d')

def all_fast(dataset_folder,n_dim,methods_set):
    path = main_path+dataset_folder
    print('reading '+dataset_folder)
    adata = load(path)
    if n_dim < adata.shape[0]:
        adata = preprocess(adata, n_dim)
        if n_dim < adata.shape[1]:
            path_out = path.replace('csv/', 'h5/')+'_'+str(n_dim)+'dims_'
            hub_paga_stab(adata, methods_set, path_out)
        else:
            print('remove '+dataset_folder)
    else:
        print('remove '+dataset_folder)

if __name__ == '__main__':
    main_path = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/dynverse/"
    tools = ['dyngen', 'dyntoy', 'prosstt', 'splatter', 'gold'] # removed 'silver'
    dataset_choice = [tool+'/csv/'+f.name for tool in tools for f in os.scandir(main_path+tool+'/csv/') if f.name.endswith("_data.csv")]
    choice1 = input("Do you want to test all datasets? y/n\n")
    if choice1 == 'y':
        dataset_set = dataset_choice
    else:
        dataset_set = input("Which datasets? \n")
    dataset_tmp = [set.split("_") for set in dataset_set]
    for i in range(len(dataset_set)):
        tmp = dataset_tmp[i]
        empty = tmp[0]
        for j in range(1,len(tmp)-1):
            empty=empty+'_'+tmp[j]
        dataset_set[i]=empty
    methods_choice = ['nothing','mp_normal','ls','ls_nicdm','dsl']
    choice2 = input("Do you mant to test all methods? y/n\n")
    if choice2 == 'y':
        methods_set = methods_choice
    else:
        methods_set = input("Which methods? \n")
    n_dim = int(input("How many PCs do you mant? \n"))
    for set in dataset_set:
        all_fast(set,n_dim,methods_set)



