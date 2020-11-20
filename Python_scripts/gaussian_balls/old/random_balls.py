import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.decomposition import PCA
import umap
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from skhubness.reduction import MutualProximity, LocalScaling, DisSimLocal
from skhubness.neighbors import NearestNeighbors
import anndata
import scanpy as sc
import gc

# parameters
seed=0

def reformat_array(arr):
    n_balls, n_points, original_dims = arr.shape
    input_points = []
    for n_ball in range(n_balls):
        for n_point in range(n_points):
            input_points.append(arr[n_ball, n_point, :])
    input_points = np.array(input_points)
    return input_points

def deformat_array(arr):
    n_points_tot, original_dims = arr.shape
    input_points = np.empty((n_balls, int(n_points_tot/n_balls), original_dims))
    for n_ball in range(n_balls):
        input_points[n_ball, :, :] = arr[range(n_ball*1000, (n_ball+1)*1000), :]
    return input_points

def get_pca(arr, n_components = 2):
    #print('in pca', arr.shape)
    n_balls, n_points, original_dims = arr.shape
    pca = PCA(n_components=n_components, random_state=seed)
    # input for pca should be (n_points, dims)
    # output of pca will be (n_points, n_components)
    input_points = reformat_array(arr)
    pca.fit(input_points)
    transformed_points = []
    for n_ball in range(n_balls):
        transformed_points.append(pca.transform(arr[n_ball, :, :], ))
    transformed = np.array(transformed_points)
    #print('out pca', transformed.shape)
    return transformed

def get_umap(arr, n_components = 2):
    #print('in umap', arr.shape)
    n_balls, n_points, original_dims = arr.shape
    reducer = umap.UMAP(n_components=n_components, random_state=seed)
    # input for umap should be (n_points, dims)
    # output of umap will be (n_points, n_components)
    input_points = reformat_array(arr)
    reducer.fit(input_points)
    transformed_points = []
    for n_ball in range(n_balls):
        transformed_points.append(reducer.transform(arr[n_ball, :, :], ))
    transformed = np.array(transformed_points)
    #print('out umap', transformed.shape)
    return transformed

def get_umap_dist(dist_arr, n_components = 2):
    #print('in umap', arr.shape)
    reducer = umap.UMAP(n_components=n_components, random_state=seed, metric='precomputed')
    # input for umap should be (n_points, n_points)
    # output of umap will be (n_points, n_components)
    #print('out umap', transformed.shape)
    transformed = reducer.fit_transform(dist_arr)
    return transformed

def plot_points(arr, idx0=0, idx1=1, title=None):
    #print('in plot', arr.shape)
    colors = ['dodgerblue', 'mediumblue', 'deeppink', 'orangered']
    for i in range(n_balls):
        plt.scatter(arr[i, :, idx0], arr[i, :, idx1], c=colors[i])
    plt.title(title)
    plt.show()

def add_dropout(arr, dropout_values=[0, 50, 80, 90]):
    dropout_percents = [elt*1./100 for elt in dropout_values]
    samples_dropout = np.array([arr for _ in range(len(dropout_percents))])
    for idx in range(len(dropout_percents)):
        samples_exp = np.copy(arr)
        samples_exp[np.random.uniform(0, 1, size=arr.shape) < dropout_percents[idx]] = 0
        samples_dropout[idx, :, :, :] = samples_exp
    return samples_dropout

def hub_reduction(arr, methods = {'nothing':None,'mp_normal':'normal','ls':'standard','ls_nicdm':'nicdm'}):
    input_points = reformat_array(arr)
    n_points_tot, original_dims = input_points.shape
    neigh_dist, neigh_ind = NearestNeighbors(n_neighbors=n_points_tot-1).fit(input_points).kcandidates()
    samples_exp = dict()
    for method_name, hubness_param in methods.items():
        if method_name == "nothing":
            samples_exp[method_name] = pairwise_distances(input_points)
        if method_name == "mp_normal":
            tmp1, tmp2 = MutualProximity(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[idx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "ls":
            tmp1, tmp2 = LocalScaling(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[idx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "ls_nicdm":
            tmp1, tmp2 = LocalScaling(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[idx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "dsl":
            tmp1, tmp2 = DisSimLocal(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[idx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
    return samples_exp

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

def hub_paga_adata(adata, methods = {'nothing':(None,None),
                                          'mp_normal':('mp',{'method': 'normal'}),
                                          'ls':('ls',None),
                                          'ls_nicdm':('ls',{'method': 'nicdm'})}, seed=seed):
    all_adata = dict()
    for method_name, (hubness, hubness_params) in methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsm[method_name] = kneighbors_graph(adata.X,
                                                            n_neighbors=10,
                                                            hubness=hubness,
                                                            hubness_params=hubness_params)
    for method_name, (hubness, hubness_params) in methods.items():
        sc.pp.neighbors(all_adata[method_name], n_neighbors=10, use_rep=method_name, random_state=seed)
        G, weights = generate_clustering_inputs(X=adata.X,
                                     metric='euclidean',
                                     n_neighbors=10,
                                     weighted=True,
                                     seed=seed,
                                     hubness=hubness,
                                     hubness_params=hubness_params)
        resol, weighted = getNclusters(all_adata[method_name], G, n_clusters=len(np.unique(adata.uns['Order'])), seed=seed, cluster_func='leiden', flavor='scanpy', weights=weights)
        sc.tl.leiden(all_adata[method_name], resolution=resol, use_weights=weighted, random_state=seed)
        sc.tl.paga(all_adata[method_name], groups="leiden")
    for method_name in all_adata.keys():
        del(all_adata[method_name].obsm[method_name])
    return all_adata

def plot_adata_leiden(adata, idx0=0, idx1=1, title=None, pca=True, seed=seed):
    if pca:
        arr = PCA(n_components=2, random_state=seed).fit_transform(adata.X)
    else:
        arr = adata.X
    cluster_id = adata.obs["leiden"]
    colors = ['dodgerblue', 'mediumblue', 'deeppink', 'orangered']
    for group in np.unique(cluster_id):
        mask = cluster_id==group
        plt.scatter(arr[mask, idx0], arr[mask, idx1], c=colors[int(group)])
    plt.title(title)
    plt.show()


# generate balls in 100D
n_balls = 3
original_dims = 100
n_points = 500
epsilon = 1
random_b = np.random.randint(10, 20, size=(original_dims,)).astype(np.float32)
random_vect = np.random.randint(0, 5, size=(original_dims,)).astype(np.float32)
random_a = (random_b + epsilon*random_vect)/np.linalg.norm(random_b + epsilon*random_vect)*np.linalg.norm(random_b)
random_c = (random_b - epsilon*random_vect)/np.linalg.norm(random_b - epsilon*random_vect)*np.linalg.norm(random_b)
print('A has '+str(np.sum(random_a<0))+' negative coordinates')
print('B has '+str(np.sum(random_b<0))+' negative coordinates')
print('C has '+str(np.sum(random_c<0))+' negative coordinates')
random_centers = np.transpose(np.stack((random_a, random_b, random_c), axis=1))
repeated_centers = np.array([random_centers for _ in range(n_points)])
#print('repeated_centers.shape', repeated_centers.shape)
repeated_centers = np.transpose(repeated_centers, (1, 0, 2))
#print('repeated_centers.shape', repeated_centers.shape)
samples = repeated_centers + np.random.normal(0, 0.1, size=repeated_centers.shape)
cluster_id = np.repeat(range(n_balls), n_points)
plot_points(samples, 1, 2, title='random axes, gaussian balls')
pca_samples = get_pca(samples)
umap_samples = get_umap(samples)
plot_points(pca_samples, title='pca of the gaussian balls')
plot_points(umap_samples, title='umap of the gaussian balls')

# add dropout
dropout_percents=[0, 10, 20, 30]
samples_dropout = add_dropout(samples, dropout_percents)
for idx in range(len(dropout_percents)):
    print(dropout_percents[idx])
    #print(np.linalg.norm(samples_exp - samples))
    pca_tmp = get_pca(samples_dropout[idx, :, :, :])
    plot_points(pca_tmp, title='pca of the gaussian balls, dropout '+str(dropout_percents[idx]))
    umap_tmp = get_umap(samples_dropout[idx, :, :, :])
    plot_points(umap_tmp, title='umap of the gaussian balls, dropout '+str(dropout_percents[idx]))

# correct graph
samples_hubred = dict()
for dropout_idx in tqdm(range(len(dropout_percents))):
    samples_hubred[dropout_percents[dropout_idx]] = hub_reduction(samples_dropout[dropout_idx])

# apply MDS
mds = MDS(n_components=2, dissimilarity='precomputed')
mds_samples_hubred = dict()
for dropout_idx in dropout_percents:
    sample_dropout = samples_hubred[dropout_idx]
    mds_samples_hubred[dropout_idx] = dict()
    for method_idx in tqdm(sample_dropout.keys()):
        sample_dropout_hubred = sample_dropout[method_idx]
        mat_upp = np.triu(sample_dropout_hubred)
        mat_low = np.tril(sample_dropout_hubred)
        mat_symm = (mat_upp + mat_low.T)/2 + (mat_upp.T + mat_low)/2
        mds_samples_hubred[dropout_idx][method_idx] = mds.fit_transform(mat_symm)
for key1 in mds_samples_hubred.keys():
    for key2 in mds_samples_hubred[key1].keys():
        plot_points(deformat_array(mds_samples_hubred[key1][key2]), title='MDS, dropout '+str(key1)+' and hub reduction '+key2)

# apply UMAP
umap_samples_hubred = dict()
for dropout_idx in dropout_percents:
    sample_dropout = samples_hubred[dropout_idx]
    umap_samples_hubred[dropout_idx] = dict()
    for method_idx in tqdm(sample_dropout.keys()):
        sample_dropout_hubred = sample_dropout[method_idx]
        mat_upp = np.triu(sample_dropout_hubred)
        mat_low = np.tril(sample_dropout_hubred)
        mat_symm = (mat_upp + mat_low.T)/2 + (mat_upp.T + mat_low)/2
        umap_samples_hubred[dropout_idx][method_idx] = get_umap_dist(mat_symm)
for key1 in umap_samples_hubred.keys():
    for key2 in umap_samples_hubred[key1].keys():
        plot_points(deformat_array(umap_samples_hubred[key1][key2]), title='UMAP, dropout '+str(key1)+' and hub reduction '+key2)


# apply PAGA
uns = {'Order':np.repeat(range(n_balls), repeats=n_points)}
samples_adata = dict()
for dropout_idx in range(len(dropout_percents)):
    dropout_tmp = reformat_array(samples_dropout[dropout_idx, :, :, :])
    samples_adata[dropout_percents[dropout_idx]] = dict()
    samples_adata[dropout_percents[dropout_idx]]['nothing'] = anndata.AnnData(X=dropout_tmp, uns=uns)
for dropout_idx in tqdm(samples_adata.keys()):
    samples_adata[dropout_idx] = hub_paga_adata(samples_adata[dropout_idx]['nothing'])

# evaluate clustering and order visually
for key1 in samples_adata.keys():
    for key2 in samples_adata[key1].keys():
        plot_adata_leiden(samples_adata[key1][key2], title='dropout '+str(key1)+' and hub reduction '+key2, pca=False)
for key1 in samples_adata.keys():
    result_over_method = []
    for key2 in samples_adata[key1].keys():
        result_over_method.append((samples_adata[key1][key2].uns['paga']['connectivities_tree'].toarray()==0).all())
    print(result_over_method)
# ???
