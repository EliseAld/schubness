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
import scanpy.external as sce
import gc
import time
import warnings
from sklearn import metrics
import sklearn
import math
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from sklearn.neighbors import NearestNeighbors
import umap
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

def _binary_search_perplexity(sqdistances, desired_perplexity, verbose):
    """Binary search for sigmas of conditional Gaussians.
    This approximation reduces the computational complexity from O(N^2) to
    O(uN).
    Parameters
    ----------
    sqdistances : array-like, shape (n_samples, n_neighbors)
        Distances between training samples and their k nearest neighbors.
        When using the exact method, this is a square (n_samples, n_samples)
        distance matrix. The TSNE default metric is "euclidean" which is
        interpreted as squared euclidean distance.
    desired_perplexity : float
        Desired perplexity (2^entropy) of the conditional Gaussians.
    verbose : int
        Verbosity level.
    Returns
    -------
    P : array, shape (n_samples, n_samples)
        Probabilities of conditional Gaussian distributions p_i|j.
    """
    # Maximum number of binary search steps
    n_steps = 100
    n_samples = sqdistances.shape[0]
    n_neighbors = sqdistances.shape[1]
    using_neighbors = n_neighbors < n_samples
    # Precisions of conditional Gaussian distributions
    beta_sum = 0.0
    # Use log scale
    desired_entropy = np.log(desired_perplexity)
    # This array is later used as a 32bit array. It has multiple intermediate
    # floating point additions that benefit from the extra precision
    P = np.zeros((n_samples, n_neighbors), dtype=np.float64)
    for i in range(n_samples):
        beta_min = -math.inf
        beta_max = math.inf
        beta = 1.0
        # Binary search of precision for i-th conditional distribution
        for l in range(n_steps):
            # Compute current entropy and corresponding probabilities
            # computed just over the nearest neighbors or over all data
            # if we're not using neighbors
            sum_Pi = 0.0
            for j in range(n_neighbors):
                if j != i or using_neighbors:
                    P[i, j] = np.exp(-sqdistances[i, j] * beta)
                    sum_Pi += P[i, j]
            if sum_Pi == 0.0:
                sum_Pi = EPSILON_DBL
            sum_disti_Pi = 0.0
            for j in range(n_neighbors):
                P[i, j] /= sum_Pi
                sum_disti_Pi += sqdistances[i, j] * P[i, j]
            entropy = np.log(sum_Pi) + beta * sum_disti_Pi
            entropy_diff = entropy - desired_entropy
            if np.abs(entropy_diff) <= PERPLEXITY_TOLERANCE:
                break
            if entropy_diff > 0.0:
                beta_min = beta
                if beta_max == math.inf:
                    beta *= 2.0
                else:
                    beta = (beta + beta_max) / 2.0
            else:
                beta_max = beta
                if beta_min == -math.inf:
                    beta /= 2.0
                else:
                    beta = (beta + beta_min) / 2.0
        beta_sum += beta
        if verbose and ((i + 1) % 1000 == 0 or i + 1 == n_samples):
            print("[t-SNE] Computed conditional probabilities for sample "
                  "%d / %d" % (i + 1, n_samples))
    if verbose:
        print("[t-SNE] Mean sigma: %f"
              % np.mean(math.sqrt(n_samples / beta_sum)))
    return P

def _joint_probabilities(distances, desired_perplexity, verbose=False):
    """Compute joint probabilities p_ij from distances.
    Parameters
    ----------
    distances : ndarray of shape (n_samples * n_neighbors,)
    desired_perplexity : float
        Desired perplexity of the joint probability distributions.
    verbose : int
        Verbosity level.
    Returns
    -------
    P : ndarray of shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    """
    # Compute conditional probabilities such that they approximately match
    # the desired perplexity
    conditional_P = _binary_search_perplexity(
        distances, desired_perplexity, verbose)
    P = conditional_P + conditional_P.T
    sum_P = np.maximum(np.sum(P), np.finfo(np.double).eps)
    P = np.maximum(squareform(P) / sum_P, np.finfo(np.double).eps)
    return P

def _kl_divergence(X_embedded, P, degrees_of_freedom):
    """t-SNE objective function: gradient of the KL divergence
    of p_ijs and q_ijs and the absolute error.
    Parameters
    ----------
    P : ndarray of shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    degrees_of_freedom : int
        Degrees of freedom of the Student's-t distribution.
    Returns
    -------
    kl_divergence : float
        Kullback-Leibler divergence of p_ij and q_ij.
    """
    #X_embedded = params.reshape(n_samples, n_components)

    # Q is a heavy-tailed distribution: Student's t-distribution
    dist = pdist(X_embedded, "sqeuclidean")
    dist /= degrees_of_freedom
    dist += 1.
    dist **= (degrees_of_freedom + 1.0) / -2.0
    Q = np.maximum(dist / (2.0 * np.sum(dist)), np.finfo(np.double).eps)

    # Optimization trick below: np.dot(x, y) is faster than
    # np.sum(x * y) because it calls BLAS

    # Objective: C (Kullback-Leibler divergence of P and Q)
    kl_divergence = 2.0 * np.dot(
            P, np.log(np.maximum(P, np.finfo(np.double).eps) / Q))

    return kl_divergence

def trustworthiness(X, X_embedded, metric='euclidean'):
    """Expresses to what extent the local structure is retained.
    The trustworthiness is within [0, 1]. It is defined as
    .. math::
        T(k) = 1 - \frac{2}{nk (2n - 3k - 1)} \sum^n_{i=1}
            \sum_{j \in \mathcal{N}_{i}^{k}} \max(0, (r(i, j) - k))
    where for each sample i, :math:`\mathcal{N}_{i}^{k}` are its k nearest
    neighbors in the output space, and every sample j is its :math:`r(i, j)`-th
    nearest neighbor in the input space. In other words, any unexpected nearest
    neighbors in the output space are penalised in proportion to their rank in
    the input space.
    * "Neighborhood Preservation in Nonlinear Projection Methods: An
      Experimental Study"
      J. Venna, S. Kaski
    * "Learning a Parametric Embedding by Preserving Local Structure"
      L.J.P. van der Maaten
    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features) or (n_samples, n_samples)
        If the metric is 'precomputed' X must be a square distance
        matrix. Otherwise it contains a sample per row.
    X_embedded : ndarray of shape (n_samples, n_components)
        Embedding of the training data in low-dimensional space.
    n_neighbors : int, default=5
        Number of neighbors k that will be considered.
    metric : str or callable, default='euclidean'
        Which metric to use for computing pairwise distances between samples
        from the original input space. If metric is 'precomputed', X must be a
        matrix of pairwise distances or squared distances. Otherwise, see the
        documentation of argument metric in sklearn.pairwise.pairwise_distances
        for a list of available metrics.
        .. versionadded:: 0.20
    Returns
    -------
    trustworthiness : float
        Trustworthiness of the low-dimensional embedding.
    """
    #n_neighbors = int(np.sqrt(X.shape[0]))
    n_neighbors = 5
    dist_X = metrics.pairwise_distances(X, metric=metric)
    # we set the diagonal to np.inf to exclude the points themselves from
    # their own neighborhood
    np.fill_diagonal(dist_X, np.inf)
    ind_X = np.argsort(dist_X, axis=1)
    # `ind_X[i]` is the index of sorted distances between i and other samples
    ind_X_embedded = NearestNeighbors(n_neighbors=n_neighbors).fit(
            X_embedded).kneighbors(return_distance=False)

    # We build an inverted index of neighbors in the input space: For sample i,
    # we define `inverted_index[i]` as the inverted index of sorted distances:
    # inverted_index[i][ind_X[i]] = np.arange(1, n_sample + 1)
    n_samples = X.shape[0]
    inverted_index = np.zeros((n_samples, n_samples), dtype=int)
    ordered_indices = np.arange(n_samples + 1)
    inverted_index[ordered_indices[:-1, np.newaxis],
                   ind_X] = ordered_indices[1:]
    ranks = inverted_index[ordered_indices[:-1, np.newaxis],
                           ind_X_embedded] - n_neighbors
    t = np.sum(ranks[ranks > 0])
    t = 1.0 - t * (2.0 / (n_samples * n_neighbors *
                          (2.0 * n_samples - 3.0 * n_neighbors - 1.0)))
    return t

def fuzzy_simplicial_set(
    X,
    n_neighbors=None,
    metric='euclidean',
    verbose=False):
    """Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        The data to be modelled as a fuzzy simplicial set.
    n_neighbors: int
        The number of neighbors to use to approximate geodesic distance.
    random_state: numpy RandomState or equivalent
        A state capable being used as a numpy random state.
    metric: string or function (optional, default 'euclidean')
        The metric to use to compute distances in high dimensional space.
    verbose: bool (optional, default False)
        Whether to report information on the current progress of the algorithm.
    return_dists: bool or None (optional, default None)
        Whether to return the pairwise distance associated with each edge.
    Returns
    -------
    fuzzy_simplicial_set: coo_matrix
        A fuzzy simplicial set represented as a sparse matrix. The (i,
        j) entry of the matrix represents the membership strength of the
        1-simplex between the ith and jth sample points.
    """
    local_connectivity = 1.0
    set_op_mix_ratio = 1.0
    if n_neighbors is None:
        n_neighbors = int(np.sqrt(X.shape[0]))
    knn_indices, knn_dists, _ = umap.umap_.nearest_neighbors(X, n_neighbors, metric, random_state=sklearn.utils.check_random_state(seed), verbose=verbose, angular=False, metric_kwds={})
    knn_dists = knn_dists.astype(np.float32)
    sigmas, rhos = umap.umap_.smooth_knn_dist(knn_dists, float(n_neighbors), local_connectivity=float(local_connectivity))
    rows, cols, vals = umap.umap_.compute_membership_strengths(knn_indices, knn_dists, sigmas, rhos)
    result = scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(X.shape[0], X.shape[0]))
    result.eliminate_zeros()
    transpose = result.transpose()
    prod_matrix = result.multiply(transpose)
    result = (set_op_mix_ratio * (result + transpose - prod_matrix) + (1.0 - set_op_mix_ratio) * prod_matrix)
    result.eliminate_zeros()
    return result

def cross_entropy(input, embedding, a=1.577, b=0.8951):
    V = fuzzy_simplicial_set(input).toarray()
    W = 1/(1+a*metrics.pairwise_distances(embedding)**(2*b))
    ce = np.sum(V*np.maximum(np.log(np.divide(np.maximum(V, EPSILON_DBL),
                                   np.maximum(W, EPSILON_DBL))),
                             EPSILON_DBL) +
                (1-V)*np.maximum(np.log(np.divide((1-V), np.maximum(1-W, EPSILON_DBL))),
                             EPSILON_DBL))
    return ce

def ti_analysis(adata,true_labels,do_norm,norm_scale, do_log,do_pca,
                n_clusters,metric,weighted,  #weighted adjmat for louvain/leiden clustering ?
                seed,n_comps,clustering_algo,n_iter,bootstrap_size):
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
        use_rep = 'X_pca'
        sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, n_comps)))
        X = adata.obsm['X_pca']
    else:
        print('pca not done!')
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
        sc.pl.paga(all_adata[kernel], show=False, random_state=seed, plot=False)
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
        sc.pl.paga(all_adata[method_name], show=False, random_state=seed, plot=False)
        sc.tl.umap(all_adata[method_name], init_pos="paga", random_state=seed)
    print('\t\t\tHubness and PAGA full pipeline:', round((time.time()-start)/60, 2), 'mn')
    start = time.time()
    ### Evaluate goodness of fit ###
    original = adata.X
    #P = _joint_probabilities(distances=metrics.pairwise_distances(original), desired_perplexity=15, verbose=False)
    cost = []
    #trustworth = []
    for method_name in all_adata.keys():
        X_paga = all_adata[method_name].obsm['X_umap']
        #ent, geo_d, geo_o = sce.cross_entropy_neighbors_in_rep(all_adata[method_name], use_rep='X_umap')
        #cost.append(_kl_divergence(X_paga, P, 1))
        cost.append(cross_entropy(original, X_paga))
        #trustworth.append(trustworthiness(original, X_paga))
    print('\t\t\tCompute match between original and PAGA:', round((time.time() - start) / 60, 2), 'mn')
    np.savetxt(get_res_path(fname)+"_costumap4.csv", cost, delimiter=',')
    #np.savetxt(get_res_path(fname)+"_trust5.csv", trustworth, delimiter=',')
#trust3 original = pca with 2 PCs
#trust4 original = pca with the same nb of PCs used to compute PAGA
#trust5 original = raw data after preprocess
#costkl with R code to compute tSNE cost function, original = raw data after preprocess
#cost2kl with python code to compute tSNE cost function, ddl=2, original = raw data after preprocess
#cost3kl with python code, ddl=1, original = raw data after preprocess
#costumap4 with umap cost function, original = raw data after preprocess

##### USAGE ######

#fixed params
resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned5/"
#resdir = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned5/"
#func to save results with proper name
get_res_path = lambda fname: f'{resdir}{fname}_norm{do_norm}_scale{norm_scale}_ncomps{n_comps}_{metric}_{clustering_algo}'
EPSILON_DBL = 1e-8
PERPLEXITY_TOLERANCE = 1e-5
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
#fnames = ['GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9', 'GSE60783_10', 'GSE67602_11']
#fnames = ['GSE70245_12', 'GSE52529_15', 'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33']


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
        if np.array([get_res_path(fname).split('/')[-1]+"_costumap" in elt for elt in os.listdir(resdir)]).any(): #dataset already computed for these settings
            continue
        print('\t\t'+fname)
        raw = readRDS(path_rds+fname+'.rds')
        adata_dict = {k:v for k, v in raw.items()}
        adata = anndata.AnnData(X=adata_dict['counts'])
        adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
        del raw, adata_dict
        #adata = anndata.read_h5ad(path_h5ad+fname+'.h5ad')

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


