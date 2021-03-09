import scanpy as sc
import anndata
import rpy2.robjects as ro
import anndata2ri
import louvain
import leidenalg
import numpy as np
import pandas as pd
import scipy
from scipy.stats import skew
import sklearn.metrics.cluster
import sklearn.cluster
import skhubness
import gc
import time
import array
import csv
import warnings

warnings.filterwarnings("ignore")
anndata2ri.activate()


class IncrementalCSRMatrix(object):
    """modified from https://maciejkula.github.io/2015/02/22/incremental-construction-of-sparse-matrices/"""

    def __init__(self, shape, dtype):

        if dtype is np.int32:
            type_flag = "i"
        elif dtype is np.int64:
            type_flag = "l"
        elif dtype is np.float32:
            type_flag = "f"
        elif dtype is np.float64:
            type_flag = "d"
        else:
            raise Exception("Dtype not supported.")

        self.dtype = dtype
        self.shape = shape

        self.rows = array.array("i")
        self.cols = array.array("i")
        self.data = array.array(type_flag)

    def append(self, i, j, v):

        m, n = self.shape

        if i >= m or j >= n:
            raise Exception("Index out of bounds")

        self.rows.append(i)
        self.cols.append(j)
        self.data.append(v)

    def tocsr(self):

        rows = np.frombuffer(self.rows, dtype=np.int32)
        cols = np.frombuffer(self.cols, dtype=np.int32)
        data = np.frombuffer(self.data, dtype=self.dtype)

        return scipy.sparse.csr_matrix((data, (rows, cols)), shape=self.shape)


def get_csr(f, shape):
    # shape = (n_dims, n_samples)
    mat = IncrementalCSRMatrix(shape, np.float32)

    for row_ind, line in enumerate(f, start=-1):
        print(row_ind, end="\r")
        if row_ind == -1:
            continue

        data = np.array(line[1:], dtype=float)
        idx = data.nonzero()[0]
        for col_ind in idx:
            mat.append(row_ind, col_ind, data[col_ind])

    return mat


def color_sign_df(value):
    """
    Colors elements in a dateframe
    green if positive and red if
    negative. Does not color NaN
    values.
    """

    if value < 0:
        color = "red"
    elif value > 0:
        color = "green"
    else:
        color = "black"

    return "color: %s" % color


def angular_dist(u, v, w=None):
    cosine_dist = scipy.spatial.distance.cosine(u, v, w)
    cosine_sim = -cosine_dist + 1
    return np.arccos(cosine_sim) / np.pi


def get_complexity(data, labels):
    """ Calculate dataset complexity as defined in Abdelaal et al."""

    unique_types = np.unique(labels)

    # calculate avg expression of each cell type
    avg_exprs = np.zeros((len(unique_types), data.shape[1]))
    for i, cell_type in enumerate(unique_types):
        avg_exprs[i, :] = np.mean(data[np.where(labels == cell_type)], axis=0)

    # correlation between avg expression of each cell type
    coeffs = np.corrcoef(avg_exprs)
    np.fill_diagonal(coeffs, 0)
    # mean of max correlation
    complexity = np.mean(np.max(coeffs, axis=1))
    return complexity


def getNclusters(
    adata,
    G,
    n_clusters,
    seed,
    clustering_algo,
    flavor,
    weights,
    range_min=0,
    range_max=3,
    max_steps=20,
):
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    weighted = weights is None
    while this_step < max_steps:
        this_resolution = this_min + ((this_max - this_min) / 2)
        if clustering_algo == "louvain":
            if flavor == "scanpy":
                sc.tl.louvain(
                    adata,
                    resolution=this_resolution,
                    random_state=seed,
                    use_weights=weighted,
                )
                clus = np.array(adata.obs["louvain"]).astype(int)
                this_clusters = adata.obs["louvain"].nunique()
            elif flavor == "base":
                part_louvain = louvain.find_partition(
                    graph=G,
                    partition_type=louvain.RBConfigurationVertexPartition,
                    weights=weights,
                    resolution_parameter=this_resolution,
                    seed=seed,
                )
                clus = np.array(part_louvain.membership)
                this_clusters = len(np.unique(clus))
        elif clustering_algo == "leiden":
            if flavor == "scanpy":
                sc.tl.leiden(
                    adata,
                    resolution=this_resolution,
                    random_state=seed,
                    use_weights=weighted,
                )
                clus = np.array(adata.obs["leiden"]).astype(int)
                this_clusters = adata.obs["leiden"].nunique()
            elif flavor == "base":
                part_leiden = leidenalg.find_partition(
                    graph=G,
                    partition_type=leidenalg.RBConfigurationVertexPartition,
                    weights=weights,
                    resolution_parameter=this_resolution,
                    seed=seed,
                )
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


def getNclusters2(
    adata,
    G,
    n_clusters,
    seed,
    cluster_func,
    flavor,
    weights,
    range_min=0,
    range_max=3,
    max_steps=20,
):
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)

    weighted = weights is None

    while this_step < max_steps:
        #         print('step ' + str(this_step))
        this_resolution = this_min + ((this_max - this_min) / 2)

        if cluster_func == "louvain":

            if flavor == "scanpy":
                sc.tl.louvain(
                    adata,
                    resolution=this_resolution,
                    random_state=seed,
                    use_weights=weighted,
                )
                clus = np.array(adata.obs["louvain"]).astype(int)
                this_clusters = adata.obs["louvain"].nunique()

            elif flavor == "base":
                louvain.set_rng_seed(seed)
                part_louvain = louvain.find_partition(
                    graph=G,
                    partition_type=louvain.RBConfigurationVertexPartition,
                    weights=weights,
                    resolution_parameter=this_resolution,
                )
                clus = np.array(part_louvain.membership)
                this_clusters = len(np.unique(clus))

        elif cluster_func == "leiden":

            if flavor == "scanpy":
                sc.tl.leiden(
                    adata,
                    resolution=this_resolution,
                    random_state=seed,
                    use_weights=weighted,
                )
                clus = np.array(adata.obs["leiden"]).astype(int)
                this_clusters = adata.obs["leiden"].nunique()

            elif flavor == "base":
                part_leiden = leidenalg.find_partition(
                    graph=G,
                    partition_type=leidenalg.RBConfigurationVertexPartition,
                    weights=weights,
                    resolution_parameter=this_resolution,
                    seed=seed,
                )
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
            return clus, dict(resolution=this_resolution, succeeded=True)
        this_step += 1

    print("Cannot find the number of clusters")
    print(
        "Clustering solution from last iteration is used:"
        + str(this_clusters)
        + " at resolution "
        + str(this_resolution)
    )

    return clus, dict(resolution=this_resolution, succeeded=False)


def generate_clustering_inputs(
    X, metric, n_neighbors, weighted, seed, hubness, hubness_params
):
    hub = skhubness.Hubness(
        k=n_neighbors,
        metric=metric,
        hubness=hubness,
        hubness_params=hubness_params,
        random_state=seed,
        store_k_occurrence=True,
        return_value="all",
    ).fit(X)
    del hub.X_train_
    gc.collect()
    knn = hub.nn_index_
    if weighted:
        adjmat = knn.kneighbors_graph(mode="distance")
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        try:
            weights = np.array(G.es["weight"]).astype(np.float64)
        except:
            weights = None
        if weights is not None:
            if not np.isfinite(np.sum(weights)):  # weights failed
                adjmat = knn.kneighbors_graph(mode="connectivity")
                G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
                weights = None
    else:
        adjmat = knn.kneighbors_graph(mode="connectivity")
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        weights = None
    del hub
    gc.collect()
    return G, weights


def generate_clustering_inputs2(
    X, metric, n_neighbors, weighted, seed, hubness, hubness_params
):

    hub = skhubness.Hubness(
        k=n_neighbors,
        metric=metric,
        hubness=hubness,
        hubness_params=hubness_params,
        random_state=seed,
        store_k_occurrence=True,
        return_value="all",
    ).fit(X)
    scores = hub.score()
    del hub.X_train_
    gc.collect()

    knn = hub.nn_index_

    if weighted:
        adjmat = knn.kneighbors_graph(mode="distance")

        affinity_matrix = 0.5 * (adjmat + adjmat.T)
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        try:
            weights = np.array(G.es["weight"]).astype(np.float64)
        except:
            weights = None

        if weights is not None:
            if not np.isfinite(np.sum(weights)):  # weights failed
                adjmat = knn.kneighbors_graph(mode="connectivity")
                affinity_matrix = 0.5 * (adjmat + adjmat.T)
                G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
                weights = None

    else:
        adjmat = knn.kneighbors_graph(mode="connectivity")

        affinity_matrix = 0.5 * (adjmat + adjmat.T)
        G = sc._utils.get_igraph_from_adjacency(adjmat, directed=True)
        weights = None

    del hub
    gc.collect()

    return adjmat, affinity_matrix, G, weights, scores


def get_scores(true_labels, res_dict, retained_cells_idx, X, metric, seed):
    keys_list = list(res_dict.keys())

    df = pd.concat(
        (
            # supervised
            pd.Series(
                data=[
                    sklearn.metrics.cluster.adjusted_rand_score(
                        true_labels, res_dict[k][retained_cells_idx]
                    )
                    for k in keys_list
                ],
                index=keys_list,
            ),
            pd.Series(
                data=[
                    sklearn.metrics.cluster.adjusted_mutual_info_score(
                        true_labels,
                        res_dict[k][retained_cells_idx],
                        average_method="arithmetic",
                    )
                    for k in keys_list
                ],
                index=keys_list,
            ),
            pd.Series(
                data=[
                    sklearn.metrics.cluster.homogeneity_score(
                        true_labels, res_dict[k][retained_cells_idx]
                    )
                    for k in keys_list
                ],
                index=keys_list,
            ),
            pd.Series(
                data=[
                    sklearn.metrics.cluster.completeness_score(
                        true_labels, res_dict[k][retained_cells_idx]
                    )
                    for k in keys_list
                ],
                index=keys_list,
            ),
            pd.Series(
                data=[
                    sklearn.metrics.cluster.v_measure_score(
                        true_labels, res_dict[k][retained_cells_idx]
                    )
                    for k in keys_list
                ],
                index=keys_list,
            ),
            # unsupervised
            # pd.Series(data=[sklearn.metrics.cluster.davies_bouldin_score(X[retained_cells_idx],res_dict[k][retained_cells_idx]) for k in keys_list],index=keys_list),
            # pd.Series(data=[sklearn.metrics.cluster.calinski_harabasz_score(X[retained_cells_idx],res_dict[k][retained_cells_idx]) for k in keys_list],index=keys_list),
            # pd.Series(data=[sklearn.metrics.cluster.silhouette_score(X,res_dict[k][retained_cells_idx],metric=metric,random_state=seed) for k in keys_list],index=keys_list)
        ),
        axis=1,
    )

    df.columns = [
        "ARI",
        "AMI",
        "Homogeneity",
        "Completeness",
        "V_measure",
    ]  # ,'Davies-Bouldin','Calinski-Harabasz']#,'Silhouette']

    return df


def recipe_duo(adata, do_log, renorm):
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)

    # keep 5000 genes with highest average expression
    exprsn = np.array(adata.X.mean(axis=0)).reshape(-1)
    keep = np.argsort(exprsn)[::-1][:5000]
    adata._inplace_subset_var(keep)
    if renorm:
        sc.pp.normalize_total(adata, target_sum=1e4)


def recipe_seurat(adata, do_log, norm_scale):
    # same as scanpy clustering tutorial except initial cells/genes prefiltering /!\ sc.pp.recipe_seurat filters non log data ?
    sc.pp.normalize_total(adata, target_sum=1e4)
    if do_log:
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata._inplace_subset_var(adata.var.highly_variable)
    if norm_scale:
        sc.pp.scale(adata, max_value=10)


def recipe_zheng17(adata, do_log, norm_scale, n_top_genes=1000):
    # same as sc.pp.recipe_zheng17 except initial cells/genes prefiltering
    sc.pp.normalize_total(adata, key_added="n_counts_all")
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X,
        flavor="cell_ranger",
        n_top_genes=min(adata.X.shape[1], n_top_genes),
        log=False,
    )
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    sc.pp.normalize_total(adata)  # renormalize after filtering
    if do_log:
        sc.pp.log1p(adata)  # log transform: X = log(X + 1)
    if norm_scale:
        sc.pp.scale(adata)


def clustering_analysis(
    adata,
    true_labels,
    do_norm,
    norm_scale,
    do_log,
    do_pca,
    n_neighbors,
    n_clusters,
    metric,
    weighted,  # weighted adjmat for louvain/leiden clustering ?
    seed,
    n_comps=50,
    hubness_methods={
        "mp_normal": ("mp", {"method": "normal"}),
        "ls": ("ls", None),
        "ls_nicdm": ("ls", {"method": "nicdm"}),
        "dsl": ("dsl", None),
    },
    retained_cells_idx=None,
):

    results_dict = {}
    results_dict["params"] = dict(
        do_norm=do_norm,
        norm_scale=norm_scale,
        do_log=do_log,
        do_pca=do_pca,
        n_neighbors=n_neighbors,
        n_clusters=n_clusters,
        metric=metric,
        weighted=weighted,
        seed=seed,
        n_comps=n_comps,
    )

    start = time.time()

    ### preprocess, prepare clustering input ###
    if retained_cells_idx is None:
        retained_cells_idx = range(len(adata.X))

    if type(do_norm) is str:
        adata.X = scipy.sparse.csr_matrix(adata.X)

        if do_norm == "seurat":
            recipe_seurat(adata, do_log, norm_scale)
            print(f"\t\tseurat norm retained {adata.X.shape[1]} genes")
        elif do_norm == "zheng17":
            recipe_zheng17(adata, do_log, norm_scale, n_top_genes=5000)
            print(f"\t\tzheng norm retained {adata.X.shape[1]} genes")
        elif do_norm == "duo":
            recipe_duo(adata, do_log, renorm=norm_scale)
            print(f"\t\tduo norm retained {adata.X.shape[1]} genes")
        else:
            raise ValueError("do_norm not in 'duo', seurat', 'zheng17'")

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    if do_log and not (type(do_norm) is str):
        print("\t\tlog_transformed data")
        sc.pp.log1p(adata)

    if do_pca:
        use_rep = "X_pca"
        sc.pp.pca(
            adata, n_comps=min(adata.X.shape[1] - 1, min(len(adata.X) - 1, n_comps))
        )
        X = adata.obsm["X_pca"]
        res_key = results_dict["X_pca"] = {}
    else:
        # already computed pca
        use_rep = "X_pca"
        X = adata.obsm["X_pca"]
        res_key = results_dict["X_pca"] = {}

    print("\t\t\tPreprocessing done:", round((time.time() - start) / 60, 2), "mn")
    start = time.time()

    adjmat, affinity_matrix, G, weights, scores = generate_clustering_inputs2(
        X,
        metric=metric,
        n_neighbors=n_neighbors,
        weighted=weighted,
        seed=seed,
        hubness=None,
        hubness_params=None,
    )

    print("\t\t\tInputs generated:", round((time.time() - start) / 60, 2), "mn")
    start = time.time()

    res_key["hubness_df"] = pd.DataFrame.from_dict(
        data=scores, orient="index", columns=["base"]
    )
    # import pdb;pdb.set_trace()

    ### base and scanpy clustering ###
    # main dictionaries
    res_key["clus"] = {}
    res_key["clus_info"] = {}
    res_key["clus_scores"] = {}

    # sub dictionaries
    clus_methods_keys = [
        "scanpy_default_umap",
        "scanpy_default_gauss",
        "base_default",
        "scanpy_umap",
        "scanpy_gauss",
        "base",
    ]
    for k in clus_methods_keys:
        res_key["clus"][k] = {}
        res_key["clus_info"][k] = {}

    # cluster with default params
    # scanpy
    for method in ["umap", "gauss"]:
        # compute neighbors
        try:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
            )
        except:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
                knn=False,
            )

        # cluster
        sc.tl.louvain(adata, resolution=1.0, random_state=seed, use_weights=weighted)
        res_key["clus"]["scanpy_default_" + method]["louvain"] = np.array(
            adata.obs["louvain"]
        ).astype(int)
        res_key["clus_info"]["scanpy_default_" + method]["louvain"] = dict(
            n_clusters=len(np.unique(np.array(adata.obs["louvain"]).astype(int)))
        )

        sc.tl.leiden(adata, resolution=1.0, random_state=seed, use_weights=weighted)
        res_key["clus"]["scanpy_default_" + method]["leiden"] = np.array(
            adata.obs["leiden"]
        ).astype(int)
        res_key["clus_info"]["scanpy_default_" + method]["leiden"] = dict(
            n_clusters=len(np.unique(np.array(adata.obs["leiden"]).astype(int)))
        )

    print(
        "\t\t\tScanpy louvain/leiden clus:", round((time.time() - start) / 60, 2), "mn"
    )
    start = time.time()

    # base
    louvain.set_rng_seed(seed)
    part_louvain = louvain.find_partition(
        graph=G,
        partition_type=louvain.RBConfigurationVertexPartition,
        weights=weights,
        resolution_parameter=1.0,
    )
    res_key["clus"]["base_default"]["louvain"] = np.array(part_louvain.membership)
    res_key["clus_info"]["base_default"]["louvain"] = dict(
        n_clusters=len(np.unique(np.array(part_louvain.membership)))
    )

    part_leiden = leidenalg.find_partition(
        graph=G,
        partition_type=leidenalg.RBConfigurationVertexPartition,
        weights=weights,
        resolution_parameter=1.0,
        seed=seed,
    )
    res_key["clus"]["base_default"]["leiden"] = np.array(part_leiden.membership)
    res_key["clus_info"]["base_default"]["leiden"] = dict(
        n_clusters=len(np.unique(np.array(part_leiden.membership)))
    )

    print("\t\t\tLouvain/leiden clus:", round((time.time() - start) / 60, 2), "mn")
    start = time.time()

    # cluster with ground truth number of clusters
    # scanpy
    for method in ["umap", "gauss"]:
        # compute neighbors
        try:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
            )
        except:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
                knn=False,
            )
        # cluster
        (
            res_key["clus"]["scanpy_" + method]["louvain"],
            res_key["clus_info"]["scanpy_" + method]["louvain"],
        ) = getNclusters2(
            adata,
            G,
            n_clusters=n_clusters,
            seed=seed,
            cluster_func="louvain",
            flavor="scanpy",
            weights=weights,
        )
        (
            res_key["clus"]["scanpy_" + method]["leiden"],
            res_key["clus_info"]["scanpy_" + method]["leiden"],
        ) = getNclusters2(
            adata,
            G,
            n_clusters=n_clusters,
            seed=seed,
            cluster_func="leiden",
            flavor="scanpy",
            weights=weights,
        )

    print(
        "\t\t\tScanpy louvain/leiden clus, searching ground truth:",
        round((time.time() - start) / 60, 2),
        "mn",
    )
    start = time.time()

    (
        res_key["clus"]["base"]["louvain"],
        res_key["clus_info"]["base"]["louvain"],
    ) = getNclusters2(
        adata,
        G,
        n_clusters=n_clusters,
        seed=seed,
        cluster_func="louvain",
        flavor="base",
        weights=weights,
    )
    (
        res_key["clus"]["base"]["leiden"],
        res_key["clus_info"]["base"]["leiden"],
    ) = getNclusters2(
        adata,
        G,
        n_clusters=n_clusters,
        seed=seed,
        cluster_func="leiden",
        flavor="base",
        weights=weights,
    )
    print(
        "\t\t\tBase louvain/leiden clus, searching ground truth:",
        round((time.time() - start) / 60, 2),
        "mn",
    )
    start = time.time()

    # store scores, info
    for k in clus_methods_keys:
        res_key["clus_scores"][k] = get_scores(
            true_labels, res_key["clus"][k], retained_cells_idx, X, metric, seed
        )
        res_key["clus_info"][k] = pd.DataFrame(res_key["clus_info"][k])

    print("\t\t\tScoring:", round((time.time() - start) / 60, 2), "mn")
    start = time.time()

    del adjmat, affinity_matrix, G, weights, scores

    ### hubness-aware clustering ###

    for method_name, (hubness, hubness_params) in hubness_methods.items():
        res_key["clus"][method_name] = {}
        res_key["clus"][method_name + "_default"] = {}
        res_key["clus_info"][method_name] = {}
        res_key["clus_info"][method_name + "_default"] = {}

        (
            hub_adjmat,
            hub_affinity_matrix,
            hub_G,
            hub_weights,
            hub_scores,
        ) = generate_clustering_inputs(
            X,
            metric=metric,
            n_neighbors=n_neighbors,
            weighted=weighted,
            seed=seed,
            hubness=hubness,
            hubness_params=hubness_params,
        )

        # store hubness information
        res_key["hubness_df"] = pd.concat(
            (
                res_key["hubness_df"],
                pd.DataFrame.from_dict(
                    data=hub_scores, orient="index", columns=[method_name]
                ),
            ),
            axis=1,
        )

        # cluster with default params
        louvain.set_rng_seed(seed)
        part_louvain = louvain.find_partition(
            graph=hub_G,
            partition_type=louvain.RBConfigurationVertexPartition,
            weights=hub_weights,
            resolution_parameter=1.0,
        )
        res_key["clus"][method_name + "_default"]["louvain"] = np.array(
            part_louvain.membership
        )
        res_key["clus_info"][method_name + "_default"]["louvain"] = dict(
            n_clusters=len(np.unique(np.array(part_louvain.membership)))
        )

        part_leiden = leidenalg.find_partition(
            graph=hub_G,
            partition_type=leidenalg.RBConfigurationVertexPartition,
            weights=hub_weights,
            resolution_parameter=1.0,
            seed=seed,
        )
        res_key["clus"][method_name + "_default"]["leiden"] = np.array(
            part_leiden.membership
        )
        res_key["clus_info"][method_name + "_default"]["leiden"] = dict(
            n_clusters=len(np.unique(np.array(part_leiden.membership)))
        )

        # cluster with ground truth number of clusters
        (
            res_key["clus"][method_name]["louvain"],
            res_key["clus_info"][method_name]["louvain"],
        ) = getNclusters2(
            adata,
            hub_G,
            n_clusters=n_clusters,
            seed=seed,
            cluster_func="louvain",
            flavor="base",
            weights=hub_weights,
        )
        (
            res_key["clus"][method_name]["leiden"],
            res_key["clus_info"][method_name]["leiden"],
        ) = getNclusters2(
            adata,
            hub_G,
            n_clusters=n_clusters,
            seed=seed,
            cluster_func="leiden",
            flavor="base",
            weights=hub_weights,
        )

        # store clustering scores, info
        res_key["clus_scores"][method_name + "_default"] = get_scores(
            true_labels,
            res_key["clus"][method_name + "_default"],
            retained_cells_idx,
            X,
            metric,
            seed,
        )
        res_key["clus_scores"][method_name] = get_scores(
            true_labels,
            res_key["clus"][method_name],
            retained_cells_idx,
            X,
            metric,
            seed,
        )

        res_key["clus_info"][method_name + "_default"] = pd.DataFrame(
            res_key["clus_info"][method_name + "_default"]
        )
        res_key["clus_info"][method_name] = pd.DataFrame(
            res_key["clus_info"][method_name]
        )

    print(
        "\t\t\tHubness methods full pipeline:",
        round((time.time() - start) / 60, 2),
        "mn",
    )

    return results_dict


def scanpy_hubness_analysis(
    adata,
    do_norm,
    norm_scale,
    do_log,
    do_pca,
    n_neighbors,
    metric,
    weighted,  # weighted adjmat for louvain/leiden clustering ?
    seed,
    n_comps=50,
    retained_cells_idx=None,
):

    results_dict = {}
    results_dict["params"] = dict(
        do_norm=do_norm,
        norm_scale=norm_scale,
        do_log=do_log,
        do_pca=do_pca,
        n_neighbors=n_neighbors,
        metric=metric,
        weighted=weighted,
        seed=seed,
        n_comps=n_comps,
    )

    start = time.time()

    ### preprocess, prepare clustering input ###
    if retained_cells_idx is None:
        retained_cells_idx = range(len(adata.X))

    if type(do_norm) is str:
        adata.X = scipy.sparse.csr_matrix(adata.X)

        if do_norm == "seurat":
            recipe_seurat(adata, do_log, norm_scale)
            print(f"\t\tseurat norm retained {adata.X.shape[1]} genes")
        elif do_norm == "zheng17":
            recipe_zheng17(adata, do_log, norm_scale, n_top_genes=5000)
            print(f"\t\tzheng norm retained {adata.X.shape[1]} genes")
        elif do_norm == "duo":
            recipe_duo(adata, do_log, renorm=norm_scale)
            print(f"\t\tduo norm retained {adata.X.shape[1]} genes")
        else:
            raise ValueError("do_norm not in 'duo', seurat', 'zheng17'")

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    if do_log and not (type(do_norm) is str):
        print("\t\tlog_transformed data")
        sc.pp.log1p(adata)

    if do_pca:
        use_rep = "X_pca"
        sc.pp.pca(
            adata, n_comps=min(adata.X.shape[1] - 1, min(len(adata.X) - 1, n_comps))
        )
        X = adata.obsm["X_pca"]
        res_key = results_dict["X_pca"] = {}
    else:
        # already computed pca
        use_rep = "X_pca"
        X = adata.obsm["X_pca"]
        res_key = results_dict["X_pca"] = {}

    skews = {}
    # scanpy
    for method in ["umap", "gauss"]:
        # compute neighbors
        try:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
            )
        except:
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors + 1,
                metric=metric,
                use_rep=use_rep,
                method=method,
                knn=False,
            )

        skews[method] = skew(adata.obsp["connectivities"].sum(axis=0).flat)
    print("\t\t\tScoring:", round((time.time() - start) / 60, 2), "mn")

    return skews


get_key_df = lambda l, key: {k: v.loc[key] for k, v in l.items()}


def parse_results_dict(results_dict):
    # extract main info
    hubness_dfs = {}
    clus = {}
    clus_scores = {}
    clus_info = {}
    for k, v in results_dict.items():
        hubness_dfs[k] = v["X_pca"]["hubness_df"]
        clus[k] = v["X_pca"]["clus"]
        clus_scores[k] = v["X_pca"]["clus_scores"]
        clus_info[k] = v["X_pca"]["clus_info"]

    # invert dicts of dicts to regroup by clustering or hubness method : {datasets:{methods}} -> {methods:{datasets}}
    get_key = lambda l, key: {k: v[key] for k, v in l.items()}

    hub_methods_keys = hubness_dfs[list(hubness_dfs.keys())[0]].columns
    clus_methods_keys = clus[list(clus.keys())[0]].keys()

    methods_clus = {}
    methods_clus_scores = {}
    methods_clus_info = {}

    for key in clus_methods_keys:
        methods_clus[key] = get_key(clus, key)
        methods_clus_scores[key] = get_key(clus_scores, key)
        methods_clus_info[key] = get_key(clus_info, key)

    methods_hubness_dfs = {key: get_key(hubness_dfs, key) for key in hub_methods_keys}

    return methods_hubness_dfs, methods_clus, methods_clus_scores, methods_clus_info


def load_rds(path, fname):
    readRDS = ro.r["readRDS"]
    raw = readRDS(path + fname + ".rds")
    adata_dict = {k: v for k, v in raw.items()}
    adata = anndata.AnnData(X=adata_dict["counts"])
    adata.uns["Order"] = adata_dict["groups_id"].iloc[:, 1].values
    del raw, adata_dict
    return adata


def load_h5ad(path, fname):
    return anndata.read_h5ad(path + fname + ".h5ad")

