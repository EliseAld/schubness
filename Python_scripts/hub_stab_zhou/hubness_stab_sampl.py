from tqdm import tqdm
from skhubness import Hubness
import matplotlib.pyplot as plt
import rpy2.robjects as ro
import anndata2ri
import anndata
import scipy
import numpy as np
import scanpy as sc
import time
import warnings

warnings.filterwarnings("ignore")
anndata2ri.activate()
readRDS = ro.r['readRDS']
path_rds = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/rds/"
path_res = "/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_zhou/"
get_res_path = lambda path_res: f'{path_res}hub_stab_sampling_norm{do_norm}_scale{norm_scale}_{metric}_{clustering_algo}.png'

#fixed params
seed = 0
n_iter = 10
do_log = True
weighted = True
n_neighbors = 10
norm_scale = True
do_norm = 'duo'
bootstrap_size = 0.9
metric = 'euclidean'
clustering_algo = 'leiden'
fnames = ['gold_hematopoiesis-gates_olsson', 'gold_germline-human-female-weeks_li', 'gold_stimulated-dendritic-cells-LPS_shalek',
          'gold_germline-human-female_guo', 'gold_mESC-differentiation_hayashi', 'gold_developing-dendritic-cells_schlitzer',
          'gold_germline-human-male-weeks_li', 'gold_pancreatic-beta-cell-maturation_zhang', 'gold_human-embryos_petropoulos',
          'gold_germline-human-male_guo', 'gold_myoblast-differentiation_trapnell', 'gold_pancreatic-alpha-cell-maturation_zhang',
          'gold_aging-hsc-young_kowalczyk', 'gold_aging-hsc-old_kowalczyk']
colors = ['red', 'gold', 'chartreuse', 'deepskyblue']

#vary euclidean
n_comps = (25, 50, 100, 500)


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

def load_data(path):
    raw = readRDS(path)
    adata_dict = {k:v for k, v in raw.items()}
    adata = anndata.AnnData(X=adata_dict['counts'])
    adata.uns['Order'] = adata_dict['groups_id'].iloc[:, 1].values
    del raw, adata_dict
    return adata

def resampling(adata, bootstrap_size=bootstrap_size,n_iter=n_iter):
    adata_sampled = dict()
    cell_iter = np.zeros((n_iter,adata.n_obs))
    for iter in tqdm(range(n_iter)):
        cell_bootstrap = np.random.uniform(0, 1, size=adata.n_obs)
        cell_bootstrap[cell_bootstrap <= bootstrap_size] = 2
        cell_bootstrap = cell_bootstrap == 2
        cell_iter[iter, :] = cell_bootstrap
        obsm = {'X_pca': adata.obsm['X_pca'][cell_bootstrap]}
        adata_sampled[iter] = anndata.AnnData(adata.X[cell_bootstrap],
                                        obsm=obsm)
        #sc.pp.neighbors(adata_sampled[iter], random_state=seed)
    return adata_sampled, cell_iter

def hub_retrieval(X, metric, k=n_neighbors):
    n_obs = X.shape[0]
    hub = Hubness(k=k, metric=metric, random_state=seed, return_value='hubs')
    hub.fit(X)
    hub_id = np.repeat('norm', n_obs)
    hub_id[hub.score()] = 'hubs'
    return hub_id

def overlap_hubs(hubs_ref, hubs_sampled, cell_iter):
    hubs_ref_match = hubs_ref[cell_iter == 1]
    size_overlap = [100*np.sum([(hubs_ref_match == 'hubs')[loop] and
                                       (hubs_sampled == 'hubs')[loop] for loop in range(len(hubs_ref_match))])/
                            np.sum(hubs_ref == 'hubs')][0]
    return size_overlap


# run functions
# get the number of hubs per dataset
size_overlap_tot = np.zeros(shape=(len(fnames), len(n_comps)))
for dim in n_comps:
    size_overlap = np.zeros(len(fnames))
    print(metric, dim, clustering_algo)
    #if np.array([get_res_path(path_res).split('/')[-1] in elt for elt in os.listdir(path_res)]).any(): #setting already computed
    #        continue
    for fname in fnames:
        adata = load_data(path_rds+fname+'.rds')
        if dim < adata.shape[0]:
            start=time.time()
            ### preprocess###
            adata.X = scipy.sparse.csr_matrix(adata.X)
            recipe_duo(adata, do_log, renorm=norm_scale)
            if scipy.sparse.issparse(adata.X):
                adata.X = adata.X.toarray()
            sc.tl.pca(adata, n_comps=min(adata.X.shape[1]-1, min(len(adata.X)-1, dim)))
            X = adata.obsm['X_pca']
            print('\t\t\tPreprocessing done:', round((time.time()-start)/60, 2), 'mn')
            start = time.time()
            ### get hubs ref ###
            if dim < adata.shape[1]:
                # Make the lists of ref hubs
                hubs_ref = hub_retrieval(X, metric=metric)
                print('\t\t\tRef hubs done:', round((time.time()-start)/60, 2), 'mn')
                start = time.time()
                ### get sampled data###
                # Make the lists of sampled hubs
                adata_sampled, cell_iter = resampling(adata)
                X_sampled = dict()
                for key in adata_sampled.keys():
                    adata_sampled[key].X = scipy.sparse.csr_matrix(adata_sampled[key].X)
                    recipe_duo(adata_sampled[key], do_log, renorm=norm_scale)
                    if scipy.sparse.issparse(adata.X):
                        adata_sampled[key].X = adata_sampled[key].X.toarray()
                    sc.tl.pca(adata_sampled[key], n_comps=min(adata_sampled[key].X.shape[1]-1, min(adata_sampled[key].X.shape[0]-1, dim)))
                    X_sampled[key] = adata_sampled[key].obsm['X_pca']
                hubs_s = dict()
                for iter in X_sampled.keys():
                    hubs_s[iter] = hub_retrieval(X_sampled[iter], metric=metric)
                print('\t\t\tSampled hubs done:', round((time.time()-start)/60, 2), 'mn')
                start = time.time()
                ### get overlap###
                size_overlap[np.argwhere(np.array(fnames) == fname)] = np.mean([overlap_hubs(hubs_ref, hubs_s[i], cell_iter[i]) for i in range(n_iter)])
                print('\t\t\tOverlap calculated:', round((time.time()-start)/60, 2), 'mn')
        else:
            size_overlap[np.argwhere(np.array(fnames) == fname)] = None
    size_overlap_tot[:, np.argwhere(np.array(n_comps) == dim)[0][0]] = size_overlap

boxplot = plt.boxplot(size_overlap_tot, labels=['25', '50', '100', '500'], patch_artist=True)
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)
for dim in range(len(n_comps)):
    plt.scatter(np.repeat(dim+1, size_overlap_tot.shape[0]), size_overlap_tot[:, dim], c=colors[::-1][dim])
plt.savefig(get_res_path(path_res))
#plt.show()
