import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from skhubness.neighbors import kneighbors_graph
import anndata
import scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.manifold import MDS
from skhubness import Hubness

# parameters
seed = 0
n_balls = 3
original_dims = 100
n_points = 500
epsilon = 0.5
n_neighbors = 50
print("dim="+str(original_dims)+"\nk="+str(n_neighbors))

def reformat_array(arr):
    n_balls, n_points, original_dims = arr.shape
    input_points = []
    for n_ball in range(n_balls):
        for n_point in range(n_points):
            input_points.append(arr[n_ball, n_point, :])
    input_points = np.array(input_points)
    return input_points

def add_dropout(adata, dropout_values):
    arr = adata.X
    cluster_id = {'Cluster':adata.uns['Cluster']}
    samples_dropout = dict()
    for dropout in tqdm(dropout_values):
        tmp = np.copy(arr)
        for sample in range(adata.n_obs):
            tmp[sample, np.random.choice(range(adata.n_vars), int(dropout/100*adata.n_vars))] = 0
        samples_dropout[dropout] = anndata.AnnData(tmp, uns=cluster_id)
    return samples_dropout

def hub_reduction(adata, methods = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'}),
                   'ls':('ls',None),
                   'ls_nicdm':('ls',{'method': 'nicdm'}),
                   'dsl':('dsl',None)
                   }):
    input_points = adata.X
    samples_exp = dict()
    for method_name, (hubness, hubness_params) in tqdm(methods.items()):
        samples_exp[method_name] = adata.copy()
        samples_exp[method_name].obsm[method_name] = kneighbors_graph(input_points,
                                                            n_neighbors=n_neighbors,
                                                            hubness=hubness,
                                                            hubness_params=hubness_params)
        sc.pp.neighbors(samples_exp[method_name], n_neighbors=n_neighbors, use_rep=method_name)
        sc.tl.umap(samples_exp[method_name], n_components=2, random_state=seed)
    return samples_exp

def plot_adata_id(adata, group_id, idx0=0, idx1=1, title=None, pca=True, umap=False, mds=False):
    if pca:
        arr = adata.obsm['X_pca']
    else:
        if umap:
            arr = adata.obsm['X_umap']
        else:
            if mds:
                arr = adata.obsm['X_mds']
            else:
                arr = adata.X
    colors = ['dodgerblue', 'mediumblue', 'deeppink', 'orangered']
    for group in np.unique(group_id):
        mask = group_id==group
        plt.scatter(arr[mask, idx0], arr[mask, idx1], c=colors[int(group)])
    plt.title(title)
    #plt.show()


# generate balls in anndata
random_b = np.random.randint(5, 10, size=(original_dims,)).astype(np.float32)
random_vect = np.random.randint(0, 10, size=(original_dims,)).astype(np.float32)
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
cluster_id = {'Cluster':np.repeat(range(n_balls), n_points)}
adata = anndata.AnnData(X=reformat_array(samples), uns=cluster_id)
sc.tl.pca(adata, svd_solver='arpack', n_comps=2, random_state=seed)
sc.pp.neighbors(adata, random_state=seed)
sc.tl.umap(adata, n_components=2, random_state=seed)

plot_adata_id(adata, adata.uns['Cluster'], pca=False, title='random axes, gaussian balls')
plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step0_randomaxes.png")
#plt.show()
plt.clf()
plot_adata_id(adata, adata.uns['Cluster'], pca=True, title='pca of the gaussian balls')
plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step0_pca.png")
#plt.show()
plt.clf()
plot_adata_id(adata, adata.uns['Cluster'], pca=False, umap=True, title='umap of the gaussian balls')
plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step0_umap.png")
#plt.show()
plt.clf()

# add dropout
dropout_percents=[0, 10, 50]
samples_dropout = add_dropout(adata, dropout_percents)
for dropout in tqdm(samples_dropout.keys()):
    sc.tl.pca(samples_dropout[dropout], svd_solver='arpack', n_comps=2, random_state=seed)
    sc.pp.neighbors(samples_dropout[dropout], random_state=seed)
    sc.tl.umap(samples_dropout[dropout], n_components=2, random_state=seed)
for dropout in dropout_percents:
    print(dropout)
    plot_adata_id(samples_dropout[dropout], samples_dropout[dropout].uns['Cluster'],
                  pca=False, title='random axes, gaussian balls, dropout '+str(dropout))
    #plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step1_randomaxes_dropout"+str(dropout)+".png")
    #plt.show()
    plt.clf()
    plot_adata_id(samples_dropout[dropout], samples_dropout[dropout].uns['Cluster'],
                  title='pca of the gaussian balls, dropout '+str(dropout))
    plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step1_pca_dropout"+str(dropout)+".png")
    #plt.show()
    plt.clf()
    plot_adata_id(samples_dropout[dropout], samples_dropout[dropout].uns['Cluster'],
                  title='umap of the gaussian balls, dropout '+str(dropout),
                  pca=False, umap=True)
    plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step1_umap_dropout"+str(dropout)+".png")
    #plt.show()
    plt.clf()

# correct graph
#skew_before = []
#for key in tqdm(samples_dropout.keys()):
#    skew_before.append(Hubness(k=10, metric='cosine').fit(samples_dropout[key]).score())
#print("Skewness before correction="+str(skew_before))
samples_hubred = dict()
for dropout in tqdm(dropout_percents):
    samples_hubred[dropout] = hub_reduction(samples_dropout[dropout])
#skew_after = []
#for key in samples_dropout.keys():
#    for method in ["mp","ls","dsl"]:
#    skew_before.append(Hubness(k=10, metric='euclidean'.fit(samples_dropout[key]).score()))
#print("Skewness before correction="+str(skew_before))

# get silhouette score
for key1 in samples_hubred.keys():
    for key2 in samples_hubred[key1].keys():
        samples_hubred[key1][key2].uns['silhouette_score_umap']=silhouette_score(samples_hubred[key1][key2].obsm["X_umap"], samples_hubred[key1][key2].uns['Cluster'],
                               random_state=seed)

# apply UMAP
for key1 in samples_hubred.keys():
    for key2 in samples_hubred[key1].keys():
        plot_adata_id(samples_hubred[key1][key2], samples_hubred[key1][key2].uns['Cluster'],
                      title='Random axes, dropout '+str(key1)+' and hub reduction '+key2+'\nsilhouette score='+str(samples_hubred[key1][key2].uns['silhouette_score_umap']),
                      pca=False)
        #plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step2_randomaxes_dropout"+str(key1)+'_hub'+key2+".png")
        #plt.show()
        plt.clf()
        plot_adata_id(samples_hubred[key1][key2], samples_hubred[key1][key2].uns['Cluster'],
                      title='UMAP, dropout '+str(key1)+' and hub reduction '+key2+'\nsilhouette score='+str(samples_hubred[key1][key2].uns['silhouette_score_umap']),
                      pca=False, umap=True)
        plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step2_umap_dropout"+str(key1)+'_hub'+key2+".png")
        #plt.show()
        plt.clf()

# apply MDS
#mds = MDS(n_components=2, dissimilarity='precomputed')
#for key1 in samples_hubred:
#    for key2 in tqdm(samples_hubred[key1].keys()):
#        sample_dropout_hubred = samples_hubred[key1][key2].obsp['distances']
#        mat_upp = np.triu(sample_dropout_hubred.toarray())
#        mat_low = np.tril(sample_dropout_hubred.toarray())
#        mat_symm = (mat_upp + mat_low.T)/2 + (mat_upp.T + mat_low)/2
#        samples_hubred[key1][key2].obsm['X_mds'] = mds.fit_transform(mat_symm)
#for key1 in samples_hubred.keys():
#    for key2 in samples_hubred[key1].keys():
#        samples_hubred[key1][key2].uns['silhouette_score_mds']=silhouette_score(samples_hubred[key1][key2].obsm['X_mds'], samples_hubred[key1][key2].uns['Cluster'],
#                               random_state=seed)
#for key1 in samples_hubred.keys():
#    for key2 in samples_hubred[key1].keys():
#        plot_adata_id(samples_hubred[key1][key2], samples_hubred[key1][key2].uns['Cluster'],
#                      title='MDS, dropout '+str(key1)+' and hub reduction '+key2+'\nsilhouette score='+str(samples_hubred[key1][key2].uns['silhouette_score_mds']),
#                      pca=False)
#        plt.savefig("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/gaussian_balls/figs_random_dim"+str(original_dims)+"_k"+str(n_neighbors)+"/RandomGaussian_step2_mds_dropout"+str(key1)+'_hub'+key2+".png")
#        #plt.show()
#        plt.clf()
