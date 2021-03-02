import numpy as np
import skdim
import matplotlib.pyplot as plt
from tqdm import tqdm
import skhubness
from sklearn.manifold import TSNE
from umap import UMAP
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
import networkx as nx
#from mpl_toolkits import mplot3d
#from mpl_toolkits.mplot3d import Axes3D

# parameters
seed = 0
EPSILON_DBL = 1e-8
PERPLEXITY_TOLERANCE = 1e-5
colors = ['tab:blue', 'tab:red']
colors2 = ['tab:red', 'gold', 'limegreen', 'deepskyblue']
n_points = 1000
n_draws = 5
n_satellites = [2, 5, 10, 20]
original_dims = int(input('Which dims? \n'))
k = int(input('Which k? \n'))
algo_proj = input('Which projection? \n')

def proj_raw(raw, reduction):
    data = pairwise_distances(raw)
    if (reduction not in ['tsne', 'umap']):
        raise ValueError("reductio has to be one of ['tsne', 'umap']")
    if reduction=='tsne':
        tsne = TSNE(n_components=2, random_state=seed, perplexity=30.0)
        embedding = tsne.fit_transform(data)
    if reduction=='umap':
        umap = UMAP(n_components=2, random_state=seed)
        embedding = umap.fit_transform(data)
    return embedding

def natural_pca(data, metric='euclidean', nPC=1000):
    pairwise_dist = pairwise_distances(data, metric=metric)
    nPC = min(nPC, pairwise_dist.shape[0])
    natural_PC = np.zeros(shape=(2, nPC), dtype=int)
    for i in tqdm(range(nPC)):
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

def QDM(input, embedding, metric='euclidean',nPC=1000):
    useful_points = natural_pca(input, metric, nPC)
    print('nPCA done')
    dis_i = pairwise_distances(input, metric=metric)
    dis_o = pairwise_distances(embedding, metric=metric)
    ravel_i = []
    ravel_o = []
    for i in range(useful_points.shape[1]):
        ravel_i.append(dis_i[useful_points[0, i], useful_points[1, i]])
        ravel_o.append(dis_o[useful_points[0, i], useful_points[1, i]])
    corr = np.corrcoef(ravel_i, ravel_o)
    return corr[0, 1]

def QNP(input, embedding, metric='euclidean', k=k):
    NN = NearestNeighbors(n_neighbors=k, metric=metric)
    neigh_i = NN.fit(input).kneighbors(return_distance=False)
    neigh_o = NN.fit(embedding).kneighbors(return_distance=False)
    qnp = np.sum([len(np.intersect1d(neigh_i[i, :], neigh_o[i, :])) for i in range(neigh_i.shape[0])])/k/input.shape[0]
    return qnp

def plot_kamada_graph(input, ref, k=k):
    adjmat = skhubness.Hubness(k=k, metric='euclidean',
                                    hubness=None, hubness_params=None,
                                    random_state=seed, store_k_occurrence=True, return_value='all').fit(
        input).nn_index_.kneighbors_graph(mode='distance')
    adjmat_ref = skhubness.Hubness(k=k, metric='euclidean',
                                    hubness=None, hubness_params=None,
                                    random_state=seed, store_k_occurrence=True, return_value='all').fit(
        ref).nn_index_.kneighbors_graph(mode='distance')
    G = nx.from_numpy_matrix(adjmat.todense(), create_using=nx.Graph)
    DG = nx.from_numpy_matrix(adjmat_ref.todense(), create_using=nx.DiGraph)
    pos = nx.kamada_kawai_layout(G, scale=2)
    norm = [np.linalg.norm(ref[i, :]) for i in range(n_points)]
    nx.draw_networkx(G, pos=pos,
                     arrowsize=5,
                     nodelist=list(np.arange(adjmat.shape[0])),
                     node_size=[0.1 + DG.in_degree()[i] * 5 for i in np.arange(adjmat.shape[0])],
                     width=.05,
                     with_labels=False,
                     node_color=norm)


# generate the modified n-sphere
sphere_init = skdim.datasets.hyperSphere(n=n_points, d=2)
sphere_init = np.concatenate((sphere_init, np.zeros((n_points, original_dims-2))), axis=1)
dist_carac = np.min(pairwise_distances(sphere_init)[pairwise_distances(sphere_init)!=0])
sphere = dict()
for j in n_satellites:
    sphere[j] = np.zeros((n_draws, n_points*(j+1), original_dims))
    for i in range(n_draws):
        gaussian_noise = np.zeros((n_points*j, original_dims))
        for k in range(n_points):
            gaussian_noise[k * j:j * (k + 1), :] = np.random.normal(sphere_init[k, :], scale=1e-2*dist_carac, size=(j, original_dims))
            gaussian_noise[k * j:j * (k + 1), :2] = sphere_init[k, :2]
        sphere[j][i, :, :] = np.concatenate((sphere_init, gaussian_noise))
print('n-spheres generated')
plt.clf()
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])
for i, j in enumerate(n_satellites):
    ax.scatter3D(sphere[j][0, :, 0], sphere[j][0, :, 1], sphere[j][0, :, 2], c=colors2[i])
plt.show()

# Evaluate skewness
sphere_skew = np.array([[skhubness.Hubness(k=k, metric='euclidean').fit(sphere[i][j, :, :]).score() for i in n_satellites] for j in range(n_draws)])
plt.clf()
box_skew = plt.boxplot(sphere_skew, patch_artist=True, labels=["2-satellites", "5-satellites", "10-satellites", "20-satellites"])
for patch, color in zip(box_skew['boxes'], colors2):
    patch.set_facecolor(color)
plt.ylabel("k-Skewness")
plt.show()
plt.savefig('/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/lowD_embedding/tibo1_skewness_k'+str(k)+'_dim'+str(original_dims)+'_'+algo_proj+'.png', bbox_inches='tight', pad_inches=0)
print('skewness calculated')

# Make 2D projection
cube_dim_proj = dict()
for i in cube_dim.keys():
    cube_dim_proj[i] = [proj_raw(cube_dim[i][j, :, :], algo_proj) for j in range(n_draws)]
print('projections done')

# Evaluate QNP & QDM
cube_dim_qdm = np.array([[QDM(cube_dim[i][j, :, :], cube_dim_proj[i][j]) for i in cube_dim.keys()] for j in range(n_draws)])
cube_dim_qnp = np.array([[QNP(cube_dim[i][j, :, :], cube_dim_proj[i][j]) for i in cube_dim.keys()] for j in range(n_draws)])
plt.clf()
box_qdm = plt.boxplot(cube_dim_qdm, patch_artist=True, labels=[10, 50, 100, 500])
for patch, color in zip(box_qdm['boxes'], colors2):
    patch.set_facecolor(color)
plt.ylabel("QDM")
plt.savefig('/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/lowD_embedding/tibo1_QMP_k'+str(k)+'_dim'+str(original_dims)+'_'+algo_proj+'.png', bbox_inches='tight', pad_inches=0)
plt.clf()
box_qnp = plt.boxplot(cube_dim_qnp, patch_artist=True, labels=[10, 50, 100, 500])
for patch, color in zip(box_qnp['boxes'], colors2):
    patch.set_facecolor(color)
plt.ylabel("QNP")
plt.savefig('/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/lowD_embedding/tibo1_QNP_k'+str(k)+'_dim'+str(original_dims)+'_'+algo_proj+'.png', bbox_inches='tight', pad_inches=0)

# Correlation between skewness and QNP&QDM
plt.clf()
plt.plot(np.mean(cube_dim_skew, axis=0), 1/np.mean(cube_dim_qdm, axis=0), c='red')
plt.plot(np.mean(cube_dim_skew, axis=0), 1/np.mean(cube_dim_qnp, axis=0), c='blue')
plt.show()
