import numpy as np
import skdim
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from umap import UMAP
from sklearn.metrics import pairwise_distances
from visualisation_common_functions import QDM, QNP, natural_pca,\
    _binary_search_perplexity, _joint_probabilities, _kl_divergence,\
    fuzzy_simplicial_set, cross_entropy
import skhubness

# parameters
seed = 0
colors = ['tab:blue', 'tab:red']
n_points = 5000
n_draws = 10
perplexity = 30
original_dims = int(input('Which dims? \n'))
k = int(input('Which k? \n'))
algo_proj = input('Which projection? umap or tsne \n')


def proj_raw(raw, reduction):
    data = pairwise_distances(raw)
    if reduction not in ['tsne', 'umap']:
        raise ValueError("reductio has to be one of ['tsne', 'umap']")
    if reduction == 'tsne':
        tsne = TSNE(n_components=2, random_state=seed, perplexity=perplexity)
        embedding = tsne.fit_transform(data)
    if reduction == 'umap':
        umap = UMAP(n_components=2, random_state=seed)
        embedding = umap.fit_transform(data)
    return embedding


# generate the distributions
cube = np.random.uniform(-0.0001, 0.0001, size=(n_draws, n_points, original_dims))
sphere = np.zeros((n_draws, n_points, original_dims))
for i in range(n_draws):
    sphere[i, :, :] = skdim.datasets.hyperSphere(n=n_points, d=original_dims)
# print('distributions generated')

# Evaluate skewness
cube_skew = [skhubness.Hubness(k=k, metric='euclidean').fit(cube[i, :, :]).score() for i in range(n_draws)]
sphere_skew = [skhubness.Hubness(k=k, metric='euclidean').fit(sphere[i, :, :]).score() for i in range(n_draws)]
skew = np.stack((cube_skew, sphere_skew), axis=1)
plt.clf()
box_skew = plt.boxplot(skew, patch_artist=True, labels=["n-Cube", "n-Sphere"])
for patch, color in zip(box_skew['boxes'], colors):
    patch.set_facecolor(color)
plt.ylabel("k-Skewness")
plt.show()
# print('skewness calculated')

# Make 2D projection
cube_proj = [proj_raw(cube[i, :, :], algo_proj) for i in range(n_draws)]
sphere_proj = [proj_raw(sphere[i, :, :], algo_proj) for i in range(n_draws)]
print('projections done')

# Evaluate QDM and QNP
# cube_qdm = Parallel(n_jobs=4, verbose=10)(QDM(cube[k, :, :], cube_proj[k]) for k in range(n_draws))
# sphere_qdm = Parallel(n_jobs=4, verbose=10)(QDM(sphere[k, :, :], sphere_proj[k]) for k in range(n_draws))
# cube_qnp = Parallel(n_jobs=4, verbose=10)(QNP(cube[k, :, :], cube_proj[k]) for k in range(n_draws))
# sphere_qnp = Parallel(n_jobs=4, verbose=10)(QNP(sphere[k, :, :], sphere_proj[k]) for k in range(n_draws))
cube_qdm = [QDM(cube[i, :, :], cube_proj[i], nPC=1000) for i in range(n_draws)]
sphere_qdm = [QDM(sphere[i, :, :], sphere_proj[i], nPC=1000) for i in range(n_draws)]
cube_qnp = [QNP(cube[i, :, :], cube_proj[i], k=k) for i in range(n_draws)]
sphere_qnp = [QNP(sphere[i, :, :], sphere_proj[i], k=k) for i in range(n_draws)]
print('metrics calculation done')
corr_tot = np.stack((cube_qdm, sphere_qdm, cube_qnp, sphere_qnp), axis=1)

# Evaluate cost functions
cube_kl = [_kl_divergence(_joint_probabilities(distances=pairwise_distances(cube[i, :, :], squared=True),
                                               desired_perplexity=perplexity),
                          cube_proj[i]) for i in range(n_draws)]
sphere_kl = [_kl_divergence(_joint_probabilities(distances=pairwise_distances(sphere[i, :, :], squared=True),
                                                 desired_perplexity=perplexity),
                            sphere_proj[i]) for i in range(n_draws)]
cube_ce = [cross_entropy(cube[i, :, :],
                         cube_proj[i]) for i in range(n_draws)]
sphere_ce = [cross_entropy(sphere[i, :, :],
                           sphere_proj[i]) for i in range(n_draws)]
print('cost function calculation done')
cost_tot = np.stack((cube_kl, sphere_kl, cube_ce, sphere_ce), axis=1)
