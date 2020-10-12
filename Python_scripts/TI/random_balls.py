import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from skhubness.reduction import MutualProximity, LocalScaling, DisSimLocal
from skhubness.neighbors import NearestNeighbors


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
    pca = PCA(n_components=n_components)
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

def plot_points(arr, idx0=0, idx1=1, title=None):
    #print('in plot', arr.shape)
    for i in range(n_balls):
        plt.scatter(arr[i, :, idx0], arr[i, :, idx1])
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
    neigh_dist, neigh_ind = NearestNeighbors(n_neighbors=n_points*n_balls-1).fit(input_points).kcandidates()
    samples_exp = dict()
    for method_name, hubness_param in methods.items():
        if method_name == "nothing":
            samples_exp[method_name] = pairwise_distances(input_points)
        if method_name == "mp_normal":
            tmp1, tmp2 = MutualProximity(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[indx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "ls":
            tmp1, tmp2 = LocalScaling(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[indx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "ls_nicdm":
            tmp1, tmp2 = LocalScaling(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[indx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
        if method_name == "dsl":
            tmp1, tmp2 = DisSimLocal(method=hubness_param).fit_transform(neigh_dist, neigh_ind, X=None)
            samples_exp[method_name] = np.empty((input_points.shape[0],input_points.shape[0]))
            for idx_rw in range(tmp1.shape[0]):
                for idx_cl in range(len(tmp2[idx_rw, :])):
                    samples_exp[method_name][idx_rw, tmp2[indx_rw, idx_cl]] = tmp1[idx_rw, idx_cl]
    return samples_exp


# generate balls in 100D
n_balls = 3
original_dims = 100
n_points = 1000
random_direction = np.random.randint(0, 10, size=(original_dims,)).astype(np.float32)
random_centers = np.array([random_direction * i for i in range(1, n_balls + 1)])
repeated_centers = np.array([random_centers for _ in range(n_points)])
#print('repeated_centers.shape', repeated_centers.shape)
repeated_centers = np.transpose(repeated_centers, (1, 0, 2))
#print('repeated_centers.shape', repeated_centers.shape)
samples = repeated_centers + np.random.normal(0, 1, size=repeated_centers.shape)
plot_points(samples, 1, 2)
pca_samples = get_pca(samples)
plot_points(pca_samples)

# add dropout
dropout_percents=[0, 50, 80, 90]
samples_dropout = add_dropout(samples, dropout_percents)
pca_dropout = np.array([pca_samples for _ in range(len(dropout_percents))])
for idx in range(len(dropout_percents)):
    print(dropout_percents[idx])
    #print(np.linalg.norm(samples_exp - samples))
    pca_tmp = get_pca(samples_dropout[idx, :, :, :])
    plot_points(pca_tmp)
    pca_dropout[idx, :, :, :] = pca_tmp

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
        plot_points(deformat_array(mds_samples_hubred[key1][key2]), title='dropout '+str(key1)+' and hub reduction '+key2)


