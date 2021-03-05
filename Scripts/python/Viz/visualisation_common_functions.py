from sklearn import metrics
import sklearn
import math
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from sklearn.neighbors import NearestNeighbors
import umap
import numpy as np
import scipy

# fixed params
EPSILON_DBL = 1e-8
PERPLEXITY_TOLERANCE = 1e-5
seed = 0


def _binary_search_perplexity(sqdistances, desired_perplexity):
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
    Returns
    -------
    P : array, shape (n_samples, n_samples)
        Probabilities of conditional Gaussian distributions p_i|j.
    """
    n_steps = 100
    n_samples = sqdistances.shape[0]
    n_neighbors = sqdistances.shape[1]
    using_neighbors = n_neighbors < n_samples
    desired_entropy = np.log(desired_perplexity)
    P = np.zeros((n_samples, n_neighbors), dtype=np.float64)
    for i in range(n_samples):
        beta_min = -math.inf
        beta_max = math.inf
        beta = 1.0
        for l in range(n_steps):
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
    return P


def _joint_probabilities(distances, desired_perplexity):
    """Compute joint probabilities p_ij from distances.
    Parameters
    ----------
    distances : ndarray of shape (n_samples * n_neighbors,)
    desired_perplexity : float
        Desired perplexity of the joint probability distributions.
    Returns
    -------
    P : ndarray of shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    """
    conditional_P = _binary_search_perplexity(distances, desired_perplexity)
    P = conditional_P + conditional_P.T
    sum_P = np.maximum(np.sum(P), np.finfo(np.double).eps)
    P = np.maximum(squareform(P) / sum_P, np.finfo(np.double).eps)
    return P


def _kl_divergence(X_embedded, P, degrees_of_freedom=1):
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
    dist = pdist(X_embedded, "sqeuclidean")
    dist /= degrees_of_freedom
    dist += 1.
    dist **= (degrees_of_freedom + 1.0) / -2.0
    Q = np.maximum(dist / (2.0 * np.sum(dist)), np.finfo(np.double).eps)
    kl_divergence = 2.0 * np.dot(
            P, np.log(np.maximum(P, np.finfo(np.double).eps) / Q))
    return kl_divergence


def fuzzy_simplicial_set(X, n_neighbors=None, metric='euclidean', random_state=seed):
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
    knn_indices, knn_dists, _ = umap.umap_.nearest_neighbors(X, n_neighbors, metric, random_state=sklearn.utils.check_random_state(random_state), verbose=False, angular=False, metric_kwds={})
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
    n_neighbors = int(np.sqrt(X.shape[0]))
    dist_X = metrics.pairwise_distances(X, metric=metric)
    np.fill_diagonal(dist_X, np.inf)
    ind_X = np.argsort(dist_X, axis=1)
    ind_X_embedded = NearestNeighbors(n_neighbors=n_neighbors).fit(
            X_embedded).kneighbors(return_distance=False)
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
