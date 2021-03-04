import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from skhubness.neighbors import kneighbors_graph
import random
from scipy import stats, special
import math
from skhubness import neighbors
from skhubness import Hubness
import matplotlib.patches as patches
import networkx as nx

# fixed parameters
seed = 0
colors = ['dodgerblue', 'deeppink', 'mediumblue', 'gold', "black"]
labels = ['Nothing', 'MP', 'LS', 'LS_nicdm', 'DSL']
n_balls = 1
original_dims = 10
n_points = 10000
density_ratio = 50
n_points_density_target = 60
k = 10
cm = plt.get_cmap('plasma')
n_gaussian = 5
verbose = True
#print('\n\nk='+str(k))

def hub_reduction(input_points,
                  methods={'nothing':(None, None),
                           'mp_normal':('mp', {'method': 'normal'}),
                           'ls': ('ls', None),
                           'ls_nicdm': ('ls', {'method': 'nicdm'}),
                           'dsl': ('dsl', None)},
                  k=k):
    samples_reduced = dict()
    for method_name, (hubness, hubness_params) in tqdm(methods.items()):
        samples_reduced[method_name] = kneighbors_graph(input_points,
                                                        n_neighbors=k,
                                                        hubness=hubness,
                                                        hubness_params=hubness_params)
    return samples_reduced

def hub_eval(input_points,
             methods={'nothing': (None, None),
                      'mp_normal': ('mp', {'method': 'normal'}),
                      'ls': ('ls', None),
                      'ls_nicdm': ('ls', {'method': 'nicdm'}),
                      'dsl': ('dsl', None)},
             k=k):
    skewness = []
    for method_name, (hubness, hubness_params) in tqdm(methods.items()):
        hub = Hubness(k=k, hubness=hubness, hubness_params=hubness_params)
        hub.fit(input_points)
        skewness.append(hub.score())
    return skewness

def plot_id(arr, group_id, ax0=0, ax1=1, title=None, colors = colors):
    for group in np.unique(group_id):
        mask = group_id == group
        plt.scatter(arr[mask, ax0], arr[mask, ax1], c=colors[int(group)])
    plt.title(title)
    #plt.show()

def autolabel(rects):
    for rect in rects:
        height = int(rect.get_height())
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x()+rect.get_width()/2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def Dijkstra(graph, source, target):
    unvisited_nodes = graph.copy()
    shortest_distance = {}
    route = []
    predecessor = {}
    for nodes in unvisited_nodes:
        shortest_distance[nodes] = math.inf
    shortest_distance[source] = 0
    while (unvisited_nodes):
        min_Node = None
        for current_node in unvisited_nodes:
            if min_Node is None:
                min_Node = current_node
            elif shortest_distance[min_Node] > shortest_distance[current_node]:
                min_Node = current_node
        for child_node, value in unvisited_nodes[min_Node].items():
            if value + shortest_distance[min_Node] < shortest_distance[child_node]:
                shortest_distance[child_node] = value + shortest_distance[min_Node]
                predecessor[child_node] = min_Node
        unvisited_nodes.pop(min_Node)
    node = target
    while node != source:
        try:
            route.insert(0, node)
            node = predecessor[node]
        except Exception:
            print('Path not reachable')
            return [source, target]
            break
    route.insert(0, source)
    if shortest_distance[target] != math.inf:
        return route

def intersection_volume(n_dim):
    return np.pi**(original_dims/2)/special.gamma(original_dims/2+1) * special.betainc((n_dim+1)/2, 1/2, np.sin(np.pi/3)**2)


# generate ball
samples = np.concatenate((np.random.normal(0, 0.0001, size=(n_points, n_gaussian)),
                          np.random.uniform(-0.0001, 0.0001, size=(n_points, original_dims-n_gaussian))),
                         axis=1)
points_right_hyperplan = samples[:, 0] > 0
idx_right_hyperplan = [i for i, e in enumerate(points_right_hyperplan) if e==True]
idx_left_hyperplan = [i for i, e in enumerate(points_right_hyperplan) if e==False]
idx_right_hyperplan = random.choices(idx_right_hyperplan,
                                     k=int(n_points/2/density_ratio))
samples_skimmed = samples[idx_right_hyperplan+idx_left_hyperplan, :]
samples_skimmed = np.concatenate((samples_skimmed, np.repeat(0, original_dims).reshape(1, original_dims)), axis=0)
points_right_hyperplan = samples_skimmed[:, 0] > 0
idx_right_hyperplan = [i for i, e in enumerate(points_right_hyperplan) if e==True]
idx_left_hyperplan = [i for i, e in enumerate(points_right_hyperplan) if e==False]
n_obs = samples_skimmed.shape[0]
plt.clf()
plot_id(samples_skimmed, points_right_hyperplan, title='random axes, dense ball')
plt.show()
if verbose:
    print("sample generated")

# correct graph
samples_hubred = hub_reduction(samples_skimmed)
if verbose:
    print("sample hub-corrected")

# eval skewness
skewness = hub_eval(samples_skimmed)
if verbose:
    print("skenwness done")

# Plot the kNN neighbors
adj_matrix = {k:v for k, v in zip(samples_hubred.keys(), [samples_hubred[method].toarray() for method in samples_hubred.keys()])}
neighbors_id = dict()
for key in adj_matrix.keys():
    N_nN = np.sum(adj_matrix[key][points_right_hyperplan, :], axis=0) > 0
    neighbors_id[key] = []
    for i in range(len(N_nN)):
        if N_nN[i]:
            if points_right_hyperplan[i]:
                neighbors_id[key].append(1)
            else:
                neighbors_id[key].append(2)
        else:
            if points_right_hyperplan[i]:
                neighbors_id[key].append(3)
            else:
                neighbors_id[key].append(0)
for key in samples_hubred.keys():
    plt.clf()
    plot_id(samples_skimmed, neighbors_id[key],
            title='Random axes with neighbors, hub reduction '+key)
    plt.show()
percentage_status = dict()
for key in neighbors_id.keys():
    percentage_status[key] = [int(100*np.sum([neighbors_id[key][loop]==0 for loop in range(n_obs)])/len(idx_left_hyperplan)),
                              int(100*np.sum([neighbors_id[key][loop]==1 for loop in range(n_obs)])/len(idx_right_hyperplan)),
                              int(100*np.sum([neighbors_id[key][loop]==2 for loop in range(n_obs)])/len(idx_left_hyperplan)),
                              int(100*np.sum([neighbors_id[key][loop]==3 for loop in range(n_obs)])/len(idx_right_hyperplan))]
L_nN = [percentage_status[key][0] for key in percentage_status.keys()]
R_N = [percentage_status[key][1] for key in percentage_status.keys()]
L_N = [percentage_status[key][2] for key in percentage_status.keys()]
R_nN = [percentage_status[key][3] for key in percentage_status.keys()]
x = np.arange(len(labels))
width = 0.2
plt.clf()
fig, ax = plt.subplots()
rects1 = ax.bar(x - 1.5*width, L_nN, width, label='L_nN', color=colors[0])
rects2 = ax.bar(x - 0.5*width, L_N, width, label='L_N', color=colors[2])
rects3 = ax.bar(x + 0.5*width, R_N, width, label='R_N', color=colors[1])
rects4 = ax.bar(x + 1.5*width, R_nN, width, label='R_nN', color=colors[3])
ax.set_ylabel('Percentage among the (L, R) groups')
ax.set_title('Repartition of the different populations in the kNN graph')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
fig.tight_layout()
plt.show()
if verbose:
    print("kNN graph repartition visualized")

### kNN-graph visualisation
# get the in and out degrees of all points
out_degree = dict()
in_degree = dict()
for key in adj_matrix.keys():
    out_degree[key] = dict()
    in_degree[key] = dict()
    for idx in tqdm(range(n_obs)):
        out_degree[key][idx] = [i for i, e in enumerate(adj_matrix[key][idx, :]) if e==True]
        in_degree[key][idx] = [i for i, e in enumerate(adj_matrix[key][:, idx]) if e==True]
# make the graph
for key in adj_matrix.keys():
    rows, cols = np.where(adj_matrix[key] == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    in_degree_weight = []
    for idx in range(n_obs):
        in_degree_weight.append(len(in_degree[key][idx]))
    plt.clf()
    nx.draw_networkx(gr, node_size=in_degree_weight[0:1000],
                     node_color=[int(points_right_hyperplan[loop]) for loop in range(n_obs)][0:1000],
                     cmap='plasma', width=0.1, edge_color='gray', with_labels=False, nodelist=range(1000))
    plt.show()

# Density
# Pick n random points in each region
# Pick source in the middle
norm = np.linalg.norm(samples_skimmed, axis=1)
idx_middle = np.argmin(norm)
idx_left_test = random.choices(np.array(idx_left_hyperplan), k=n_points_density_target)
idx_right_test = random.choices(np.array(idx_right_hyperplan), k=n_points_density_target)
# Check that start and end are different
while (np.array([idx_right_test[loop]==idx_middle for loop in range(n_points_density_target)]).any()==True
    or np.array([idx_left_test[loop]==idx_middle for loop in range(n_points_density_target)]).any()==True):
    if verbose:
        print("Pick new points")
    idx_left_test = random.choices(np.array(idx_left_hyperplan), k=n_points_density_target)
    idx_right_test = random.choices(np.array(idx_right_hyperplan), k=n_points_density_target)
# Plot the position of these points
point_test_status = np.repeat("Neutral", n_obs)
point_test_status[idx_left_hyperplan] = 0
point_test_status[idx_right_hyperplan] = 1
point_test_status[np.array(idx_middle)] = 4
point_test_status[np.array(idx_left_test)] = 2
point_test_status[np.unique(idx_right_test)] = 3
plt.clf()
plot_id(samples_skimmed, point_test_status,
        title='random axes, dense ball')
plt.show()

# Get reference density value
gaussian_kde = stats.gaussian_kde(samples_skimmed.T)
density_left = np.mean(gaussian_kde.logpdf(samples_skimmed[idx_left_test, :].T))
density_right = np.mean(gaussian_kde.logpdf(samples_skimmed[idx_right_test, :].T))
if verbose:
    print("true log density ratio is:" + str(density_left/density_right))
eta_d = np.pi**(original_dims/2)/special.gamma(original_dims/2+1)
v_d = intersection_volume(original_dims)
value_base = (density_left - density_right)/(eta_d/k/v_d)
if verbose:
    print("true log-corrected density difference is:" + str(value_base))

# Density estimate
# Plot the distribution of in degree for each key
plt.clf()
for key in adj_matrix.keys():
    plt.hist([len(in_degree[key][loop]) for loop in range(n_obs)])
    plt.title('Histogram of in degree, hub correction '+key)
    plt.show()
if verbose:
    print("in and out degree done")
# Write the graph in the dict form
dict_knn = dict()
for key in adj_matrix.keys():
    dict_knn[key] = dict()
    for idx in tqdm(range(n_obs)):
        val = [i for i, e in enumerate(adj_matrix[key][idx, :]) if e==True]
        dict_knn[key][idx] = {k:v for k, v in zip(val, np.repeat(1, len(val)))}
# Find the shortest path from start to end (Dijkstra)
path_collection_left = dict()
path_collection_right = dict()
for key in dict_knn.keys():
    path_collection_left[key] = []
    path_collection_right[key] = []
    for iter in tqdm(range(n_points_density_target)):
        path_collection_left[key].append(Dijkstra(dict_knn[key], idx_middle, idx_left_test[iter]))
        path_collection_right[key].append(Dijkstra(dict_knn[key], idx_middle, idx_right_test[iter]))
# remove unreachable path
idx_rm_left = dict()
idx_rm_right = dict()
for key in dict_knn.keys():
    idx_rm_left[key] = []
    idx_rm_right[key] = []
    for i in range(n_points_density_target):
        if len(path_collection_left[key][i]) < 3:
            idx_rm_left[key].append(i)
        if len(path_collection_right[key][i]) < 3:
            idx_rm_right[key].append(i)
idx_left_test = {key:idx_left_test for key in dict_knn.keys()}
idx_right_test = {key:idx_right_test for key in dict_knn.keys()}
for key in dict_knn.keys():
    path_collection_left[key] = [element for i, element in enumerate(path_collection_left[key]) if i not in idx_rm_left[key]]
    path_collection_right[key] = [element for i, element in enumerate(path_collection_right[key]) if i not in idx_rm_right[key]]
    idx_left_test[key] = [element for i, element in enumerate(idx_left_test[key]) if i not in idx_rm_left[key]]
    idx_right_test[key] = [element for i, element in enumerate(idx_right_test[key]) if i not in idx_rm_right[key]]
if verbose:
    print("shortest path finding done")
# Plot the paths
plt.clf()
for key in dict_knn.keys():
    plt.scatter(samples_skimmed[:, 0], samples_skimmed[:, 1], color="black")
    for path in path_collection_right[key]:
        plt.plot(samples_skimmed[path, 0], samples_skimmed[path, 1], color="red")
    for path in path_collection_left[key]:
        plt.plot(samples_skimmed[path, 0], samples_skimmed[path, 1], color="gold")
    plt.show()
# get the left and right ensembles for each path + the density
density_left = dict()
density_right = dict()
for key in dict_knn.keys():
    density_left[key] = []
    density_right[key] = []
    for path in tqdm(path_collection_left[key]):
        if not path is None:
            len_path = len(path)
            if len_path > 2:
                left = []
                right = []
                for step in range(1, (len_path-1)):
                    x = path[step]
                    x_l = path[step-1]
                    x_r = path[step+1]
                    out_x = out_degree[key][x]
                    in_x_l = in_degree[key][x_l]
                    in_x_r = in_degree[key][x_r]
                    left.append(len(np.intersect1d(out_x, in_x_l)))
                    right.append(len(np.intersect1d(out_x, in_x_r)))
                density_left[key].append(np.sum(np.array(right)-np.array(left)))
    for path in tqdm(path_collection_right[key]):
        if not path is None:
            len_path = len(path)
            if len_path > 2:
                left = []
                right = []
                for step in range(1, (len_path-1)):
                    x = path[step]
                    x_l = path[step-1]
                    x_r = path[step+1]
                    out_x = out_degree[key][x]
                    in_x_l = in_degree[key][x_l]
                    in_x_r = in_degree[key][x_r]
                    left.append(len(np.intersect1d(out_x, in_x_l)))
                    right.append(len(np.intersect1d(out_x, in_x_r)))
                density_right[key].append(np.sum(np.array(right)-np.array(left)))
# we get target - source -> negative !! instead make L-R => source - target (hence diff in logG - logD = (L-R)d - (L-R)g
density_mean_left = []
density_mean_right = []
for key in dict_knn.keys():
    density_mean_left.append(np.mean(density_left[key]))
    density_mean_right.append(np.mean(density_right[key]))
value_estimate = density_mean_left[0] - density_mean_right[0]
if verbose:
    print("true density difference is :" + str(value_base))
    print("estimated density difference is:" + str(value_estimate))
    print(density_mean_left)
    print(density_mean_right)
# Visualize the different densities
methods = [key for key in dict_knn.keys()]
range_density_left = [item for sublist in [v for v in density_left.values()] for item in sublist]
range_density_right = [item for sublist in [v for v in density_right.values()] for item in sublist]
range_density_order = np.argsort(range_density_right+range_density_left)
length_l = {key:len(density_left[key]) for key in methods}
length_r = {key:len(density_right[key]) for key in methods}
range_density_order_l = dict()
range_density_order_r = dict()
idx_l = len(range_density_right)
idx_r = 0
for key in methods:
    range_density_order_l[key] = range_density_order[idx_l:(idx_l+length_l[key])]
    range_density_order_r[key] = range_density_order[idx_r:(idx_r+length_r[key])]
    idx_l+=length_l[key]
    idx_r+=length_r[key]
for key in dict_knn.keys():
    plt.clf()
    for i in range(length_r[key]):
        plt.plot(samples_skimmed[[idx_middle, idx_left_test[key][i]], 0],
                 samples_skimmed[[idx_middle, idx_left_test[key][i]], 1], color=cm(1-range_density_order_l[key][i]/len(range_density_order)))
    for i in range(length_r[key]):
        plt.plot(samples_skimmed[[idx_middle, idx_right_test[key][i]], 0],
                 samples_skimmed[[idx_middle, idx_right_test[key][i]], 1], color=cm(1-range_density_order_r[key][i]/len(range_density_order)))
    plt.axvline(x=0, color="black")
    plt.title('Random axes, hub correction '+key)
    plt.show()

# Viz the difference
xy = [[x, y] for x, y in zip(range(1, len(methods)+1), np.array(density_mean_right))]
height = np.array(density_mean_left)-np.array(density_mean_right)
col = np.repeat("green", len(methods))
col[height < 0] = "red"
x = 1+np.arange(len(labels))
width = 0.8
plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', adjustable='box')
for i in range(len(methods)):
    ax.add_patch(
    patches.Rectangle(xy=xy[i], height=height[i], width=width, color=col[i]))
plt.xlim((0, 6))
plt.ylim((min(min(density_mean_right), np.min(np.array(density_mean_right)+height))-2, np.max(np.array(density_mean_right)+height)+2))
ax.set_ylabel('Density')
ax.set_title('Density difference between the two half spheres')
ax.set_xticks(x)
ax.set_xticklabels(labels, Rotation=90)
plt.show()
