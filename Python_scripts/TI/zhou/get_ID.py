import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
import seaborn as sns


# params
scale = True
n_comps = (25, 50, 100, 500)
do_norm = ('seurat', 'duo')
params_list = list(product(n_comps, do_norm))

gid = {}
lid = {}
for dim, norm in params_list:
    with open("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/ID/"+f'all_id_TI_{dim}_{norm}_True.pkl','rb') as f:
        gid[(dim, norm, scale)], lid[(dim, norm, scale)] = pickle.load(f)
gid={k:pd.DataFrame(v[k]).assign(metric=str(k)) for k,v in gid.items()}
lid={k:pd.DataFrame(v[k]).assign(metric=str(k)) for k,v in lid.items()}

n_datasets = gid[(25, 'duo', True)].shape[0]
methods = [name for name in gid[(25, 'duo', True)].columns][:4]

# Plot the PCA GID estimate
gid_sns = dict()
for method in methods:
    gid_sns[method] = pd.DataFrame(np.zeros((n_datasets*len(n_comps)*len(do_norm), 3)), columns=['GID', 'Recipe', 'Dimension'])
    gid_sns[method]['Dimension'] = np.tile(np.repeat(n_comps, n_datasets), len(do_norm))
    gid_sns[method]['Recipe'] = np.repeat(do_norm, len(n_comps)*n_datasets)
    for key in gid.keys():
        idx_range = np.intersect1d(np.argwhere(np.array(gid_sns[method]['Dimension']==key[0])),
                                       np.argwhere(np.array(gid_sns[method]['Recipe']==key[1])))
        gid_sns[method]['GID'][idx_range] = [gid[key][method][key3] for key3 in gid[key][method].keys()]
for method in methods:
    sns.boxplot(x="Dimension", y="GID", hue="Recipe", data=gid_sns[method], palette="Set1", fliersize=0)
    sns.stripplot(x="Dimension", y="GID", hue="Recipe", data=gid_sns[method], jitter=True, dodge=True)
    plt.ylabel('GID '+method)
    plt.hlines(y=25, xmin=-.5, xmax=0.5, colors="red")
    plt.hlines(y=50, xmin=.5, xmax=1.5, colors="red")
    plt.hlines(y=100, xmin=1.5, xmax=2.5, colors="red")
    plt.hlines(y=500, xmin=2.5, xmax=3.5, colors="red")
    plt.semilogy()
    plt.show()

# Retrieve name of datasets that have GID = 25 with seurat and 25
gid[(25, 'seurat', True)]['pca']>=25

