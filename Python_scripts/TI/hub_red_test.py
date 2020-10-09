# https://pypi.org/project/scikit-hubness/
# https://scikit-hubness.readthedocs.io/en/latest/documentation/concepts.html
# https://scikit-hubness.readthedocs.io/en/latest/getting_started/example.html

import pandas as pd
from skhubness import Hubness
from skhubness.neighbors import kneighbors_graph

# load data
data = pd.read_csv("/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/GSE67602_Joost_et_al_expression.txt", index_col=0, sep="\t")
print(f'data = {data.shape}') # 36280 x 8617

# evaluate hubness before anything
skew = Hubness(k=10, metric='cosine').fit(data).score()
print(f'Skewness = {skew:.3f}')

# Apply hub reduction
hubness_methods = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'})
                   #'ls':('ls',None),
                   #'ls_nicdm':('ls',{'method': 'nicdm'}),
                   #'dsl':('dsl',None)
                   }
hubness_reduction_data = dict()
for method_name, (hubness, hubness_params) in hubness_methods.items():
    hub_red = Hubness(k=10, hubness=hubness, hubness_params=hubness_params, metric='euclidean')
    hub_red.fit(data)
    hubness_reduction_data[method_name] = hub_red
    skew = hub_red.score()
    print(method_name+' completed')
    print(f'Skewness: {skew:.3f}')

# Export the corrected kNN graph
hubness_reduction_graph = dict()
for method_name, (hubness, hubness_params) in hubness_methods.items():
    hubness_reduction_data[method_name] = kneighbors_graph(data,
                                      n_neighbors=10,
                                      hubness=hubness,
                                      hubness_params=hubness_params)
    print(method_name+' completed')

