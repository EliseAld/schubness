# https://pypi.org/project/scikit-hubness/
# https://scikit-hubness.readthedocs.io/en/latest/documentation/concepts.html
import numpy as np
import pandas as pd
import skhubness
from skhubness import Hubness
from sklearn.model_selection import cross_val_score
from skhubness.neighbors import KNeighborsClassifier
from skhubness.neighbors import kneighbors_graph

# load data
satija = pd.read_csv("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv", index_col=0)
print(f'satija = {satija.shape}') # 36280 x 8617

# evaluate hubness
hub = Hubness(k=10, metric='cosine')
hub.fit(satija)
k_skew = hub.score()
print(f'Skewness = {k_skew:.3f}')
print(f'Robin hood index: {hub.robinhood_index:.3f}')
print(f'Antihub occurrence: {hub.antihub_occurrence:.3f}')
print(f'Hub occurrence: {hub.hub_occurrence:.3f}')

# There is considerable hubness
# Let's see, whether hubness reduction works
# Did it actually reduce hubness?
hub_mp = Hubness(k=10, metric='cosine',
                 hubness='mp')
hub_ls = Hubness(k=10, metric='cosine',
                 hubness='ls')
hub_dsl = Hubness(k=10, metric='euclidean',
                 hubness='dsl')
hub_mp.fit(satija)
hub_ls.fit(satija)
hub_dsl.fit(satija)
k_skew_mp = hub_mp.score()
k_skew_ls = hub_ls.score()
k_skew_dsl = hub_dsl.score()

print(f'Skewness after MP: {k_skew_mp:.3f} '
      f'(reduction of {k_skew - k_skew_mp:.3f})')
print(f'Robin hood: {hub_mp.robinhood_index:.3f} '
      f'(reduction of {hub.robinhood_index - hub_mp.robinhood_index:.3f})')
print(f'Skewness after LS: {k_skew_ls:.3f} '
      f'(reduction of {k_skew - k_skew_ls:.3f})')
print(f'Robin hood: {hub_ls.robinhood_index:.3f} '
      f'(reduction of {hub.robinhood_index - hub_ls.robinhood_index:.3f})')
print(f'Skewness after DSL: {k_skew_dsl:.3f} '
      f'(reduction of {k_skew - k_skew_dsl:.3f})')
print(f'Robin hood: {hub_dsl.robinhood_index:.3f} '
      f'(reduction of {hub.robinhood_index - hub_dsl.robinhood_index:.3f})')

