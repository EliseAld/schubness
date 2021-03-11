# import required packages
import scEntropy.scEntropy as scEntropy  # install scEntropy as in https://github.com/jzlei/scEntropy
import pandas as pd
import numpy as np

# fixed param
dataset = ['dataset1.csv', 'dataset2.csv', ...]
# list of datasets

# DATA PATH TO CHANGE AND EXPORT PATH
data_file_path = '/path/to/data/'  # PCA data
export_entropy_path = '/path/to/entropy/results/'

# load data
data = {k:v for k, v in zip(dataset, [pd.read_csv(data_file_path+data,
                                                  index_col=0, header=0, sep=" ") for data in dataset])}
# [data[key].shape for key in data.keys()]
# data is a dictionary of PCs x cell matrices

# Run single cell entropy for all cells
scentropy = {k:v for k, v in zip(dataset,
                                 [scEntropy.scEntropy(data[key], ref_vec=None, option='RCSA')
                                  for key in data.keys()])}

# Export scentropy
for data in dataset:
    np.savetxt(export_entropy_path+data, scentropy[data].tolist(), delimiter=',')
