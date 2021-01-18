# import required packages
import scEntropy.scEntropy as scEntropy
import pandas as pd
import numpy as np

# fixed param
dataset = ["Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq"]

# DATA PATH TO CHANGE AND EXPORT PATH
data_file_path = '/Users/elise/Desktop/Github/Hubness_sc/Data/Duo_10kHVG/pca/'
export_entropy_path = '/Users/elise/Desktop/Github/Hubness_sc/Data/Duo_10kHVG/entropy/'

# load data
data = {k:v for k, v in zip(dataset, [pd.read_csv(data_file_path+data+'_pca_readyforhubness.csv',
                                                  index_col=0, header=0, sep=" ") for data in dataset])}
[data[key].shape for key in data.keys()]
# data is a dict of PCs by cell matrices

# Run single cell entropy for all cells
scentropy = {k:v for k, v in zip(dataset,
                                 [scEntropy.scEntropy(data[key], ref_vec=None, option='RCSA')
                                  for key in data.keys()])}

# Export scentropy
for data in dataset:
    np.savetxt(export_entropy_path+data+'.csv', scentropy[data].tolist(), delimiter=',')
