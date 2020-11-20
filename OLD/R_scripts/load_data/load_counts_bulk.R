#library(rhdf5)
#h5ls("/Users/elise/Downloads/human_matrix.h5")
#data <- h5read("/Users/elise/Downloads/human_matrix.h5", "/data")

#python
#import h5py
#import numpy as np
#import pandas as pd
#f1 = h5py.File('human_matrix.h5','r+')  
#data=f1['data/expression']
#genes=f1['meta/genes']
#df=np.array(data[:1000,:])
#df=pd.DataFrame(df)
#df.to_csv('data.csv')
#df = np.array(genes)
#df=pd.DataFrame(df)
#df.to_csv('genes.csv')

#data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/data1.csv",header=T, row.names = 1)
#genes <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/genes.csv", row.names = 1)
#library(stringr)
#genes <- sapply(genes,
                function(x) str_remove(x, 'b'))
#colnames(data) <- genes
#write.csv(data,file="/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/data.csv")
data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/bulkandrei/data.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/bulkandrei/bulk_simul0_pca_readyforhubness.txt")
#entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/entropy_lei/baron_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_bulk.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_bulk.rds"
