data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/data/Pancreatic/data.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Pancreatic/afterPCA/pancreatic_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Pancreatic/entropy_lei/pancreatic_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Pancreatic.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Pancreatic.rds"
