data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/data.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/afterPCA/baron_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/entropy_lei/baron_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Baron.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Baron.rds"
