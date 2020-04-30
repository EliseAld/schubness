data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/data/Darmanis/data.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Darmanis/afterPCA/darmanis_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Darmanis/entropy_lei/darmanis_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Darmanis.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Darmanis.rds"
