data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/data/Koh/data.csv",header=T, row.names = 1, sep=";", dec = ",")
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Koh/afterPCA/koh_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Koh/entropy_lei/koh_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Koh.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Koh.rds"
