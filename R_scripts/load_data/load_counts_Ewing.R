data <- read.table("/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/Andrei_Ewing/ASP14_TS_u.txt",header=T)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Ewing/afterPCA/ewing_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Ewing/entropy_lei/scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Ewing.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Ewing.rds"
