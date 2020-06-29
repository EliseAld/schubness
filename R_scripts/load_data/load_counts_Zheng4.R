data <- read.csv("/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/DimRedPaper/Zheng4/data.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Zheng4/afterPCA/zheng4_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Zheng4/entropy_lei/zheng4_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_Zheng4.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_Zheng4.rds"
