data <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX3/pdx3_merged.txt",header=T)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX3/afterPCA/pdx3_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX3/entropy_lei/pdx3_raw_merged_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pdx3.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_pdx3.rds"
