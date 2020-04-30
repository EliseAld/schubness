data <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX1/pdx1_raw_merged.txt",header=T)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX1/afterPCA/pdx1_pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/PDX1/entropy_lei/pdx1_raw_merged_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pdx1.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_pdx1.rds"
