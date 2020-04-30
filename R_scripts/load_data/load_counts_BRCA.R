data <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/brca.csv",header=T, row.names = 1)
data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/bulkRNAseq/brca_pca_readyforhubness.txt")
#entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/entropy_lei/baron_scentropy_v1.csv")
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_BRCA.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_BRCA.rds"
