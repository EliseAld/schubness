data <- read.table("/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/Guo_lung_2018/GSE99254_NSCLC.TCell.S9055.count.labeled.txt",header=T)
#data_pca <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Guo/afterPCA/GSE99254_NSCLC.TCell.S9047.count.labeled.pca_readyforhubness.txt")
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Guo/entropy_lei/GSE99254_NSCLC.TCell.S9055_scentropy_v1.csv")
# Remove cells that look abnormal
nFeature <- apply(data,2,function(x) sum(x!=0))
nCount <- apply(data,2,sum)
data <- data[,nFeature<7500 & nCount<35e+05] # 9047 cells for Guo
entropy <- entropy[nFeature<7500 & nCount<35e+05,]
rm(nCount,nFeature)
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_guo.rds"
path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_guo.rds"
