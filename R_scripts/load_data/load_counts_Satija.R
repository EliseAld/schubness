data <- read.csv("/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/10x/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv", row.names = 1)
entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Satija/entropy_lei/GSE100866_scentropy_v1.csv")
# Remove cells that look abnormal
data <- data[grep("HUMAN_",rownames(data)),] # 20400 genes
rownames(data) <- gsub("HUMAN_","",rownames(data))
nFeature <- apply(data,2,function(x) sum(x!=0))
nCount <- apply(data,2,sum)
data <- data[,nFeature<4000 & nCount<15e+03] # 8604 cells for Satija
entropy <- entropy[nFeature<4000 & nCount<15e+03,]
rm(nCount,nFeature)
path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_satija.rds"
