dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC")
data <- pblapply(dataset,
                 function(x) read.csv(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/duo/counts/",x,".csv"), sep=""))
data_pca <- pblapply(dataset,
                     function(x) read.csv(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/duo/pca/",x,"_pca_readyforhubness.csv"), sep=""))
#entropy <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/data/Baron/entropy_lei/baron_scentropy_v1.csv")
#path_neighbors <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_bulk.rds"
#path_neighbors_pca <- "/Users/elise/Desktop/GitHub/Hubness_sc/R_scripts/hub_network/neighbors_pca_bulk.rds"
