library(future.apply)
library(pbapply)
plan("multiprocess",workers=30)

data_brca <- read.csv("/shared/projects/mne/elise/bulk/TCGA_BRCA/brca.csv", row.names=1)
data_archs4 <- read.csv("/shared/projects/mne/elise/bulk/ARCHS4/data.csv", row.names=1)
load("/shared/projects/mne/elise/bulk/TCGA_KIRC/TCGA_KIRC_RPKM_withHGNCsymbol.Rdata")
data_kirc <- matrix
rm(matrix)
data <- list(data_brca, data_kirc, data_archs4)
names(data) <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")


generate_dropout <- function(array, dropout_percent) {
  if (dropout_percent==0) {
    return(array)
  }
  else {
  n_genes <- nrow(array)
  zeros <- sapply(seq(ncol(array)),
                  function(x) sample(seq(n_genes)[array[,x]!=0], 1+round(sum(array[,x]!=0)*dropout_percent/100)))
  simul_data <- array
  for (i in seq(ncol(array))) {
    simul_data[zeros[[i]],i] <- 0
  }
  return(simul_data)
  }
}
zero_percent <- function(data) { # gene x cell matrix
  return(sum(data==0)/ncol(data)/nrow(data)*100)
}

dropout_percent <- c(0,10,20,30,40,50,60,70,80,90,95)
data_simul <- future_lapply(data,
                     function(y) future_lapply(dropout_percent,function(x) generate_dropout(y,x)))
sapply(data_simul,
       function(x) sapply(x,zero_percent))
mapply(x=data_simul, y=names(data_simul),
       function(x,y) pbmapply(t=x, z=dropout_percent,
                            function(t,z) write.table(t, file=paste0("/shared/projects/mne/elise/bulk/", y, "/simulated/simul_dropout", z,".csv"))))

print("SUCCESS SIMU")
