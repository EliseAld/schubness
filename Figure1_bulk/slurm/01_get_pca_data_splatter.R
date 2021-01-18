# path and functions
library(future.apply)
library(pbapply)
library(ggplot2)
library(ggpubr)
plan("multiprocess",workers=30)

dropout_percent <- c(0,10,20,30,40,50,60,70,80,90,95)
datasets <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")
data_in <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/shared/projects/mne/elise/bulk/", x, "/splatter_simu/simul_dropout", y,".csv")))
data_out <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/shared/projects/mne/elise/bulk/", x, "/splatter_pca/simul_dropout", y,"_pca.csv")))
sdev_out <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/shared/projects/mne/elise/bulk/", x, "/splatter_pca/simul_dropout", y,"_sdev.csv")))

data_prep <- function(data) {
  data <- log10(data+1)
  return(data)
}


# load data
data <- future_lapply(X=data_in, FUN=function(x) future_lapply(x,
                                                   function(y) read.table(file=y)))
# Keep 10k HVG
data <- future_lapply(data, function(x) future_lapply(x,
                                          function(y) {hvg<-names(sort(apply(y,1,var),decreasing=T)[1:1e4]); return(y[hvg,])}))
# Log transfo
data <- future_lapply(data, function(y) future_lapply(y, function(x) data_prep(x)))
# PCA
pca <- future_lapply(data, function(y) future_lapply(y, function(x) prcomp(x, center=T, scale.=F)))
data_proj <- lapply(pca, function(y) pblapply(y, function(x) t(x$rotation)))

mapply(z=data_proj, t=data_out, function(z,t) pbmapply(x=z, y=t, function(x,y) write.table(x, file=y)))
mapply(z=pca, t=sdev_out, function(z,t) pbmapply(x=z, y=t, function(x,y) write.table(x$sdev, file=y)))
print("SUCCESS PCA SPLATTER")
