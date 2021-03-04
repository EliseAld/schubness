#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param n_hvg number of Highly Variable Genes to retain for the PCA (default set to 10,000)

#' @return PCA projection of the data (PCs in rows, cell in columns)

#' @example pca <- data_PCA(data_prePCA(data))

### Functions
source("~/evaluation_common_functions.R")
data_prep <- function(data) {
  data <- log10(data+1)
  return(data)
}

data_prePCA <- function(data, n_hvg=1e4) {
  hvg<-names(sort(apply(data,1,var),decreasing=T)[1:n_hvg])
  data_filter = data[hvg,]
  return(data_prep(data_filter))
}

data_PCA <- function(data) {
  return(t(prcomp(data, center=T, scale.=F)$rotation))
}