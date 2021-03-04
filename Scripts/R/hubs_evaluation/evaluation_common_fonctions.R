#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param k size of the neighborhood to consider (default to the square root of the cardinality)

#' @return PCA projection of the data (PCs in rows, cell in columns)

#' @example pca <- data_PCA(data_prePCA(data))

### Functions
# Make the kNN graph
kNN_ <- function(data, k) {
  k.nn <- RANN::nn2(t(data[1:n_dim,]), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  return(k.nn)
}


# Count the sparsity
sparsity <- function(data) { # gene x cell matrix
  return(sum(data==0)/ncol(data)/nrow(data)*100)
}