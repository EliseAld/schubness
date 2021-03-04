#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param k size of the neighborhood to consider (default to the square root of the cardinality)
#' @param n_dim number of PCs to retain to build the kNN graph

#' @return in-degree scores of the data samples

#' @example in_degree <- get_in_degree(data, k=10, n_dim=50)

### Load librairies
library(RANN)

#### Functions
# Make the kNN graph
kNN_ <- function(data, k) {
  k.nn <- RANN::nn2(t(data[1:n_dim,]), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  return(k.nn)
}

# Get the in-degree from the kNN graph
k.occurence_ <- function(k.nn) {
  occurence <- sapply(seq(nrow(k.nn)), function(x) sum(k.nn==x))
  return(occurence)
}

# Get the in-degree from the data
get_in_degree <- function(data, k, n_dim) {
  data_pc <- data[1:n_dim,]
  score <- k.occurence_(kNN_(data_pc, k))
  return(score)
}