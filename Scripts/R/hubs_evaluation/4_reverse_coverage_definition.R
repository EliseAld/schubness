#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param k size of the neighborhood to consider (default to the square root of the cardinality)
#' @param score vector of in-degree scores
#' @param n number of putative hubs to consider
#' @param n_val range of values for the parameter n
#' @param thd percentage of increment in reverse coverage size below which we consider we reached a plateau (default set to 1)

#' @return hubness magnitude evaluator

#' @example hubs <- names(sort(score, decreasing = T)[1:get_hub_nb(data, k, score, n_val, thd)])

### Load librairies
library(RANN)

### Functions to define hubs via the reverse coverage method
# Make the kNN graph
kNN_ <- function(data, k) {
  k.nn <- RANN::nn2(t(data[1:n_dim,]), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  return(k.nn)
}

# Calculate the size of the reverse coverage of specific cells
rev_cov_hub <- function(data, k, score, n) {
  knn_graph = kNN_(data, k)
  putative_hubs <- which(order(score, decreasing = T) %in% seq(n))
  size_cov <- sum(apply(knn_graph, 1, function(x)
    any(x %in% putative_hubs)))
  return(size_cov/nrow(knn_graph)*100)
}

rev_cov_ctrl <- function(data, k, n) {
  knn_graph = kNN_(data, k)
  putative_hubs <- sample(seq(nrow(knn_graph)), n)
  size_cov <- sum(apply(knn_graph, 1, function(x)
    any(x %in% putative_hubs)))
  return(size_cov/nrow(knn_graph)*100)
}

# Compute the increase in size
seq_size_hub <- function(data, k, score, n_val) {
  knn_graph <- kNN_(data, k)
  idx=1
  cov_size=rev_cov_hub(knn_graph, score, n_val[idx])
  while (cov_size[idx]<100 & idx < length(n_val) & n_val[idx] <= ncol(data)) {
    idx=idx+1
    cov_size=c(cov_size, rev_cov_hub(knn_graph, score, n_val[idx]))
  }
  return(data.frame("N_percent"=n_val[seq(idx)]/ncol(data),
                    "N"=n_val[seq(idx)],
                    "Size"=cov_size))
}

seq_size_ctrl <- function(data, k, score, n_val) {
  knn_graph <- kNN_(data, k)
  idx=1
  cov_size=rev_cov_ctrl(knn_graph, n_val[idx])
  while (cov_size[idx]<100 & idx < length(n_val) & n_val[idx] <= ncol(data)) {
    idx=idx+1
    cov_size=c(cov_size, rev_cov_ctrl(knn_graph, n_val[idx]))
  }
  return(data.frame("N_percent"=n_val[seq(idx)]/ncol(data),
                    "N"=n_val[seq(idx)],
                    "Size"=cov_size))
}

# Detect the plateau
get_plateau <- function(size_df, n_val, thd=1) {
  increment = diff(size_df[,"Size"])
  return(min(ifelse(all(increment>=thd),
                                  length(increment)+1,
                                  which(increment<thd)[1]+1),
                           length(n.val)))
}

# Retrieve the number of hubs
get_hub_nb <- function(data, k, score, n_val, thd) {
  seq_size_df = seq_size_hub(data, k, score, n_val)
  return(seq_size_df[,"N"][get_plateau(seq_size_df, n_val, thd)])
}

