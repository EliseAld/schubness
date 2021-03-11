#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param k size of the neighborhood to consider (default to the square root of the cardinality)
#' @param n_dim number of PCs to retain to build the kNN graph
#' @param score vector of in-degree scores

#' @return hubness magnitude evaluator

#' @example skewness <- evaluator_skewness(score)

### Load librairies
library(knn.covertree)

### Functions to compute the various evaluators of hubness
evaluator_max <- function(score) { # Maximum hubness score
   return(max(score)/length(score))
}

evaluator_mu_sd <- function(score) { # Counting nb of cells with an in-degree above mean + 3*sd (Tomasev et al.)
   thd = mean(score)+3*sd(score)
   return(sum(score>=thd)/length(score)*100)
}

evaluator_2k <- function(score, k, log=F) { # Counting nb of cells with an in-degree above 2k
   if (log==T) {
      hub_nb = sum(score>=2*log(k))
   }
   else {
      hub_nb = sum(score>=2*k)
   }
   return(hub_nb/length(score)*100)
}

evaluator_skewness <- function(score) { # skewness
   return(mean((score-mean(score))^3)/sd(score)^3)
}

evaluator_antihubs <- function(score) { # Counting antihubs
   return(sum(score==0)/length(score)*100)
}

count_uniD_edges_ <- function(knn_graph) {
   increment = 0
   for (i in 1:nrow(knn_graph)) {
      for (j in 1:ncol(knn_graph)) {
         if (i %in% knn_graph[knn_graph[i,j],]) {
            increment <- increment + 1 
         }
      }
   }
   return(100*(1-increment/(nrow(knn_graph)*ncol(knn_graph))))
}

evaluator_asymmetry <- function(data, k, n_dim) { # asymmetry
  knn_graph <- knn.covertree::find_knn(t(data[1:n_dim,]),k)$index
  asym <- count_uniD_edges_(knn_graph)
  return(asym)
}