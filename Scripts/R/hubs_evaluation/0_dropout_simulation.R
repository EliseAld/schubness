#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param dropout_percent percentage of dropout added to the expression matrix

#' @return expression matrix with additional simulated dropout

#' @example simulated_data <- generate_dropout(data, 10)
#' @example splatter_data <- splatSimDropout(data, 10)

### Functions to generate dropout
source("~/evaluation_common_functions.R")
# Generate random dropout
generate_dropout <- function(data, dropout_percent) {
  if (dropout_percent==0) {
    return(data)
  }
  else {
    n_genes <- nrow(data)
    zeros <- sapply(seq(ncol(data)),
                    function(x) sample(seq(n_genes)[data[,x]!=0], 1+round(sum(data[,x]!=0)*dropout_percent/100)))
    simul_data <- data
    for (i in seq(ncol(data))) {
      simul_data[zeros[[i]],i] <- 0
    }
    return(simul_data)
  }
}

# Generate splatter dropout
splatSimDropout <- function(data, dropout_percentage) {
    dropout.mid <- -10
    dropout.shape <- -1
    nCells <- ncol(data)
    nGenes <- nrow(data)
    drop.prob <- vapply(seq_len(nCells), function(idx) {
        eta <- log(data[,idx])
        return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
    }, as.numeric(seq_len(nGenes)))
    keep <- matrix(rbinom(nCells * nGenes, 1, 1-drop.prob),
                   nrow = nGenes, ncol = nCells)
    counts <- list(data*keep)
    #print(paste0("Adding ", round(mean(keep[data!=0]==0)*100), "% of dropout"))
    for (idx in seq(2,length(dropout_percentage))) {
      while (mean(keep[data!=0]==0) < dropout_percentage[idx]/100) {
        dropout.mid <- dropout.mid + 0.5
        drop.prob <- vapply(seq_len(nCells), function(idx) {
        eta <- log(data[,idx])
        return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        }, as.numeric(seq_len(nGenes)))
        keep <- matrix(rbinom(nCells * nGenes, 1, 1-drop.prob),
                   nrow = nGenes, ncol = nCells)
      }
      counts <-  c(counts, list(data*keep))
      #print(paste0("Adding ", round(mean(keep[data!=0]==0)*100), "% of dropout"))
    }
    names(counts) <-  dropout_percentage
    return(counts)
}

logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}
