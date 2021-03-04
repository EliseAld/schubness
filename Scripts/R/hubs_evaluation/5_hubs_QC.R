#' @param data matrix or dataframe object of the expression matrix (gene in rows, cell in columns)
#' @param k size of the neighborhood to consider (default to the square root of the cardinality)
#' @param entropy scentropy cellular vector calculated with python
#' @param score vector of in-degree scores
#' @param n_val range of values for the parameter n
#' @param thd percentage of increment in reverse coverage size below which we consider we reached a plateau (default set to 1)
#' @param n_dims number of dimensions to consider for the UMAP (default set to 50)

#' @return Seurat object with the QC in the metadata slot, and UMAP and PCA coordinates

#' @example seurat <- QC_computation(data, k, score, n_val, thd, entropy)

### Load libraires
library(Seurat)
source("~/3_evaluators.R")
source("~/5_hubs_QC.R")

### Functions to perform QC
entropy

QC_computation <- function(data, k, score, n_val, thd, entropy, n_dims=50) {
  names(score) <- colnames(data)
  hub_nb <- get_hub_nb(data, k, score, n_val, thd)
  hubs <- names(sort(score, decreasing = T)[1:hub_nb])
  antihubs <- names(score)[score==0]
  seurat <- Seurat::CreateSeuratObject(data)
  seurat$Antihubs <- c("normal","antihub")[factor(colnames(seurat) %in% antihubs,
                                                  levels = c("FALSE","TRUE"))]
  seurat$Hubs <- c("normal","hub")[factor(colnames(seurat) %in% hubs,
                                                  levels = c("FALSE","TRUE"))]
  seurat$percent.mt <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat$percent.ribo <- PercentageFeatureSet(seurat, pattern = c("^RPS","^RPL"))
  seurat$Dropout_rate <- apply(data, 2, function(z) mean(z==0)*100)
  seurat$Entropy <- entropy
  seurat$Hub_status="Normal"
  for (y in seq(ncol(seurat))) {
    seurat$Hub_status[y] <- ifelse(seurat$Hubs[y]=="hub",
                                   "Hub",
                                   ifelse(seurat$Antihubs[y]=="antihub",
                                          "Antihub",
                                          "Normal"))
    }
  seurat$Hub_status <- factor(seurat$Hub_status,
                              levels=c("Antihub","Normal","Hub"))
  seurat$Hubness_score <- score
  seurat <- NormalizeData(sseurat)
  seurat <- ScaleData(seurat)
  seurat <- FindVariableFeatures(seurat)
  seurat <- RunPCA(seurat, verbose=F)
  seurat <- RunUMAP(seurat, dims=1:n_dims)
  return(seurat)
}
