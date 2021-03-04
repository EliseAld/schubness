#' @param path path to .h5ad python PAGA output
#' @param dataset_path path to .h5seurat R PAGA object
#' @param metric_choice name of quality metrics used (default to c("correlation","F1_branches","featureimp_wcor"))
#' @param clustering clustering algorithm used to perform PAGA
#' @param truth_path path to the .rds object if available

#' @return quality metrics

#' @example scores <- make_traj_and_score(dataset_path, clustering, metric_choice, truth_path)

### Load librairies
library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)
library(Xmisc)

### Functions to compute the three quality scores (Corr, F1, Wcor)
# Convert the python PAGA output .h5ad to a R-readable format
convert_PtoR <- function(path) {
  if (file.exists(paste0(path))) {
    SeuratDisk::Convert(path, dest = "h5seurat", overwrite = TRUE, verbose = F)
    #file.remove(paste0(path))
  }
  else {
    warning('No file')
  }
}

# Dynverse trajectory object
prepare_dynverse_traj <- function(dataset_path, clustering, metric_choice, truth_path) {
  if (file.exists(dataset_path)) {
    if (Xmisc::endsWith(truth_path, ".rds")) {
      truth = readRDS(file=truth_path)
    }
    seurat <- LoadH5Seurat(dataset_path, verbose=F)
    paga <- seurat@misc$paga$connectivities_tree
    n_col <- sapply(seq(ncol(paga)),
                    function(x) paga@p[x+1] - paga@p[x])
    from <- unlist(sapply(seq(n_col), function(x) rep(as.character(x), n_col[x])))
    to <- as.character(paga@i+1)
    missing_cluster_ids <- setdiff(as.character(as.numeric(seurat@meta.data[, clustering])),
                                   c(from,to))
    from <- c(from,missing_cluster_ids)
    to <- c(to,missing_cluster_ids)
    length <- c(paga@x,rep(0,length(missing_cluster_ids)))
    directed <- rep(T, length(length))
    milestone_network <- tibble::tibble("from"=from,
                                        "to"=to,
                                        "length"=length,
                                        "directed"=directed)
    dataset <- wrap_expression(counts = t(seurat@assays$RNA@counts),
                               expression = t(seurat@assays$RNA@data))
    traj <- dataset %>%
      add_grouping(as.character(as.numeric(seurat@meta.data[, clustering]))) %>%
      add_cluster_graph2(milestone_network) %>%
      add_dimred(dyndimred::dimred_mds,
                 expression_source = dataset$expression)
    if (Xmisc::endsWith(truth_path, ".rds")) {
      traj_truth <- dataset %>% add_trajectory(milestone_ids = truth$milestone_ids,
                                             milestone_network = truth$milestone_network,
                                             milestone_percentages = truth$milestone_percentages,
                                             divergence_regions = truth$divergence_regions)
      return(list("seurat"=seurat,
                  "trajectory"=traj,
                  "milestone_network"=milestone_network,
                  "ground_truth"=traj_truth,
                  "dataset"=dataset))
    }
    else {
      return(list("seurat"=seurat,
                  "trajectory"=traj,
                  "milestone_network"=milestone_network,
                  "ground_truth"=NA,
                  "dataset"=dataset))
    }
  }
  else {
    stop('There is no file at the given path')
  }
}

# Compute the scores
score_traj <- function(dynverse_traj) {
  dataset <- dynverse_traj$dataset
  seurat <- dynverse_traj$seurat
  if (is.na(dynverse_traj$ground_truth)) {
    if ("misc" %in% slotNames(seurat) & "Order" %in% colnames(seurat@misc)) {
      pseudotime <- as.numeric(seurat@misc$Order)
      names(pseudotime) <- colnames(dynverse_traj$seurat)
      GroundTruth <- dataset %>%
        add_linear_trajectory(pseudotime) %>%
        add_dimred(dyndimred::dimred_mds,
                 expression_source = dataset$expression) %>%
        add_cell_waypoints()
    }
    else {
      stop('No Ground Truth')
    }
  }
  else {
    GroundTruth <- dynverse_traj$ground_truth %>%
      add_cell_waypoints()
  }
  model <- dynverse_traj$trajectory%>%
    add_cell_waypoints()
  metric_ids <- dyneval::metrics %>%
    filter(category != "average") %>%
    filter(perfect==1) %>%
    filter(metric_id %in% metric_choice)%>%
    pull(metric_id)
  metrics <- calculate_metrics2(GroundTruth, model, metric_ids)
  metrics <- metrics[, metric_ids]
  colnames(metrics) <- metric_ids
  return(metrics)
}

# Single function
make_traj_and_score <- function(dataset_path, clustering, metric_choice, truth_path) {
  #print(paste0("analyzing dataset ",dataset_path))
  traj_output <- prepare_dynverse_traj(dataset_path, clustering, metric_choice, truth_path)
  score <- data.frame(score_traj(traj_output))
  return(score)
}