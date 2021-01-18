library(SeuratDisk)
library(Seurat)
library(pbapply)
library(ggpubr)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)
library(furrr)
library(future.apply)
plan("multiprocess",workers=50)
options(future.globals.maxSize=1048576000)

#https://zenodo.org/record/1443566
# ====
path <- "/shared/projects/hubness/TI/cytotrace/h5_jo/"

dataset_choice = sort(c('GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9',
                        'GSE60783_10', 'GSE67602_11', 'GSE70245_12', 'GSE52529_15',
                        'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33'))

do_log = T #already done with do_norm
do_pca = T
weighted=T
norm_scale = "True"
n_neighbors = 10
seed = 0
n_iter = 10
bootstrap_size = 0.95
metric='cosine'
n_comps=c(25,50,100,500)
do_norm = c('duo','seurat')
clustering_algo = c('louvain', 'leiden')

datasets <- c()
for (met in metric) {
  for (comp in n_comps) {
    for (norm in do_norm) {
      for (algo in clustering_algo) {
        for (set in dataset_choice) {
          datasets <- c(datasets, paste0(set,"_norm",norm,
                                        "_scale",norm_scale,
                                        "_ncomps",comp,
                                        "_",met,
                                        "_",algo))
        }
      }
    }
  }
}

source("/home/externe/curie/eamblard/dev/TI/dynverse_test_functions.R")
methods <- c("nothing","mp_normal","ls","ls_nicdm","dsl")

choice="cytotrace_id"
metric_choice <- c(#"him",
  "correlation","F1_branches","featureimp_wcor")
settings <- c('normduo_scaleTrue_ncomps25_cosine_louvain', 'normduo_scaleTrue_ncomps25_cosine_leiden', 
              'normseurat_scaleTrue_ncomps25_cosine_louvain', 'normseurat_scaleTrue_ncomps25_cosine_leiden', 
              'normduo_scaleTrue_ncomps50_cosine_louvain', 'normduo_scaleTrue_ncomps50_cosine_leiden',
              'normseurat_scaleTrue_ncomps50_cosine_louvain', 'normseurat_scaleTrue_ncomps50_cosine_leiden',
              'normduo_scaleTrue_ncomps100_cosine_louvain', 'normduo_scaleTrue_ncomps100_cosine_leiden',
              'normseurat_scaleTrue_ncomps100_cosine_louvain', 'normseurat_scaleTrue_ncomps100_cosine_leiden',
              'normduo_scaleTrue_ncomps500_cosine_louvain', 'normduo_scaleTrue_ncomps500_cosine_leiden',
              'normseurat_scaleTrue_ncomps500_cosine_louvain', 'normseurat_scaleTrue_ncomps500_cosine_leiden')
# ====
make_paga_traj <- function(dataset_id) {
  if (file.exists(paste0(path,dataset_id, "_",methods[1],".h5seurat"))) {
    clustering_algo <- ifelse(length(grep("leiden",dataset_id))>0,"leiden","louvain")
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"_", x,".h5seurat"),
                                              verbose=F))
    seurat <- lapply(seurat,
                     function(x) {x <- subset(x, cells=colnames(x)[!is.na(x@misc$Order)]);
                       x@misc$Order <- x@misc$Order[!is.na(x@misc$Order)];
                       return(x)})
    seurat <- lapply(seq(length(methods)),
                     function(x) {seurat[[x]]@misc$paga <- seurat[[x]]@misc$paga$connectivities_tree;
                     return(seurat[[x]])})
    paga <- lapply(seurat,
                   function(x) x@misc$paga)
    names(seurat) <- methods
    names(paga) <- methods
    n_col <- lapply(paga,
                    function(pg) sapply(seq(ncol(pg)),
                                        function(x) pg@p[x+1] - pg@p[x]))
    from <- lapply(n_col,
                   function(y) unlist(sapply(seq(y), function(x) rep(as.character(x),y[x]))))
    to <- lapply(paga,
                 function(pg) as.character(pg@i+1))
    missing_cluster_ids <- lapply(seq(length(seurat)),
                                  function(x) setdiff(as.character(as.numeric(seurat[[x]]@meta.data[,clustering_algo])),
                                                      c(from[[x]],to[[x]])))
    from <- lapply(seq(length(from)),
                   function(x) c(from[[x]],missing_cluster_ids[[x]]))
    to <- lapply(seq(length(to)),
                 function(x) c(to[[x]],missing_cluster_ids[[x]]))
    length <- lapply(seq(length(methods)),
                     function(pg) c(paga[[pg]]@x,rep(0,length(missing_cluster_ids[[pg]]))))
    directed <- lapply(length,
                       function(x) rep(T, length(x)))
    milestone_network <- lapply(seq(length(methods)),
                                function(x) tibble("from"=from[[x]],
                                                   "to"=to[[x]],
                                                   "length"=length[[x]],
                                                   "directed"=directed[[x]]))
    dataset <- wrap_expression(counts = t(seurat$nothing@assays$RNA@counts),
                               expression = t(seurat$nothing@assays$RNA@data))
    traj <- lapply(seq(length(methods)),
                   function(x) {dataset %>% add_grouping(as.character(as.numeric(seurat[[x]]@meta.data[,clustering_algo]))) %>% add_cluster_graph2(milestone_network[[x]])})
    traj <- lapply(seq(length(methods)),
                   function(x) traj[[x]] %>% add_dimred(dyndimred::dimred_mds,
                                                        expression_source = dataset$expression))
    return(list("seurat"=seurat,
                "trajectory"=traj,
                "milestone_network"=milestone_network,
                "dataset_id"=dataset_id))
  }
  else {
  }
}
score_paga_traj <- function(paga_traj_output) {
  dataset_id <- paga_traj_output$dataset_id
  dataset <- wrap_expression(counts = t(paga_traj_output$seurat$nothing@assays$RNA@counts),
                             expression = t(paga_traj_output$seurat$nothing@assays$RNA@data))
  pseudotime <- as.numeric(paga_traj_output$seurat$nothing@misc$Order)
  names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
                                                                              expression_source = dataset$expression) %>% add_cell_waypoints()
  model <- paga_traj_output$trajectory
  names(model) <- methods
  model <- lapply(model,
                  function(x) add_cell_waypoints(x))
  # take same metrics as in paper https://www-nature-com.proxy.insermbiblio.inist.fr/articles/s41587-019-0071-9#Sec9
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% metric_choice) %>% pull(metric_id)
  metrics <- lapply(model,
                    function(x) calculate_metrics2(GroundTruth,x,metric_ids))
  metrics <- lapply(metrics,
                    function(x) x[,metric_ids])
  names(metrics) <- methods
  return(metrics)
}
make_traj_and_score <- function(dataset_id) {
  if (file.exists(paste0(path,dataset_id, "_",methods[1],".h5seurat"))) {
    print(paste0("analyzing dataset ",dataset_id))
    paga_traj_output <- make_paga_traj(dataset_id)
    score <- score_paga_traj(paga_traj_output)
    score <- data.frame(do.call(rbind,score))
    rownames(score)<-methods
    return(score)
  }
}
# ====
#GT_labz <- lapply(dataset_choice, show_ground_truth_annotations)
future_map(settings,
         function(setting) {
           if (!file.exists(paste0(
             "/shared/projects/hubness/TI/cytotrace/scores/",
             setting,
             ".rds"))) {
             dataset_choice <- datasets[grep(setting,datasets)];
             print(paste0("doing setting ",setting))
             scores_result <- future_map(dataset_choice,
                                       make_traj_and_score);
             scores_result <- lapply(scores_result,
                                     function(x1) {idx = rep(T,length(x1));
                                     idx <- sapply(seq(length(x1)),
                                                   function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
                                     x1<-x1[idx]}); # remove empty spots
             saveRDS(scores_result,
                     file=paste0("/shared/projects/hubness/TI/cytotrace/scores/",
                                 setting,
                                 ".rds"))
           }
           else {print(paste0(setting," done"))} })
