library(SeuratDisk)
library(Seurat)
library(pbapply)
library(ggpubr)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)
library(future.apply)
library(furrr)

#https://zenodo.org/record/1443566
# ====
path <- "/shared/projects/hubness/TI/zhou/h5_jo_ktuned/"

dataset_choice = sort(c('gold_hematopoiesis-gates_olsson', 'gold_germline-human-female-weeks_li',
                        'gold_stimulated-dendritic-cells-LPS_shalek',
                        'gold_germline-human-female_guo', 'gold_mESC-differentiation_hayashi',
                        'gold_developing-dendritic-cells_schlitzer',
                        'gold_germline-human-male-weeks_li',
                        'gold_pancreatic-beta-cell-maturation_zhang',
                        'gold_human-embryos_petropoulos', 
                        'gold_germline-human-male_guo', 'gold_myoblast-differentiation_trapnell',
                        'gold_pancreatic-alpha-cell-maturation_zhang',
                        'gold_aging-hsc-young_kowalczyk','gold_aging-hsc-old_kowalczyk'
))

do_log = T #already done with do_norm
do_pca = T
weighted=T
norm_scale = "True"
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

choice="zhou"
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
settings <- c(settings,gsub('cosine','euclidean',settings))
# ====
make_paga_traj <- function(dataset_id) {
  if (file.exists(paste0(path,dataset_id, "_",methods[1],".h5seurat"))) {
    rds_truth = readRDS(file=paste0(gsub("h5_jo_ktuned","rds2",path),
                                    strsplit(dataset_id, "_norm")[[1]][1],
                                    ".rds"))
    clustering_algo <- ifelse(length(grep("leiden",dataset_id))>0,"leiden","louvain")
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"_", x,".h5seurat"),
                                              verbose=F))
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
    dataset <- wrap_expression(counts = rds_truth$counts,
                               expression = rds_truth$expression)
    traj <- lapply(seq(length(methods)),
                   function(x) {dataset %>% add_grouping(as.character(as.numeric(seurat[[x]]@meta.data[,clustering_algo]))) %>% add_cluster_graph2(milestone_network[[x]])})
    traj <- lapply(seq(length(methods)),
                   function(x) traj[[x]] %>% add_dimred(dyndimred::dimred_mds,
                                                        expression_source = dataset$expression))
    traj_truth <- dataset %>% add_trajectory(milestone_ids = rds_truth$milestone_ids,
                                             milestone_network = rds_truth$milestone_network,
                                             milestone_percentages = rds_truth$milestone_percentages,
                                             divergence_regions = rds_truth$divergence_regions )
    return(list("seurat"=seurat,
                "trajectory"=traj,
                "milestone_network"=milestone_network,
                "dataset_id"=dataset_id,
                "ground_truth"=traj_truth,
                "dataset"=dataset))
  }
  else {
  }
}
score_paga_traj <- function(paga_traj_output) {
  dataset_id <- paga_traj_output$dataset_id
  dataset <- paga_traj_output$dataset
  #pseudotime <- paga_traj_output$seurat$nothing@misc$Order
  #if (length(grep("kowalczyk",dataset_id))>0) {
  #  pseudotime <- as.numeric(gsub("LT-HSC",1,
  #                     gsub("ST-HSC",2,
  #                          gsub("MPP",3,pseudotime))))
  #}
  #if (length(grep("schlitzer",dataset_id))>0) {
  #  pseudotime <- as.numeric(gsub("MDP",1,
  #                     gsub("CDP",2,
  #                          gsub("PreDC",3,pseudotime))))
  #}
  #if (length(grep("guo",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                       function(x) substr(x,2,nchar(x)-1)))
  #}
  #if (length(grep("weeks_li",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                                  function(x) substr(x,3,nchar(x)-4)))
  #}
  #if (length(grep("olsson",dataset_id))>0) {
  #  pseudotime <- as.numeric(gsub("Lsk",1,
  #                     gsub("Cmp",2,
  #                          gsub("Gmp",3,pseudotime))))
  #}
  #if (length(grep("petropoulos",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                                  function(x) strsplit(x," ")[[1]][3]))
  #}
  #if (length(grep("hayashi",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                                  function(x) substr(x,1,2)))
  #}
  #if (length(grep("trapnell",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                                  function(x) substr(x,2,nchar(x))))
  #}
  #if (length(grep("pancreatic-alpha",dataset_id))>0) {
  #  pseudotime <- sapply(pseudotime,
  #                       function(x) strsplit(x," ")[[1]][2])
  #  pseudotime <- as.numeric(gsub("E17.5",0,
  #                     gsub("P0",1,
  #                          gsub("P9",10,
  #                               gsub("P15",16,
  #                                    gsub("P18",19,
  #                                         gsub("P60",61,pseudotime)))))))
  #}
  #if (length(grep("pancreatic-beta",dataset_id))>0) {
  #  pseudotime <- sapply(pseudotime,
  #                       function(x) strsplit(x," ")[[1]][2])
  #  pseudotime <- as.numeric(gsub("E17.5",0,
  #                                gsub("P0",1,
  #                                     gsub("P3",4,
  #                                          gsub("P9",10,
  #                                               gsub("P15",16,
  #                                                    gsub("P18",19,
  #                                                         gsub("P60",61,pseudotime))))))))
  #}
  #if (length(grep("shalek",dataset_id))>0) {
  #  pseudotime <- as.numeric(sapply(pseudotime,
  #                                  function(x) substr(x,nchar(x)-1,nchar(x)-1)))
  #}
  #names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  #GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
  #                                                                            expression_source = dataset$expression) %>% add_cell_waypoints()
  GroundTruth <- paga_traj_output$ground_truth %>% add_cell_waypoints()
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
pblapply(settings,
         function(setting) {
           if (!file.exists(paste0(
             "/shared/projects/hubness/TI/zhou/scores_ktuned/",
             setting,
             ".rds"))) {
             dataset_choice <- datasets[grep(setting,datasets)];
             print(paste0("doing setting ",setting))
             scores_result <- pblapply(dataset_choice,
                                       make_traj_and_score);
             scores_result <- lapply(scores_result,
                                     function(x1) {idx = rep(T,length(x1));
                                     idx <- sapply(seq(length(x1)),
                                                   function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
                                     x1<-x1[idx]}); # remove empty spots
             saveRDS(scores_result,
                     file=paste0("/shared/projects/hubness/TI/zhou/scores_ktuned/",
                                 setting,
                                 ".rds"))
           }
           else {print(paste0(setting," done"))} })
