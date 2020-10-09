library(SeuratDisk)
library(Seurat)
library(dynwrap)
library(dyno)
library(dyneval)
library(tidyverse)
library(assertthat)
library(ggpubr)
library(pbapply)
library(BBmisc)
library(Matrix)
source("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/dynverse_test_functions.R")
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/"
methods <- c("nothing","mp_normal","ls","ls_nicdm","dsl")
dataset_choice = c('GSE36552','GSE45719','GSE52529','GSE52583','GSE52583b','GSE59114','GSE60783','GSE64447',
                  'GSE67123','GSE69761','GSE74767b','GSE75330b','GSE75748','GSE76408','GSE85066', # 1st lot
                  'GSE86146','GSE87375','GSE87375b','GSE90047','GSE90860','GSE92332','GSE92332b','GSE93421',
                  'GSE94641','GSE95753b','GSE98451','GSE98664','GSE102066','GSE106587','GSE107122','GSE107910'
                  #'GSE67602','GSE70245', 'GSE75330','GSE74767','GSE95753','GSE97391','GSE97391b','GSE99933','GSE103633','GSE109774','GSE109774b'   # excluded because  not all cells are ordered
                  )
n_dim_choice <- c(50,100)
metric_choice <- c(#"him",
  "correlation","F1_branches","featureimp_wcor")
n_iter=10
#pblapply(dataset_choice,
#         function(dts) lapply(n_dim_choice,
#                              function(dim) lapply(methods,
#                                                   function(mtd) {if (file.exists(paste0(path,dts, "/",dim,"_dims_",mtd[1],".h5ad"))) {
#                                                                     Convert(paste0(path,dts, "/",dim,"_dims_",mtd,".h5ad"),
#                                                                             dest = "h5seurat", overwrite = TRUE, verbose = F)}})))
#sapply(dataset_choice,function(x) {seurat <- LoadH5Seurat(paste0(path,x,"/50_dims_nothing.h5seurat"),verbose=F);
#tmp<-ifelse(ncol(seurat)==length(seurat@misc$Order),"",paste0("remove",x));
#tmp<-tmp[tmp!=""];
#print(unname(tmp))}
#)
make_paga_traj <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id, "/",n_dim,"_dims_",methods[1],".h5ad"))) {
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"/",n_dim,"_dims_", x,".h5seurat"),
                                              verbose=F))
    iters <- lapply(methods,
                    function(x) read.csv(paste0(path,dataset_id,"/",n_dim,"_dims_", x,".csv"), header=FALSE, row.names=1))
    iters <- lapply(seq(length(methods)),
                    function(x) {lapply(seq(n_iter),
                                        function(i) {tmp<-strsplit(iters[[x]][i,],"\t")[[1]];
                                        tmp<-unlist(strsplit(tmp,"\n"));
                                        values<-as.numeric(tmp[-grep(",",tmp)]);
                                        idx<-sapply(tmp[grep(",",tmp)], function(j) 1+as.numeric(strsplit(j,
                                                                                                        split="")[[1]][c(4,7)]));
                                        M=matrix(0,
                                                 ncol=ncol(seurat[[x]]@misc$paga$connectivities_tree),
                                                 nrow=nrow(seurat[[x]]@misc$paga$connectivities_tree))
                                        for (z in seq(length(values))) {
                                          M[idx[1,z],idx[2,z]]=as.numeric(values[z])
                                        }
                                        M=as(M,"sparseMatrix")
                    })})
    seurat <- lapply(seq(length(methods)),
                     function(x) {seurat[[x]]@misc$paga0 <- seurat[[x]]@misc$paga$connectivities_tree;
                     seurat[[x]]@misc$paga <- NULL
                     for (i in seq(n_iter)) {
                       seurat[[x]]@misc[[paste0("paga",i)]] <- iters[[x]][[i]]
                     }
                     return(seurat[[x]])})
    paga <- lapply(seurat,
                   function(x) x@misc[grep("paga",names(x@misc))])
    names(seurat) <- methods
    names(paga) <- methods
    n_col <- lapply(paga,
                    function(pg1) lapply(pg1,
                                        function(pg) sapply(seq(ncol(pg)),
                                                            function(x) pg@p[x+1] - pg@p[x])))
    from <- lapply(n_col,
                   function(y1) lapply(y1,
                                      function(y) unlist(sapply(seq(y), function(x) rep(as.character(x),y[x])))))
    to <- lapply(paga,
                 function(pg1) lapply(pg1,
                                      function(pg) as.character(pg@i+1)))
    missing_cluster_ids <- lapply(seq(length(seurat)),
                                  function(x1) lapply(seq(n_iter+1),
                                                     function(x) setdiff(as.character(as.numeric(seurat[[x1]]$leiden)),
                                                                         c(from[[x1]][[x]],to[[x1]][[x]]))))
    from <- lapply(seq(length(from)),
                   function(x1) lapply(seq(n_iter+1),
                                       function(x) c(from[[x1]][[x]],missing_cluster_ids[[x1]][[x]])))
    to <- lapply(seq(length(to)),
                 function(x1) lapply(seq(n_iter+1),
                                     function(x) c(to[[x1]][[x]],missing_cluster_ids[[x1]][[x]])))
    length <- lapply(seq(length(methods)),
                     function(pg1) lapply(seq(n_iter+1),
                                          function(pg) c(paga[[pg1]][[pg]]@x,rep(0,length(missing_cluster_ids[[pg1]][[pg]])))))
    directed <- lapply(length,
                       function(x1) lapply(x1,
                                           function(x) rep(T, length(x))))
    milestone_network <- lapply(seq(length(methods)),
                                function(x1) lapply(seq(n_iter+1),
                                                    function(x) tibble("from"=from[[x1]][[x]],
                                                                       "to"=to[[x1]][[x]],
                                                                       "length"=length[[x1]][[x]],
                                                                       "directed"=directed[[x1]][[x]])))
    dataset <- wrap_expression(counts = t(seurat$nothing@assays$RNA@counts),
                               expression = t(seurat$nothing@assays$RNA@data))
    traj <- lapply(seq(length(methods)),
                   function(x1) lapply(seq(n_iter+1),
                                       function(x) {dataset %>% add_grouping(as.character(as.numeric(seurat[[x1]]$leiden))) %>% add_cluster_graph2(milestone_network[[x1]][[x]])}))
    traj <- lapply(seq(length(methods)),
                   function(x1) lapply(seq(n_iter+1),
                                       function(x) traj[[x1]][[x]] %>% add_dimred(dyndimred::dimred_mds,
                                                                      expression_source = dataset$expression)))
    return(list("seurat"=seurat,
                "trajectory"=traj,
                "milestone_network"=milestone_network,
                "dataset_id"=dataset_id,
                "n_dim"=n_dim))
  }
  else {
  }
}
plot_paga_traj <- function(paga_traj_output) {
  idx_to_plot=seq(length(methods))
  if (paga_traj_output$dataset_id == 'GSE45719' & paga_traj_output$n_dim==50) {
    idx_to_plot=c(1,2,3,5)
  }
  pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/",
                  paga_traj_output$dataset_id, "/",
                  paga_traj_output$n_dim,"_dims_paga_graph.pdf" ))
  ggarrange(plotlist=lapply(idx_to_plot,
                              function(x) plot_dimred(paga_traj_output$trajectory[[x]], 
                                                      expression_source = t(paga_traj_output$seurat[[x]]@assays$RNA@data))),
              labels=methods)
  dev.off()
}
plot_compare_paga_traj <- function(paga_traj_output) {
  dataset <- wrap_expression(counts = t(paga_traj_output$seurat$nothing@assays$RNA@counts),
                             expression = t(paga_traj_output$seurat$nothing@assays$RNA@data))
  pseudotime <- paga_traj_output$seurat$nothing@misc$Order
  names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
                                                                                   expression_source = dataset$expression) %>% add_cell_waypoints()
  model <- paga_traj_output$trajectory
  names(model) <- methods
  model <- map(model, add_cell_waypoints)
  # take same metrics as in paper https://www-nature-com.proxy.insermbiblio.inist.fr/articles/s41587-019-0071-9#Sec9
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% c(#"him",
                                                                                                                      "correlation","F1_branches","featureimp_wcor")) %>% pull(metric_id)
  metrics <- lapply(model,
                    function(x) calculate_metrics2(GroundTruth,x,metric_ids))
  metrics_df <- data.frame("metric_value"=c(sapply(metrics,function(x) return(t(x[,metric_ids])))),
                        "metric_id"=metric_ids,
                        "method_id"=rep(methods,each=length(metric_ids)))
  metrics_df <- bind_rows(metrics_df, data.frame("metric_value"=colMeans(sapply(metrics,function(x) return(t(x[,metric_ids])))),
                                                 "metric_id"="Mean",
                                                 "method_id"=methods))
  metrics_df$metric_id <- factor(metrics_df$metric_id, levels = c(metric_ids,"Mean"))
  #pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/",
  #                paga_traj_output$dataset_id, "/",
  #                paga_traj_output$n_dim,"_dims_paga_heatmap_comp.pdf" ))
  ggplot(metrics_df, aes(metric_id, method_id, fill = metric_value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 90))
  #dev.off()
}
score_paga_traj <- function(paga_traj_output) {
  dataset <- wrap_expression(counts = t(paga_traj_output$seurat$nothing@assays$RNA@counts),
                             expression = t(paga_traj_output$seurat$nothing@assays$RNA@data))
  pseudotime <- paga_traj_output$seurat$nothing@misc$Order
  names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
                                                                              expression_source = dataset$expression) %>% add_cell_waypoints()
  model <- paga_traj_output$trajectory
  names(model) <- methods
  model <- lapply(model,
                  function(x) map(x, add_cell_waypoints))
  # take same metrics as in paper https://www-nature-com.proxy.insermbiblio.inist.fr/articles/s41587-019-0071-9#Sec9
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% metric_choice) %>% pull(metric_id)
  metrics <- lapply(model,
                    function(x1) lapply(x1,
                                        function(x) calculate_metrics2(GroundTruth,x,metric_ids)))
  metrics <- lapply(metrics,
                    function(x1) lapply(x1,
                                        function(x) x[,metric_ids]))
  metrics <- lapply(metrics,
                    function(x) colMeans(do.call(rbind,x)))
  return(metrics)
}
compare_byScore <- function(scores) {
  # first normalize scores across methods for each dataset
  score_norm <- lapply(scores,
                       function(x) lapply(x,
                                          function(y) normalize(as.numeric(y))))
  # apply unit proba density function
  score_norm <- lapply(score_norm,
                       function(x) lapply(x,
                                          function(y) dnorm(y)))
  # arithmetic mean over datasets
  avg_score <- sapply(methods,
                      function(x) {tmp<-sapply(score_norm,
                                         function(y) y[[x]]);
                      rowMeans(matrix(as.numeric(tmp), nrow=length(metric_choice)))})
  rownames(avg_score) <- metric_choice
  data.frame(t(avg_score))
}
compare_allScore <- function(scores) {
  avg_score <- compare_byScore(scores)
  # geometric mean
  data.frame("OverallScore"=apply(avg_score,
        1,
        function(x) exp(mean(log(x[x!=0])))))
}
make_traj_and_score <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id, "/",n_dim,"_dims_",methods[1],".h5ad"))) {
    paga_traj_output <- make_paga_traj(dataset_id, n_dim)
    score <- score_paga_traj(paga_traj_output)
    return(score)
  }
}
cluster_nb <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id, "/",n_dim,"_dims_",methods[1],".h5ad"))) {
   seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"/",n_dim,"_dims_", x,".h5seurat"),
                                              verbose=F))
   return("cluster_nb"=sapply(seurat,function(x) length(unique(x$leiden))))
  }
  else {
  }
}


#paga_traj_output <- pblapply(n_dim_choice,
#                       function(dims) pblapply(dataset_choice[1:15],
#                                            function(set) make_paga_traj(set,dims)))
#paga_traj_output <- lapply(paga_traj_output,
#                           function(x1) {idx = rep(T,length(x1));
#                           idx <- sapply(seq(length(x1)),
#                                    function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
#                           x1<-x1[idx]}) # remove empty spots
#scores_result <- pblapply(paga_traj_output,
#                          function(x) pblapply(x, score_paga_traj))
cluster_nb <- pblapply(n_dim_choice,
                       function(dims) pblapply(dataset_choice,
                                               function(set) cluster_nb(set,dims)))
scores_result <- pblapply(n_dim_choice,
                          function(dims) pblapply(dataset_choice,
                                                  function(set) make_traj_and_score(set,dims)))
scores_result <- lapply(scores_result,
                        function(x1) {idx = rep(T,length(x1));
                        idx <- sapply(seq(length(x1)),
                                      function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
                        x1<-x1[idx]}) # remove empty spots
scores <- lapply(scores_result,
                 function(x) {tmp<-cbind(compare_byScore(x),compare_allScore(x));
                 data.frame("Score"=unlist(tmp),
                            "Metric"=rep(colnames(tmp), each= nrow(tmp)),
                            "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)))})
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/bubbleplot_score_AllDatasets_",
                n_iter,"iter_unsup.pdf"))
ggarrange(plotlist=lapply(scores,
                          function(x) ggplot(x, aes(Metric, Hub_reduction)) +
                            geom_point(aes(size=Score, color = Score)) +
                            geom_text(aes(label=round(Score,3))) +
                            scale_size(range = c(0,20)) +
                            scale_color_gradient2(low="navyblue",
                                                  mid="mediumvioletred",
                                                  high="yellow",
                                                  midpoint=0.3) +
                            theme(axis.text.x = element_text(angle = 90))),
          labels=paste0(n_dim_choice," PCs"),
          common.legend = T)
dev.off()


