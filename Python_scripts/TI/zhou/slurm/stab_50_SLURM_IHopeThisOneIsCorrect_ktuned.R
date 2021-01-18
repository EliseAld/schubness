library(SeuratDisk)
library(Seurat)
library(pbapply)
library(furrr)
library(ggpubr)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)
library(tictoc)
library(future.apply)


path <- "/shared/projects/hubness/TI/zhou/h5_jo_ktuned/"
path_res <- "/shared/projects/hubness/TI/zhou/scores_ktuned/stab_"

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
n_iter = 10
seed = 0
metric=c('cosine','euclidean')
n_comps=50
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

methods <- c("nothing","mp_normal","ls","ls_nicdm","dsl")
metric_choice <- c(#"him",
  "correlation","F1_branches","featureimp_wcor")
choice="zhou"
settings <- c('normduo_scaleTrue_ncomps50_cosine_louvain', 'normduo_scaleTrue_ncomps50_cosine_leiden',
              'normseurat_scaleTrue_ncomps50_cosine_louvain', 'normseurat_scaleTrue_ncomps50_cosine_leiden',
	      'normduo_scaleTrue_ncomps50_euclidean_louvain', 'normduo_scaleTrue_ncomps50_euclidean_leiden',
              'normseurat_scaleTrue_ncomps50_euclidean_louvain', 'normseurat_scaleTrue_ncomps50_euclidean_leiden')
# ====
source("/home/externe/curie/eamblard/dev/TI/dynverse_test_functions.R")
make_stab_score <- function(dataset_id) {
  print(paste0("Starting ",dataset_id))
  if (file.exists(paste0(path,dataset_id, "_",methods[1],".h5seurat"))) {
    clustering_algo <- ifelse(length(grep("leiden",dataset_id))>0,"leiden","louvain")
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id, "_",methods[1],".h5seurat"),
                                              verbose=F))
    iters <- lapply(methods,
                    function(x) read.csv(paste0(path,dataset_id,"_",x,"_stab.csv"), header=FALSE, row.names=1))
    cell_mask <- lapply(methods,
                        function(x) read.csv(paste0(path,dataset_id,"_",x,"_stab_cell.csv"), header=FALSE))
    feat_mask <- lapply(methods,
                        function(x) read.csv(paste0(path,dataset_id,"_",x,"_stab_feat.csv"), header=FALSE))
    mask <- list(cell_mask,feat_mask)
    mask <- lapply(mask,
                   function(msk1) {tmp=lapply(msk1,
                                              function(msk2) {mtp<-lapply(seq(n_iter-1),
                                                                          function(iter) {iter1=as.logical(gsub(F,NA,as.logical(msk2[iter,])));
                                                                          iter2=as.logical(gsub(F,NA,as.logical(msk2[(iter+1),])));
                                                                          iter1[intersect(which(is.na(iter2)),
                                                                                          which(iter1==T))]<-F;
                                                                          iter2[intersect(which(is.na(iter1)),
                                                                                          which(iter2==T))]<-F;
                                                                          return(list(iter1[!is.na(iter1)],
                                                                                      iter2[!is.na(iter2)]))});
                                              names(mtp)<-paste(seq(n_iter-1),seq(2,n_iter));
                                              return(mtp)});
                   names(tmp)=methods;
                   return(tmp)})
    names(mask) <- c("cell","feat")
    mask2 <- list(cell_mask,feat_mask)
    mask2 <- lapply(mask2,
                    function(msk1) {tmp=lapply(msk1,
                                               function(msk2) {mtp<-lapply(seq(n_iter-1),
                                                                           function(iter) {iter1=as.logical(msk2[iter,]);
                                                                           iter2=as.logical(msk2[(iter+1),]);
                                                                           return(iter1 & iter2)});
                                               names(mtp)<-paste(seq(n_iter-1),seq(2,n_iter));
                                               return(mtp)});
                    names(tmp)=methods;
                    return(tmp)})
    names(mask2) <- c("cell","feat")
    iters <- lapply(seq(length(methods)),
                    function(x) {pblapply(seq(n_iter),
                                          function(i) {tmp<-strsplit(iters[[x]][i,],"\t")[[1]];
                                          if (length(tmp)>0) {
                                            tmp<-unlist(strsplit(tmp,"\n"));
                                            values<-as.numeric(tmp[-grep(",",tmp)]);
                                            idx<-sapply(tmp[grep(",",tmp)], function(j) {splitted=strsplit(j,
                                                                                                           split="")[[1]];
                                            if (length(splitted)==8) {
                                              index=1+as.numeric(splitted[c(4,7)])
                                            }
                                            else {
                                              before1=4
                                              after1=grep(",",splitted)[1]-1
                                              before2=after1+3
                                              after2=length(splitted)-1
                                              if (before1==after1) {
                                                nb1=as.numeric(splitted[before1])+1
                                              }
                                              else {
                                                nb1=as.numeric(paste(splitted[before1:after1], collapse=""))+1
                                              }
                                              if (before2==after2) {
                                                nb2=as.numeric(splitted[before2])+1
                                              }
                                              else {
                                                nb2=as.numeric(paste(splitted[before2:after2], collapse=""))+1
                                              }
                                              index=c(nb1,nb2)
                                            }
                                            return(index)});
                                            M=matrix(0,
                                                     ncol=max(idx),
                                                     nrow=max(idx))
                                            for (z in seq(length(values))) {
                                              M[idx[1,z],idx[2,z]]=as.numeric(values[z])
                                            }
                                            M=as(M,"sparseMatrix")
                                          }
                                          else {
                                            M=as(matrix(0),"sparseMatrix")
                                          }
                                          return(M)
                                          })})
    seurat <- lapply(seq(length(methods)),
                     function(x) {seurat[[x]]@misc$paga <- NULL;
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
                                  function(x1) lapply(seq(n_iter),
                                                      function(x) setdiff(as.character(as.numeric(seurat[[x1]]@meta.data[,clustering_algo])),
                                                                          c(from[[x1]][[x]],to[[x1]][[x]]))))
    from <- lapply(seq(length(from)),
                   function(x1) lapply(seq(n_iter),
                                       function(x) c(from[[x1]][[x]],missing_cluster_ids[[x1]][[x]])))
    to <- lapply(seq(length(to)),
                 function(x1) lapply(seq(n_iter),
                                     function(x) c(to[[x1]][[x]],missing_cluster_ids[[x1]][[x]])))
    length <- lapply(seq(length(methods)),
                     function(pg1) lapply(seq(n_iter),
                                          function(pg) c(paga[[pg1]][[pg]]@x,rep(0,length(missing_cluster_ids[[pg1]][[pg]])))))
    directed <- lapply(length,
                       function(x1) lapply(x1,
                                           function(x) rep(T, length(x))))
    milestone_network <- lapply(seq(length(methods)),
                                function(x1) lapply(seq(n_iter),
                                                    function(x) tibble("from"=from[[x1]][[x]],
                                                                       "to"=to[[x1]][[x]],
                                                                       "length"=length[[x1]][[x]],
                                                                       "directed"=directed[[x1]][[x]])))
    names(milestone_network)<-methods
    dataset <- lapply(methods,
                      function(method) {tmp<-lapply(seq(n_iter-1),
                                                    function(iter) {wrap_expression(counts = t(seurat[[method]]@assays$RNA@counts[mask2$feat[[method]][[iter]],
                                                                                                                                  mask2$cell[[method]][[iter]]]),
                                                                                    expression = t(seurat[[method]]@assays$RNA@data[mask2$feat[[method]][[iter]],
                                                                                                                                    mask2$cell[[method]][[iter]]]))});
                      names(tmp) <- paste(seq(n_iter-1),seq(2,n_iter));
                      return(tmp)})
    names(dataset)<-methods
    traj <- lapply(methods,
                   function(method) {tmp<-lapply(seq(n_iter-1),
                                                 function(iter) lapply(0:1,
                                                                       function(comp_idx) {dataset[[method]][[paste(iter,iter+1)]] %>% 
                                                                           add_grouping(as.character(as.numeric(seurat[[method]]@meta.data[,clustering_algo][mask2$cell[[method]][[iter]]]))) %>% 
                                                                           add_cluster_graph2(milestone_network[[method]][[iter+comp_idx]])}));
                   names(tmp) <- paste(seq(n_iter-1),seq(2,n_iter));
                   return(tmp)})
    names(traj)<-methods
    traj <- lapply(methods,
                     function(method) lapply(seq(n_iter-1),
                                             function(iter) lapply(seq(2),
                                                                   function(comp_idx) traj[[method]][[iter]][[comp_idx]] %>% add_dimred(dyndimred::dimred_mds,
                                                                                                                                        expression_source = dataset[[method]][[paste(iter,iter+1)]]$expression))))
    
    names(traj)<-methods
    traj <- lapply(traj,
                   function(traj1) lapply(traj1,
                                          function(traj2) lapply(traj2,
                                                                 function(traj3) traj3 %>% add_cell_waypoints())))
    metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% metric_choice) %>% pull(metric_id)
    metrics <- lapply(metric_ids,
                        function(metric) lapply(methods,
                                                function(method) lapply(seq(n_iter-1),
                                                                        function(iter) {
                                                                          calculate_metrics2(traj[[method]][[iter]][[1]],
                                                                                             traj[[method]][[iter]][[2]],
                                                                                             metric,
                                                                                             dataset[[method]][[paste(iter,iter+1)]]$expression)[,metric]})))
    metrics <- lapply(metrics,
                      function(score) lapply(score,
                                             function(iter) mean(unlist(iter))))
    metrics <- data.frame(sapply(metrics,
                                 unlist))
    colnames(metrics) <- metric_ids
    rownames(metrics) <- methods
    return(metrics)
  }
  else {
  }
}
# ====

plan("multiprocess",workers=30)
future_lapply(settings,
         function(setting) {
           if (!file.exists(paste0(
             path_res,
             setting,
             ".rds"))) {
             dataset_choice <- datasets[grep(setting,datasets)];
             print(paste0("doing setting ",setting))
             tic("setting took:")
             scores_result <- furrr::future_map(dataset_choice, make_stab_score, .progress = T)
             toc()
             saveRDS(scores_result,
                     file=paste0(path_res,
                                 setting,
                                 ".rds"))
           }
           else {print(paste0(setting," done"))} })
