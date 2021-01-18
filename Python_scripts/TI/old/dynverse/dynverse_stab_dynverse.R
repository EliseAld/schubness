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
library(FedData)
#https://zenodo.org/record/1443566
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/dynverse/"
source("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/dynverse_test_functions.R")
methods <- c("nothing","mp_normal","ls","ls_nicdm","dsl")
datasets_choice = unname(sapply(c('dyntoy_linear_2.rds', 'gold_psc-astrocyte-maturation-glia_sloan.rds',
                                   'gold_germline-human-female-weeks_li.rds', 'dyngen_1.rds',
                                   'dyngen_10.rds', 'gold_mESC-differentiation_hayashi.rds',
                                   'prosstt_bifurcating_1.rds', 'dyngen_60.rds', 'prosstt_tree_1.rds',
                                   'prosstt_tree_3.rds', 'dyntoy_cyclic_2.rds', 'prosstt_linear_3.rds',
                                   'gold_human-embryos_petropoulos.rds', 'prosstt_binary_tree_1.rds',
                                   'dyntoy_tree_3.rds', 'dyngen_50.rds', 'splatter_binary_tree_2.rds',
                                   'dyngen_22.rds', 'splatter_multifurcating_3.rds',
                                   'dyntoy_bifurcating_3.rds', 'dyntoy_diverging_converging_2.rds',
                                   'splatter_multifurcating_5.rds', 'gold_myoblast-differentiation_trapnell.rds',
                                   'splatter_tree_2.rds', 'splatter_linear_1.rds'),
                                function(x) strsplit(x,".rds")[[1]]))
choice="selection"
dataset_choice <- unname(sapply(sort(datasets_choice),
                                function(x) paste0(strsplit(x,"_")[[1]][1],"/h5/",x)))
n_dim_choice <- c(25,50,100,500)
metric_choice <- c(#"him",
  "correlation","F1_branches","featureimp_wcor")

make_stab_score <- function(dataset_id, n_dim, n_iter=10) {
  print(paste0("Starting ",dataset_id," with ",n_dim," dims"))
  if (file.exists(paste0(path,dataset_id, "_",n_dim,"dims_",methods[1],"2.h5seurat"))) {
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"_",n_dim,"dims_", x,"2.h5seurat"),
                                              verbose=F))
    iters <- lapply(methods,
                    function(x) read.csv(paste0(path,dataset_id,"_",n_dim,"dims_", x,"_stab.csv"), header=FALSE, row.names=1))
    cell_mask <- lapply(methods,
                    function(x) read.csv(paste0(path,dataset_id,"_",n_dim,"dims_", x,"_stab_cell.csv"), header=FALSE))
    feat_mask <- lapply(methods,
                    function(x) read.csv(paste0(path,dataset_id,"_",n_dim,"dims_", x,"_stab_feat.csv"), header=FALSE))
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
                                                      function(x) setdiff(as.character(as.numeric(seurat[[x1]]$leiden)),
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
                                           add_grouping(as.character(as.numeric(seurat[[method]]$leiden[mask2$cell[[method]][[iter]]]))) %>% 
                                           add_cluster_graph2(milestone_network[[method]][[iter+comp_idx]])}));
                   names(tmp) <- paste(seq(n_iter-1),seq(2,n_iter));
                   return(tmp)})
    names(traj)<-methods
    traj <- pblapply(methods,
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
    metrics <- pblapply(metric_ids,
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
compare_byScore <- function(scores) {
  # first normalize scores across methods for each dataset + apply unit proba density function
  score_norm <- lapply(scores,
                       function(x) {tmp<-normalize(x);
                       tmp<-apply(tmp,2,dnorm);
                       return(tmp)})
  # arithmetic mean over datasets
  avg_score <- sapply(methods,
                      function(x) {tmp<-sapply(score_norm,
                                               function(y) y[x,]);
                      rowMeans(tmp)})
  data.frame(t(avg_score))
}
compare_allScore <- function(scores) {
  avg_score <- compare_byScore(scores)
  # geometric mean
  data.frame("OverallScore"=apply(avg_score,
                                  1,
                                  function(x) exp(mean(log(x[x!=0])))))
}


scores25 <- pblapply(25,
                          function(dims) {tmp<-pblapply(dataset_choice,
                                                        function(set) make_stab_score(set,dims));
                          names(tmp) <- dataset_choice[1:length(tmp)];
                          return(tmp)})
scores50 <- pblapply(50,
                     function(dims) {tmp<-pblapply(dataset_choice,
                                                   function(set) make_stab_score(set,dims));
                     names(tmp) <- dataset_choice[1:length(tmp)];
                     return(tmp)})
scores100 <- pblapply(100,
                     function(dims) {tmp<-pblapply(dataset_choice,
                                                   function(set) make_stab_score(set,dims));
                     names(tmp) <- dataset_choice[1:length(tmp)];
                     return(tmp)})
scores500 <- pblapply(500,
                     function(dims) {tmp<-pblapply(dataset_choice,
                                                   function(set) make_stab_score(set,dims));
                     names(tmp) <- dataset_choice[1:length(tmp)];
                     return(tmp)})
names(scores_result) <- n_dim_choice
scores_result <- lapply(scores_result,
                        function(x1) {idx = rep(T,length(x1));
                        idx <- sapply(seq(length(x1)),
                                      function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
                        x1<-x1[idx]}) # remove empty spots
scores_result <- scores_result[lengths(scores_result)!=0]
scores <- lapply(scores_result,
                 function(x) {tmp<-cbind(compare_byScore(x),compare_allScore(x));
                 data.frame("Score"=unlist(tmp),
                            "Metric"=rep(colnames(tmp), each= nrow(tmp)),
                            "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)))})
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/bubbleplot_score_",choice,"_sup.pdf"),
    width=12,height=8)#sup refers to the nb of clusters chosen or not
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
saveRDS(scores_result,
        file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/score_",choice,"_sup.rds"))

