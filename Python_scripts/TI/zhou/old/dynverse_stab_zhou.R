library(SeuratDisk)
library(Seurat)
library(pbapply)
library(parallel)
library(ggpubr)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)
library(rstatix)

#https://zenodo.org/record/1443566
# ====
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo/"

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
n_neighbors = 10
seed = 0
n_iter = 10
metric=c('cosine','euclidean')
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

source("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/dynverse_test_functions.R")
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
              'normseurat_scaleTrue_ncomps500_cosine_louvain', 'normseurat_scaleTrue_ncomps500_cosine_leiden',
              'normduo_scaleTrue_ncomps25_euclidean_louvain', 'normduo_scaleTrue_ncomps25_euclidean_leiden',
              'normseurat_scaleTrue_ncomps25_euclidean_louvain', 'normseurat_scaleTrue_ncomps25_euclidean_leiden',
              'normduo_scaleTrue_ncomps50_euclidean_louvain', 'normduo_scaleTrue_ncomps50_euclidean_leiden',
              'normseurat_scaleTrue_ncomps50_euclidean_louvain', 'normseurat_scaleTrue_ncomps50_euclidean_leiden',
              'normduo_scaleTrue_ncomps100_euclidean_louvain', 'normduo_scaleTrue_ncomps100_euclidean_leiden',
              'normseurat_scaleTrue_ncomps100_euclidean_louvain', 'normseurat_scaleTrue_ncomps100_euclidean_leiden',
              'normduo_scaleTrue_ncomps500_euclidean_louvain', 'normduo_scaleTrue_ncomps500_euclidean_leiden',
              'normseurat_scaleTrue_ncomps500_euclidean_louvain', 'normseurat_scaleTrue_ncomps500_euclidean_leiden')
# ====
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
compare_byScore <- function(scores) {
  # first normalize scores across methods for each dataset + apply unit proba density function
  score_norm <- lapply(scores,
                       function(x) {tmp<-BBmisc::normalize(x);
                       tmp<-apply(tmp,2,dnorm);
                       rownames(tmp)<-methods;
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
# ====

#GT_labz <- lapply(dataset_choice, show_ground_truth_annotations)
pblapply(settings,
         function(setting) {
           if (!file.exists(paste0(
             "/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores/stab_",
             setting,
             ".rds"))) {
             dataset_choice <- datasets[grep(setting,datasets)];
             print(paste0("doing setting ",setting))
             scores_result <- pblapply(dataset_choice, make_stab_score)
             saveRDS(scores_result,
                     file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores/stab_",
                                 setting,
                                 ".rds"))
           }
           else {print(paste0(setting," done"))} })
# ====
setting_combination <- c(outer(X=c(outer(X=do_norm,
                                         Y=metric,
                                         FUN="paste")),
                               Y=clustering_algo,
                               FUN="paste"))
settings_across_comps <- lapply(setting_combination,
                                function(combi) {
                                  setting_combi <- strsplit(combi," ")[[1]];
                                  sapply(n_comps,
                                         function(comp) paste0("norm",setting_combi[1],
                                                               "_scale",norm_scale,
                                                               "_ncomps",comp,
                                                               "_",setting_combi[2],
                                                               "_",setting_combi[3]))
                                })


scores_result <- lapply(settings_across_comps,
                        function(x) {tmp=lapply(x,
                                                function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores/",
                                                                                y,
                                                                                ".rds")));
                        names(tmp) <- n_comps;
                        return(tmp)})
names(scores_result) <- setting_combination
scores_result <- lapply(scores_result,
                        function(x) x[lengths(x)!=0])
scores_result <- scores_result[lengths(scores_result)!=0]
scores_result <- lapply(scores_result,
                        function(x1) lapply(x1,
                                            function(x2) {idx = rep(T,length(x2));
                                            idx <- sapply(seq(length(x2)),
                                                          function(x3) idx[x3]=ifelse(is.null(x2[[x3]]),F,T));
                                            x2<-x2[idx]}))
scores <- lapply(scores_result,
                 function(x1) lapply(x1,
                                     function(y) {tmp<-cbind(compare_byScore(y),compare_allScore(y));
                                     data.frame("Score"=unlist(tmp),
                                                "Metric"=rep(colnames(tmp), each=nrow(tmp)),
                                                "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)))}))
scores2 <- lapply(scores_result,
                  function(x1) {tmp <- lapply(seq(length(x1)),
                                              function(idx) {tmp<-cbind(compare_byScore(x1[[idx]]),compare_allScore(x1[[idx]]));
                                              data.frame("Score"=unlist(tmp),
                                                         "Metric"=rep(colnames(tmp), each=nrow(tmp)),
                                                         "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)),
                                                         "Dimension"=factor(n_comps, levels=n_comps)[idx])});
                  do.call(rbind,tmp)})
scores2 <- lapply(scores2,
                  function(x) {tmp<-x;
                  tmp$X<-gsub(25,1,
                              gsub(50,2,
                                   gsub(100,3,
                                        gsub(500,4,tmp$Dimension))));
                  tmp$Y<-gsub("dsl",1,
                              gsub("ls",2,
                                   gsub("ls_nicdm",3,
                                        gsub("mp_normal",4,
                                             gsub("nothing",5,tmp$Hub_reduction)))));
                  tmp$diff_score<-tmp$Score;
                  for (idx in seq(nrow(tmp))) {
                    tmp$diff_score[idx]<-tmp$diff_score[idx]-tmp$Score[5*(floor((idx-1)/5))+1]};
                  return(tmp)})

# ====
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_stab_sup.pdf"),
    width=18,height=15) #sup refers to the nb of clusters chosen or not
lapply(seq(length(scores)),
       function(x1) {fig=ggarrange(plotlist=lapply(scores[[x1]],
                                                   function(x) ggplot(x, aes(Metric, Hub_reduction)) +
                                                     geom_point(aes(size=Score, color = Score)) +
                                                     geom_text(aes(label=round(Score,3)), size=2) +
                                                     scale_size(range = c(0,20)) +
                                                     scale_color_gradient2(low="navyblue",
                                                                           mid="mediumvioletred",
                                                                           high="yellow",
                                                                           midpoint=0.3) +
                                                     theme(axis.text.x = element_text(angle = 90),
                                                           text = element_text(size=12),
                                                           axis.ticks.y = element_blank(),
                                                           axis.ticks.x = element_blank())),
                                   labels=paste0(n_comps," PCs"),
                                   common.legend=T);
       annotate_figure(fig, top=paste(names(scores)[x1],collapse="_"))}
)
dev.off()
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_stab_sup2.pdf"),
    width=18,height=15) #sup refers to the nb of clusters chosen or not
lapply(seq(length(scores2)),
       function(x) {fig=ggplot(scores2[[x]], aes(x=Dimension, y=Hub_reduction)) +
         facet_wrap(~Metric) +
         geom_point(aes(size=Score, color = Score)) +
         geom_text(aes(label=round(Score,3)), size=2) +
         scale_size(range = c(0,15)) +
         scale_color_viridis_c(option="magma") +
         theme(axis.text.x = element_text(angle = 90),
               text = element_text(size=12),
               axis.ticks.y = element_blank(),
               axis.ticks.x = element_blank());
       annotate_figure(fig, top=names(scores2)[x])}
)
dev.off()
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_stab_sup3.pdf"),
    width=18,height=15) #sup refers to the nb of clusters chosen or not
lapply(seq(length(scores2)),
       function(x) {fig=ggplot(scores2[[x]], aes(x=X, y=Y)) +
         facet_wrap(~Metric) +
         geom_rect(aes(ymin=as.numeric(Y), ymax=as.numeric(Y)+(diff_score/0.13*1), xmin=as.numeric(X)-0.4, xmax=as.numeric(X)+0.4, fill = diff_score)) +
         scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
         theme(text = element_text(size=12),
               axis.ticks.y = element_blank()) +
         scale_y_discrete(labels=c("1"=methods[5],
                                   "2"=methods[4],
                                   "3"=methods[3],
                                   "4"=methods[2],
                                   "5"=methods[1])) +
         ylab("Hub_reduction") +
         scale_x_discrete(labels=c("1"="25",
                                   "2"="50",
                                   "3"="100",
                                   "4"="500")) +
         xlab("Dimension");
       annotate_figure(fig, top=names(scores2)[x])}
)
dev.off()
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_stab_sup4.pdf"),
    width=18,height=15) #sup refers to the nb of clusters chosen or not
lapply(seq(length(scores2)),
       function(x) {fig=ggplot(scores2[[x]], aes(x=X, y=Y)) +
         facet_wrap(~Metric) +
         geom_point(aes(color = diff_score, size=abs(diff_score), )) +
         scale_size(range=c(5,20)) +
         scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
         theme(text = element_text(size=12),
               axis.ticks.y = element_blank()) +
         scale_y_discrete(labels=c("1"=methods[5],
                                   "2"=methods[4],
                                   "3"=methods[3],
                                   "4"=methods[2],
                                   "5"=methods[1])) +
         ylab("Hub_reduction") +
         scale_x_discrete(labels=c("1"="25",
                                   "2"="50",
                                   "3"="100",
                                   "4"="500")) +
         xlab("Dimension");
       annotate_figure(fig, top=names(scores2)[x])}
)
dev.off()
# ====

# Make a figure with overallscore and all preprocess merged together
tmp <- lapply(seq(length(scores2)),
              function(idx) {scores2[[idx]]$Preprocess <- names(scores2)[idx];
              return(scores2[[idx]][scores2[[idx]]$Metric=="OverallScore",])})
overallscore_single <- data.frame("Score"=sapply(seq(nrow(tmp[[1]])),
                                                 function(cond) mean(sapply(tmp,
                                                                            function(x) x$Score[cond]))),
                                  "Metric"="OverallScore",
                                  'Hub_reduction'=tmp[[1]]$Hub_reduction,
                                  "Dimension"=tmp[[1]]$Dimension,
                                  "X"=tmp[[1]]$X)
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/overallscore_single_stab.pdf"),
    width=18,height=15)
ggplot(overallscore_single, aes(x=Dimension, y=Score, fill=Hub_reduction)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(c(0,0.4)) +
  scale_fill_viridis_d() +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
dev.off()
# ====

# Make seurat leiden euclidean with high ID datasets (from get_ID.py)
keep_id_datasets <- c(2,3,6,8,14)
scores2_id <- lapply(scores_result,
                     function(x1) {tmp <- lapply(seq(length(x1)),
                                                 function(idx) {tmp<-cbind(compare_byScore(x1[[idx]][keep_id_datasets]),compare_allScore(x1[[idx]][keep_id_datasets]));
                                                 data.frame("Score"=unlist(tmp),
                                                            "Metric"=rep(colnames(tmp), each=nrow(tmp)),
                                                            "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)),
                                                            "Dimension"=factor(n_comps, levels=n_comps)[idx])});
                     do.call(rbind,tmp)})
scores2_id <- lapply(scores2_id,
                     function(x) {tmp<-x;
                     tmp$X<-gsub(25,1,
                                 gsub(50,2,
                                      gsub(100,3,
                                           gsub(500,4,tmp$Dimension))));
                     tmp$Y<-gsub("dsl",1,
                                 gsub("ls",2,
                                      gsub("ls_nicdm",3,
                                           gsub("mp_normal",4,
                                                gsub("nothing",5,tmp$Hub_reduction)))));
                     tmp$diff_score<-tmp$Score;
                     for (idx in seq(nrow(tmp))) {
                       tmp$diff_score[idx]<-tmp$diff_score[idx]-tmp$Score[5*(floor((idx-1)/5))+1]};
                     return(tmp)})
lapply(seq(length(scores2_id)),
       function(x) {fig=ggplot(scores2_id[[x]], aes(x=X, y=Y)) +
         facet_wrap(~Metric) +
         geom_point(aes(color = diff_score, size=abs(diff_score), )) +
         scale_size(range=c(5,20)) +
         scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
         theme(text = element_text(size=12),
               axis.ticks.y = element_blank()) +
         scale_y_discrete(labels=c("1"=methods[5],
                                   "2"=methods[4],
                                   "3"=methods[3],
                                   "4"=methods[2],
                                   "5"=methods[1])) +
         ylab("Hub_reduction") +
         scale_x_discrete(labels=c("1"="25",
                                   "2"="50",
                                   "3"="100",
                                   "4"="500")) +
         xlab("Dimension");
       annotate_figure(fig, top=names(scores2_id)[x])}
)
tmp_id <- lapply(seq(length(scores2_id)),
              function(idx) {scores2_id[[idx]]$Preprocess <- names(scores2_id)[idx];
              return(scores2_id[[idx]][scores2_id[[idx]]$Metric=="OverallScore",])})
overallscore_single_id <- data.frame("Score"=sapply(seq(nrow(tmp_id[[1]])),
                                                    function(cond) mean(sapply(tmp_id,
                                                                               function(x) x$Score[cond]))),
                                     "Metric"="OverallScore",
                                     'Hub_reduction'=tmp_id[[1]]$Hub_reduction,
                                     "Dimension"=tmp_id[[1]]$Dimension,
                                     "X"=tmp_id[[1]]$X)
# ====
# Figure paper TI
# Stab and stab high ID
overall_ci <- do.call(rbind,tmp)
overall_ci$Metric <- unname(sapply(overall_ci$Preprocess, function(x)
  strsplit(x, " ")[[1]][2]))
overall_ci_id <- do.call(rbind,tmp_id)
overall_ci_id$Metric <- unname(sapply(overall_ci_id$Preprocess, function(x)
  strsplit(x, " ")[[1]][2]))
# remove rows with cosine dsl
overall_ci <- overall_ci[!(overall_ci$Metric=="cosine" & 
                             overall_ci$Hub_reduction=="dsl"),]
overall_ci_id <- overall_ci_id[!(overall_ci_id$Metric=="cosine" & 
                                   overall_ci_id$Hub_reduction=="dsl"),]

ggplot(overall_ci, aes(x=Dimension, y=Score, fill=Hub_reduction)) + # 3x5
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),
             aes(color=Hub_reduction, shape=Metric), fill="black") +
  ylim(c(0.2,0.35)) +
  ylab("PAGA stability score") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
ggplot(overall_ci_id, aes(x=Dimension, y=Score, fill=Hub_reduction)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),
             aes(color=Hub_reduction, shape=Metric), fill="black") +
  ylim(c(0.2,0.35)) +
  ylab("PAGA stability score") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))


stat.test = overall_ci %>%
  group_by(Dimension) %>%
  t_test(Score ~ Hub_reduction, ref.group = "nothing")
stat.test.id = overall_ci_id %>%
  group_by(Dimension) %>%
  t_test(Score ~ Hub_reduction, ref.group = "nothing")



