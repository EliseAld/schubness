library(SeuratDisk)
library(Seurat)
library(pbapply)
library(ggpubr)
library(tidyverse)
library(dyno)
library(assertthat)
library(dyneval)
library(BBmisc)

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
n_iter = 10-1
bootstrap_size = 0.95
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
pblapply(datasets,
         function(dts) {print(dts);
           lapply(methods,
                  function(mtd) {if (file.exists(paste0(path,dts,"_",mtd,".h5ad"))) {
                    Convert(paste0(path,dts,"_",mtd,".h5ad"),
                            dest = "h5seurat", overwrite = TRUE, verbose = F);
                    file.remove(paste0(path,dts,"_",mtd,".h5ad"))}})})
# ====
make_paga_traj <- function(dataset_id) {
  if (file.exists(paste0(path,dataset_id, "_",methods[1],".h5seurat"))) {
    rds_truth = readRDS(file=paste0(gsub("h5_jo","rds2",path),
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
compare_byScore <- function(scores) {
  # first normalize scores across methods for each dataset + apply unit proba density function
  score_norm <- lapply(scores,
                       function(x) {tmp<-normalize(x);
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
show_ground_truth_annotations <- function(dataset_id) {
  if (file.exists(paste0(path,"csv/",dataset_id, "_order.csv")))
    return(unique(read_csv(file = paste0(path,"csv/",dataset_id, "_order.csv"),
                           )[,2]))
}
# ====

#GT_labz <- lapply(dataset_choice, show_ground_truth_annotations)
pblapply(settings,
         function(setting) {
           if (!file.exists(paste0(
             "/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores/",
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
                     file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores/",
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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_sup.pdf"),
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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_sup2.pdf"),
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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_sup3.pdf"),
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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_sup4.pdf"),
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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/overallscore_single.pdf"),
    width=18,height=15)
ggplot(overallscore_single, aes(x=Dimension, y=Score, fill=Hub_reduction)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(c(0,0.4)) +
  scale_fill_viridis_d()
dev.off()

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
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/bubbleplot_score_",choice,"_sup_id4.pdf"),
    width=18,height=15) #sup refers to the nb of clusters chosen or not
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
dev.off()
tmp <- lapply(seq(length(scores2_id)),
              function(idx) {scores2_id[[idx]]$Preprocess <- names(scores2_id)[idx];
              return(scores2_id[[idx]][scores2_id[[idx]]$Metric=="OverallScore",])})
overallscore_single_id <- data.frame("Score"=sapply(seq(nrow(tmp[[1]])),
                                                 function(cond) mean(sapply(tmp,
                                                                            function(x) x$Score[cond]))),
                                  "Metric"="OverallScore",
                                  'Hub_reduction'=tmp[[1]]$Hub_reduction,
                                  "Dimension"=tmp[[1]]$Dimension,
                                  "X"=tmp[[1]]$X)
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/overallscore_single_id.pdf"),
    width=18,height=15)
ggplot(overallscore_single_id, aes(x=Dimension, y=Score, fill=Hub_reduction)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(c(0,0.4)) +
  scale_fill_viridis_d()
dev.off()

# Reproduce Jo's figure with Over all score and merge all preprocess to choose the best hub reduction
scores_jo <- lapply(scores_result,
                    function(preprocess) {datafr <- do.call(rbind,lapply(preprocess,
                                                function(dim) {df <- do.call(rbind,lapply(dim,
                                                                      function(score) data.frame("OverallScore"=apply(score,1,
                                                                                                                      function(row) exp(mean(log(row[row!=0])))),
                                                                                                 "Hub_reduction"=rownames(score))));
                                                df$Dataset <- rep(dataset_choice, each=length(methods));
                                                return(df)}));
                    datafr$Dimension <- rep(n_comps, each=length(methods)*length(dataset_choice));
                    return(datafr)})
length_nrow <- nrow(scores_jo$`duo cosine louvain`)
names_preprocess <- names(scores_jo)
scores_jo <- do.call(rbind,scores_jo)
scores_jo$Preprocess <- rep(names_preprocess, each=length_nrow)
scores_jo$Preprocess2 <- paste(scores_jo$Preprocess,scores_jo$Dimension)
scores_jo$Dimension <- factor(scores_jo$Dimension, levels=n_comps)
ggplot(scores_jo[scores_jo$Preprocess==names_preprocess[1],], aes(y=OverallScore, x=Dimension, fill=Hub_reduction)) +
  geom_boxplot()
