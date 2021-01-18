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
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned/"

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
score_norm <- function(scores) {
  tmp = scores
  for (col in grep("Met_",colnames(tmp), value = T)) {
    col_norm = paste0(col,"_norm")
    tmp = tmp %>% group_by(Dataset_id) %>% mutate(!!col_norm:=pnorm(normalize(get(col))))
  }
  return(tmp)
}
score_aggregation <- function(scores) {
  tmp = score_norm(scores)
  tmp = tmp %>% mutate(comp1=paste(Traj_type,Dataset_source,Method_id))
  for (col in grep("_norm",colnames(tmp), value = T)) {
    col_norm = paste0(col,"2")
    tmp = tmp %>%
      group_by(comp1) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp1)) %>% mutate(comp2=paste(Traj_type,Method_id))
  for (col in grep("_norm2",colnames(tmp), value = T)) {
    col_norm = gsub("2","3",col)
    tmp = tmp %>%
      group_by(comp2) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp2)) %>% group_by(Method_id)
  for (col in grep("_norm3",colnames(tmp), value = T)) {
    col_norm = gsub("3","4",col)
    tmp = tmp %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(Method_id))
  tmp$Overall_score = apply(data.frame(tmp), 1, function(x) {
    cols = grep("_norm4",colnames(tmp));
    vals = as.numeric(unname(unlist(x[cols])));
    mean(vals)})
  tmp = tmp[,c("Method_id",grep("_norm4",colnames(tmp),value=T),"Overall_score")]
  return(tmp)
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
traj <- sapply(dataset_choice, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$trajectory)
source <- sapply(dataset_choice, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$source)
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
                                                function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores_ktuned/stab_",
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

# ====
scores_df <- do.call(rbind,lapply(seq(scores_result), function(preprocess)
  do.call(rbind,lapply(seq(scores_result[[preprocess]]), function(dim)
    do.call(rbind,lapply(seq(scores_result[[preprocess]][[dim]]), function(set) {
      dataset = scores_result[[preprocess]][[dim]][[set]];
      dataset$Dataset_id = dataset_choice[set];
      dataset$Traj_type = traj[set];
      dataset$Dataset_source = "real/gold";
      dataset$Method_id = rownames(dataset);
      dataset$Dimension = names(scores_result[[preprocess]])[dim];
      dataset$Preprocess = names(scores_result)[preprocess];
      dataset$Met_corr = dataset$correlation;
      dataset$Met_F1 = dataset$F1_branches;
      dataset$Met_feat = dataset$featureimp_wcor;
      dataset = dataset[,-c(1:3)];
      return(dataset)}))))))
# Retain high ID only (from get_ID.py)
keep_id_datasets <- dataset_choice[c(2,3,6,8,14)]
scores_df <- scores_df[scores_df$Dataset_id %in% keep_id_datasets,]

# ====
# Figure TI paper 2nd try better
# Plot all datasets for the different methods (ie preprocess)
scores_normalized <- score_norm(scores_df)
scores_normalized$Metric <- sapply(scores_normalized$Preprocess, function(x)
  strsplit(x, " ")[[1]][2])
scores_normalized <- scores_normalized[!(scores_normalized$Metric=="cosine" &
                                         scores_normalized$Method_id=="dsl"),]
scores_normalized <- scores_normalized[,-grep("_feat",colnames(scores_normalized))]

scores_normalized$Overall_score <- apply(data.frame(scores_normalized),1,function(x)
  mean(as.numeric(unlist(x[grep("_norm", colnames(scores_normalized))]))))

scores_normalized$Dimension <- factor(scores_normalized$Dimension, levels=c(25,50,100,500))

ggplot(data.frame(scores_normalized), aes(x=Dimension, y=Overall_score, fill=Method_id)) + #3x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  #geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
  ylim(c(0,1)) +
  ylab("PAGA stability score") +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(~Metric) +
  #stat_compare_means(ref.group = "nothing", label = "p.signif") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))


# ====
# Figure TI paper 1st try
# ====
#scores_perdim_perpreprocess = lapply(setting_combination, function(preprocess) {
  liste = lapply(n_comps, function(dim)
    scores_df[scores_df$Dimension==dim &
                scores_df$Preprocess==preprocess,]);
  names(liste) = n_comps;
  return(liste)})
#names(scores_perdim_perpreprocess) = setting_combination
#scores_perdim_perpreprocess_avg <- lapply(scores_perdim_perpreprocess, function(x)
#  lapply(x, function(y) score_aggregation(y)))
## Stab for high ID eucli or cosine
## remove rows with cosine dsl
#scores_perdim_perpreprocess_avg <- lapply(seq(scores_perdim_perpreprocess_avg), function(x) {
#  if (strsplit(names(scores_perdim_perpreprocess_avg)[x]," ")[[1]][2]=="cosine") {
#    lapply(scores_perdim_perpreprocess_avg[[x]], function(y)
#      y[y$Method_id!="dsl",])
#  }
#  else {
#    scores_perdim_perpreprocess_avg[[x]]
#  }
#  })
#names(scores_perdim_perpreprocess_avg) = setting_combination
#
#scores_perdim_perpreprocess_avg <- do.call(rbind,lapply(seq(scores_perdim_perpreprocess_avg), function(x) {
# do.call(rbind,lapply(seq(scores_perdim_perpreprocess_avg[[x]]), function(y) {
#   data = scores_perdim_perpreprocess_avg[[x]][[y]]
#   data$Dimension =  n_comps[y]
#   data$Metric = strsplit(setting_combination[x]," ")[[1]][2]
#   data
# }))
#}))
#
#scores_perdim_perpreprocess_avg$Dimension <- factor(scores_perdim_perpreprocess_avg$Dimension)
#p2=ggplot(data.frame(scores_perdim_perpreprocess_avg), aes(x=Dimension, y=Overall_score, fill=Method_id)) + # 3x5
#  geom_boxplot(outlier.shape=NA, alpha=0.6) +
#  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),
#             aes(color=Method_id), fill="black") +
#  ylim(c(0,1)) +
#  ylab("PAGA stability score") +
#  scale_fill_brewer(palette="Set1") +
#  scale_color_brewer(palette="Set1") +
#  facet_wrap(~Metric) +
#  theme(panel.background = element_rect(fill="grey98"),
#        panel.grid = element_line(colour = "grey80"),
#        axis.text = element_text(size=15))
#
#ggarrange(p1,p2,nrow=2)

#stat.test = overall_ci %>%
#  group_by(Dimension) %>%
#  t_test(Score ~ Hub_reduction, ref.group = "nothing")
#stat.test.id = overall_ci_id %>%
#  group_by(Dimension) %>%
#  t_test(Score ~ Hub_reduction, ref.group = "nothing")



