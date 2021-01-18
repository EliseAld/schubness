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
# Set the different params
# ====
path_c <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned/"
path_z <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned/"


dataset_choice_c = sort(c('GSE90860_3', 'GSE67123_6', 'GSE98451_7', 'GSE94641_9',
                        'GSE60783_10', 'GSE67602_11', 'GSE70245_12', 'GSE52529_15',
                        'GSE52583_21', 'GSE52583_29', 'GSE69761_30', 'GSE64447_33'))
dataset_choice_z = sort(c('gold_hematopoiesis-gates_olsson', 'gold_germline-human-female-weeks_li',
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

norm_scale = "True"
metric=c('cosine','euclidean')
n_comps=c(25,50,100,500)
do_norm = c('duo','seurat')
clustering_algo = c('louvain', 'leiden')

datasets_c <- c()
for (met in metric) {
  for (comp in n_comps) {
    for (norm in do_norm) {
      for (algo in clustering_algo) {
        for (set in dataset_choice_c) {
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
datasets_z <- c()
for (met in metric) {
  for (comp in n_comps) {
    for (norm in do_norm) {
      for (algo in clustering_algo) {
        for (set in dataset_choice_z) {
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

choice="highid"
metric_choice <- c("correlation","F1_branches","featureimp_wcor")
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
# Functions
# ====
get_weighted_relative_diff <- function(unnorm_scores) {
  tmp=unnorm_scores
  for (col in grep("Met_",colnames(tmp), value = T)) {
    col_ref = paste0(col,"_ref")
    col_diff = paste0(col,"_diff")
    tmp[,col_ref]=sapply(seq(nrow(tmp)),function(x) {
      if (tmp$Method_id[x]=="nothing") {
        unlist(tmp[,col][x])}
      else {
        idx = max(which(tmp$Method_id[seq(x)]=="nothing"));
        unlist(tmp[,col][idx])
      }
    }
    )
    tmp=tmp%>%mutate(!!col_diff:=(get(col)-get(col_ref))*get(col_ref))
  }
  return(tmp)
}
score_norm <- function(scores) {
  tmp = scores %>% group_by(Dataset_id)
  for (col in grep("Met_",colnames(tmp), value = T)) {
    col_norm = paste0(col,"_norm")
    tmp = tmp %>% mutate(!!col_norm:=pnorm(normalize(get(col))))
  }
  return(tmp)
}
score_aggregation <- function(scores) {
  #tmp = score_norm(scores)
  tmp=scores
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
score_aggregation_w <- function(scores) {
  #tmp = score_norm(scores)
  tmp=scores
  tmp = tmp %>% mutate(comp1=paste(Traj_type,Dataset_source,Method_id))
  for (col in grep("_diff",colnames(tmp), value = T)) {
    col_norm = paste0(col,"2")
    tmp = tmp %>%
      group_by(comp1) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp1)) %>% mutate(comp2=paste(Traj_type,Method_id))
  for (col in grep("_diff2",colnames(tmp), value = T)) {
    col_norm = gsub("2","3",col)
    tmp = tmp %>%
      group_by(comp2) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp2)) %>% group_by(Method_id)
  for (col in grep("_diff3",colnames(tmp), value = T)) {
    col_norm = gsub("3","4",col)
    tmp = tmp %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(Method_id))
  tmp$Overall_score = apply(data.frame(tmp), 1, function(x) {
    cols = grep("_diff4",colnames(tmp));
    vals = as.numeric(unname(unlist(x[cols])));
    mean(vals)})
  tmp = tmp[,c("Method_id",grep("_diff4",colnames(tmp),value=T),"Overall_score")]
  return(tmp)
}

# ====
# Retrieve results
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


scores_result_c <- lapply(settings_across_comps,
                          function(x) {tmp=lapply(x,
                                             function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/cytotrace/scores_ktuned/",
                                                                             y,
                                                                             ".rds")));
                          names(tmp) <- n_comps;
                          return(tmp)})
scores_result_z <- lapply(settings_across_comps,
                          function(x) {tmp=lapply(x,
                                                  function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores_ktuned/",
                                                                                  y,
                                                                                  ".rds")));
                          names(tmp) <- n_comps;
                          return(tmp)})
keep_id_datasets <- c(2,3,6,8,14)
dataset_choice_z <- dataset_choice_z[keep_id_datasets]
scores_result_z <- lapply(scores_result_z, function(x)
  lapply(x, function(y) y[keep_id_datasets]))

# ====
# Get the traj and source info for both sets
# ====
#sapply(dataset_choice_c, function(x)
#  SeuratDisk::Convert(paste0(gsub("h5_jo_ktuned","raw",path_c),"raw_",x,".h5ad"),
#                      dest = "h5seurat", overwrite = TRUE, verbose = F))
#test=lapply(dataset_choice_c, function(x)
#  LoadH5Seurat(paste0(gsub("h5_jo_ktuned","raw",path_c),"raw_",x,".h5seurat"),
#               verbose=F)) # info not displayed, lets assume it's linear as well
traj_c <- rep("linear",length(dataset_choice_c))
source_c <- rep("real/gold",length(dataset_choice_c))
source_c[c(3,7,8)] <- "real/silver" # from the supp table
traj_z <- sapply(dataset_choice_z, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$trajectory) #all zhou are linear
source_z <- sapply(dataset_choice_z, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$source) #all zhou are gold
# ====
# Put all scores together
# ====
scores_result <- lapply(seq(setting_combination), function(x) {
  liste=lapply(seq(n_comps), function(y) {
    liste2=c(scores_result_c[[x]][[y]],scores_result_z[[x]][[y]]);
    names(liste2)=c(dataset_choice_c,dataset_choice_z);
    liste2})
  names(liste)=n_comps;
  liste})
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
traj<-unname(c(traj_c,traj_z))
source<-unname(c(source_c,source_z))
dataset_choice<-c(dataset_choice_c,dataset_choice_z)

# ====
# Get the averaged scores and everything
# ====
scores_df <- do.call(rbind,lapply(seq(scores_result), function(preprocess)
  do.call(rbind,lapply(seq(scores_result[[preprocess]]), function(dim)
    do.call(rbind,lapply(seq(scores_result[[preprocess]][[dim]]), function(set) {
      dataset = scores_result[[preprocess]][[dim]][[set]];
      dataset$Dataset_id = dataset_choice[set];
      dataset$Traj_type = traj[set];
      dataset$Dataset_source = source[set];
      dataset$Method_id = rownames(dataset);
      dataset$Dimension = names(scores_result[[preprocess]])[dim];
      dataset$Preprocess = names(scores_result)[preprocess];
      dataset$Metric = sapply(names(scores_result)[preprocess], function(x) strsplit(x," ")[[1]][2]);
      dataset$Met_corr = dataset$correlation;
      dataset$Met_F1 = dataset$F1_branches;
      dataset$Met_feat = dataset$featureimp_wcor;
      dataset = dataset[,-c(1:3)];
      return(dataset)}))))))
# Remove cosine dsl
scores_df <- scores_df[!(scores_df$Metric=="cosine" &
                           scores_df$Method_id=="dsl"),]

# ====
# Figure TI paper 2nd try better
# Plot all datasets for the different methods (ie preprocess)
# =====
scores_normalized <- score_norm(scores_df)

scores_normalized$Overall_score <- apply(data.frame(scores_normalized),1,function(x)
  mean(as.numeric(unlist(x[grep("_norm", colnames(scores_normalized))]))))

scores_normalized$Dimension <- factor(scores_normalized$Dimension, levels=c(25,50,100,500))

ggplot(data.frame(scores_normalized), aes(x=Dimension, y=Overall_score, fill=Method_id)) + #3x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  #geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
  ylim(c(0,1)) +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(~Metric) +
  #stat_compare_means(ref.group = "nothing", label = "p.signif") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))

# Make aggregation
scores_perdim_perpreprocess = lapply(setting_combination, function(preprocess) {
  liste = lapply(n_comps, function(dim)
    scores_normalized[scores_normalized$Dimension==dim &
                        scores_normalized$Preprocess==preprocess,]);
  names(liste) = n_comps;
  return(liste)})
names(scores_perdim_perpreprocess) = setting_combination
scores_perdim_perpreprocess_avg <- lapply(scores_perdim_perpreprocess, function(x)
  lapply(x, function(y) score_aggregation(y)))

# ====
# Figure TI paper A & B
# ====
df <- do.call(rbind,scores_perdim_perpreprocess_avg$`seurat euclidean leiden`)[,-c(1)]
df <- data.frame("Score"=unlist(df),
                 "Metric"=c(rep(c(metric_choice,"Overall_score"),each=nrow(df))),
                 "Dimension"=factor(rep(rep(n_comps,each=length(methods)),ncol(df)),levels = n_comps),
                 "Hub_red"=rep(methods,ncol(df)*length(n_comps)))
df$Diff <- df$Score
for (idx in seq(nrow(df))) {
  df$Diff[idx]<-df$Diff[idx]-df$Score[5*(floor((idx-1)/5))+1]}
df$X <- gsub(25,1,
             gsub(50,2,
                  gsub(100,3,
                       gsub(500,4,df$Dimension))))
df$Y <- gsub("dsl",1,
             gsub("ls",2,
                  gsub("ls_nicdm",3,
                       gsub("mp_normal",4,
                            gsub("nothing",5,df$Hub_red)))))
ggplot(df, aes(x=Dimension, y=Hub_red)) + #6x7.5
  facet_wrap(~Metric) +
  geom_point(aes(size=Score, color=Score)) +
  geom_text(aes(label=round(Score,3)), size=2) +
  scale_size(range = c(0,15)) +
  scale_color_viridis_c(option="inferno") +
  scale_y_discrete(labels=c("DSL","LS_nicdm","LS","MP","Nothing")) +
  theme(text = element_text(size=12),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
# Diff score
ggplot(df[df$Hub_red!="nothing",], aes(x=X, y=Y)) +
  facet_wrap(~Metric) +
  geom_rect(aes(ymin=as.numeric(Y), ymax=as.numeric(Y)+(Diff/0.55*1), xmin=as.numeric(X)-0.4, xmax=as.numeric(X)+0.4, fill = Diff)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label=round(Diff,3), y=as.numeric(Y)+(Diff/0.55*1)/2), size=3) +
  scale_y_discrete(labels=c("DSL","LS_nicdm","LS","MP")) +
  theme(text = element_text(size=12),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15)) +
  ylab("Hub_reduction") +
  scale_x_discrete(labels=c("1"="25",
                            "2"="50",
                            "3"="100",
                            "4"="500")) +
  xlab("Dimension")


# ====
# Weighted diff (maybe not..!)
# ====
weighted_diff_score <- get_weighted_relative_diff(scores_df)
# Make aggregation
weighted_perdim_perpreprocess = lapply(setting_combination, function(preprocess) {
  liste = lapply(n_comps, function(dim)
    weighted_diff_score[weighted_diff_score$Dimension==dim &
                          weighted_diff_score$Preprocess==preprocess,]);
  names(liste) = n_comps;
  return(liste)})
names(weighted_perdim_perpreprocess) = setting_combination

weighted_perdim_perpreprocess_avg <- lapply(weighted_perdim_perpreprocess, function(x)
  lapply(x, function(y) score_aggregation_w(y)))

df <- do.call(rbind,weighted_perdim_perpreprocess_avg$`seurat euclidean leiden`)[,-c(1)]
df <- data.frame("Diff"=as.numeric(unname(unlist(df))),
                 "Metric"=c(rep(c(metric_choice,"Overall_score"),each=nrow(df))),
                 "Dimension"=factor(rep(rep(n_comps,each=length(methods)),ncol(df)),levels = n_comps),
                 "Hub_red"=rep(methods,ncol(df)*length(n_comps)))
df$X <- gsub(25,1,
             gsub(50,2,
                  gsub(100,3,
                       gsub(500,4,df$Dimension))))
df$Y <- gsub("dsl",1,
             gsub("ls",2,
                  gsub("ls_nicdm",3,
                       gsub("mp_normal",4,
                            gsub("nothing",5,df$Hub_red)))))
ggplot(df[df$Hub_red!="nothing",], aes(x=X, y=Y)) +
  facet_wrap(~Metric) +
  geom_rect(aes(ymin=as.numeric(Y), ymax=as.numeric(Y)+(Diff/0.55*1), xmin=as.numeric(X)-0.4, xmax=as.numeric(X)+0.4, fill = Diff)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label=round(Diff,3), y=as.numeric(Y)+(Diff/0.55*1)/2), size=3) +
  scale_y_discrete(labels=c("DSL","LS_nicdm","LS","MP")) +
  theme(text = element_text(size=12),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15)) +
  ylab("Hub_reduction") +
  scale_x_discrete(labels=c("1"="25",
                            "2"="50",
                            "3"="100",
                            "4"="500")) +
  xlab("Dimension")
