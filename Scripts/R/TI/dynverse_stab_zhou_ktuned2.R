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
# Set the different params
# ====
path_c <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned4/"
path_z <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned4/"


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
          datasets_c <- c(datasets_c, paste0(set,"_norm",norm,
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
          datasets_z <- c(datasets_z, paste0(set,"_norm",norm,
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
methods <- c("nothing","umap","gauss","mp_normal","ls","ls_nicdm","dsl")

choice="highid"
metric_choice <- c("correlation","F1_branches","featureimp_wcor"#,"F1_milestones"
)
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
score_norm <- function(scores) {
  tmp = scores %>% group_by(Dataset_id)
  for (col in grep("Met_",colnames(tmp), value = T)) {
    col_norm = paste0(col,"_normed")
    tmp = tmp %>% mutate(!!col_norm:=pnorm(normalize(get(col))))
  }
  return(tmp)
}
score_aggregation <- function(scores,to_grep) {
  #tmp = score_norm(scores)
  tmp=scores
  tmp = tmp %>% mutate(comp1=paste(Traj_type,Dataset_source,Method_id))
  for (col in grep(to_grep,colnames(tmp), value = T)) {
    col_norm = paste0(col,"_norm2")
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
# Get the traj and source info for both sets
# ====

traj <- sapply(dataset_choice, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$trajectory)
source <- sapply(dataset_choice, function(x)
  readRDS(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/rds2/",
                 x,
                 ".rds"))$source)
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


resc <- lapply(settings_across_comps,
               function(x) {tmp=lapply(x,
                                       function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/cytotrace/scores_ktuned4/stab_",
                                                                       y,
                                                                       "2.rds")));
               names(tmp) <- n_comps;
               return(tmp)})
resz <- lapply(settings_across_comps,
               function(x) {tmp=lapply(x,
                                       function(y) readRDS(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/zhou/scores_ktuned4/stab_",
                                                                       y,
                                                                       "2.rds")));
               names(tmp) <- n_comps;
               return(tmp)})

# ====
# Put all scores together
# ====
scores_result <- lapply(seq(setting_combination), function(x) {
  liste=lapply(seq(n_comps), function(y) {
    liste2=c(resc[[x]][[y]],resz[[x]][[y]]);
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
#traj=traj_z
source<-unname(c(source_c,source_z))
#source=source_z
dataset_choice<-c(dataset_choice_c,dataset_choice_z)
#dataset_choice=dataset_choice_z
rm(resc,resz,traj_c,traj_z,source_c,source_z)

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
      dataset$Met_corr = dataset$correlation;
      dataset$Met_F1 = dataset$F1_branches;
      #dataset$Met_feat = dataset$featureimp_wcor;
      #dataset$Met_F1m = dataset$F1_milestones;
      dataset$Metric = sapply(names(scores_result)[preprocess], function(x) strsplit(x," ")[[1]][2]);
      dataset = dataset[,-c(1:4)];
      return(dataset)}))))))
# Retain high ID only (from get_ID.py)
keep_id_datasets <- dataset_choice[c(2,3,6,8,14)]
# Remove cosine dsl
scores_df <- scores_df[!(scores_df$Metric=="cosine" &
                           scores_df$Method_id=="dsl"),]

# ====
# Figure TI paper 2nd try better
# Plot all datasets for the different methods (ie preprocess)
scores_aggreg <- scores_df[scores_df$Preprocess %in% c("seurat euclidean leiden",
                                                       "seurat cosine leiden"),]
scores_aggreg1 <- scores_df[scores_df$Preprocess %in% c("duo euclidean louvain",
                                                       "duo cosine louvain"),]
scores_aggreg2 <- scores_df[scores_df$Preprocess %in% c("duo euclidean leiden",
                                                       "duo cosine leiden"),]
scores_aggreg3 <- scores_df[scores_df$Preprocess %in% c("seurat euclidean louvain",
                                                       "seurat cosine louvain"),]
scores_aggreg4 <- scores_df[scores_df$Preprocess %in% c("seurat euclidean leiden",
                                                       "seurat cosine leiden"),]
#scores_aggreg <- score_norm(scores_aggreg)
scores_aggreg$Overall_score <- apply(data.frame(scores_aggreg),1,function(x)
  mean(as.numeric(unlist(x[grep("Met_", colnames(scores_aggreg))]))))
scores_aggreg1$Overall_score <- apply(data.frame(scores_aggreg1),1,function(x)
  mean(as.numeric(unlist(x[grep("Met_", colnames(scores_aggreg1))]))))
scores_aggreg2$Overall_score <- apply(data.frame(scores_aggreg2),1,function(x)
  mean(as.numeric(unlist(x[grep("Met_", colnames(scores_aggreg2))]))))
scores_aggreg3$Overall_score <- apply(data.frame(scores_aggreg3),1,function(x)
  mean(as.numeric(unlist(x[grep("Met_", colnames(scores_aggreg3))]))))
scores_aggreg4$Overall_score <- apply(data.frame(scores_aggreg4),1,function(x)
  mean(as.numeric(unlist(x[grep("Met_", colnames(scores_aggreg4))]))))

scores_aggreg$Dimension <- factor(scores_aggreg$Dimension, levels=n_comps)
scores_aggreg1$Dimension <- factor(scores_aggreg1$Dimension, levels=rev(n_comps))
scores_aggreg2$Dimension <- factor(scores_aggreg2$Dimension, levels=rev(n_comps))
scores_aggreg3$Dimension <- factor(scores_aggreg3$Dimension, levels=rev(n_comps))
scores_aggreg4$Dimension <- factor(scores_aggreg4$Dimension, levels=rev(n_comps))

scores_aggreg$Method_id <- factor(scores_aggreg$Method_id, levels=rev(methods))
scores_aggreg1$Method_id <- factor(scores_aggreg1$Method_id, levels=rev(methods))
scores_aggreg2$Method_id <- factor(scores_aggreg2$Method_id, levels=rev(methods))
scores_aggreg3$Method_id <- factor(scores_aggreg3$Method_id, levels=rev(methods))
scores_aggreg4$Method_id <- factor(scores_aggreg4$Method_id, levels=rev(methods))

scores_aggreg <- scores_aggreg[scores_aggreg$Dataset_id %in% keep_id_datasets,]
scores_aggreg1b <- scores_aggreg1[!(scores_aggreg1$Dataset_id %in% keep_id_datasets),]
scores_aggreg2b <- scores_aggreg2[!(scores_aggreg2$Dataset_id %in% keep_id_datasets),]
scores_aggreg3b <- scores_aggreg3[!(scores_aggreg3$Dataset_id %in% keep_id_datasets),]
scores_aggreg4b <- scores_aggreg4[!(scores_aggreg4$Dataset_id %in% keep_id_datasets),]
scores_aggreg1 <- scores_aggreg1[scores_aggreg1$Dataset_id %in% keep_id_datasets,]
scores_aggreg2 <- scores_aggreg2[scores_aggreg2$Dataset_id %in% keep_id_datasets,]
scores_aggreg3 <- scores_aggreg3[scores_aggreg3$Dataset_id %in% keep_id_datasets,]
scores_aggreg4 <- scores_aggreg4[scores_aggreg4$Dataset_id %in% keep_id_datasets,]

ggplot(scores_aggreg, aes(x=Dimension, y=Met_corr, fill=Method_id)) + #3x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
  ylim(c(0,1)) +
  scale_colour_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                               rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                               rgb(50,116,161,max=255))) +
  scale_fill_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                             rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                             rgb(50,116,161,max=255))) +
  facet_wrap(~Metric) +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
ggsave(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure2_TI/panelD1.pdf"),
       width=8,height=3)
ggplot(scores_aggreg, aes(x=Dimension, y=Met_F1, fill=Method_id)) + #3x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
  ylim(c(0,1)) +
  scale_colour_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                               rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                               rgb(50,116,161,max=255))) +
  scale_fill_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                             rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                             rgb(50,116,161,max=255))) +
  facet_wrap(~Metric) +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
ggsave(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure2_TI/panelD2.pdf"),
       width=8,height=3)
ggplot(scores_aggreg, aes(x=Dimension, y=Overall_score, fill=Method_id)) + #3x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
  ylim(c(0,1)) +
  scale_colour_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                               rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                               rgb(50,116,161,max=255))) +
  scale_fill_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                             rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                             rgb(50,116,161,max=255))) +
  facet_wrap(~Metric) +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
ggsave(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure2_TI/panelD3.pdf"),
       width=8,height=3)


scores_aggreg_liste <- list(scores_aggreg1,scores_aggreg2,scores_aggreg3,scores_aggreg4,
                            scores_aggreg1b,scores_aggreg2b,scores_aggreg3b,scores_aggreg4b)
names(scores_aggreg_liste) <- c("duo louvain highid","duo leiden highid","seurat louvain highid","seurat leiden highid",
                                "duo louvain lowid", "duo leiden lowid", "seurat louvain lowid", "seurat leiden lowid")
for (i in seq(scores_aggreg_liste)) {
  for (metric in c("Met_corr","Met_F1","Overall_score")) {
    title_id = strsplit(names(scores_aggreg_liste)[i]," ")[[1]][3]
    title_preprocess = paste(strsplit(names(scores_aggreg_liste)[i]," ")[[1]][1:2], collapse="_")
    ggplot(scores_aggreg_liste[[i]], aes(x=Dimension, y=Overall_score, fill=Method_id)) + #3x8
      geom_boxplot(outlier.shape=NA, alpha=0.6) +
      geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
      stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), shape=24, size=2, fill="white", color="black", aes(group=Method_id)) +
      ylim(c(0,1)) +
      scale_colour_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                                   rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                                   rgb(50,116,161,max=255))) +
      scale_fill_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                                 rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                                 rgb(50,116,161,max=255))) +
      facet_wrap(~Metric, nrow=2) +
      coord_flip() +
      theme(panel.background = element_rect(fill="grey98"),
            panel.grid = element_line(colour = "grey80"),
            axis.text = element_text(size=15))
    ggsave(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure5_supp/supp_stabaggreg_all_",title_id,"_",metric,"_",title_preprocess,".pdf"),
           width=8,height=6)
  }
}
