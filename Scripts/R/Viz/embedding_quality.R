library(ggplot2)
library(cluster)
library(dplyr)

# ====
# Load data
# ====
path_c <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/cytotrace/h5_jo_ktuned6/"
path_z <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/h5_jo_ktuned6/"


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
                          'gold_aging-hsc-young_kowalczyk','gold_aging-hsc-old_kowalczyk'))

norm_scale = "True"
metric=c('cosine','euclidean')
n_comps=c(50,500)
do_norm = 'seurat'
clustering_algo = 'leiden'

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

methods <- c("gauss","umap","nothing","mp_normal","ls","ls_nicdm","dsl")
reference = c("Raw_kl","500D_kl","2D_kl",
              "Raw_ce","500D_ce","2D_ce",
              "Raw_st","500D_st","2D_st",
              "Raw_cor","500D_cor","2D_cor")
embedding = c("tSNE","UMAP","PAGA+UMAP")

# ====
# Functions
# ====
prep_df <- function(scores, dataset, name, origin) {
  scores2 <- do.call(rbind,lapply(seq(scores), function(x) {
    tmp <- lapply(seq(scores[[x]]), function(y) {
      scores[[x]][[y]] <- sapply(scores[[x]][[y]], as.numeric);
      scores[[x]][[y]] <- c(scores[[x]][[y]])
      scores[[x]][[y]] <- data.frame("Score" = scores[[x]][[y]],
                                     "Method_id" = rep(methods, each=ifelse(name=="Cost",length(reference),length(reference[1:3]))),
                                     "info" = dataset[y])
      scores[[x]][[y]]});
    tmp2 <- do.call(rbind,tmp);
    tmp2$Embedding <- embedding[x];
    tmp2
  }))
  if (name=="Cost") {
    scores2$CostFunction <- sapply(reference, function(x) {
      strsplit(x, "_")[[1]][2]})
    scores2$Reference <-  sapply(reference, function(x) {
      strsplit(x, "_")[[1]][1]})
  }
  else {
    scores2$CostFunction <- "trust"
    scores2$Reference <- sapply(reference[1:3], function(x) {
      strsplit(x, "_")[[1]][1]})
  }
  if (origin=="cyto") {
    scores2$Metric <- sapply(scores2$info, function(x) strsplit(x, "_")[[1]][6])
    scores2$Dimension <- sapply(sapply(scores2$info, function(x) strsplit(x, "_")[[1]][5]), function(y) strsplit(y,"ncomps")[[1]][2])
    scores2$Recipe <- sapply(sapply(scores2$info, function(x) strsplit(x, "_")[[1]][3]), function(y) strsplit(y,"norm")[[1]][2])
    scores2$Clustering <- sapply(scores2$info, function(x) strsplit(x, "_")[[1]][7])
    scores2$Dataset_id <- sapply(scores2$info, function(x) strsplit(x, "_norm")[[1]][1])
  }
  else {
    scores2$Metric <- sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[2])
    scores2$Dimension <- sapply(sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[3]), function(y) strsplit(y,"ncomps")[[1]][2])
    scores2$Recipe <- sapply(sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[5]), function(y) strsplit(y,"norm")[[1]][2])
    scores2$Clustering <- sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[1])
    scores2$Dataset_id <- sapply(scores2$info, function(x) strsplit(x, "_norm")[[1]][1])
  }
  scores2 <- scores2[,-3]
  return(scores2)
}
prep_viz <- function(coord, dataset, origin, clust) {
  scores2 <- do.call(rbind,lapply(seq(coord), function(x) {
    tmp <- lapply(seq(coord[[x]]), function(y) {
      do.call(rbind,lapply(seq(coord[[x]][[y]]), function(z) {
      coord[[x]][[y]][[z]] <- data.frame(sapply(coord[[x]][[y]][[z]], as.numeric));
      colnames(coord[[x]][[y]][[z]]) <- c('X1','X2');
      coord[[x]][[y]][[z]]$info = dataset[y];
      coord[[x]][[y]][[z]]$Cluster = factor(c(clust[[y]][[z]])[[1]]);
      coord[[x]][[y]][[z]]$Method_id = methods[z];
      coord[[x]][[y]][[z]]}))});
    tmp2 <- do.call(rbind,tmp);
    tmp2$Embedding <- strsplit(strsplit(viz_extension[x],".csv")[[1]][1],"_")[[1]][2];
    tmp2
  }))
  if (origin=="cyto") {
    scores2$Metric <- sapply(scores2$info, function(x) strsplit(x, "_")[[1]][6])
    scores2$Dimension <- sapply(sapply(scores2$info, function(x) strsplit(x, "_")[[1]][5]), function(y) strsplit(y,"ncomps")[[1]][2])
    scores2$Recipe <- sapply(sapply(scores2$info, function(x) strsplit(x, "_")[[1]][3]), function(y) strsplit(y,"norm")[[1]][2])
    scores2$Clustering <- sapply(scores2$info, function(x) strsplit(x, "_")[[1]][7])
    scores2$Dataset_id <- sapply(scores2$info, function(x) strsplit(x, "_norm")[[1]][1])
  }
  else {
    scores2$Metric <- sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[2])
    scores2$Dimension <- sapply(sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[3]), function(y) strsplit(y,"ncomps")[[1]][2])
    scores2$Recipe <- sapply(sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[5]), function(y) strsplit(y,"norm")[[1]][2])
    scores2$Clustering <- sapply(scores2$info, function(x) rev(strsplit(x, "_")[[1]])[1])
    scores2$Dataset_id <- sapply(scores2$info, function(x) strsplit(x, "_norm")[[1]][1])
  }
  scores2 <- scores2[,-3]
  return(scores2)
}
prep_aggreg <- function(df1,df2,keep_id=TRUE) {
  res <- rbind(df1,df2)
  res$Dimension <- factor(res$Dimension, levels=n_comps)
  res$Method_id <- factor(res$Method_id, levels=c("dsl","ls_nicdm","ls","mp_normal","gauss","umap","nothing"))
  if (keep_id) {
    res <- res[res$Dataset_id %in% keep_id_datasets,]
  }
  res <- res[!(res$Metric=="cosine" &
                 res$Method_id=="dsl"),]
  return(res)
}

# ====
# Retrieve results
# ====
embed_extension <- c("_tsne_c.csv","_umap_c.csv","_paga_c.csv")
trust_extension <- c("_tsne_t.csv","_umap_t.csv","_paga_t.csv")
viz_extension <- c("_tsne.csv","_umap.csv","_paga.csv")

cost_c <- lapply(embed_extension, function(y)
  lapply(datasets_c, function(x)
    read.csv(file=paste0(path_c,x,y), header=FALSE, dec=",")))
cost_z <- lapply(embed_extension, function(y)
  lapply(datasets_z, function(x)
    read.csv(file=paste0(path_z,x,y), header=FALSE, dec=",")))
trust_c <- lapply(trust_extension, function(y)
  lapply(datasets_c, function(x)
    read.csv(file=paste0(path_c,x,y), header=FALSE, dec=",")))
trust_z <- lapply(trust_extension, function(y)
  lapply(datasets_z, function(x)
    read.csv(file=paste0(path_z,x,y), header=FALSE, dec=",")))
viz_c <- lapply(viz_extension, function(y)
  lapply(datasets_c, function(x) lapply(methods, function(z)
    read.csv(file=paste0(path_c,x,'_',z,y), header=FALSE, dec=","))))
viz_z <- lapply(viz_extension, function(y)
  lapply(datasets_z, function(x) lapply(methods, function(z)
    read.csv(file=paste0(path_z,x,'_',z,y), header=FALSE, dec=","))))
#clust_c <- lapply(datasets_c, function(x) lapply(methods, function(z)
#    read.csv(file=paste0(path_c,x,'_cluster.csv'), header=T, row.names=1)))
#clust_z <- lapply(datasets_z, function(x) lapply(methods, function(z)
#    read.csv(file=paste0(path_z,x,'_',z,'_cluster.csv'), header=T, row.names=1)))
label_c <- lapply(datasets_c, function(x) lapply(methods, function(z)
  read.table(file=paste0(path_c,x,'_',z,'_label.csv'), quote="\"", comment.char="")))
label_z <- lapply(datasets_z, function(x) lapply(methods, function(z)
  read.table(file=paste0(path_z,x,'_',z,'_label.csv'), quote="\"", comment.char="")))

keep_id_datasets <- c(dataset_choice_c,dataset_choice_z[c(2,3,6,8,14)])

# ====
# Put all scores together
# ====
resc <- prep_df(cost_c, datasets_c, "Cost", "cyto")
resz <- prep_df(cost_z, datasets_z, "Cost", "zhou")

trustc <- prep_df(trust_c, datasets_c, "Trust", "cyto")
trustz <- prep_df(trust_z, datasets_z, "Trust", "zhou")

vizc <- prep_viz(viz_c, datasets_c, "cyto", label_c)
vizz <- prep_viz(viz_z, datasets_z, "zhou", label_z)

# ====
# Plot
# ====
res <- prep_aggreg(resc, resz)
trust <- prep_aggreg(trustc, trustz)
viz <- prep_aggreg(vizc, vizz)
score <- prep_aggreg(res,trust)
for (recipe in do_norm) {
  for (reference in c("Raw","500D","2D")) {
    for (embeddings in embedding) {
      for (costfn in c("kl","ce","st","cor","trust")) {
        pdf(file = paste0("/Users/elise/Desktop/Github/Hubness_sc/Figure3_embedding/fig_",
                          recipe,embeddings,costfn,reference,
                          ".pdf"), width = 8, height = 3)
        print(
        ggplot(score[score$Recipe==recipe &
                       score$Clustering=="leiden" &
                       score$Reference==reference &
                       score$Embedding==embeddings &
                       score$CostFunction==costfn &
                       score$Method_id!="nothing",], aes(x=Dimension, y=Score, fill=Method_id)) +
          geom_boxplot(outlier.shape=NA, alpha=0.6) +
          facet_wrap(~Metric, scales = "free") +
          geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
          scale_y_log10() +
          scale_fill_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                                     rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                                     rgb(50,116,161,max=255))) +
          scale_color_manual(values=c(rgb(214,132,189,max=255),rgb(132,91,83,max=255),rgb(147,114,178,max=255),
                                      rgb(192,61,62,max=255),rgb(58,146,58,max=255),rgb(225,129,44,max=255),
                                      rgb(50,116,161,max=255))) +
          theme(panel.background = element_rect(fill="grey98"),
                panel.grid = element_line(colour = "grey80"),
                axis.text = element_text(size=15)))
        dev.off()
      }
    }
  }
}
for (recipe in do_norm) {
  for (clustering in clustering_algo) {
    for (embeddings in c("tsne",'umap','paga')) {
      for (metric in c('euclidean','cosine')) {
        for (dim in n_comps) {
          for (dataset in unique(viz$Dataset_id)) {
            if (!file.exists( paste("/Users/elise/Desktop/Github/Hubness_sc/Figure3_embedding/proj",
                                    dataset,recipe,dim,metric,clustering,embeddings,
                                    ".pdf", sep = "_"))) {
            tmp <- viz[viz$Recipe==recipe &
                         viz$Clustering==clustering &
                         viz$Embedding==embeddings &
                         viz$Metric==metric &
                         viz$Dimension==dim &
                         viz$Dataset_id==dataset,]
            sil=c()
            for (met in c('nothing','umap','gauss','mp_normal','ls','ls_nicdm','dsl')) {
              sm = summary(silhouette(as.numeric(tmp$Cluster[tmp$Method_id==met]),
                                      dist(t(matrix(c(tmp$X1[tmp$Method_id==met],
                                                      tmp$X2[tmp$Method_id==met]),
                                                    nrow=2,byrow=T)))))
              if (is.na(silhouette(as.numeric(tmp$Cluster[tmp$Method_id==met]),
                                   dist(t(matrix(c(tmp$X1[tmp$Method_id==met],
                                                   tmp$X2[tmp$Method_id==met]),
                                                 nrow=2,byrow=T)))))) {
                sil = c(sil,NA)
              }
              else {
                sil = c(sil,
                      mean(sm$clus.avg.widths))
              }
            }
            sil = data.frame("Method_id"=c('nothing','umap','gauss','mp_normal','ls','ls_nicdm','dsl'),
                             "Silhouette"=paste0("silhouette=",round(sil,2)))
            pdf(file = paste("/Users/elise/Desktop/Github/Hubness_sc/Figure3_embedding/proj",
                              dataset,recipe,dim,metric,clustering,embeddings,
                              ".pdf", sep = "_"), width = 8, height = 8)
            print(
            ggplot(viz[viz$Recipe==recipe &
                         viz$Clustering==clustering &
                         viz$Embedding==embeddings &
                         viz$Metric==metric &
                         viz$Dimension==dim &
                         viz$Dataset_id==dataset,], aes(x=X1, y=X2, color=Cluster)) +
              geom_point() +
              facet_wrap(~factor(Method_id, levels=c('nothing','umap','gauss','mp_normal','ls','ls_nicdm','dsl')), scales = "free") +
              theme(panel.background = element_rect(fill="grey98"),
                    panel.grid = element_line(colour = "grey80"),
                    axis.text = element_text(size=15)) +
              geom_text(data=sil, aes(x=-Inf, y=Inf, label=Silhouette),
                        hjust = -0.1, vjust = 1.1,
                        colour="black", inherit.aes=FALSE, parse=FALSE)
            )
            dev.off()
            }
          }
        }
      }
    }
  }
}
