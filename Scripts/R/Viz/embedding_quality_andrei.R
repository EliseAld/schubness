library(ggplot2)
library(cluster)
library(dplyr)

# ====
# Load data
# ====
path_c <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/cytotrace/h5_jo_ktuned7/"
path_z <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/forTI/zhouDR/datasets/h5_jo_ktuned7/"


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
reference = c("Raw_qdm","500D_qdm","2D_qdm",
              "Raw_qnp","500D_qnp","2D_qnp")
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
                                     "Method_id" = rep(methods, each=length(reference)),
                                     "info" = dataset[y])
      scores[[x]][[y]]});
    tmp2 <- do.call(rbind,tmp);
    tmp2$Embedding <- embedding[x];
    tmp2
  }))
  scores2$CostFunction <- sapply(reference, function(x) {
      strsplit(x, "_")[[1]][2]})
  scores2$Reference <-  sapply(reference, function(x) {
      strsplit(x, "_")[[1]][1]})
  
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
prep_aggreg <- function(df1,df2,keep_id=TRUE, keep_low=FALSE) {
  res <- rbind(df1,df2)
  res$Dimension <- factor(res$Dimension, levels=n_comps)
  res$Method_id <- factor(res$Method_id, levels=c("dsl","ls_nicdm","ls","mp_normal","gauss","umap","nothing"))
  if (keep_id) {
    res <- res[res$Dataset_id %in% keep_id_datasets,]
  }
  if (keep_low) {
    res <- res[!(res$Dataset_id %in% keep_id_datasets),]
  }
  res <- res[!(res$Metric=="cosine" &
                 res$Method_id=="dsl"),]
  return(res)
}

# ====
# Retrieve results
# ====
embed_extension <- c("_tsne_q.csv","_umap_q.csv","_paga_q.csv")

cost_c <- lapply(embed_extension, function(y)
  lapply(datasets_c, function(x)
    read.csv(file=paste0(path_c,x,y), header=FALSE, dec=",")))
cost_z <- lapply(embed_extension, function(y)
  lapply(datasets_z, function(x)
    read.csv(file=paste0(path_z,x,y), header=FALSE, dec=",")))

keep_id_datasets <- c(dataset_choice_c,dataset_choice_z[c(2,3,6,8,14)])

# ====
# Put all scores together
# ====
resc <- prep_df(cost_c, datasets_c, "Cost", "cyto")
resz <- prep_df(cost_z, datasets_z, "Cost", "zhou")

# ====
# Plot
# ====
res <- prep_aggreg(resc, resz)
res2 <- prep_aggreg(resc, resz, keep_id = F, keep_low = T)

for (recipe in do_norm) {
  for (reference in c("Raw","500D","2D")) {
    for (embeddings in embedding) {
      for (costfn in c("qdm","qnp")) {
        pdf(file = paste0("/Users/elise/Desktop/Github/Hubness_sc/Figure3_embedding/andrei_low_",
                          recipe,embeddings,costfn,reference,
                          ".pdf"), width = 8, height = 3)
        print(
        ggplot(res2[res2$Recipe==recipe &
                      res2$Clustering=="leiden" &
                      res2$Reference==reference &
                      res2$Embedding==embeddings &
                      res2$CostFunction==costfn &
                      res2$Method_id!="nothing",], aes(x=Dimension, y=Score, fill=Method_id)) +
          geom_boxplot(outlier.shape=NA, alpha=0.6) +
          facet_wrap(~Metric, scales = "free") +
          geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
          #scale_y_log10() +
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
