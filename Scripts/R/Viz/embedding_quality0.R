library(ggplot2)

# ====
# Load data
# ====
path_c <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/cytotrace/h5_jo_ktuned5/"
path_z <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/h5_jo_ktuned5/"


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

methods <- c("gauss","umap","nothing","mp_normal","ls","ls_nicdm","dsl")

# ====
# Functions
# ====
prep_df <- function(scores, dataset, name, origin) {
  scores2 <- do.call(rbind,lapply(seq(scores), function(x) {
    tmp <- lapply(seq(scores[[x]]), function(y) {
      scores[[x]][[y]] <- sapply(scores[[x]][[y]], as.numeric);
      scores[[x]][[y]] <- c(scores[[x]][[y]])
      scores[[x]][[y]] <- data.frame("Score" = scores[[x]][[y]],
                                     "Method_id" = rep(methods),
                                     "info" = dataset[y],
                                     "CostFn" = ifelse(name=="Cost", extension[x], extension2[x]))
      scores[[x]][[y]]});
    tmp2 <- do.call(rbind,tmp);
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
  return(scores2)
}
prep_aggreg <- function(df1,df2,keep_id=TRUE) {
  res <- rbind(df1,df2)
  res$Dimension <- factor(res$Dimension, levels=n_comps)
  res$Method_id <- factor(res$Method_id, levels=c("dsl","ls_nicdm","ls","mp_normal","umap","gauss","nothing"))
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
extension <- c("_costkl.csv","_costkl2.csv","_costkl3.csv","_costumap4.csv")
extension2 <- c("_trust2.csv","_trust3.csv","_trust4.csv","_trust5.csv")

cost_c <- lapply(extension, function(y)
  lapply(datasets_c, function(x)
    read.table(file=paste0(path_c,x,y))))
cost_z <- lapply(extension, function(y)
  lapply(datasets_z, function(x)
    read.table(file=paste0(path_z,x,y))))
trust_c <- lapply(extension2, function(y)
  lapply(datasets_c, function(x)
    read.table(file=paste0(path_c,x,y))))
trust_z <- lapply(extension2, function(y)
  lapply(datasets_z, function(x)
    read.table(file=paste0(path_z,x,y))))

keep_id_datasets <- c(dataset_choice_c,dataset_choice_z[c(2,3,6,8,14)])

# ====
# Put all scores together
# ====
resc <- prep_df(cost_c, datasets_c, "Cost", "cyto")
resz <- prep_df(cost_z, datasets_z, "Cost", "zhou")

trustc <- prep_df(trust_c, datasets_c, "Trust", "cyto")
trustz <- prep_df(trust_z, datasets_z, "Trust", "zhou")

# ====
# Plot
# ====
res <- prep_aggreg(resc, resz)
trust <- prep_aggreg(trustc, trustz)

ggplot(res[res$Recipe=="seurat" &
             res$Clustering=="leiden" &
             res$CostFn==extension[4],], aes(x=Dimension, y=Score, fill=Method_id)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  facet_wrap(~Metric, scales = "free") +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  scale_y_log10() +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
ggplot(trust[trust$Recipe=="seurat" &
               trust$Clustering=="leiden" &
               trust$CostFn==extension2[4],], aes(x=Dimension, y=Score, fill=Method_id)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  facet_wrap(~Metric, scales = "free") +
  geom_point(size=0.8,position=position_jitterdodge(jitter.width = 0.1),aes(color=Method_id), fill="black") +
  #scale_y_log10() +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
