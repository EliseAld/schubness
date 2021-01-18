library(future.apply)
library(furrr)
library(pbapply)
library(Seurat)
library(ggpubr)
library(RANN)
plan("multiprocess", workers=30)
options(future.globals.maxSize=2097152000)

data_load <- function(path) {
   source(file=path)
   return(hubness)
}
get_hub3 <- function(df,kval,pcval,log=F) { # Using the method from [1] >= 2k
   if (log==T) {
      hub_nb = c()
      ymin = c()
      ymax = c()
      y_increment = 1
      for (d in pcval) {
         for (k in kval) {
               val = df$score[df$Dimension==d & df$k==k]
               hub_nb = c(hub_nb, sum(val>=2*log(k)))
               ymin = c(ymin,y_increment)
               ymax = c(ymax,y_increment+1)
         }
      y_increment <- y_increment+1
      }
   }
   else {
      hub_nb = c()
      ymin = c()
      ymax = c()
      y_increment = 1
      for (d in pcval) {
         for (k in kval) {
               val = df$score[df$Dimension==d & df$k==k]
               hub_nb = c(hub_nb, sum(val>=2*k))
               ymin = c(ymin,y_increment)
               ymax = c(ymax,y_increment+1)
         }
      y_increment <- y_increment+1
      }
   }
   hub_val <- data.frame("k"=rep(kval,times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=hub_nb,
                         "ymin"=ymin,
                         "ymax"=ymax)
   hub_val$Threshold <- 2*hub_val$k
   return(hub_val)
}
get_hub6 <- function(df,kval,pcval) { # anti hubs
   anti_hub = c()
   for (d in pcval) {
      for (k in kval) {
           anti_hub = c(anti_hub,unique(df$Antihubs[df$Dimension==d & df$k==k]))
      }
   }
   hub_val <- data.frame("k"=rep(kval,times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=anti_hub)
   hub_val$Dimension <- factor(hub_val$Dimension, levels = pcval)
   return(hub_val)
}
rev_cov_hub <- function(knn_graph, hub_scores, n) {
  attempted_hubs <- which(order(hub_scores, decreasing = T) %in% seq(n))
  size_cov <- sum(apply(knn_graph,1,function(x)
    any(x %in% attempted_hubs)))
  return(size_cov/length(hub_scores)*100)
}
get_hub_rev_cov <- function(data, hub_scores, k=10, thd=0.1) {
  k.nn <- nn2(t(data), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  n.val=seq(from=10, to=ncol(data), by=5)
  idx=1
  cov_size <- rev_cov_hub(k.nn, hub_scores, n.val[idx])
  while (cov_size[idx]<100 & idx < length(n.val)) {
                         idx=idx+1
                         cov_size=c(cov_size,rev_cov_hub(k.nn, hub_scores, n.val[idx]))
  }
  increment <- diff(cov_size)
  idx_plateau <- min(ifelse(all(increment>=thd),length(increment)+1,which(increment<thd)[1]+1), length(n.val))
  hub_nb <- n.val[idx_plateau]
  return(colnames(data)[which(order(hub_scores, decreasing = T) %in% seq(hub_nb))])
}

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
data_hubness <- data_load("/home/externe/cnrs/aregimbeau/dev/elise/slurm/hubness_load_write_data_duo.R")
rm(hubness,i)
data_in <- "/shared/projects/mne/elise/sc/data/"
data_counts <- pblapply(dataset,
                 function(x) read.table(paste0(data_in,x,".csv")))

data_hubness <- lapply(seq(length(data_hubness)),
                    function(x) {data_hubness[[x]]$cellID <- names(hubness_scores_data[[x]][[1]][[1]]);
                    return(data_hubness[[x]])})
rm(hubness_scores_data)

entropy <- pblapply(dataset, function(x) 
  c(read.table(paste0("/shared/projects/mne/elise/sc/entropy/",x,".csv"), quote="\"", comment.char=""))[[1]])
print("data loaded!")

params <- data.frame(10,50)
names(params) <- c("k","Dimension")

hub_scores <- pblapply(seq(length(data_hubness)),
                       function(x) data_hubness[[x]]$score[data_hubness[[x]]$Dimension==params$Dimension & data_hubness[[x]]$k==params$k])

hub_scores <- pblapply(seq(length(hub_scores)),
                       function(x) {names(hub_scores[[x]]) <- data_hubness[[x]]$cellID[seq(length(hub_scores[[x]]))];
                       return(hub_scores[[x]])})


hub3 <- future_lapply(seq(length(dataset)),
                 function(x) colnames(data_counts[[x]])[hub_scores[[x]] >= 2*params$k])

hub6 <- future_lapply(seq(length(dataset)),
                 function(x) colnames(data_counts[[x]])[hub_scores[[x]]==0])
hub_rev <- future_mapply(SIMPLIFY = F,
                  x=data_counts, y=hub_scores, function(x,y)
                    get_hub_rev_cov(x, y))
print("hub done!")

seurat <- future_lapply(data_counts,CreateSeuratObject)
seurat <- future_lapply(seq(length(seurat)),
                   function(x) {seurat[[x]]$KThd <- c("normal","hub")[factor(colnames(seurat[[x]]) %in% hub3[[x]], levels = c("FALSE","TRUE"))];
                   seurat[[x]]$Antihubs <- c("normal","antihub")[factor(colnames(seurat[[x]]) %in% hub6[[x]], levels = c("FALSE","TRUE"))];
                   seurat[[x]][["percent.mt"]] <- PercentageFeatureSet(seurat[[x]], pattern = "^MT-");
                   seurat[[x]][["percent.ribo"]] <- PercentageFeatureSet(seurat[[x]], pattern = c("^RPS","^RPL"));
                   seurat[[x]][["Dropout_rate"]] <- apply(data_counts[[x]],2,function(z) mean(z==0)*100);
                   seurat[[x]][["Rev_cov"]] <- c("normal","hub")[factor(colnames(seurat[[x]]) %in% hub_rev[[x]], levels = c("FALSE","TRUE"))];
                   seurat[[x]]$Entropy <- entropy[[x]];
                   return(seurat[[x]])})
seurat <- future_lapply(seurat,
                 function(x) {x$Hub_status_k="Normal";
                 x$Hub_status_rev="Normal";
                 for (y in seq(ncol(x))) {
                      x$Hub_status_k[y] <- ifelse(x$KThd[y]=="hub",
                                                "Hub",
                                                ifelse(x$Antihubs[y]=="antihub",
                                                       "Antihub",
                                                       "Normal"))
                      x$Hub_status_rev[y] <- ifelse(x$Rev_cov[y]=="hub",
                                                "Hub",
                                                ifelse(x$Antihubs[y]=="antihub",
                                                       "Antihub",
                                                       "Normal"))};
                 x$Hub_status_rev <- factor(x$Hub_status_rev, levels=c("Antihub","Normal","Hub"))
                 x$Hub_status_k <- factor(x$Hub_status_k, levels=c("Antihub","Normal","Hub"))
                 return(x)})
seurat <- future_mapply(SIMPLIFY = F,
                 x=seurat, y=data_hubness, function(x,y)
                   {x[["Hubness_score"]] <- y$score[y$k==params$k & y$Dimension==params$Dimension];
                   return(x)})
# Make UMAP
seurat <- future_map(seurat,NormalizeData)
seurat <- future_map(seurat,ScaleData)
seurat <- future_map(seurat,FindVariableFeatures)
seurat <- future_map(seurat,function(x) RunPCA(x, verbose = F))
seurat <- future_map(seurat,function(x) RunUMAP(x, dims=1:50))
saveRDS(seurat, file="/shared/projects/mne/elise/sc/qc_seurat.rds")
#seurat <- readRDS(file="/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/qc_seurat.rds")

#for (idx in seq(length(seurat))) {
#print(DimPlot(seurat[[idx]],
#        reduction=reduction,
#        group.by="Hub_status_k",
#        pt.size = 1) +
#  scale_color_brewer(palette = "Set1") +
#    xlab("UMAP1") + ylab("UMAP2"))
#print(DimPlot(seurat[[idx]],
#        reduction="pca",
#        group.by="Hub_status_k",
#        pt.size = 1) +
#  scale_color_brewer(palette = "Set1") +
#    xlab("PC1") + ylab("PC2"))
#print(DimPlot(seurat[[idx]],
#        reduction=reduction,
#        group.by="Hub_status_rev",
#        pt.size = 1) +
#  scale_color_brewer(palette = "Set1") +
#    xlab("UMAP1") + ylab("UMAP2"))
#print(DimPlot(seurat[[idx]],
#        reduction="pca",
#        group.by="Hub_status_rev",
#        pt.size = 1) +
#  scale_color_brewer(palette = "Set1") +
#    xlab("PC1") + ylab("PC2"))
#}


#for (idx in seq(length(seurat))) {
#qc_df <- seurat[[idx]]@meta.data

# Dropout rate
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=Dropout_rate, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#    xlab("Hub status") +
#    ylab("Dropout rate") +
#    labs(fill = "Hub status") + 
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25)))

# nfeature rate
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=nFeature_RNA, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#    xlab("Hub status") +
#    ylab("Number of total features") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25)))

# ncount rate
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=nCount_RNA, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#    xlab("Hub status") +
#    ylab("Number of unique genes") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25)))

# entropy
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=Entropy, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#    xlab("Hub status") +
#    ylab("scEntropy") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25)))

# percent ribo rate
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=percent.ribo, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15)))

# percent mito rate
#print(ggplot(qc_df, aes(x=Hub_status_rev,
#                  y=percent.mt, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group=".all.") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15)))

# Dropout rate histo
#ggplot(qc_df, aes(Dropout_rate, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# nfeature histo
#ggplot(qc_df, aes(nFeature_RNA, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# ncount histo
#ggplot(qc_df, aes(nCount_RNA, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# percent mito histo
#ggplot(qc_df, aes(percent.mt, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# percent ribo histo
#ggplot(qc_df, aes(percent.ribo, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()
#}


#for (idx in seq(length(seurat))) {
#qc_df <- seurat[[idx]]@meta.data

# Dropout rate
#print(ggplot(qc_df, aes(x=Hub_status_k,
#                  y=Dropout_rate, fill=Hub_status_k)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal"))

# nfeature rate
#print(ggplot(qc_df, aes(x=Hub_status_k,
#                  y=nFeature_RNA, fill=Hub_status_k)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal"))

# ncount rate
#print(ggplot(qc_df, aes(x=Hub_status_k,
#                  y=nCount_RNA, fill=Hub_status_k)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal"))

# percent ribo rate
#print(ggplot(qc_df, aes(x=Hub_status_k,
#                  y=percent.ribo, fill=Hub_status_k)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal"))

# percent mito rate
#print(ggplot(qc_df, aes(x=Hub_status_k,
#                  y=percent.mt, fill=Hub_status_k)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  stat_compare_means(method="t.test", label="p.signif", ref.group=".all."))

# Dropout rate histo
#ggplot(qc_df, aes(Dropout_rate, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# nfeature histo
#ggplot(qc_df, aes(nFeature_RNA, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# ncount histo
#ggplot(qc_df, aes(nCount_RNA, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# percent mito histo
#ggplot(qc_df, aes(percent.mt, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()

# percent ribo histo
#ggplot(qc_df, aes(percent.ribo, fill=Hub_status)) +
#  geom_bar(position="dodge") +
#  scale_fill_viridis_d()
#}

# Dropout x Hubness
#for (idx in seq(length(seurat))) {
#df = data.frame("UMAP1"=seurat[[idx]]@reductions$umap@cell.embeddings[,1],
#                "UMAP2"=seurat[[idx]]@reductions$umap@cell.embeddings[,2],
#                "PC1"=seurat[[idx]]@reductions$pca@cell.embeddings[,1],
#                "PC2"=seurat[[idx]]@reductions$pca@cell.embeddings[,2],
#                "Dropout"=seurat[[idx]]$Dropout_rate,
#                "Hubness"=seurat[[idx]]$Hubness_score)
#print(ggplot(df, aes(x=UMAP1, y=UMAP2, color=Hubness, alpha=Dropout)) +
#  geom_point() +
#  scale_color_viridis_c(option='plasma', direction = (-1)))
#print(ggplot(df, aes(x=PC1, y=PC2, color=Hubness, alpha=Dropout)) +
#  geom_point() +
#  scale_color_viridis_c(option='plasma', direction = (-1)))
#}

#qc_df_supp <- do.call(rbind, lapply(seurat, function(x) x@meta.data))
#qc_df_supp$Dataset <- unlist(mapply(x=dataset, y=seurat, function(x,y)
#  rep(x, ncol(y))))


# Dropout rate
#ggplot(qc_df_supp, aes(x=Hub_status_rev, y=Dropout_rate, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  xlab("Hub status") +
#  ylab("") +
#  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#  theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25),
#        axis.text.x = element_text(angle=90),
#        axis.title.x = element_blank(),
#        strip.text.x = element_text(size=15),
#        legend.position ="none")

# nfeature rate
#print(ggplot(qc_df_supp, aes(x=Hub_status_rev, y=nFeature_RNA, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  xlab("Hub status") +
#  ylab("") +
#  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#  theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25),
#        axis.text.x = element_text(angle=90),
#        axis.title.x = element_blank(),
#        strip.text.x = element_text(size=15),
#        legend.position ="none"))

# ncount rate
#print(ggplot(qc_df_supp, aes(x=Hub_status_rev,
#                  y=nCount_RNA, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#  xlab("Hub status") +
#  ylab("") +
#  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#  theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25),
#        axis.text.x = element_text(angle=90),
#        axis.title.x = element_blank(),
#        strip.text.x = element_text(size=15),
#        legend.position ="none"))

# entropy
#print(ggplot(qc_df_supp, aes(x=Hub_status_rev,
#                  y=Entropy, fill=Hub_status_rev)) +
#  geom_violin() +
#  scale_fill_brewer(palette = "Set1") +
#    xlab("Hub status") +
#    ylab("") +
#  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
#  stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=25),
#        axis.text.x = element_text(angle=90),
#        axis.title.x = element_blank(),
#        strip.text.x = element_text(size=15),
#        legend.position ="none"))

#ggarrange(plotlist=lapply(seurat, function(z)
#  DimPlot(z,
#        reduction="umap",
#        group.by="Hub_status_rev",
#        pt.size = 1) +
#    scale_color_brewer(palette = "Set1") +
#    xlab("UMAP1") +
#    ylab("UMAP2")),
#  common.legend = T)
  
#ggarrange(plotlist=lapply(seurat, function(z)
#  DimPlot(z,
#        reduction="pca",
#        group.by="Hub_status_rev",
#        pt.size = 1) +
#    scale_color_brewer(palette = "Set1") +
#    xlab("PC1") +
#    ylab("PC2")),
#  common.legend = T,
#  nrow = 2, ncol=6)

#euclid <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_duo/hub_stab_sampling_normduo_scaleTrue_euclidean.csv")
#cosine <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_duo/hub_stab_sampling_normduo_scaleTrue_cosine.csv")

#stab_df <- rbind(data.frame("Hub_overlap"= unlist(c(euclid)),
#                        "Dataset"=rep(dataset, times=ncol(euclid)),
#                        "Dimension"=factor(rep(c(25,50,100,500), each=nrow(euclid))),
#                        "Metric"="Euclidean"),
#                 data.frame("Hub_overlap"= unlist(c(cosine)),
#                        "Dataset"=rep(dataset, times=ncol(cosine)),
#                        "Dimension"=factor(rep(c(25,50,100,500), each=nrow(cosine))),
#                        "Metric"="Cosine"))
#ggplot(stab_df, aes(x=Dimension, y=Hub_overlap, fill=Dimension)) +
#  geom_boxplot(outlier.shape=NA, alpha=0.6) +
#  geom_jitter(aes(color=Dimension), fill="black") +
#  ylim(c(0,90)) +
#  facet_wrap(~Metric) +
#  theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80")) +
#  ylab("Hubs overlap (%)")
print("SUCCESS")
