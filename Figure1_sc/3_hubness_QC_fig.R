library(future.apply)
library(furrr)
library(pbapply)
library(Seurat)
library(ggpubr)
library(RANN)

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")

params <- data.frame(10,50)
names(params) <- c("k","Dimension")

# Figure 1 - one Example with idx=10
seurat <- readRDS(file="/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/slurm/qc_seurat.rds")
idx=10
qc_df <- seurat[[idx]]@meta.data
# UMAP
tplot=DimPlot(seurat[[idx]], #5x7
              reduction="umap",
              group.by="Hub_status_rev",
              pt.size = c(0.5,0.5,1)[seurat[[idx]]$Hub_status_rev]) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text = element_text(size=15))
tplot[[1]]$layers[[1]]$aes_params$alpha = ifelse(seurat[[idx]]$Hub_status_rev=="Normal", .1, 1)
tplot
# PCA
tplot=DimPlot(seurat[[idx]],
        reduction="pca",
        group.by="Hub_status_rev",
        pt.size = c(0.5,0.5,1)[seurat[[idx]]$Hub_status_rev]) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text = element_text(size=15))
tplot[[1]]$layers[[1]]$aes_params$alpha = ifelse(seurat[[idx]]$Hub_status_rev=="Normal", .1, 1)
tplot
# Dropout rate
print(ggplot(qc_df, aes(x=Hub_status_rev, #4x7
                  y=Dropout_rate, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
    xlab("Hub status") +
    ylab("Dropout rate") +
    labs(fill = "Hub status") + 
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25)))
# nfeature rate
print(ggplot(qc_df, aes(x=Hub_status_rev,
                  y=nFeature_RNA, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
    xlab("Hub status") +
    ylab("Number of total features") +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25)))

# ncount rate
print(ggplot(qc_df, aes(x=Hub_status_rev,
                  y=nCount_RNA, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
    xlab("Hub status") +
    ylab("Number of unique genes") +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25)))

# entropy
print(ggplot(qc_df, aes(x=Hub_status_rev,
                  y=Entropy, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
    xlab("Hub status") +
    ylab("scEntropy") +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25)))

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


# Figure 1 supp
qc_df_supp <- do.call(rbind, lapply(seurat, function(x) x@meta.data))
qc_df_supp$Dataset <- unlist(mapply(x=dataset, y=seurat, function(x,y)
  rep(x, ncol(y))))
# Dropout rate
ggplot(qc_df_supp, aes(x=Hub_status_rev, y=Dropout_rate, fill=Hub_status_rev)) +
  geom_violin() + #3,5x28
  scale_fill_brewer(palette = "Set1") +
  xlab("Hub status") +
  ylab("") +
  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
  theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25),
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=15),
        legend.position ="none")

# nfeature rate
print(ggplot(qc_df_supp, aes(x=Hub_status_rev, y=nFeature_RNA, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Hub status") +
  ylab("") +
  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
  theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25),
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=15),
        legend.position ="none"))

# ncount rate
print(ggplot(qc_df_supp, aes(x=Hub_status_rev,
                  y=nCount_RNA, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Hub status") +
  ylab("") +
  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
  theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25),
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=15),
        legend.position ="none"))

# entropy
print(ggplot(qc_df_supp, aes(x=Hub_status_rev,
                  y=Entropy, fill=Hub_status_rev)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set1") +
    xlab("Hub status") +
    ylab("") +
  facet_wrap(~Dataset, scales="free", ncol=length(dataset)) +
  #stat_compare_means(method="t.test", label="p.signif", ref.group="Normal") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15),
         axis.title = element_text(size=25),
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=15),
        legend.position ="none"))

# PCA
ggarrange(plotlist=lapply(seurat, function(z) #4x16
  {tplot=DimPlot(z,
          reduction="pca",
          group.by="Hub_status_rev",
          pt.size = c(0.5,0.5,1)[z$Hub_status_rev]) +
    scale_color_brewer(palette = "Set1") +
    xlab("PC1") +
    ylab("PC2") +
    theme(axis.text = element_text(size=15));
  tplot[[1]]$layers[[1]]$aes_params$alpha = ifelse(z$Hub_status_rev=="Normal", .2, 1)
  tplot}),
  common.legend = T,
  nrow = 2, ncol=6)
# UMAP
ggarrange(plotlist=lapply(seurat, function(z)
  DimPlot(z,
        reduction="umap",
        group.by="Hub_status_rev",
        pt.size = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("UMAP1") +
    ylab("UMAP2")),
  common.legend = T,
  nrow = 2, ncol=6)

euclid <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_duo/hub_stab_sampling_normduo_scaleTrue_euclidean.csv")
cosine <- read.table("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/hub_stab_duo/hub_stab_sampling_normduo_scaleTrue_cosine.csv")

stab_df <- rbind(data.frame("Hub_overlap"= unlist(c(euclid)),
                        "Dataset"=rep(dataset, times=ncol(euclid)),
                        "Dimension"=factor(rep(c(25,50,100,500), each=nrow(euclid))),
                        "Metric"="Euclidean"),
                 data.frame("Hub_overlap"= unlist(c(cosine)),
                        "Dataset"=rep(dataset, times=ncol(cosine)),
                        "Dimension"=factor(rep(c(25,50,100,500), each=nrow(cosine))),
                        "Metric"="Cosine"))
ggplot(stab_df, aes(x=Dimension, y=Hub_overlap, fill=Dimension)) + #3.5x8
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_jitter(aes(color=Dimension), fill="black") +
  ylim(c(0,90)) +
  facet_wrap(~Metric) +
  theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15))
