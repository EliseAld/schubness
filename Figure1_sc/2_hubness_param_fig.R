library(ggplot2, quietly = T)
library(viridis, quietly = T)
library(ggridges, quietly = T)
library(scales, quietly = T)
library(ggpubr, quietly = T)
library(pbapply)
library(knn.covertree)
library(ggh4x)
library(future.apply)

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC", "Zhengmix4eq", "Zhengmix4uneq","Zhengmix8eq")

# Figure 1
df <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/slurm/hubness_param_ggplot1.rds")
df2 <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/slurm/hubness_param_ggplot1bis.rds")

labz_method = rep(c("2k (% of hubs)", "Antihubs (% of hubs)", "Asymmetry (% of edges)", "Skewness"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("2k","Antihubs","Asymmetry","Skewness")))
scales_y <- list(scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,22)))

ggplot(df, aes(x=Dimension, y=GHubness)) + # 3x14
   geom_line(aes(group=Dataset, color=Dataset)) +
   geom_point(size=0.1) +
   expand_limits(x=c(0,6500)) +
   facet_wrap(~Method, scales="free_y", labeller = labeller(facet = labz_method), nrow=1) +
   facetted_pos_scales(y=scales_y, x=list(scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10())) +
   ylab(NULL) +
   scale_color_viridis_d(option="plasma") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15))

# Figure 1 supp
df_supp <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/slurm/hubness_param_ggplot2.rds")
df2_supp <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/slurm/hubness_param_ggplot2bis.rds")

labz_method = rep(c("Mean (% of hubs)", "Maximum hub score"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("Mean+Sd","Max")))
scales_y <- list(scale_y_continuous(limits = c(0,1)),
                 scale_y_continuous(limits = c(0,4.5)))

ggplot(df2_supp, aes(x=Dimension, y=GHubness)) + #3x7
   geom_line(aes(group=Dataset, color=Dataset)) +
   geom_point(size=0.1) +
   facet_wrap(~Method, scales = "free_y", labeller = labeller(facet = labz_method), nrow = 1) +
   facetted_pos_scales(y=scales_y) +
   ylab(NULL) +
   scale_color_viridis_d(option="plasma") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"),
         axis.text = element_text(size=15))

# estimate ID of each dataset
# Screeplot
#sdev <- data.frame("Sdev_"=unlist(lapply(data_sdev, function(x)
#   unlist(x)/unlist(x)[1])),
#   "Dataset"=unname(unlist(mapply(x=dataset, y=data_sdev, function(x,y)
#      rep(x, nrow(y)), SIMPLIFY = T))),
#   "PC"=unlist(sapply(data_sdev, function(x)
#      seq(nrow(x)))))
#ggplot(sdev[sdev$PC<30,], aes(x=PC, y=Sdev_, color=Dataset)) +
#   geom_point() +
#   geom_line(aes(group=Dataset)) +
#   geom_hline(yintercept=0.1) +
#   scale_color_viridis_d(option="plasma") +
#   ylab("Sd/Sd1") +
#   xlab("Principal Components")

#linear_dim_df <- data.frame("GID_pca"=unlist(lapply(data_sdev, function(x)  return(which(unlist(x)<=unlist(x)[1]/10)[1]))),
#                         "Dataset"=dataset)
#ggplot(linear_dim_df, aes(x=Dataset, y=GID_pca, color=Dataset)) +
#  scale_color_viridis_d(option="plasma") +
#  geom_point() +
#  ylab("Intrinsic Dimension")

# Get dropout
#data_in <- "/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/DimRedPaper/DuoClustering2018/sce_full/sce_full_"
#data_counts <- pblapply(dataset,
#                 function(x) {tmp <- readRDS(paste0(data_in,x,".rds"));
#                 tmp <- tmp@assays$data@listData$counts;
#                 hvg <- names(sort(apply(tmp,1,var),
#                                   decreasing = T)[1:min(1e4,nrow(tmp))]);
#                 tmp <- tmp[hvg,]
#                 return(tmp)})
#linear_dim_df$Dropout <- sapply(data_counts, function(x)
#   mean(x==0)*100)

# ID x Dropout
#ggplot(linear_dim_df[-grep("TCC",linear_dim_df$Dataset),], aes(x=Dropout, y=GID_pca, color=Dataset)) +
#  scale_color_viridis_d(option="plasma") +
#  geom_point() +
#  ylab("Intrinsic Dimension")

# Dropout x Hubness
#k=10
#Dimension=50
#hubness_dropout <- mapply(x=data_hubness, y=data_counts, function(x,y)
#   data.frame("Hubness"=x$score[x$Dimension==Dimension & x$k==k],
#              "Dropout"=apply(y,2,function(z) mean(z==0)*100)), SIMPLIFY = F)
#mapply(z=hubness_dropout, y=dataset, function(z,y)
#   ggplot(z, aes(x=Dropout, y=Hubness, color=Hubness)) +
#      geom_point() +
#      scale_color_viridis(option="plasma") +
#      scale_y_log10() +
#      scale_x_log10() +
#      ggtitle(y), SIMPLIFY = F)

