library(ggplot2, quietly = T)
library(ggpubr, quietly = T)
library(colorspace, quietly = T)
library(pbapply, quietly = T)
library(gridExtra, quietly = T)
library(RANN, quietly = T)
library(dplyr, quietly = T)

dataset_sc <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC", "Zhengmix4eq", "Zhengmix4uneq","Zhengmix8eq")
dataset_bulk <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")
dropout_percent <- c(0,10,20,30,40,50,60,70,80,90,95)

# function
fn_zarb <- function(knn_list) {
   increment = 0
   for (i in 1:nrow(knn_list)) {
      for (j in 1:ncol(knn_list)) {
         if (i %in% knn_list[knn_list[i,j],]) {
            increment <- increment + 1 
         }
      }
   }
   return(100*(1-increment/(nrow(knn_list)*ncol(knn_list))))
}
get_hub7 <- function(data,kval,pcval) { #asymmetry
   knn_list <- lapply(pcval,
                             function(x) lapply(kval,
                                                       function(y) knn.covertree::find_knn(t(data[1:x,]),y)$index))
   result <- lapply(knn_list,
                           function(x) sapply(x,
                                              function(y) fn_zarb(y)))
   asy_val <- data.frame("k"=rep(kval, times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=unlist(result))
   asy_val$Dimension <- factor(asy_val$Dimension, levels = pcval)
   return(asy_val)
}

# load data
data_in <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/Duo_10kHVG/data/"
data_sc <- pblapply(dataset_sc,
                        function(x) read.table(paste0(data_in,x,".csv")))
data_bulk <- pblapply(dataset_bulk,
       function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/simulated/simul_dropout", y,".csv"))))
data_bulk_splatter <- pblapply(dataset_bulk,
                      function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_simu/simul_dropout", y,".csv"))))
sdev_sc <- pblapply(dataset_sc, function(x) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/Duo_10kHVG/pca/sdev_",x,".csv")))
sdev_bulk <- pblapply(dataset_bulk, function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/pca/simul_dropout", y,"_sdev.csv"))))
sdev_bulk_splatter <- pblapply(dataset_bulk, function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_pca/simul_dropout", y,"_sdev.csv"))))
pca_sc <- pblapply(dataset_sc,
                   function(x) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data/Duo_10kHVG/pca/",x,"_pca_readyforhubness.csv")))
pca_bulk <- pblapply(dataset_bulk, function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/pca/simul_dropout", y,"_pca.csv"))))
pca_bulk_splatter <- pblapply(dataset_bulk, function(x) lapply(dropout_percent, function(y) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_pca/simul_dropout", y,"_pca.csv"))))

# Var x Sparsity
VarSpars <- data.frame("Variance"=c(unlist(sapply(data_sc, function(x) sapply(x,var))),
                                    unlist(sapply(data_bulk, function(x) sapply(x[[1]],var))),
                                    unlist(sapply(data_bulk, function(x) sapply(x[[6]],var)))),
                       "Sparsity"=c(unlist(sapply(data_sc, function(x) sapply(x, function(y) mean(y==0)))),
                                    unlist(sapply(data_bulk, function(x) sapply(x[[1]], function(y) mean(y==0)))),
                                    unlist(sapply(data_bulk, function(x) sapply(x[[6]], function(y) mean(y==0))))),
                       "Dataset"=c(unlist(mapply(x=dataset_sc,y=sapply(data_sc,ncol), function(x,y) rep(x,y))),
                                   unlist(sapply(seq(data_bulk), function(x) rep(dataset_bulk[x], ncol(data_bulk_splatter[[x]][[1]])))),
                                   unlist(sapply(seq(data_bulk), function(x) rep(paste(dataset_bulk[x],"splatter"), ncol(data_bulk_splatter[[x]][[1]]))))))
ggplot(VarSpars[-c(grep("ARCH", VarSpars$Dataset),
                  grep("TCGA", VarSpars$Dataset)),], aes(Variance, color=Dataset)) + geom_density() + scale_y_log10()
ggplot(VarSpars[-c(grep("TCC", VarSpars$Dataset),
                   grep("Sim", VarSpars$Dataset)),], aes(x=Sparsity, y=Variance, color=Dataset)) + geom_point() + scale_y_log10()
PCA_proj <- data.frame("PC1"=c(unlist(sapply(pca_sc, function(x) x[1,])),
                               unlist(sapply(pca_bulk, function(x) x[[1]][1,]))),
                       "PC2"=c(unlist(sapply(pca_sc, function(x) x[2,])),
                               unlist(sapply(pca_bulk, function(x) x[[1]][2,]))),
                       "Dataset"=c(unlist(mapply(x=dataset_sc,y=sapply(data_sc,ncol), function(x,y) rep(x,y))),
                                   unlist(sapply(seq(data_bulk), function(x) rep(dataset_bulk[x], ncol(data_bulk_splatter[[x]][[1]]))))))
ggplot(PCA_proj[-c(grep("TCC", PCA_proj$Dataset),
                   grep("Sim", PCA_proj$Dataset)),], aes(x=PC1, y=PC2, color=Dataset)) + geom_point()

# estimate ID of each dataset
# Screeplot
sdev_df_sc <- data.frame("Sdev_"=unlist(lapply(sdev_sc, function(x)
   unlist(x)/unlist(x)[1])),
   "Dataset"=unname(unlist(mapply(x=dataset_sc, y=sdev_sc, function(x,y)
      rep(x, nrow(y)), SIMPLIFY = T))),
   "PC"=unlist(sapply(sdev_sc, function(x)
      seq(nrow(x)))))
sdev_df_sc <- sdev_df_sc[sdev_df_sc$PC<30,]
sdev_df_bulk <- lapply(sdev_bulk, function(z) {data.frame("Sdev_"=unname(unlist(lapply(z, function(x)
   unlist(x)/unlist(x)[1]))),
   "Dropout"=unname(c(mapply(x=dropout_percent, y=z, function(x,y)
      rep(x, nrow(y)), SIMPLIFY = T))),
   "PC"=c(sapply(z, function(x)
      seq(nrow(x)))))})
sdev_df_bulk <- lapply(sdev_df_bulk, function(x) x[x$PC<30,])
sdev_df_bulk_splatter <- lapply(sdev_bulk_splatter, function(z) {data.frame("Sdev_"=unname(unlist(lapply(z, function(x)
   unlist(x)/unlist(x)[1]))),
   "Dropout"=unname(c(mapply(x=dropout_percent, y=z, function(x,y)
      rep(x, nrow(y)), SIMPLIFY = T))),
   "PC"=c(sapply(z, function(x)
      seq(nrow(x)))))})
sdev_df_bulk_splatter <- lapply(sdev_df_bulk_splatter, function(x) x[x$PC<30,])
nrow <- sapply(sdev_df_bulk,nrow)
sdev_df_bulk <- do.call(rbind,sdev_df_bulk)
sdev_df_bulk$Dataset <- c(sapply(seq(nrow), function(x) rep(dataset_bulk[x],nrow[[x]])))
nrow <- sapply(sdev_df_bulk_splatter,nrow)
sdev_df_bulk_splatter <- do.call(rbind,sdev_df_bulk_splatter)
sdev_df_bulk_splatter$Dataset <- c(sapply(seq(nrow), function(x) rep(dataset_bulk[x],nrow[[x]])))
ggplot(sdev_df_sc[sdev_df_sc$PC<30,], aes(x=PC, y=Sdev_, color=Dataset)) +
   geom_point() +
   geom_line(aes(group=Dataset)) +
   geom_hline(yintercept=0.1) +
   scale_color_viridis_d(option="plasma") +
   ylab("Sd/Sd1") +
   xlab("Principal Components")
ggplot(sdev_df_bulk[sdev_df_bulk$PC<30,], aes(x=PC, y=Sdev_, color=Dropout)) +
   geom_point() +
   geom_line(aes(group=Dropout)) +
   geom_hline(yintercept=0.1) +
   facet_wrap(~Dataset) +
   scale_color_viridis_c() +
   ylab("Sd/Sd1") +
   xlab("Principal Components")
ggplot(sdev_df_bulk_splatter[sdev_df_bulk_splatter$PC<30,], aes(x=PC, y=Sdev_, color=Dropout)) +
   geom_point() +
   geom_line(aes(group=Dropout)) +
   geom_hline(yintercept=0.1) +
   facet_wrap(~Dataset) +
   scale_color_viridis_c() +
   ylab("Sd/Sd1") +
   xlab("Principal Components")

# Make a subsample of Zhengmix4eq and Trapnell
idx=11
sample100 <- sort(sample(seq(ncol(data_sc[[idx]])), 100))
sample50 <- sort(sample(seq(ncol(data_sc[[idx]])), 50))
sdev100 <- prcomp(log10(data_sc[[idx]][,sample100]+1), center=T, scale.=F)$sdev
sdev50 <- prcomp(log10(data_sc[[idx]][,sample50]+1), center=T, scale.=F)$sdev
zheng_sampled <- data.frame("GID_pca"=c(which(sdev100<=sdev100[1]/10)[1],
                                        which(sdev50<=sdev50[1]/10)[1]),
                            "Dataset"=paste0("Zhengmix4uneq_sampled",c(100,50)),
                            #"Dataset"=rep("Zhengmix4uneq",2),
                            "Sparsity"=c(mean(data_sc[[idx]][,sample100]==0)*100,
                                         mean(data_sc[[idx]][,sample50]==0)*100),
                            "Cardinality"=c(100,50),
                            "Median_eigenvalue"=c(median(sdev100),median(sdev50)),
                            "SNR_proxy"=c(max(sdev100)/median(sdev100),
                                             max(sdev50)/median(sdev50)),
                            "Type"="sc")
idx=8
sample100 <- sort(sample(seq(ncol(data_sc[[idx]])), 100))
sample50 <- sort(sample(seq(ncol(data_sc[[idx]])), 50))
sdev100 <- prcomp(log10(data_sc[[idx]][,sample100]+1), center=T, scale.=F)$sdev
sdev50 <- prcomp(log10(data_sc[[idx]][,sample50]+1), center=T, scale.=F)$sdev
trapnell_sampled <- data.frame("GID_pca"=c(which(sdev100<=sdev100[1]/10)[1],
                                        which(sdev50<=sdev50[1]/10)[1]),
                            "Dataset"=paste0("Trapnell_sampled",c(100,50)),
                            #"Dataset"=rep("Trapnell",2),
                            "Sparsity"=c(mean(data_sc[[idx]][,sample100]==0)*100,
                                         mean(data_sc[[idx]][,sample50]==0)*100),
                            "Cardinality"=c(100,50),
                            "Median_eigenvalue"=c(median(sdev100),median(sdev50)),
                            "SNR_proxy"=c(max(sdev100)/median(sdev100),
                                          max(sdev50)/median(sdev50)),
                            "Type"="sc")

# ID x Sparsity with bulk and sc
# Sparsity
sparsity_sc <- sapply(data_sc, function(x)
   mean(x==0)*100)
sparsity_bulk <- pblapply(data_bulk, function(x)
   sapply(x, function(y) mean(y==0)*100))
sparsity_bulk_splatter <- pblapply(data_bulk_splatter, function(x)
   sapply(x, function(y) mean(y==0)*100))
id_sc <- data.frame("GID_pca"=unlist(lapply(sdev_sc, function(x)  return(which(unlist(x)<=unlist(x)[1]/10)[1]))),
                         "Dataset"=dataset_sc,
                    "Sparsity"=sparsity_sc,
                    "Cardinality"=sapply(data_sc,ncol),
                    "Median_eigenvalue"=sapply(sdev_sc,function(x) median(x[,1])),
                    "SNR_proxy"=sapply(sdev_sc,max)/sapply(sdev_sc,function(x) median(x[,1])),
                    "Type"="sc")
df <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/slurm/hubness_param_ggplot1.rds")
df <- df[df$Method=="Asymmetry" & df$Dimension==100,]
df_asym <- inner_join(id_sc, df, by="Dataset")
df_asym$Interplay <- df_asym$Sparsity/df_asym$SNR_proxy
df_asym$Interplay2 <- df_asym$Sparsity/df_asym$SNR_proxy/df_asym$Cardinality
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Interplay", y = "GHubness", color="Dataset",
          add = "reg.line", #3.3x5
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "GID_pca", y = "GHubness", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Interplay", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Sparsity", y = "SNR_proxy", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Cardinality", y = "SNR_proxy", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Sparsity", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "SNR_proxy", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")
ggscatter(df_asym[-c(grep("Sim",df_asym$Dataset)),], x = "Cardinality", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman")

id_bulk <- do.call(rbind,lapply(seq(sdev_bulk), function(y) data.frame("GID_pca"=unname(unlist(lapply(sdev_bulk[[y]], function(x)  return(which(unlist(x)<=unlist(x)[1]/10)[1])))),
                    "Dataset"=rep(dataset_bulk[y],length(dropout_percent)),
                    "Sparsity"=sparsity_bulk[[y]],
                    "Cardinality"=rep(ncol(data_bulk[[y]][[1]]),length(dropout_percent)),
                    "Median_eigenvalue"=sapply(sdev_bulk[[y]],function(x) median(x[,1])),
                    "SNR_proxy"=sapply(sdev_bulk[[y]],max)/sapply(sdev_bulk[[y]],function(x) median(x[,1])),
                    "Type"="bulk")))
df <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_ggplot1.rds")
df <- do.call(rbind,lapply(seq(df), function(x) {df[[x]]$Dataset<-dataset_bulk[x];return(df[[x]])}))
df <- df[df$Method=="Asymmetry" & df$Dimension==100,]
df_asym <- cbind(id_bulk, df[,-c(grep("Dataset",colnames(df)))])
df_asym$Interplay <- df_asym$Sparsity/df_asym$SNR_proxy
df_asym$Interplay2 <- df_asym$Sparsity/df_asym$SNR_proxy/df_asym$Cardinality
ggscatter(df_asym, x = "Sparsity", y = "GHubness", color="Dataset",
          add = "reg.line", # 3x5
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")
ggscatter(df_asym, x = "Sparsity", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")
ggscatter(df_asym, x = "GID_pca", y = "GHubness", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")

ggplot(id_bulk, aes(x=Sparsity, y=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point() +
   scale_y_log10() +
   scale_x_log10()
id_bulk_splatter <- do.call(rbind,lapply(seq(sdev_bulk_splatter), function(y) data.frame("GID_pca"=unname(unlist(lapply(sdev_bulk_splatter[[y]], function(x)  return(which(unlist(x)<=unlist(x)[1]/10)[1])))),
                                                                       "Dataset"=paste0(rep(dataset_bulk[y],length(dropout_percent)),"_splatter"),
                                                                       "Sparsity"=sparsity_bulk_splatter[[y]],
                                                                       "Cardinality"=rep(ncol(data_bulk[[y]][[1]]),length(dropout_percent)),
                                                                       "Median_eigenvalue"=sapply(sdev_bulk_splatter[[y]],function(x) median(x[,1])),
                                                                       "SNR_proxy"=sapply(sdev_bulk_splatter[[y]],max)/sapply(sdev_bulk_splatter[[y]],function(x) median(x[,1])),
                                                                       "Type"="bulk")))
ggplot(id_bulk_splatter, aes(x=Sparsity, y=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point()
df <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_splatter_ggplot1.rds")
df <- do.call(rbind,lapply(seq(df), function(x) {df[[x]]$Dataset<-dataset_bulk[x];return(df[[x]])}))
df <- df[df$Method=="Asymmetry" & df$Dimension==100,]
df_asym <- cbind(id_bulk_splatter, df[,-c(grep("Dataset",colnames(df)))])
df_asym$Interplay <- df_asym$Sparsity/df_asym$SNR_proxy
df_asym$Interplay2 <- df_asym$Sparsity/df_asym$SNR_proxy/df_asym$Cardinality
ggscatter(df_asym, x = "Sparsity", y = "GHubness", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")
ggscatter(df_asym, x = "Interplay", y = "GHubness", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")
ggscatter(df_asym, x = "Interplay", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")
ggscatter(df_asym, x = "GID_pca", y = "GHubness", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   stat_cor(method="spearman")


id_tot <- rbind(rbind(rbind(rbind(id_bulk,id_bulk_splatter),id_sc),zheng_sampled),trapnell_sampled)
id_tot$Dataset <- factor(id_tot$Dataset, levels=c(sort(c(dataset_bulk,paste0(dataset_bulk,"_splatter"))),sort(unique(c(id_sc$Dataset,zheng_sampled$Dataset,trapnell_sampled$Dataset)))))
ggplot(id_tot[-c(grep("TCC",id_tot$Dataset),
                 grep("Sim",id_tot$Dataset)),], aes(x=Sparsity, y=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
  geom_point(aes(size=Cardinality)) +
  ylab("Intrinsic Dimension") +
   scale_y_log10()
ggplot(id_tot[-c(grep("TCC",id_tot$Dataset),
                 grep("Sim",id_tot$Dataset)),], aes(x=Sparsity, y=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point(aes(size=Cardinality, alpha=SNR_proxy)) +
   ylab("Intrinsic Dimension") +
   scale_y_log10() +
   scale_alpha(trans="log")
ggplot(id_tot, aes(x=Sparsity, y=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point(aes(size=Cardinality, alpha=SNR_proxy)) +
   ylab("Intrinsic Dimension") +
   scale_y_log10() +
   scale_alpha(trans="log") +
   scale_color_discrete_divergingx(palette = "Spectral") +
   geom_hline(yintercept = 25, linetype='longdash') +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))
ggplot(rbind(id_tot[id_tot$Type!="bulk",][-c(grep("_sampled", id_tot$Dataset[id_tot$Type!="bulk"])),],id_tot[-c(grep("splatter", id_tot$Dataset)),] %>% filter(Type=="bulk") %>% group_by(Dataset) %>% filter(Sparsity==min(Sparsity))),
       aes(x=Sparsity, y=GID_pca, color=Dataset)) + # Same but with slightly less info (no added dropout for the bulk, remove the splatter)
   geom_point(aes(size=Cardinality, alpha=SNR_proxy)) +
   ylab("Intrinsic Dimension") +
   scale_y_log10() +
   scale_alpha(trans="log") +
   scale_color_discrete_divergingx(palette = "Spectral") +
   #geom_hline(yintercept = 25, linetype='longdash') +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))
# add side panel
resolution=10
setting_id_test=rbind(id_tot[id_tot$Type!="bulk",][-c(grep("_sampled", id_tot$Dataset[id_tot$Type!="bulk"])),],id_tot[-c(grep("splatter", id_tot$Dataset)),] %>% filter(Type=="bulk") %>% group_by(Dataset) %>% filter(Sparsity==min(Sparsity)))
brks <- sapply(with(setting_id_test, seq(min(log10(GID_pca), na.rm=T), max(log10(GID_pca), na.rm=T), length.out = resolution)), function(x)
   exp(x)^log(10))
grps <- gsub("\\(","",with(setting_id_test, cut(GID_pca, breaks = brks, include.lowest = TRUE)))
grps <- gsub("\\[","",grps)
grps <- gsub("\\]","",grps)
lower_bound <- as.numeric(unname(sapply(grps, function(x) if (is.na(x)) {return(x)}
               else {as.numeric(strsplit(x, ",")[[1]][1])})))
upper_bound <- as.numeric(unname(sapply(grps, function(x) if (is.na(x)) {return(x)}
                                        else {as.numeric(strsplit(x, ",")[[1]][2])})))
mean_snr <- sapply(sort(unique(lower_bound[!is.na(lower_bound)])), function(x)
   mean(setting_id_test$SNR_proxy[which(lower_bound==x)]))
mean_cardinality <- sapply(sort(unique(lower_bound[!is.na(lower_bound)])), function(x)
   mean(setting_id_test$Cardinality[which(lower_bound==x)]))
mean_sparsity <- sapply(sort(unique(lower_bound[!is.na(lower_bound)])), function(x)
   mean(setting_id_test$Sparsity[which(lower_bound==x)]))
side_panel <- data.frame("SNR_proxy"=mean_snr,
                         "Cardinality"=mean_cardinality,
                         "Sparsity"=mean_sparsity,
                         "GID_pca"=paste0("[",sort(unique(lower_bound[!is.na(lower_bound)])),
                                          ",",sort(unique(upper_bound[!is.na(upper_bound)])),
                                          "]"),
                         "Interval"=seq(length(unique(lower_bound[!is.na(lower_bound)]))),
                         "Interplay"=mean_snr*mean_cardinality,
                         "Xmin"=unique(lower_bound[!is.na(lower_bound)]),
                         "Xmax"=unique(upper_bound[!is.na(upper_bound)]))

 p1=ggplot(id_tot, aes(x=Sparsity, y=GID_pca, color=Dataset)) + # 7.5x10
   geom_line(aes(group=Dataset)) +
   geom_point(aes(size=Cardinality, alpha=SNR_proxy)) +
   ylab("Intrinsic Dimension") +
   scale_y_log10() +
   scale_alpha(trans="log") +
   scale_color_discrete_divergingx(palette = "Spectral") +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80")) +
   geom_hline(yintercept = brks[1], color="grey60") +
   geom_hline(yintercept = brks[2], color="grey60") +
   geom_hline(yintercept = brks[3], color="grey60") +
   geom_hline(yintercept = brks[4], color="grey60") +
   geom_hline(yintercept = brks[5], color="grey60") +
   geom_hline(yintercept = brks[6], color="grey60") +
   geom_hline(yintercept = brks[7], color="grey60") +
   geom_hline(yintercept = brks[8], color="grey60") +
   geom_hline(yintercept = brks[9], color="grey60") +
   geom_hline(yintercept = brks[10], color="grey60")
p2=ggplot(side_panel, aes(y=Interval, fill=SNR_proxy)) +
   geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=0,xmax=1,fill=SNR_proxy)) +
   scale_fill_continuous_divergingx(palette='Geyser') +
   theme_void()
p3=ggplot(side_panel, aes(y=Interval, fill=Cardinality)) +
   geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=0,xmax=1,fill=Cardinality)) +
   scale_fill_continuous_divergingx(palette='Geyser') +
   theme_void()
p4=ggplot(side_panel, aes(y=Interval, fill=Interplay)) +
   geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=0,xmax=1,fill=Interplay)) +
   scale_fill_continuous_divergingx(palette="Geyser") +
   theme_void()
grid.arrange(p2,p3,p4,nrow=1) # 7.5 x 3.5
p1  # 7.5x12

# Regimb alternative HUM!!!
#ggplot(id_tot, aes(color=Dataset, size=Cardinality)) + # 7.5x10
#   geom_point(aes(x=GID_pca, y=100-Sparsity), shape=1) +
#   geom_point(aes(x=GID_pca, y=SNR_proxy*100/80), shape=15) +
#   geom_line(aes(x=GID_pca, y=100-Sparsity, group=Dataset)) +
#   geom_line(aes(x=GID_pca, y=SNR_proxy, group=Dataset)) +
#   xlab("Intrinsic Dimension") +
#   scale_y_continuous(name="Sparsity", sec.axis=sec_axis(~./100*80, name="SNR")) +
#   scale_x_log10() +
#   scale_color_discrete_divergingx(palette = "Spectral") +
#   theme(panel.background = element_rect(fill="grey98"),
#         panel.grid = element_line(colour = "grey80"))
#ggplot(side_panel) +
#   geom_rect(aes(xmin=Xmin, xmax=Xmax, ymin=SNR_proxy-1, ymax=SNR_proxy+1, color=SNR_proxy),
#             alpha=0.5, fill="white") +
#   geom_rect(aes(xmin=Xmin, xmax=Xmax, ymin=101, ymax=103, fill=Interplay),
#             alpha=0.5) +
#   scale_x_log10() +
#   scale_fill_gradient(low="white", high="black") +
#   scale_color_gradient(low="white", high="black")

ggplot(id_tot, aes(y=SNR_proxy, x=GID_pca, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point() +
   scale_x_log10() +
   scale_y_log10() +
   geom_vline(xintercept = brks[1]) +
   geom_vline(xintercept = brks[2]) +
   geom_vline(xintercept = brks[3]) +
   geom_vline(xintercept = brks[4]) +
   geom_vline(xintercept = brks[5]) +
   geom_vline(xintercept = brks[6]) +
   geom_vline(xintercept = brks[7]) +
   geom_vline(xintercept = brks[8]) +
   geom_vline(xintercept = brks[9]) +
   geom_vline(xintercept = brks[10])
ggplot(id_tot, aes(y=SNR_proxy, x=Sparsity, color=Dataset)) +
   geom_line(aes(group=Dataset)) +
   geom_point() +
   scale_y_log10()

# Make two figures : 1 scatter plot, 1 heatmap
ggplot(side_panel, aes(y=Interval)) +
   geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=0,xmax=1,fill=Sparsity)) +
   #geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=1,xmax=2,fill=SNR_proxy)) +
   #geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=2,xmax=3,fill=Cardinality)) +
   #geom_rect(aes(ymin=Interval-0.5,ymax=Interval+0.5,xmin=3,xmax=4,fill=Interplay)) +
   scale_fill_continuous_divergingx(palette="Geyser") +
   theme_void()


# Get the total phenotypic volume for each data set (Azizi et al)
# do it at 200D
# pick randomly 200 cells for each dataset (do it 10 times)
n_iter=10
size_sample=200
dim=100
pheno_vol_sc <- lapply(pca_sc, function(x) {
   data<-x[seq(dim),]
   pheno_vol <- rep(0,n_iter)
   for (i in seq(n_iter)) {
      cells <- sample(seq(ncol(data)), size_sample)
      eigenvalue <- prcomp(t(data[,cells]))$sdev
      pheno_vol[i] <- prod(eigenvalue[eigenvalue>1e-10])
   }
   return(pheno_vol)
   })
pheno_vol_bulk <- lapply(pca_bulk, function(y) lapply(y, function(x){
   data<-x[seq(dim),]
   pheno_vol <- rep(0,n_iter)
   for (i in seq(n_iter)) {
      cells <- sample(seq(ncol(data)), size_sample)
      eigenvalue <- prcomp(t(data[,cells]))$sdev
      pheno_vol[i] <- prod(eigenvalue[eigenvalue>1e-10])
   }
   return(pheno_vol)
}))
pheno_vol_bulk_splatter <- lapply(pca_bulk_splatter, function(y) lapply(y, function(x){
   data<-x[seq(dim),]
   pheno_vol <- rep(0,n_iter)
   for (i in seq(n_iter)) {
      cells <- sample(seq(ncol(data)), size_sample)
      eigenvalue <- prcomp(t(data[,cells]))$sdev
      pheno_vol[i] <- prod(eigenvalue[eigenvalue>1e-10])
   }
   return(pheno_vol)
}))
# Add asymmetry info
asym_sc <- c(pbsapply(pca_sc, function(x) rep(get_hub7(x, 10, dim)$GHubness,n_iter)))
asym_bulk <- c(sapply(pca_bulk, function(y) c(pbsapply(y, function(x) rep(get_hub7(x, 10, dim)$GHubness,n_iter)))))
asym_bulk_splatter <- c(sapply(pca_bulk_splatter, function(y) c(pbsapply(y, function(x) rep(get_hub7(x, 10, dim)$GHubness,n_iter)))))

pheno_df <- data.frame("PhenotypicVolume"=c(unlist(pheno_vol_sc),unlist(pheno_vol_bulk),unlist(pheno_vol_bulk_splatter)),
                       "Dataset"=c(rep(dataset_sc,each=n_iter),
                                   rep(rep(dataset_bulk,each=n_iter*length(dropout_percent)),2)),
                       "Dropout"=c(rep(NA,n_iter*length(data_sc)),
                                   rep(rep(dropout_percent,each=n_iter),2*length(dataset_bulk))),
                       "Method"=c(rep("nothing",n_iter*length(data_sc)),
                                  rep("manual",n_iter*length(dataset_bulk)*length(dropout_percent)),
                                  rep("splatter",n_iter*length(dataset_bulk)*length(dropout_percent))),
                       "Asymmetry"=c(asym_sc,asym_bulk,asym_bulk_splatter))
pheno_df$Dataset <- factor(pheno_df$Dataset, levels=c(sort(dataset_bulk),dataset_sc))
ggplot(pheno_df[pheno_df$PhenotypicVolume>0 & pheno_df$Method!='splatter',], aes(x=Dataset, y=PhenotypicVolume, fill=paste(Dataset,Dropout))) +
   geom_boxplot() +
   theme(axis.text.x = element_text(angle=90)) +
   scale_y_log10()
ggplot(pheno_df[pheno_df$PhenotypicVolume>0 & pheno_df$Method!='splatter',], aes(x=Asymmetry, y=PhenotypicVolume, color=paste(Dataset,Dropout))) +
   geom_point() +
   scale_y_log10() +
   scale_x_log10()


#tot_df <- cbind(id_tot[1:78,],data.frame("PhenoVol"=c(c(sapply(pheno_vol_bulk,function(x)
#   sapply(x, mean))),
#   c(sapply(pheno_vol_bulk_splatter,function(x)
#      sapply(x, mean))),
#   sapply(pheno_vol_sc,mean))))
#ggplot(tot_df, aes(x=PhenoVol, y=Cardinality, color=Dataset)) +
#   geom_line(aes(group=Dataset)) +
#   geom_point() +
#   scale_x_log10()

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

# Alternative : 3 scatter plots
id_sc <- data.frame("GID_pca"=unlist(lapply(sdev_sc, function(x)  return(which(unlist(x)<=unlist(x)[1]/10)[1]))),
                    "Dataset"=dataset_sc,
                    "Sparsity"=sparsity_sc,
                    "Cardinality"=sapply(data_sc,ncol),
                    "Median_eigenvalue"=sapply(sdev_sc,function(x) median(x[,1])),
                    "SNR_proxy"=sapply(sdev_sc,max)/sapply(sdev_sc,function(x) median(x[,1])),
                    "Type"="sc",
                    "Sampled"="no")
id_sc <- id_sc[-c(grep("Sim",id_sc$Dataset)),]
id_tot2<-id_sc
for (idx in seq(data_sc)) {
   sample100 <- sort(sample(seq(ncol(data_sc[[idx]])), 100))
   sample50 <- sort(sample(seq(ncol(data_sc[[idx]])), 50))
   sdev100 <- prcomp(log10(data_sc[[idx]][,sample100]+1), center=T, scale.=F)$sdev
   sdev50 <- prcomp(log10(data_sc[[idx]][,sample50]+1), center=T, scale.=F)$sdev
   df_sampled <- data.frame("GID_pca"=c(which(sdev100<=sdev100[1]/10)[1],
                                        which(sdev50<=sdev50[1]/10)[1]),
                            "Dataset"=dataset_sc[idx],
                            "Sparsity"=c(mean(data_sc[[idx]][,sample100]==0)*100,
                                         mean(data_sc[[idx]][,sample50]==0)*100),
                            "Cardinality"=c(100,50),
                            "Median_eigenvalue"=c(median(sdev100),
                                                  median(sdev50)),
                            "SNR_proxy"=c(max(sdev100)/median(sdev100),
                                          max(sdev50)/median(sdev50)),
                            "Type"="sc",
                            "Sampled"="yes")
   id_tot2 <- rbind(id_tot2,df_sampled)
}
id_tot2$Interplay <- 1/id_tot2$SNR_proxy*id_tot2$Sparsity
id_tot2$Interplay2 <- 1/id_tot2$SNR_proxy*id_tot2$Sparsity/id_tot2$Cardinality

ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Sparsity", y = "SNR_proxy", color="Dataset", #3x3
          add = "reg.line",
          add.params = list(color = "gray50", fill = "gray90"),
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # independent parameters
ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Cardinality", y = "SNR_proxy", color="Dataset",
          add = "reg.line",  
          add.params = list(color = "gray50", fill = "gray90"),
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # dependent parameters
ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Sparsity", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"),
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # partial explanation
ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Cardinality", y = "GID_pca", color="Dataset",
          add = "reg.line",  
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # partial explanation but better
ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "SNR_proxy", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # partial explanation but even better

ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Interplay", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE) +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # 0.77
ggscatter(id_tot2[id_tot2$Sampled=="no",], x = "Interplay2", y = "GID_pca", color="Dataset",
          add = "reg.line", 
          add.params = list(color = "gray50", fill = "gray90"), 
          conf.int = TRUE)  +
   scale_color_viridis_d(option="plasma") +
   stat_cor(method="spearman") # not better

# Last try scatter plot ID = f(sparsity) with groups of SNR
resolution=3
id_tot3 <- id_tot2[id_tot2$Sampled=="no",]
brks <- with(id_tot3, seq(min(SNR_proxy, na.rm=T), max(SNR_proxy, na.rm=T), length.out = resolution))
grps <- gsub("\\(","",with(id_tot3, cut(SNR_proxy, breaks = brks, include.lowest = TRUE)))
grps <- gsub("\\[","",grps)
grps <- gsub("\\]","",grps)
lower_bound <- as.numeric(unname(sapply(grps, function(x) if (is.na(x)) {return(x)}
                                        else {as.numeric(strsplit(x, ",")[[1]][1])})))
upper_bound <- as.numeric(unname(sapply(grps, function(x) if (is.na(x)) {return(x)}
                                        else {as.numeric(strsplit(x, ",")[[1]][2])})))
mean_snr_group <- sapply(lower_bound, function(x) which(sort(unique(lower_bound))==x))
id_tot3$SNR_group <- as.character(mean_snr_group)
ggscatter(id_tot3[id_tot3$Sampled=="no",], x = "Sparsity", y = "GID_pca", color="Dataset",
                                    add = "reg.line", palette="uchicago",
                                    conf.int = TRUE)  +
   #scale_color_viridis_d(option="plasma") +
   stat_cor(aes(col=SNR_group), label.x = 3)
