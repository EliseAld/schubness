# libs
library(ggplot2, quietly = T)
library(viridis, quietly = T)
library(ggridges, quietly = T)
library(scales, quietly = T)
library(ggpubr, quietly = T)
library(pbapply)
library(knn.covertree)
library(ggh4x)

# functions
data_load <- function(path) {
   source(file=path)
   return(hubness)
}
dataset <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")
data_hubness <- data_load("/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/hubness_load_write_data_bulk_splatter.R")
data_pca <- lapply(dataset,
                   function(y) pblapply(dropout_percent,
                     function(x) read.table(paste0(
                        "/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/",
                        y,"/splatter_pca/",
                        "simul_dropout",x,
                        "_pca.csv"))))

get_hub1 <- function(df,kval,pcval) { # Maximum hubness score
   max = c()
   ymin = c()
   ymax = c()
   y_increment = 1
   for (d in pcval) {
      for (k in kval) {
         val = df$score[df$Dimension==d & df$k==k]
         max = c(max, max(val))
         ymin = c(ymin,y_increment)
         ymax = c(ymax,y_increment+1)
      }
     y_increment <- y_increment+1
   }
   hub_val <- data.frame("k"=rep(kval,times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=max,
                         "ymin"=ymin,
                         "ymax"=ymax)
   return(hub_val)
}
get_hub2 <- function(df,kval,pcval) { # Using the mean + 3*sd from Tomasev
   x_thd = c()
   hub_nb = c()
   ymin = c()
   ymax = c()
   y_increment = 1
   for (d in pcval) {
      for (k in kval) {
         val = df$score[df$Dimension==d & df$k==k]
         x_thd = c(x_thd, mean(val)+3*sd(val))
         hub_nb = c(hub_nb, sum(val>=(mean(val)+3*sd(val))))
         ymin = c(ymin,y_increment)
         ymax = c(ymax,y_increment+1)
      }
     y_increment <- y_increment+1
   }
   hub_val <- data.frame("k"=rep(kval,times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=hub_nb,
                         "Threshold"=x_thd,
                         "ymin"=ymin,
                         "ymax"=ymax)
   return(hub_val)
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
get_hub4 <- function(df,kval,pcval) { # skewness
   s = c()
   for (d in pcval) {
      for (k in kval) {
            val = df$score[df$Dimension==d & df$k==k]
            s = c(s, mean((val-mean(val))^3)/sd(val)^3)
      }
   }
   skewness_val <- data.frame("k"=rep(kval,times=length(pcval)),
                         "Dimension"=rep(pcval,each=length(kval)),
                         "GHubness"=s)
   skewness_val$Dimension <- factor(skewness_val$Dimension, levels = pcval)
   return(skewness_val)
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
make_merge <- function(hubness,data,kval,pcval,name,size) {
#hub_max <- get_hub1(hubness, kval=kval,pcval=pcval)
#hub_mean <- get_hub2(hubness, kval=kval,pcval=pcval)
hub_k <- get_hub3(hubness, kval=kval,pcval=pcval)
hub_skew <- get_hub4(hubness, kval=kval,pcval=pcval)
hub_anti <- get_hub6(hubness, kval=kval,pcval=pcval)
hub_asym <- get_hub7(data, kval=kval,pcval=pcval)
hub_k <- hub_k[,colnames(hub_skew)]
#hub_mean <- hub_mean[,colnames(hub_skew)]
#hub_max <- hub_max[,colnames(hub_skew)]
hub_k$GHubness <- hub_k$GHubness*100/size
#hub_mean$GHubness <- hub_mean$GHubness*100/size
hub_anti$GHubness <- hub_anti$GHubness*100/size
#hub_max$GHubness <- hub_max$GHubness/size
#hub_max$Method <- "Max"
#hub_mean$Method <- "Mean+Sd"
hub_k$Method <- "2k"
hub_skew$Method <- "Skewness"
hub_anti$Method <- "Antihubs"
hub_asym$Method <- "Asymmetry"
hub <- rbind(hub_skew,hub_k,hub_asym,hub_anti)
hub$Dropout <- name
return(hub)
}

# First figure
# Merge the dfs
pcval=lapply(dim_nb,
               function(x) return(c(2,5,10,20,30,40,50,100,200,500,x-1)))
pcval2=lapply(dim_nb,
               function(x) return(c(2,5,10,20,30,40,50,100,200,500)))
df <- lapply(seq(n_dataset),
             function(y) do.call(rbind,pblapply(seq(length(dropout_percent)),
             function(x) make_merge(data_hubness[[y]][[x]],data_pca[[y]][[x]],10,pcval[[y]],dropout_percent[x],dim_nb[y]))))
df2 <- lapply(seq(n_dataset),
             function(y) do.call(rbind,pblapply(seq(length(dropout_percent)),
             function(x) make_merge(data_hubness[[y]][[x]],data_pca[[y]][[x]],10,pcval2[[y]],dropout_percent[x],dim_nb[y]))))
for (i in seq(n_dataset)) {
   df[[i]]$Dimension <- as.numeric(as.character(df[[i]]$Dimension))
   df2[[i]]$Dimension <- as.numeric(as.character(df2[[i]]$Dimension))
}

# Viz
labz_method = c("2k (% of hubs)", "Antihubs (% of hubs)", "Asymmetry (% of edges)", "Skewness")
names(labz_method) <- c("2k","Antihubs","Asymmetry","Skewness")
ggarrange(plotlist=lapply(df,
                          function(x) ggplot(x, aes(x=Dimension, y=GHubness)) +
                             geom_line(aes(group=Dropout, color=Dropout)) +
                             geom_point(size=0.1) +
                             facet_wrap(~Method, scales = "free", labeller = labeller(Method = labz_method)) + 
                             ylab(NULL) +
                             scale_color_viridis_c()),
          labels=dataset)
ggarrange(plotlist=lapply(df2,
                          function(x) ggplot(x, aes(x=Dimension, y=GHubness)) +
                             geom_line(aes(group=Dropout, color=Dropout)) +
                             geom_point(size=0.1) +
                             facet_wrap(~Method, scales = "free", labeller = labeller(Method = labz_method)) +
                             ylab(NULL) +
                             scale_color_viridis_c()),
          labels=dataset)

df_ <- do.call(rbind, df)
df_$Dataset <- unlist(lapply(seq(dataset), function(x) rep(dataset[x], nrow(df[[x]]))))
df_$facet <- paste(df_$Dataset,df_$Method)
labz_method = rep(c("2k (% of hubs)", "Antihubs (% of hubs)", "Asymmetry (% of edges)", "Skewness"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("2k","Antihubs","Asymmetry","Skewness")))
scales_y <- list(scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,20)))

ggplot(df_, aes(x=Dimension, y=GHubness)) +
   geom_line(aes(group=Dropout, color=Dropout)) +
   geom_point(size=0.1) +
   facet_wrap(~facet, scales="free", labeller = labeller(facet = labz_method), nrow = 2) +
   facetted_pos_scales(y=scales_y) +
   ylab(NULL) +
   scale_color_viridis_c() +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))

df2_ <- do.call(rbind, df2)
df2_$Dataset <- unlist(lapply(seq(dataset), function(x) rep(dataset[x], nrow(df2[[x]]))))
df2_$facet <- paste(df2_$Dataset,df2_$Method)
labz_method = rep(c("2k (% of hubs)", "Antihubs (% of hubs)", "Asymmetry (% of edges)", "Skewness"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("2k","Antihubs","Asymmetry","Skewness")))

ggplot(df2_, aes(x=Dimension, y=GHubness)) +
   geom_line(aes(group=Dropout, color=Dropout)) +
   geom_point(size=0.1) +
   facet_wrap(~facet, scales = "free", labeller = labeller(facet = labz_method), nrow = 2) +
   ylab(NULL) +
   scale_color_viridis_c() +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))

saveRDS(df,
        file="/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/hubness_param_ggplot1_splatter.rds")
saveRDS(df2,
        file="/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/hubness_param_ggplot1bis_splatter.rds")


# Figure supp
make_merge_supp <- function(hubness,data,kval,pcval,name,size) {
hub_max <- get_hub1(hubness, kval=kval,pcval=pcval)
hub_mean <- get_hub2(hubness, kval=kval,pcval=pcval)
hub_skew <- get_hub4(hubness, kval=kval,pcval=pcval)
hub_mean <- hub_mean[,colnames(hub_skew)]
hub_max <- hub_max[,colnames(hub_skew)]
hub_mean$GHubness <- hub_mean$GHubness*100/size
hub_max$GHubness <- hub_max$GHubness/size
hub_max$Method <- "Max"
hub_mean$Method <- "Mean+Sd"
hub <- rbind(hub_max,hub_mean)
hub$Dropout <- name
return(hub)
}
df2_supp <- lapply(seq(n_dataset),
             function(y) do.call(rbind,pblapply(seq(length(dropout_percent)),
             function(x) make_merge_supp(data_hubness[[y]][[x]],data_pca[[y]][[x]],10,pcval2[[y]],dropout_percent[x],dim_nb[y]))))
for (i in seq(n_dataset)) {
   df2_supp[[i]]$Dimension <- as.numeric(as.character(df2_supp[[i]]$Dimension))
}

# Viz
labz_method_supp = c("Mean (% of hubs)", "Maximum hub score")
names(labz_method_supp) <- c("Mean+Sd","Max")
ggarrange(plotlist=lapply(df2_supp,
                          function(x) ggplot(x, aes(x=Dimension, y=GHubness)) +
                             geom_line(aes(group=Dropout, color=Dropout)) +
                             geom_point(size=0.1) +
                             facet_wrap(~Method, scales = "free", labeller = labeller(Method = labz_method_supp)) +
                             ylab(NULL) +
                             scale_color_viridis_c()),
          labels=dataset)

df2_supp_ <- do.call(rbind, df2_supp)
df2_supp_$Dataset <- unlist(lapply(seq(dataset), function(x) rep(dataset[x], nrow(df2_supp[[x]]))))
df2_supp_$facet <- paste(df2_supp_$Dataset,df2_supp_$Method)
labz_method = rep(c("Mean (% of hubs)", "Maximum hub score"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("Mean+Sd","Max")))

ggplot(df2_supp_, aes(x=Dimension, y=GHubness)) +
   geom_line(aes(group=Dropout, color=Dropout)) +
   geom_point(size=0.1) +
   facet_wrap(~facet, scales = "free", labeller = labeller(facet = labz_method), nrow = 2) +
   ylab(NULL) +
   scale_color_viridis_c() +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))

saveRDS(df2_supp,
        file="/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/hubness_param_ggplot2_splatter.rds")

