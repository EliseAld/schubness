# Libs & paths
library(pbapply)
library(ggplot2)
library(ggpubr)
library(rgl)
dropout_percent <- c(0,10,20,30,40,50,60,70,80,90,95)
datasets <- c("TCGA_BRCA","ARCHS4")
data_counts <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_simu/simul_dropout", y,".csv")))
data_pca <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_pca/simul_dropout", y,"_pca.csv")))
data_sdev <- lapply(datasets,
                  function(x) sapply(dropout_percent, function(y) paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/splatter_pca/simul_dropout", y,"_sdev.csv")))

# functions
data_load <- function(path) {
   source(file=path)
   return(hubness)
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
#hub_k <- get_hub3(hubness, kval=kval,pcval=pcval)
hub_skew <- get_hub4(hubness, kval=kval,pcval=pcval)
#hub_anti <- get_hub6(hubness, kval=kval,pcval=pcval)
hub_asym <- get_hub7(data, kval=kval,pcval=pcval)
#hub_k <- hub_k[,colnames(hub_skew)]
#hub_mean <- hub_mean[,colnames(hub_skew)]
#hub_max <- hub_max[,colnames(hub_skew)]
#hub_k$GHubness <- hub_k$GHubness*100/size
#hub_mean$GHubness <- hub_mean$GHubness*100/size
#hub_anti$GHubness <- hub_anti$GHubness*100/size
#hub_max$GHubness <- hub_max$GHubness/size
#hub_max$Method <- "Max"
#hub_mean$Method <- "Mean+Sd"
#hub_k$Method <- "2k"
hub_skew$Method <- "Skewness"
#hub_anti$Method <- "Antihubs"
hub_asym$Method <- "Asymmetry"
hub <- rbind(hub_skew,hub_asym)
hub$Dropout <- name
return(hub)
}

# load data
counts <- lapply(data_counts, function(x) pblapply(x,function(y) read.table(file=y)))
pca <- lapply(data_pca, function(x) pblapply(x,function(y) read.table(file=y)))
sdev <- lapply(data_sdev, function(x) pblapply(x,function(y) read.table(file=y)))

dataset <- c("TCGA_BRCA","ARCHS4")
data_hubness <- data_load("/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/hubness_load_write_data_bulk_splatter.R")

# plot sdev and ID=f(Dropout)
sdev <- lapply(sdev, function(y) lapply(y, function(x) return(unlist(x)/unlist(x)[1])))
sdev_df <- lapply(sdev,
                  function(y) data.frame("Sdev" = unlist(y),
                      "Dropout"=as.factor(rep(dropout_percent,each=length(y[[1]]))),
                      "PC"=1:length(y[[1]])))
pdf(file="/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_supp/splatter_screeplot.pdf")
ggarrange(ggplot(sdev_df[[1]][sdev_df[[1]]$PC<30,], aes(x=PC, y=Sdev, color=Dropout)) +
                              geom_point() +
                              geom_line(aes(group=Dropout)) +
                              geom_hline(yintercept=0.1) +
                              scale_color_viridis_d() +
                              ylab("Sd/Sd1") +
            xlab("Principal Components"),
          ggplot(sdev_df[[2]][sdev_df[[2]]$PC<30,], aes(x=PC, y=Sdev, color=Dropout)) +
                              geom_point() +
                              geom_line(aes(group=Dropout)) +
                              geom_hline(yintercept=0.1) +
                              scale_color_viridis_d() +
                              ylab("") +
            xlab("Principal Components"),
          labels = datasets,
          common.legend = T)
dev.off()

pdf(file="/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_supp/splatter_IDxDropout.pdf")
linear_dim_df <- data.frame("GID_pca"=unlist(lapply(sdev, function(y) sapply(y, function(x) return(which(x<=0.1)[1])))),
                         "Dataset"=rep(datasets, each=length(dropout_percent)),
                         "Dropout"=dropout_percent)
ggplot(linear_dim_df, aes(x=Dropout, y=GID_pca)) +
  geom_point() +
  geom_line(aes(group=Dataset, color=Dataset)) +
  ylab("Intrinsic Dimension") +
  scale_x_log10()
dev.off()

# plot dropout, hubness and ID at the same time
n_methods_shown = 2
pcval2=lapply(dim_nb,
               function(x) return(c(2,5,10,20,30,40,50,100,200)))
df <- lapply(seq(n_dataset),
             function(y) do.call(rbind, pblapply(seq(length(dropout_percent)),
             function(x) make_merge(data_hubness[[y]][[x]],pca[[y]][[x]],10,pcval2[[y]],dropout_percent[x],dim_nb[y]))))
for (i in seq(n_dataset)) {
   df[[i]]$Dimension <- as.numeric(as.character(df[[i]]$Dimension))
   df[[i]]$ID <- rep(sapply(sdev[[i]], function(x) return(which(unlist(x)/unlist(x)[1]<=0.1)[1])), each=n_methods_shown*length(pcval2[[i]]))
   df[[i]]$ID2 <- mapply(function(x,y) min(df[[i]]$Dimension[y], ifelse(is.na(x), dim_nb[i], x)), df[[i]]$ID, seq(nrow(df[[i]])))
}

# Viz
labz_method = c("Asymmetry (% of edges)", "Skewness")
names(labz_method) <- c("Asymmetry","Skewness")
#pdf(file="/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_supp/IDxHubness.pdf")
#ggarrange(plotlist=lapply(df,
#                          function(x) ggplot(x, aes(x=Dimension, y=GHubness)) +
#                             geom_line(aes(group=Dropout, color=ID2)) +
#                             geom_point(size=0.1) +
#                             facet_wrap(~Method, scales = "free", labeller = labeller(Method = labz_method)) +
#                             ylab(NULL) +
#                             #geom_vline(xintercept=x$Dropout[is.na(x$ID)][1]) +
#                             scale_color_viridis_c()),
#          labels=dataset)
#dev.off()

pdf(file="/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_supp/splatter_DropoutxHubness.pdf")
ggarrange(plotlist=lapply(df,
                          function(x) ggplot(x, aes(x=Dimension, y=GHubness)) +
                             geom_line(aes(group=Dropout, color=Dropout)) +
                             geom_point(size=0.1) +
                             facet_wrap(~Method, scales = "free", labeller = labeller(Method = labz_method)) +
                             ylab(NULL) +
                             #geom_vline(xintercept=x$Dropout[is.na(x$ID)][1]) +
                             scale_color_viridis_c()),
          labels=dataset)
dev.off()

brks <- with(df[[i]][df[[i]]$Method=="Asymmetry",], seq(min(Dropout), max(Dropout), length.out = length(dropout_percent)))
grps <- with(df[[i]][df[[i]]$Method=="Asymmetry",], cut(Dropout, breaks = brks, include.lowest = TRUE))
increment=1
idx=1
colorby = rep(0,length(grps))
while (idx <= length(grps)-1) {
  change = which(grps[idx+1:(length(grps)-1)] != grps[idx])[1]+idx
  colorby[idx:change] <- increment
  idx=change+1
  increment=increment+1
}
colorby[colorby==0] <- max(colorby)+1
for (i in seq(n_dataset)) {
   for (method in c("Asymmetry", "Skewness")) {
      d3_df = df[[i]][df[[i]]$Method==method,]
      d3_df$colorby = factor(colorby)
      open3d()
      par3d()
      with(d3_df, plot3d(Dimension, ID2, GHubness, col = colorby,
                         xlab = "Embedding dimension", ylab = "Intrinsic dimension", zlab = method))
      legend3d("topright", legend = unique(d3_df$Dropout), col = levels(d3_df$colorby), pch=18)
      for (dropout in unique(d3_df$Dropout)) {
         with(d3_df[d3_df$Dropout==dropout,], lines3d(Dimension, ID2, GHubness, col = colorby))
      }
      title3d(paste0(method, " of ", dataset[i]))
      rgl.postscript(paste0("splatter_3DHubnessxIDxDropout_",dataset[i],"_",method,".pdf"),"pdf") 
      movie3d(spin3d(), dir="/Users/elise/Desktop/GitHub/Hubness_sc/Figure3_supp/", clean=T, duration=15, movie=paste0("splatter_3DHubnessxIDxDropout_",dataset[i],"_",method))
   }
}


