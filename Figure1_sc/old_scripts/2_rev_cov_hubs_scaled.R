# libs
library(pbapply)
library(RANN)
library(ggplot2)
library(dplyr)

# functions
data_load <- function(path) {
  source(file=path)
  return(hubness)
}
kNN <- function(data,k) {
  k.nn <- nn2(t(data), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  return(k.nn)
}
rev_cov_hub <- function(knn_graph, hub_scores, n) {
  attempted_hubs <- which(order(hub_scores, decreasing = T) %in% seq(n))
  size_cov <- sum(apply(knn_graph,1,function(x)
    any(x %in% attempted_hubs)))
  return(size_cov/length(hub_scores)*100)
}
rev_cov_random <- function(knn_graph, hub_scores, n) {
  attempted_hubs <- sample(seq(nrow(knn_graph)), n)
  size_cov <- sum(apply(knn_graph,1,function(x)
    any(x %in% attempted_hubs)))
  return(size_cov/length(hub_scores)*100)
}
get_plateau <- function(df_rev_cov, thd=1) {
  increment=diff(df_rev_cov[,"Rev_cov_size"])
  return(df_rev_cov[,"N"][ifelse(all(increment>=thd),length(increment)+1,which(increment<thd)[1]+1)])
}

# load data
params <- list(10,50)
names(params) <- c("k","Dimension") 
dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
data_hubness <- data_load("/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/hubness_load_write_data_duo_scaled.R")
rm(hubness_scores_data, hubness)
data_hubness <- lapply(data_hubness, function(x)
  x$score[x$k==params$k & x$Dimension==params$Dimension])
data_in <- "/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/DimRedPaper/DuoClustering2018/sce_full/sce_full_"
data_counts <- pblapply(dataset,
                        function(x) {tmp <- readRDS(paste0(data_in,x,".rds"));
                        tmp <- tmp@assays$data@listData$counts;
                        hvg <- names(sort(apply(tmp,1,var),
                                          decreasing = T)[1:min(1e4,nrow(tmp))]);
                        tmp <- tmp[hvg,]
                        return(tmp)})

# Check the skewness
scores_df <- data.frame("InDegree"=unlist(data_hubness))
scores_df$Dataset <- unname(unlist(mapply(x=dataset, y=data_hubness, function(x,y)
  rep(x, length(y)))))
ggplot(scores_df, aes(InDegree, fill=Dataset)) +
  geom_density(alpha=0.6) +
  scale_fill_viridis_d(option="plasma")

# reverse get all the data fixing k but increasing N
# get the neighborhood
k_neighbor <- pblapply(data_counts, function(x)
  kNN(x, 5))
# get the rev cov with different N values and ctrl
n_val <- lapply(data_counts, function(x)
  seq(from=10, to=ncol(x), by=5))
rev_cov_size <- mapply(SIMPLIFY = F,
                       x=k_neighbor, y=data_hubness, z=n_val, function(x,y,z) {
                         idx=1
                       cov_size=rev_cov_hub(x,y,z[idx])
                       while (cov_size[idx]<100 & idx < length(z)) {
                         idx=idx+1
                         cov_size=c(cov_size,rev_cov_hub(x,y,z[idx]))
                       }
                       return(data.frame("N"=z[seq(idx)],
                                         "Rev_cov_size"=cov_size))
                       })
rev_cov_ctrl <- mapply(SIMPLIFY = F,
                      x=k_neighbor, y=data_hubness, z=n_val, function(x,y,z) {
                        idx=1
                        cov_size=rev_cov_random(x,y,z[idx])
                        while (cov_size[idx]<100 & idx < length(z)) {
                          idx=idx+1
                          cov_size=c(cov_size,rev_cov_random(x,y,z[idx]))
                        }
                        return(data.frame("N"=z[seq(idx)],
                                          "Rev_cov_size"=cov_size))
                      })
rev_cov <- list(rev_cov_size,rev_cov_ctrl)
names(rev_cov) <- c("Hubs", "Negative control")
rev_cov_df <- lapply(rev_cov, function(x) do.call(rbind,x))
rev_cov_df <- mapply(SIMPLIFY=F, z=seq(length(rev_cov)), function(z) {
  rev_cov_df[[z]]$Dataset <- unlist(mapply(x=rev_cov[[z]], y=dataset, function(x,y)
    rep(y, nrow(x))));
  rev_cov_df[[z]]$N_percentage <- unlist(mapply(x=rev_cov[[z]], y=data_counts, function(x,y)
    rep(ncol(y), nrow(x))));
  rev_cov_df[[z]]$N_percentage <- rev_cov_df[[z]]$N/rev_cov_df[[z]]$N_percentage*100;
  rev_cov_df[[z]] <- rev_cov_df[[z]] %>%
    group_by(Dataset) %>%
    mutate(Increment=c(diff(Rev_cov_size),0))
  return(rev_cov_df[[z]])
})
rev_df <- do.call(rbind, rev_cov_df)
rev_df$Setting <- c(rep("Hub",nrow(rev_cov_df[[1]])),
                        rep("Control",nrow(rev_cov_df[[2]])))
ggplot(rev_df, aes(x=N_percentage, y=Rev_cov_size, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  facet_wrap(~Setting) +
  scale_color_viridis_d(option="plasma") +
  xlab("Number of cells (%)") +
  ylab("Size of the reverse coverage (%)") +
  scale_x_log10() +
  scale_y_log10()
ggplot(rev_df, aes(x=N_percentage, y=Increment, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  facet_wrap(~Setting) +
  scale_color_viridis_d(option="plasma") +
  xlab("Number of cells (%)") +
  ylab("Increment in reverse coverage size (%)") +
  scale_x_log10()

# get hubs
hub_nb <- lapply(rev_cov, function(y)
  sapply(y, function(x)
  get_plateau(x, thd=0.1)))
hub_prop <- lapply(hub_nb, function(y)
  mapply(x=y, z=data_counts, function(x,z)
    x/ncol(z)*100))

# Remove hubs
data_counts_skimmed <- mapply(SIMPLIFY = F,
                              x=data_counts, y=hub_nb$Hubs, z=data_hubness, function(x,y,z)
                                x[,-which(order(z, decreasing = T) %in% seq(y))])

# Save data
pbmapply(x=data_counts_skimmed, y=dataset, function(x,y)
  write.csv(x, file=paste0("/Users/elise/Desktop/Github/Hubness_sc/Data/Duo_hub_rm_scaled/",y,".csv")))

# make PCA
data_counts_skimmed <- pblapply(data_counts_skimmed, function(x) log10(1+x))
pca <- pblapply(data_counts_skimmed, function(x) prcomp(x, center=T, scale.=F))
data_proj <- pblapply(pca, function(x) t(x$rotation))
pblapply(1:length(data_proj),
         function(x)
           write.table(data_proj[[x]],
                                 file=paste0("/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/pca_hub_rm_scaled/",dataset[x],"_pca_readyforhubness.csv")))
pblapply(1:length(data_proj),
         function(x)
           write.table(pca[[x]]$sdev,
                       file=paste0("/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/pca_hub_rm_scaled/",dataset[x],".csv")))

# Check the screeplot
sdev <- lapply(pca, function(x) x$sdev/x$sdev[1])
sdev_df <- data.frame("Sdev" = unlist(sdev),
                      "Dataset"=unlist(mapply(x=dataset,y=lengths(sdev), function(x,y)
                        rep(x, y))),
                      "PC"=unlist(sapply(lengths(sdev), seq)))
ggplot(sdev_df[sdev_df$PC<30,], aes(x=PC, y=Sdev, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  geom_hline(yintercept=0.1) +
  scale_color_viridis_d(option="plasma") +
  ylab("") +
  xlab("Principal Components")

# get the scores (from the .R file)

# Evaluate the asymmetry post removal
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
n_dim <- sapply(data_proj, ncol)
pc.val <- lapply(n_dim,function(x)
  c(2,5,10,20,30,40,50,100,x-1))
k.val <- c(5,10,50,100)
asymmetry <- pbmapply(x=data_proj, y=pc.val, function(x,y)
                      get_hub7(x, params$k, y), SIMPLIFY = F)
asymmetry_df <- do.call(rbind,asymmetry)
asymmetry_df$Dataset <- rep(dataset, each=lengths(pc.val[1]))
asymmetry_df$Dimension <- as.numeric(as.character(asymmetry_df$Dimension))
ggplot(asymmetry_df, aes(x=Dimension, y=GHubness, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  scale_color_viridis_d(option="plasma") +
  scale_x_log10() +
  scale_y_log10()
  
  