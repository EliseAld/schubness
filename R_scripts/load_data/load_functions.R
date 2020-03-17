library(ggplot2, quietly = T)
library(viridis, quietly = T)
library(ggpubr, quietly = T)
library(Seurat, quietly = T)
library(cluster, quietly = T)
library(RaceID, quietly = T)
library(BiocNeighbors, quietly = T)

get_hub2 <- function(df,kval,pcval,pval) { # Using the mean + 3*sd from Tomasev
  x_thd = c()
  hub_nb = c()
  ymin = c()
  ymax = c()
  y_increment = 1
  for (d in pcval) {
    for (k in kval) {
      for (p in pval) {
        val = df$score[df$Dimension==d & df$p==p & df$k==k]
        x_thd = c(x_thd, mean(val)+3*sd(val))
        hub_nb = c(hub_nb, sum(val>=(mean(val)+3*sd(val))))
        ymin = c(ymin,y_increment)
        ymax = c(ymax,y_increment+1)
      }
    }
    y_increment <- y_increment+1
  }
  hub_val <- data.frame("p"=rep(pval,times=length(pcval)*length(kval)),
                        "k"=rep(rep(kval,each=length(pval)),times=length(pcval)),
                        "Dimension"=rep(pcval,each=length(kval)*length(pval)),
                        "GHubness"=hub_nb,
                        "Threshold"=x_thd,
                        "ymin"=ymin,
                        "ymax"=ymax)
  return(hub_val)
}
get_hub6 <- function(df,kval,pcval,pval) { # anti hubs
  anti_hub = c()
  for (d in pcval) {
    for (k in kval) {
      for (p in pval) {
        anti_hub = c(anti_hub,unique(df$Antihubs[df$Dimension==d & df$p==p & df$k==k]))
      }
    }
  }
  hub_val <- data.frame("p"=rep(pval,times=length(pcval)*length(kval)),
                        "k"=rep(rep(kval,each=length(pval)),times=length(pcval)),
                        "Dimension"=rep(pcval,each=length(kval)*length(pval)),
                        "GHubness"=anti_hub)
  hub_val$Dimension <- factor(hub_val$Dimension, levels = pcval)
  return(hub_val)
}
get_hub7 <- function(df,kval,pcval,pval) { # median + 3MAD
  x_thd = c()
  hub_nb = c()
  ymin = c()
  ymax = c()
  y_increment = 1
  for (d in pcval) {
    for (k in kval) {
      for (p in pval) {
        val = df$score[df$Dimension==d & df$p==p & df$k==k]
        x_thd = c(x_thd, median(val)+3*mad(val,
                                           #center=ifelse(median(val)!=0,median(val),1)
        ))
        hub_nb = c(hub_nb, sum(val>=(median(val)+3*mad(val,
                                                       #center=ifelse(median(val)!=0,median(val),1)
        ))))
        ymin = c(ymin,y_increment)
        ymax = c(ymax,y_increment+1)
      }
    }
    y_increment <- y_increment+1
  }
  hub_val <- data.frame("p"=rep(pval,times=length(pcval)*length(kval)),
                        "k"=rep(rep(kval,each=length(pval)),times=length(pcval)),
                        "Dimension"=rep(pcval,each=length(kval)*length(pval)),
                        "GHubness"=hub_nb,
                        "Threshold"=x_thd,
                        "ymin"=ymin,
                        "ymax"=ymax)
  return(hub_val)
}

histoplot_hub_spec <- function(seurat, slot, slot_name, id) {
  ggplot(data.frame(slot_name=slot,
                    "hubs"=id),
         aes(x=slot_name, fill=hubs)) +
    geom_histogram()
}
vlnplot_hub_spec <- function(seurat, feature, group.by, ylim_inf, ylim_sup, ypos) {
  gplot <- VlnPlot(seurat, features = feature, group.by = group.by)
  gplot +
    geom_hline(yintercept = mean(seurat@meta.data[,feature][seurat@meta.data[,group.by] != "normal"]),
               col="#F8766D") +
    geom_hline(yintercept = mean(seurat@meta.data[,feature][seurat@meta.data[,group.by] == "normal"]),
               col="#00BFC4") +
    stat_compare_means(method = "t.test", p.adjust.method = "BH",
                       comparisons = list(c(ifelse(group.by=="hub6","antihub","hub"), "normal")),
                       label.y = ypos) +
    ylim(ylim_inf,ylim_sup)
}
vlnplot_hub_spec_connect <- function(seurat, feature, group.by, ylim_inf, ylim_sup, ypos) {
  gplot <- VlnPlot(seurat, features = feature, group.by = group.by)
  gplot +
    geom_hline(yintercept = mean(seurat@meta.data[,feature][seurat@meta.data[,group.by] != "free"]),
               col="#F8766D") +
    geom_hline(yintercept = mean(seurat@meta.data[,feature][seurat@meta.data[,group.by] == "free"]),
               col="#00BFC4") +
    stat_compare_means(method = "t.test", p.adjust.method = "BH",
                       comparisons = list(c("connected", "free")),
                       label.y = ypos) +
    ylim(ylim_inf,ylim_sup)
}

coverage_random_hub <- function(neighbor_list, random_cell_nb) {
  cell_randomized <- sample(rownames(neighbor_list),nrow(neighbor_list))
  random_cells <- cell_randomized[1:random_cell_nb]
  k <- 5
  coverage <- length(unique(unlist(as.list(neighbor_list[random_cells,1:k]))))
  while (coverage[k/5] < cell_nb*0.9 & k < k_max) {
    k <- k+5
    coverage <- c(coverage, length(unique(unlist(as.list(neighbor_list[random_cells,1:k])))))
  }
  return(list(coverage=coverage,k=k))
}
coverage_hub <- function(neighbor_list, hub_scores, hub_nb) {
  cell_scores_sorted <- names(sort(hub_scores, decreasing  =T))
  hubs <- cell_scores_sorted[1:hub_nb]
  k <- 5
  coverage <- length(unique(unlist(as.list(neighbor_list[hubs,1:k]))))
  while (coverage[k/5] < cell_nb*0.9 & k < k_max) {
    k <- k+5
    coverage <- c(coverage, length(unique(unlist(as.list(neighbor_list[hubs,1:k])))))
  }
  coverage_rdm <- coverage_random_hub(neighbor_list, hub_nb)
  plot(seq(5,k,5),coverage, ylim=c(0,length(hub_scores)),
       pch=16,
       xlab="Neighborhood size",
       main=paste0("Coverage of the data with ",hub_nb," hubs"))
  points(seq(5,coverage_rdm$k,5),coverage_rdm$coverage, col="lightblue", pch=3)
  abline(h=length(hub_scores)*0.9, col='red', lty='dashed')
  abline(v=k, col="red", lty="dashed")
  abline(v=coverage_rdm$k, col="lightblue", lty="dotted")
}
coverage_random_k <- function(neighbor_list, k) {
  random_cell_nb <- 1
  cell_randomized <- sample(rownames(neighbor_list),nrow(neighbor_list))
  random_cells <- cell_randomized[1:random_cell_nb]
  coverage <- length(unique(unlist(as.list(neighbor_list[random_cells,1:k]))))
  while (coverage[random_cell_nb] < cell_nb*0.9 & random_cell_nb < 1000) {
    random_cell_nb <- random_cell_nb+1
    random_cells <- cell_randomized[1:random_cell_nb]
    coverage <- c(coverage, length(unique(unlist(as.list(neighbor_list[random_cells,1:k])))))
  }
  return(list(coverage=coverage,random_cell_nb=random_cell_nb))
}
coverage_k <- function(neighbor_list, hub_scores, k) {
  cell_scores_sorted <- names(sort(hub_scores, decreasing  =T))
  hub_nb <- 1
  hubs <- cell_scores_sorted[1:hub_nb]
  coverage <- length(unique(unlist(as.list(neighbor_list[hubs,1:k]))))
  while (coverage[hub_nb] < cell_nb*0.9 & hub_nb < 1000) {
    hub_nb <- hub_nb+1
    hubs <- cell_scores_sorted[1:hub_nb]
    coverage <- c(coverage, length(unique(unlist(as.list(neighbor_list[hubs,1:k])))))
  }
  coverage_rdm <- coverage_random_k(neighbor_list, k)
  plot(1:hub_nb,coverage, ylim=c(0,length(hub_scores)),
       pch=16,
       xlab="Number of hubs",
       main=paste0("Coverage of the data looking at ",k," neighbors"))
  points(1:coverage_rdm$random_cell_nb,coverage_rdm$coverage, col="lightblue", pch=3)
  abline(h=length(hub_scores)*0.9, col='red', lty='dashed')
  abline(v=hub_nb, col="red", lty="dashed")
  abline(v=coverage_rdm$random_cell_nb, col="lightblue", lty="dotted")
}