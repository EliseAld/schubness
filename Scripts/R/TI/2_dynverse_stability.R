#' @param path path to .h5ad python PAGA output
#' @param dataset_path path to .h5seurat R PAGA object
#' @param stab_path path to .csv python PAGA output
#' @param metric_choice name of quality metrics used (default to c("correlation","F1_branches"))
#' @param clustering clustering algorithm used to perform PAGA

#' @return stability quality metrics

#' @example stability scores <- stab_score_traj(dataset_path, stab_path, clustering, metric_choice)

### Load librairies
library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(dyno)
library(assertthat)
library(testthat)
library(glue)
library(dyneval)
library(BBmisc)
library(Xmisc)
library(rstatix)
source(file = "~/dynverse_common_functions.R")

### Function to compute stability of the two quality scores (Corr, F1)
# Convert the python PAGA output .h5ad to a R-readable format
convert_PtoR <- function(path) {
  if (file.exists(paste0(path))) {
    SeuratDisk::Convert(path, dest = "h5seurat", overwrite = TRUE, verbose = F)
    #file.remove(paste0(path))
  }
  else {
    warning('No file')
  }
}

stab_score_traj <- function(dataset_path, stab_path, clustering, metric_choice) {
  #print(paste0("Starting ", dataset_path))
  if (file.exists(dataset_path)) {
    seurat <- LoadH5Seurat(dataset_path, verbose=F)
    iters <- read.csv(stab_path[1], header=FALSE, row.names=1)
    cell_mask <- read.csv(stab_path[2], header=FALSE)
    feat_mask <- read.csv(stab_path[3], header=FALSE)
    mask <- list(cell_mask, feat_mask)
    mask <- lapply(mask, function(msk) {
      tmp<-lapply(seq(n_iter-1), function(iter) {
                    iter1=as.logical(gsub(F,NA,as.logical(msk[iter,])))
                    iter2=as.logical(gsub(F,NA,as.logical(msk[(iter+1),])));
                    iter1[intersect(which(is.na(iter2)),
                                    which(iter1==T))]<-F;
                    iter2[intersect(which(is.na(iter1)),
                                    which(iter2==T))]<-F;
                    return(list(iter1[!is.na(iter1)],
                                iter2[!is.na(iter2)]))});
      names(tmp)<-paste(seq(n_iter-1), seq(2,n_iter));
      return(tmp)})
    names(mask) <- c("cell", "feat")
    mask2 <- list(cell_mask, feat_mask)
    mask2 <- lapply(mask2, function(msk) {
      tmp<-lapply(seq(n_iter-1), function(iter) {
        iter1=as.logical(msk[iter,]);
        iter2=as.logical(msk[(iter+1),]);
        return(iter1 & iter2)});
      names(tmp)<-paste(seq(n_iter-1), seq(2,n_iter));
      return(tmp)})
    names(mask2) <- c("cell","feat")
    iters <- lapply(seq(n_iter), function(i) {
      tmp<-strsplit(iters[i,],"\t")[[1]];
      if (length(tmp)>0) {
        tmp<-unlist(strsplit(tmp,"\n"));
        values<-as.numeric(tmp[-grep(",",tmp)]);
        idx<-sapply(tmp[grep(",",tmp)], function(j) {
          splitted=strsplit(j,
                            split="")[[1]];
          if (length(splitted)==8) {
            index=1+as.numeric(splitted[c(4,7)])
          }
          else {
            before1=4
            after1=grep(",",splitted)[1]-1
            before2=after1+3
            after2=length(splitted)-1
            if (before1==after1) {
              nb1=as.numeric(splitted[before1])+1
            }
            else {
              nb1=as.numeric(paste(splitted[before1:after1], collapse=""))+1
            }
            if (before2==after2) {
              nb2=as.numeric(splitted[before2])+1
            }
            else {
              nb2=as.numeric(paste(splitted[before2:after2], collapse=""))+1
            }
            index=c(nb1,nb2)
          }
          return(index)});
        M=matrix(0,
                 ncol=max(idx),
                 nrow=max(idx));
        for (z in seq(length(values))) {
          M[idx[1,z],idx[2,z]]=as.numeric(values[z])
        };
        M=as(M,"sparseMatrix")
        };
      else {
        M=as(matrix(0),"sparseMatrix")
      };
      return(M)})
    seurat@misc$paga <- NULL
    for (i in seq(n_iter)) {
      seurat@misc[[paste0("paga",i)]] <- iters[[i]]
    }
    paga <- seurat@misc[grep("paga",names(seurat@misc))]
    n_col <- lapply(paga, function(pg)
      sapply(seq(ncol(pg)), function(x)
        pg@p[x+1] - pg@p[x]))
    from <- lapply(n_col, function(y)
      unlist(sapply(seq(y), function(x)
        rep(as.character(x),y[x]))))
    to <- lapply(paga, function(pg)
      as.character(pg@i+1))
    missing_cluster_ids <- lapply(seq(n_iter), function(x)
      setdiff(as.character(as.numeric(seurat@meta.data[, clustering])),
              c(from[[x]],to[[x]])))
    from <- lapply(seq(n_iter), function(x)
      c(from[[x]], missing_cluster_ids[[x]]))
    to <- lapply(seq(n_iter), function(x)
      c(to[[x]], missing_cluster_ids[[x]]))
    length <- lapply(seq(n_iter), function(y)
      c(paga[[y]]@x, rep(0,length(missing_cluster_ids[[y]]))))
    directed <- lapply(length, function(x)
      rep(T, length(x)))
    milestone_network <- lapply(seq(n_iter), function(x)
      tibble("from"=from[[x]],
             "to"=to[[x]],
             "length"=length[[x]],
             "directed"=directed[[x]]))
    dataset <- lapply(seq(n_iter-1), function(iter) {
      wrap_expression(counts = t(seurat@assays$RNA@counts[mask2$feat[[iter]],
                                                          mask2$cell[[iter]]]),
                      expression = t(seurat@assays$RNA@data[mask2$feat[[iter]],
                                                            mask2$cell[[iter]]]))})
    names(dataset) <- paste(seq(n_iter-1), seq(2,n_iter))
    traj <- lapply(seq(n_iter-1), function(iter)
      lapply(0:1, function(comp_idx) {
        dataset[[paste(iter,iter+1)]] %>%
          add_grouping(as.character(as.numeric(seurat@meta.data[, clustering][mask2$cell[[iter]]]))) %>%
          add_cluster_graph2(milestone_network[[iter+comp_idx]])}))
    names(traj) <- paste(seq(n_iter-1),seq(2,n_iter))
    traj <- lapply(seq(n_iter-1), function(iter)
      lapply(seq(2), function(comp_idx)
        traj[[iter]][[comp_idx]] %>%
          add_dimred(dyndimred::dimred_mds,
                     expression_source = dataset[[paste(iter,iter+1)]]$expression)))
    traj <- lapply(traj, function(traj1)
      lapply(traj1, function(traj2)
        traj2 %>% add_cell_waypoints()))
    metric_ids <- dyneval::metrics %>%
      filter(category != "average") %>%
      filter(perfect==1) %>%
      ilter(metric_id %in% metric_choice) %>%
      pull(metric_id)
    metrics <- lapply(metric_ids,
                      function(metric) lapply(seq(n_iter-1), function(iter) {
                        calculate_metrics2(traj[[iter]][[1]],
                                           traj[[iter]][[2]],
                                           metric,
                                           dataset[[paste(iter,iter+1)]]$expression)[,metric]}))
    metrics <- lapply(metrics, function(iter) mean(unlist(iter)))
    metrics <- data.frame(metrics)
    colnames(metrics) <- metric_ids
    return(metrics)
  }
  else {
    stop('There is no file at the given path')
  }
}
