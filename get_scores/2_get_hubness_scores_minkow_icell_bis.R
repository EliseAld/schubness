library(Seurat)
library(philentropy)
library(pbapply)

# MAKE THE DIST MATRIX WITH ALL ENTRIES EXCEPT DIAG INSTEAD OF UPPER TRI

matrix <- read.table("GSE99254_NSCLC.TCell.S4519.count.labeled.icell_pca_readyforhubness.txt")
print("data loaded!")
n_cell = ncol(matrix)

# matrix is PCs x cell

# Range of k values
k.val <- c(5,50,100,200)

# Range of p values for Lp norm
p.val <- c(0.1,0.5,1,1.5,2,4,10)

# Range of PCs nb values
pc.val <- c(2,5,10,50,100,500,1000,n_cell-1)

# Write the function for the kNN graph
distance_dim <- function(data,p) {
  mat <- t(data)
  dist.matrix <- philentropy::distance(mat,method="minkowski",p)
  diag(dist.matrix) <- NA # remove diag only
  return(dist.matrix)
}
kNN <- function(dist.matrix,k) {
  k.nn <- t(apply(dist.matrix,1,function(x) order(x,na.last=T)[1:k]))
  return(k.nn)
}
# Get the k-occurences from the kNN graph
k.occurence <- function(k.nn,data) {
  occurence <- sapply(X=1:nrow(k.nn), FUN=function(x,y) {sum(y==x)}, y=k.nn)
  names(occurence) <- colnames(data)
  return(occurence)
}
# Get the scores
get_scores <- function(data,k.val,p.val,pc.val) {
  for (pc in pc.val) {
    data_pc <- data[1:pc,]
    for (p in p.val) {
      dist.matrix <- distance_dim(data=data_pc,p=p)
      minkow <- pbapply::pblapply(X=k.val, FUN=function(x,y,z) {k.occurence(kNN(y,x),z)}, y=dist.matrix, z=data)
      saveRDS(minkow, file=paste0("kNN_occurence_",p,"_pca",pc,"_minkow_icell_bis.rds"))
      print(paste0(p, " pval done for PC",pc))
    }
  }
}
get_scores(matrix, k.val=k.val, p.val=p.val, pc.val=pc.val)
