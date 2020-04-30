library(Seurat)
library(philentropy)
library(pbapply)

matrix <- read.table("GSE100866_pca_readyforhubness.csv")
print("data loaded!")

# matrix is PCs x cell

# Range of k values
k.val <- c(5,50,100,200)

# Range of p values for Lp norm
p.val <- c(0.1,0.5,1,1.5,2,4,10)

# Range of PCs nb values
pc.val <- c(20,30,40)

# Write the function for the kNN graph
distance_dim <- function(data,p) {
  mat <- t(data)
  dist.matrix <- philentropy::distance(mat,method="minkowski",p)
  dist.matrix[lower.tri(dist.matrix, diag=T)] <- NA # make the matrix upper triangular
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
      saveRDS(minkow, file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/results/Satija/bis/kNN_occurence_",p,"_pca",pc,"_minkow.rds"))
      print(paste0(p, " pval done for PC",pc))
    }
  }
}
get_scores(matrix, k.val=k.val, p.val=p.val, pc.val=pc.val)
