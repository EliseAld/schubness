library(Seurat)
library(philentropy)

matrix <- read.table("GSE99254_NSCLC.TCell.S9055.count.labeled.pca.txt")

# matrix is PCs x cell

# Range of p values for Lp norm
p.val <- c(0.1,0.5,1,1.5,2,4)

# Range of k values
k.val <- c(5,10,20,50,100,200)
k.val <- k.val[5]

# Nb of PCs
pc.nb <- 5
matrix <- matrix[1:pc.nb,]

# Write the function for the kNN graph
kNN <- function(data,method="minkowski",k,p) {
  dist.matrix <- philentropy::distance(t(data),method,p)
  dist.matrix[lower.tri(dist.matrix, diag=T)] <- NA # make the matrix upper triangular
  k.nn <- t(apply(dist.matrix,1,function(x) order(x,na.last=T)[1:k]))
  return(k.nn)
}

# Get the k-occurences from the kNN graph
k.occurence <- function(k.nn) {
  occurence <- rep(0,nrow(k.nn))
  names(occurence) <- 1:length(occurence)
  for (i in 1:length(occurence)) {
    occurence[i] <- sum(k.nn==i)
  }
  return(occurence)
}

# test for Minkow
minkow <- list()
for (p in p.val) {
  k.nn <- kNN(data=matrix,k=k.val,p=p)
  occu <- k.occurence(k.nn=k.nn)
   minkow[[which(p.val==p)]] <- occu
}
saveRDS(minkow, file="kNN_occurence_100_pca5_minkow.rds")
