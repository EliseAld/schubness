library(pbapply)

# MAKE THE DIST MATRIX WITH ALL ENTRIES EXCEPT DIAG INSTEAD OF UPPER TRI

Dataset1 = "Zheng4"
dataset1="zheng4"
Dataset2 = "Guo"
dataset2="guo"

matrix1 <- read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/forTAC/results/",Dataset1,"/",dataset1,"_pca_readyforhubness.txt"))
matrix2 <- read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/forTAC/results/",Dataset2,"/",dataset2,"_pca_readyforhubness.txt"))
print("data loaded!")
n_dim1 = nrow(matrix1)
n_dim2 = nrow(matrix2)

# matrix is PCs x cell

# Range of PCs values
pc.val <- c(2,5,10,20,30,40,50,100,500,min(n_dim1,n_dim2)-1)

# Write the function for the kNN graph
get_mnn <- function(data1,data2) {
  distance1to2 <- pbapply(data1,
                          2,
                          function(x) {tmp<-apply(data2,
                                                  2,
                                                  function(y) dist(rbind(x,y)));
                          return(names(which.min(tmp)))})
  distance2to1 <- pbapply(data2,
                          2,
                          function(x) {tmp<-apply(data1,
                                                  2,
                                                  function(y) dist(rbind(x,y)));
                          return(names(which.min(tmp)))})
  df1 <- data.frame("from"=colnames(data1),
                    "to"=distance1to2)
  df2 <- data.frame("from"=colnames(data2),
                    "to"=distance2to1)
  df <- data.frame("from"=df1$from,
                   "to"=df1$to,
                   "back"=df2[df1$to,]$to)
  df <- df[as.character(df$from)==as.character(df$back),]
  return(data.frame("MNN1"=df$from,
                    "MNN2"=df$to))
}

# Get the k-occurences from the kNN graph
k.occurence <- function(k.nn,data) {
  occurence <- sapply(X=1:nrow(k.nn), FUN=function(x,y) {sum(y==x)}, y=k.nn)
  names(occurence) <- colnames(data)
  return(occurence)
}
# Get the scores
get_scores <- function(data,k.val,p.val,pc.val,name) {
  for (pc in pc.val) {
    data_pc <- data[1:pc,]
    for (p in p.val) {
      dist.matrix <- distance_dim(data=data_pc,p=p)
      minkow <- pbapply::pblapply(X=k.val, FUN=function(x,y,z) {k.occurence(kNN(y,x),z)}, y=dist.matrix, z=data)
      saveRDS(minkow, file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/forTAC/results/",name,"/kNN_occurence_",p,"_pca",pc,"_minkow_bis.rds"))
      print(paste0(p, " pval done for PC",pc))
    }
  }
}
get_scores(matrix, k.val=k.val, p.val=p.val, pc.val=pc.val, Dataset)

