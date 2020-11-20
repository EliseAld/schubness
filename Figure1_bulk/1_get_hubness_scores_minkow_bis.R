library(pbapply)
library(RANN)

# MAKE THE DIST MATRIX WITH ALL ENTRIES EXCEPT DIAG INSTEAD OF UPPER TRI

dataset <- c("TCGA_BRCA", "ARCHS4")
dropout_percent <- c(0,5,10,20,30,40,50,60,70,75,80,85,90,95)
matrix <- lapply(dataset,
                   function(x) pblapply(dropout_percent,
                                        function(y)
                                        read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/", x, "/pca/simul_dropout", y,"_pca.csv"))))
names(matrix) <- dataset
print("data loaded!")
n_dim = sapply(matrix, function(x) unique(sapply(x,nrow)))

# matrix is PCs x cell

# Range of k values
k.val <- c(5,10,50,100)

# Range of PCs nb values
pc.val <- lapply(n_dim,
                   function(x) return(c(2,5,10,20,30,40,50,100,200,x-1)))

# Write the function for the kNN graph
kNN <- function(data,k) {
  k.nn <- nn2(t(data), k=k+1)$nn.idx
  k.nn <- k.nn[, 2:(k+1)]
  return(k.nn)
}
# Get the k-occurences from the kNN graph
k.occurence <- function(k.nn,data) {
  occurence <- sapply(seq(nrow(k.nn)), function(x) sum(k.nn==x))
  names(occurence) <- colnames(data)
  return(occurence)
}
# Get the scores
get_scores <- function(data,k.val,pc.val,name,dropout) {
  print(paste0("Doing ",name))
  for (pc in pc.val) {
    data_pc <- data[1:pc,]
    score <- lapply(k.val,
                               function(x) {k.occurence(kNN(data_pc,x), data_pc)})
    saveRDS(score, file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/scores/",name,"/Dropout",dropout,"_kNN_occurence_",pc,"pca_minkow_bis.rds"))
  }
}
lapply(dataset,
         function(x) pblapply(seq(length(dropout_percent)),
                              function(y) get_scores(matrix[[x]][[y]], k.val=k.val, pc.val=pc.val[[x]], name=x, dropout=dropout_percent[y])))

