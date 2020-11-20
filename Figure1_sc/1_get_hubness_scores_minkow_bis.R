library(philentropy)
library(pbapply)

# MAKE THE DIST MATRIX WITH ALL ENTRIES EXCEPT DIAG INSTEAD OF UPPER TRI

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy",
             "SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC",
             "Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
matrix <- pblapply(dataset,
                   function(x) read.table(paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/pca_csv/",x,"_pca_readyforhubness.csv")))
print("data loaded!")
n_dim = pbsapply(matrix,
                 function(x) nrow(x))

# matrix is PCs x cell

# Range of k values
k.val <- c(5,10,50,100)

# Range of PCs nb values
pc.val <- c(2,5,10,20,30,40,50,100)
pc.val <- pblapply(n_dim,
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
get_scores <- function(data,k.val,pc.val,name) {
  print(paste0("Doing ",name))
  for (pc in pc.val) {
    data_pc <- data[1:pc,]
    score <- pbapply::pblapply(k.val,
                               function(x) {k.occurence(kNN(data_pc,x), data_pc)})
    saveRDS(score, file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/scores/",name,"_kNN_occurence_",pc,"pca_minkow_bis.rds"))
  }
}
pblapply(1:length(matrix),
         function(x) get_scores(matrix[[x]], k.val=k.val, pc.val=pc.val[[x]], name=dataset[x]))

