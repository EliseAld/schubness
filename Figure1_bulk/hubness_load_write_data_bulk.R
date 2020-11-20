# ====
# Data loading
# ====
n_dataset = 2
dim_nb = c(1213,1670)
cell_nb = c(1213,1670)
dropout_percent = c(0,5,10,20,30,40,50,60,70,75,80,85,90,95)
pc.val <- lapply(dim_nb,
                   function(x) return(c(2,5,10,20,30,40,50,100,200,x-1)))
k.val <- c(5,10,50,100)

dataset <- c("TCGA_BRCA","ARCHS4")
localization = "/Users/elise/Desktop/Github/Hubness_sc/Figure1_bulk/scores/"
all_paths=lapply(seq(n_dataset), function(z) lapply(dropout_percent,
                                               function(y) paste0(localization,"",
                                                                  dataset[z],"/Dropout",
                                                                  y,
                                                                  "_kNN_occurence_",
                                                                  pc.val[[z]],
                                                                  "pca_minkow_bis.rds")))
hubness_scores_data <- lapply(all_paths,
                              function(z) pblapply(z, function(x) lapply(x,function(y) {readRDS(y)}))) # First 7 lists are pca2, then pca5 etc

# ====
# Data writing
# ====
hubness <- lapply(seq(n_dataset),
                  function(x) {tmp=lapply(seq(length(dropout_percent)),
                                     function(y) data.frame("score"=unlist(hubness_scores_data[[x]][[y]]),
                                                            "Dimension"=rep(pc.val[[x]], each=length(k.val)*cell_nb[x]),
                                                            "k"=rep(k.val, each=cell_nb[x], times=length(pc.val[[x]])),
                                                            "Dropout"=dropout_percent[y]));
                     names(tmp)=dropout_percent;
                     return(tmp)})
names(hubness)=dataset
# Reorder factor
#hubness$Dimension <- factor(hubness$Dimension, levels = pcval)
# Removing cells with hubness score = 0 DO NOT DO IT
find_antihub_nb_per_condition <- function(hubness_df, cell_nb) {
   pc_nb = length(unique(hubness_df$Dimension))
   k_nb = length(unique(hubness_df$k))
   antihb = c()
   for (i in 1:(pc_nb*k_nb)) {
      borne_inf=(cell_nb*(i-1)+1)
      borne_sup=cell_nb*i
      antihb = c(antihb,sum(hubness_df$score[borne_inf:borne_sup]==0))
   }
   return(antihb)
}
for (i in seq(n_dataset)) {
   for (d in seq(length(dropout_percent))) {
   hubness[[i]][[d]]$Antihubs <- rep(find_antihub_nb_per_condition(hubness[[i]][[d]], cell_nb = cell_nb[i]), each=cell_nb[i])
}
}
# Log to zoom on differences BAD IDEA
#hubness_log <- hubness
#hubness_log$score_log <- log(hubness_log$score + 1)
rm(localization,find_antihub_nb_per_condition,all_paths)
