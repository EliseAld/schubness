# ====
# Data loading
# ====
n_dataset = 12
dim_nb = c(531,531,246,263,500,499,499,222,227,3994,6498,3994)
cell_nb = c(531,531,246,263,500,499,499,222,227,3994,6498,3994)
pc.val <- lapply(dim_nb,
                   function(x) return(c(2,5,10,20,30,40,50,100,200,x-1)))
k.val <- c(5,10,50,100)

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC", "Zhengmix4eq", "Zhengmix4uneq","Zhengmix8eq")
localization = "/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/scores_scaled/"
all_paths=lapply(1:length(dataset), function(z) paste0(localization,
                                                    dataset[z],
                                                    "_kNN_occurence_",
                                                    pc.val[[z]],
                                                    "pca_minkow_bis.rds"))
hubness_scores_data <- pblapply(all_paths, function(x) lapply(x,function(y) {readRDS(y)})) # First 7 lists are pca2, then pca5 etc

# ====
# Data writing
# ====
hubness <- lapply(seq(n_dataset), function(x) data.frame("score"=unlist(hubness_scores_data[[x]]),
                                                         "Dimension"=rep(pc.val[[x]], each=length(k.val)*cell_nb[x]),
                                                         "k"=rep(k.val, each=cell_nb[x], times=length(pc.val[[x]]))))
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
for (i in 1:length(hubness)) {
   hubness[[i]]$Antihubs <- rep(find_antihub_nb_per_condition(hubness[[i]], cell_nb = cell_nb[i]),each=cell_nb[i])
}
# Log to zoom on differences BAD IDEA
#hubness_log <- hubness
#hubness_log$score_log <- log(hubness_log$score + 1)
rm(localization,find_antihub_nb_per_condition,all_paths)
