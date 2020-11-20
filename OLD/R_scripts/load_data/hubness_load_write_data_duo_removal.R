# ====
# Data loading
# ====
dim_nb = c(471,471,206,223,450,439,439,192,167)
cell_nb = c(471,471,206,223,450,439,439,192,167)
pval=2
pcval=pblapply(dim_nb,
               function(x) return(c(2,5,10,20,30,40,50,100,x-1)))
kval=c(5,50,100,200)

dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC")
localization = "/Users/elise/Desktop/Github/Hubness_sc/results/"
path = sapply(dataset, function(x) paste0("duo_removal/",x,"/"))
all_paths=lapply(1:length(path), function(z) paste0(localization,
                                                    path[z],
                                                    "kNN_occurence_",
                                                    pval,
                                                    "_pca",
                                                    pcval[[z]],
                                                    "_minkow_bis.rds"))
hubness_scores_data <- pblapply(all_paths, function(x) lapply(x,function(y) {readRDS(y)})) # First 7 lists are pca2, then pca5 etc

# ====
# Data writing
# ====

hubness <- lapply(1:length(pcval), function(x) data.frame("score"=unlist(hubness_scores_data[[x]]),
                                                              "Dimension"=rep(pcval[[x]], each=length(pval)*length(kval)*cell_nb[x]),
                                                              "p"=rep(pval, each=length(kval)*cell_nb[x], times=length(pcval[[x]])),
                                                              "k"=rep(kval, each=cell_nb[x], times=length(pcval[[x]])*length(pval))))
# Reorder factor
#hubness$Dimension <- factor(hubness$Dimension, levels = pcval)
# Removing cells with hubness score = 0 DO NOT DO IT
find_antihub_nb_per_condition <- function(hubness_df, cell_nb) {
   pc_nb = length(unique(hubness_df$Dimension))
   p_nb = length(unique(hubness_df$p))
   k_nb = length(unique(hubness_df$k))
   antihb = c()
   for (i in 1:(pc_nb*p_nb*k_nb)) {
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
rm(localization,find_antihub_nb_per_condition,path,all_paths)
