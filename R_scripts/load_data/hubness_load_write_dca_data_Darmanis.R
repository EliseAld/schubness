# ====
# Data loading
# ====
cell_nb = 419
n_dim = 32

pval=c(0.1,0.5,1,1.5,2,4,10)
pcval=c(2,5,10,n_dim)
kval=c(5,50,100,200)

localization = "/Users/elise/Desktop/Github/Hubness_sc/results/"
path = "Darmanis/dca/"
all_paths=unlist(lapply(X=pcval,
                        FUN=function(x,y) paste0(localization,
                                                 path,
                                                 "kNN_occurence_",
                                                 y,
                                                 "_pca",
                                                 x,
                                                 "_minkow_latent_bis.rds"),
                        y=pval))
hubness_scores_data <- lapply(all_paths,function(x) {readRDS(x)}) # First 7 lists are pca2, then pca5 etc

# ====
# Data writing
# ====

hubness <- data.frame("score"=unlist(hubness_scores_data),
                      "Dimension"=rep(pcval, each=length(pval)*length(kval)*cell_nb),
                      "p"=rep(pval, each=length(kval)*cell_nb, times=length(pcval)),
                      "k"=rep(kval, each=cell_nb, times=length(pcval)*length(pval)))
# Reorder factor
hubness$Dimension <- factor(hubness$Dimension, levels = pcval)
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
hubness$Antihubs <- rep(find_antihub_nb_per_condition(hubness, cell_nb = cell_nb),each=cell_nb)
# Log to zoom on differences BAD IDEA
#hubness_log <- hubness
#hubness_log$score_log <- log(hubness_log$score + 1)

rm(localization,find_antihub_nb_per_condition,path,all_paths)
