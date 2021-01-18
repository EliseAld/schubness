library(ggplot2, quietly = T)
library(viridis, quietly = T)
library(ggridges, quietly = T)
library(scales, quietly = T)
library(ggpubr, quietly = T)
library(pbapply)
library(knn.covertree)
library(ggh4x)

dataset <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")

# Figure 1
df <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_splatter_ggplot1.rds")
df2 <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_splatter_ggplot1bis.rds")
df_ <- do.call(rbind, df)
df_$Dataset <- unlist(lapply(seq(dataset), function(x) rep(dataset[x], nrow(df[[x]]))))
df_$facet <- paste(df_$Dataset,df_$Method)
labz_method = rep(c("2k (% of hubs)", "Antihubs (% of hubs)", "Asymmetry (% of edges)", "Skewness"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("2k","Antihubs","Asymmetry","Skewness")))
scales_y <- list(scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,22)),
                 scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,22)),
                 scale_y_continuous(limits = c(0,20)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,100)),
                 scale_y_continuous(limits = c(0,22)))

ggplot(df_, aes(x=Dimension, y=GHubness)) +
   geom_line(aes(group=Dropout, color=Dropout)) +
   geom_point(size=0.1) +
   facet_wrap(~facet, scales="free_y", labeller = labeller(facet = labz_method), nrow = length(dataset)) +
   facetted_pos_scales(y=scales_y, x=list(scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10(),scale_x_log10())) +
   ylab(NULL) +
   expand_limits(x=c(0,6500)) +
   scale_color_viridis_c() +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))



# Figure 1 supp
df_supp <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_splatter_ggplot2.rds")
df2_supp <- readRDS("/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_bulk/slurm/hubness_param_splatter_ggplot2bis.rds")

df_supp_ <- do.call(rbind, df_supp)
df_supp_$Dataset <- unlist(lapply(seq(dataset), function(x) rep(dataset[x], nrow(df_supp[[x]]))))
df_supp_$facet <- paste(df_supp_$Dataset,df_supp_$Method)
labz_method = rep(c("Mean (% of hubs)", "Maximum hub score"), length(dataset))
names(labz_method) <- sapply(dataset, function(x) paste(x, c("Mean+Sd","Max")))
scales_y <- list(scale_y_continuous(limits = c(0,1)),
                 scale_y_continuous(limits = c(0,4.5)),
                 scale_y_continuous(limits = c(0,1)),
                 scale_y_continuous(limits = c(0,4.5)),
                 scale_y_continuous(limits = c(0,1)),
                 scale_y_continuous(limits = c(0,4.5)))

ggplot(df_supp_, aes(x=Dimension, y=GHubness)) +
   geom_line(aes(group=Dropout, color=Dropout)) +
   geom_point(size=0.1) +
   facet_wrap(~facet, scales = "free_y", labeller = labeller(facet = labz_method), nrow = length(dataset)) +
   facetted_pos_scales(y=scales_y) +
   ylab(NULL) +
   scale_color_viridis_c() +
   theme(panel.background = element_rect(fill="grey98"),
         panel.grid = element_line(colour = "grey80"))
