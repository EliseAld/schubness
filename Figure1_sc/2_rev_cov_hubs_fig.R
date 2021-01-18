# libs
library(pbapply)
library(RANN)
library(ggplot2)
library(dplyr)

# load data
params <- list(10,50)
names(params) <- c("k","Dimension") 
dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")

# Check the skewness
#scores_df <- data.frame("InDegree"=unlist(data_hubness))
#scores_df$Dataset <- unname(unlist(mapply(x=dataset, y=data_hubness, function(x,y)
#  rep(x, length(y)))))
#ggplot(scores_df, aes(InDegree, fill=Dataset)) +
#  geom_density(alpha=0.6) +
#  scale_fill_viridis_d(option="plasma")

# reverse get all the data fixing k but increasing N
rev_df <- readRDS("/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/slurm/hubness_revcov_ggplot1.rds")
ggplot(rev_df[rev_df$Setting=="Hub",], aes(x=N_percentage, y=Rev_cov_size, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  scale_color_viridis_d(option="plasma") +
  xlab("Putative hubs (%)") +
  ylab("Size of the reverse coverage (%)") +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=25))
ggplot(rev_df, aes(x=N_percentage, y=Increment, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  facet_wrap(~Setting) +
  scale_color_viridis_d(option="plasma") +
  xlab("Putative hubs (%)") +
  ylab("Increment in reverse coverage size (%)") +
  scale_x_log10() +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=25))

# get hubs

# Remove hubs

# Save data

# make PCA

# Check the screeplot
#sdev <- lapply(pca, function(x) x$sdev/x$sdev[1])
#sdev_df <- data.frame("Sdev" = unlist(sdev),
#                      "Dataset"=unlist(mapply(x=dataset,y=lengths(sdev), function(x,y)
#                        rep(x, y))),
#                      "PC"=unlist(sapply(lengths(sdev), seq)))
#ggplot(sdev_df[sdev_df$PC<30,], aes(x=PC, y=Sdev, color=Dataset)) +
#  geom_point() +
#  geom_line(aes(group=Dataset)) +
#  geom_hline(yintercept=0.1) +
#  scale_color_viridis_d(option="plasma") +
#  ylab("") +
#  xlab("Principal Components")

# get the scores (from the .R file) USELESS

# Evaluate the asymmetry post removal
asymmetry_df <- readRDS("/Users/elise/Desktop/Github/Hubness_sc/Figure1_sc/slurm/hubness_revcov_ggplot2.rds")
ggplot(asymmetry_df, aes(x=Dimension, y=GHubness, color=Dataset)) +
  geom_point() +
  geom_line(aes(group=Dataset)) +
  scale_color_viridis_d(option="plasma") +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=25)) +
  ylab("Asymmetry (%)")
  
