library(tibble)
library(ggplot2)
df <- read.csv("~/Downloads/modularity.csv")
df <- data.frame("Hub_red"=df$Mode,
                 "Dim"=df$PCA,
                 "Mod"=df$Louvain.modularity,
                 "Dataset"=df$GSE)
df$Hub_red <- gsub("Standard KNN graph","Base k-NN graph",
                   gsub("mutual_proximity","MP k-NN graph",
                        gsub("local_scaling","LS k-NN graph",df$Hub_red)))
ggplot(df[df$Hub_red!='LS k-NN graph',], aes(x=Hub_red, y=Mod, fill=Hub_red)) +#5x10
  geom_boxplot() +
  geom_line(aes(group=Dataset), color="grey", alpha=0.5) +
  facet_wrap(~Dim) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_fill_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
ggplot(df[df$Hub_red!='MP k-NN graph',], aes(x=Hub_red, y=Mod, fill=Hub_red)) +
  geom_boxplot() +
  geom_line(aes(group=Dataset), color="grey", alpha=0.5) +
  facet_wrap(~Dim) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_fill_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
