library(tibble)
library(ggplot2)
df <- read.csv("~/Downloads/degree_distrib.csv", row.names=1)
n_iter=100
n_bins=199
dim=c(2,10,50,500)
distrib=c("Gaussian","Uniform")
mean <- sapply(df, function(x) {grp <- list();
for (i in seq(length(dim)*length(distrib))) {
  grp[[i]]=x[((i-1)*n_iter+1):(i*n_iter)]
}
sapply(grp,mean)})
df <- data.frame("Dimension"=factor(rep(rep(dim,each=n_bins),length(distrib))),
                 "Distribution"=rep(distrib,each=n_bins*length(dim)),
                 "Bin"=rep(seq(n_bins),length(dim)*length(distrib)),
                 "PDF"=c(t(mean)))

ggplot(df, aes(x=Bin, y=PDF, color=Dimension)) +#5x10
  geom_point() +
  geom_line() +
  facet_wrap(~Distribution, scales = "free") +
  xlim(c(0,30)) +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))


ggplot(df[df$Dimension!="2" & df$Dimension!="10",], aes(x=Bin, y=PDF, color=Dimension)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~Distribution, scales = "free") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
