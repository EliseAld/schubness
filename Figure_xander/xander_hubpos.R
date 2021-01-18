library(ggplot2)
df1 <- read.csv("~/Desktop/Github/Hubness_sc/Figure_xander/df_score_10e5.csv", row.names = 1)
df2 <- read.csv("~/Desktop/Github/Hubness_sc/Figure_xander/df_norm_10e5.csv", row.names = 1)
df3 <- read.csv("~/Desktop/Github/Hubness_sc/Figure_xander/df_rank_10e5.csv", row.names = 1)
dim = c(2,5,8,10,20)
distrib = c("unif","gauss")

rank_list = list()
norm_list = list()
for (j in seq(ncol(df1))) {
  d=dim[ifelse(j%%5==0,5,j%%5)]
  print(j)
  tmp_norm = c()
  tmp_rank = c()
  for (i in seq(max(df1))) {
    mask = df1[,j]>= i
    tmp_norm = c(tmp_norm,ifelse(sum(mask)<100,NA,mean(sqrt(df2[mask,j]/d))))
    tmp_rank = c(tmp_rank,ifelse(sum(mask)<100,NA,mean(df3[mask,j])))
  }
  norm_list[[j]] = tmp_norm
  rank_list[[j]] = tmp_rank
}

names(norm_list) <- colnames(df1)
names(rank_list) <- colnames(df1)

df <- data.frame("Distrib"=rep(distrib,each=length(dim)*max(df1)),
                 "Dim"=factor(rep(rep(dim,each=max(df1)),length(distrib))),
                 "Av_rank"=unlist(rank_list),
                 "Av_norm"=unlist(norm_list),
                 "x"=rep(seq(max(df1)),length(distrib)*length(dim)))
df <- df[!is.na(df$Av_rank),]

ggplot(df, aes(x=x, y=Av_rank, color=Dim)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Distrib, scales = "free") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
ggplot(df, aes(x=x, y=Av_norm, color=Dim)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Distrib, scales = "free") +
  scale_color_brewer(palette="Set1") +
  theme(panel.background = element_rect(fill="grey98"),
        panel.grid = element_line(colour = "grey80"))
