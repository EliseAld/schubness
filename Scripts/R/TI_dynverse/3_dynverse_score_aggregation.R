#' @param scores dataframe of quality metric scores with columns Dataset_id (identifier for each dataset),
#' Traj_type (type of trajectory, provided with the ground truth information),
#' Dataset_source (whether gold, silver or simulated),
#' Method_id (method used to compute the kNN graph used for the TI task)

#' @return scores ready to be plotted (normalized, aggregated, averaged)

#' @example df <- overallscore_aggregated(scores)

### Load librairies
library(dplyr)

### Function to compare & aggregate quality scores
# Score normalization
score_norm <- function(scores) {
  tmp = scores %>% group_by(Dataset_id)
  for (col in grep("Met_",colnames(tmp), value = T)) {
    col_norm = paste0(col,"_normed")
    tmp = tmp %>% mutate(!!col_norm:=pnorm(normalize(get(col))))
  }
  return(tmp)
}

# Score aggregation
score_aggregation <- function(scores) {
  tmp = score_norm(scores)
  tmp = tmp %>% mutate(comp1=paste(Traj_type,Dataset_source,Method_id))
  for (col in grep("_normed", colnames(tmp), value = T)) {
    col_norm = paste0(col,"_norm2")
    tmp = tmp %>%
      group_by(comp1) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp1)) %>% mutate(comp2=paste(Traj_type,Method_id))
  for (col in grep("_norm2",colnames(tmp), value = T)) {
    col_norm = gsub("2","3",col)
    tmp = tmp %>%
      group_by(comp2) %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(comp2)) %>% group_by(Method_id)
  for (col in grep("_norm3",colnames(tmp), value = T)) {
    col_norm = gsub("3","4",col)
    tmp = tmp %>%
      mutate(!!col_norm:=mean(get(col)))
  }
  tmp = tmp %>% filter(!duplicated(Method_id))
  tmp = tmp[,c("Method_id",grep("_norm4",colnames(tmp),value=T))]
  return(tmp)
}

# Overall aggregated score
overallscore_aggregated <- function(scores) {
  tmp = score_aggregation(scores)
  tmp$Overall_score = apply(data.frame(tmp), 1, function(x) {
    cols = grep("_norm4",colnames(tmp));
    vals = as.numeric(unname(unlist(x[cols])));
    mean(vals)})
  return(tmp)
}