library(Seurat)
library(pbapply)
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/zhouDR/datasets/rds/"
dataset <- sort(c('gold_hematopoiesis-gates_olsson.rds', 'gold_germline-human-female-weeks_li.rds',
                   'gold_stimulated-dendritic-cells-LPS_shalek.rds',
                   'gold_germline-human-female_guo.rds', 'gold_mESC-differentiation_hayashi.rds',
                   'gold_developing-dendritic-cells_schlitzer.rds',
                   'gold_germline-human-male-weeks_li.rds',
                   'gold_pancreatic-beta-cell-maturation_zhang.rds',
                   'gold_human-embryos_petropoulos.rds', 'gold_aging-hsc-old_kowalczyk.rds',
                   'gold_germline-human-male_guo.rds', 'gold_myoblast-differentiation_trapnell.rds',
                   'gold_pancreatic-alpha-cell-maturation_zhang.rds',
                   'gold_aging-hsc-young_kowalczyk.rds'))
dataRDS <- pblapply(dataset,
                    function(x) readRDS(paste0(path,x)))
dataRDS <- lapply(dataRDS[7],
                  function(x) {tmp1=data.frame(x$prior_information$groups_id);
                  rownames(tmp1)<-tmp1$cell_id;
                  tmp1<-tmp1[x$cell_ids,];
                  tmp2=x$counts;
                  return(list("counts"=tmp2,"groups_id"=tmp1))})
saveRDS(dataRDS, file = paste0(path,dataset[7]))
rm(dataRDS)


