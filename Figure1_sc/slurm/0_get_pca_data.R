library(ggplot2)
library(dplyr)
library(pbapply)
library(biomaRt)
library(future.apply)
library(furrr)

#data_in<- "/Users/elise/Desktop/Thèse/scRNAseq/Data_sc/DimRedPaper/DuoClustering2018/#sce_full/sce_full_"
dataset <- c("Koh","KohTCC","Kumar","KumarTCC","SimKumar4easy","SimKumar4hard","SimKumar8hard","Trapnell","TrapnellTCC","Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
data_out <- "/Users/elise/Desktop/GitHub/Hubness_sc/Data/Duo_10kHVG/"
pca_out <- "/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/pca_csv/"
sdev_out <- "/Users/elise/Desktop/GitHub/Hubness_sc/Figure1_sc/pca_csv/sdev_"

filter <- function(df) {
  feature <- apply(df,2,function(x) sum(x!=0, na.rm = T))
  count <- apply(df,2,function(x) sum(x, na.rm=T))
  print(plot(feature,count))
}
log_transfo <- function(data) {
  data <- log10(data+1)
  return(data)
}

plan("multiprocess",workers=30)
# load data while keeping only 10k HVG
#data <- pblapply(dataset,
#                 function(x) {tmp <- readRDS(paste0(data_in,x,".rds"));
#                 tmp <- tmp@assays$data@listData$counts;
#                 if (length(grep("TCC",x))==0 & length(grep("Sim",x))==0) {
#                 gene.names=data.frame("ensemblID"=rownames(tmp))
#                 gene.names$ensembl <- gsub('\\..+$', '', gene.names$ensembl)
#                 ensembl <- useDataset("hsapiens_gene_ensembl",
#                                       useMart("ensembl"))
#                 genemap <- getBM(attributes = c("ensembl_gene_id",
#                                                 "hgnc_symbol"),
#                                   filters = "ensembl_gene_id",
#                                   values = gene.names$ensembl,
#                                   mart = ensembl )
#                 if (x=="Kumar") {
#                 ensembl <- useDataset("mmusculus_gene_ensembl",
#                                       useMart("ensembl"))
#                 genemap <- getBM(attributes = c("ensembl_gene_id",
#                                                 "external_gene_name"),
#                                   filters = "ensembl_gene_id",
#                                   values = gene.names$ensembl,
#                                   mart = ensembl )
#                 genemap$hgnc_symbol = genemap$external_gene_name
#                 }
#                 idx <- match( gene.names$ensembl, genemap$ensembl_gene_id )
#                 gene.names$hgnc_symbol <- genemap$hgnc_symbol[idx]
#                 rownames(gene.names) <- gene.names$ensemblID
#                 gene.names <- gene.names[rownames(tmp),]
#                 bool_rm <- gene.names$hgnc_symbol=="" |
#                   is.na(gene.names$hgnc_symbol) |
#                   duplicated(gene.names$hgnc_symbol) |
#                   is.null(gene.names$hgnc_symbol)
#                 tmp <- tmp[!bool_rm,]
#                 rownames(tmp) <- gene.names$hgnc_symbol[!bool_rm]
#                 }
#                 hvg <- names(sort(apply(tmp,1,var),
#                                   decreasing = T)[1:min(1e4,nrow(tmp))]);
#                 tmp <- tmp[hvg,]
#                 return(tmp)})
#pblapply(data,filter)
#pbsapply(data,dim)
#future_lapply(1:length(data),
#         function(x) write.table(data[[x]], file=paste0(data_out,dataset[x],".csv")))
data <- future_lapply(seq(dataset), function(x)
         read.table(paste0(data_out,dataset[x],".csv")))

# Log transfo
data <- future_lapply(data, function(x) log_transfo(x))

# PCA
pca <- furrr::future_map(data, function(x) prcomp(x, center=T, scale.=F), .progress = T)
data_proj <- future_lapply(pca, function(x) t(x$rotation))

future_lapply(1:length(data_proj),
         function(x) write.table(data_proj[[x]], file=paste0(pca_out,dataset[x],"_pca_readyforhubness.csv")))
future_lapply(1:length(data_proj),
          function(x) write.table(pca[[x]]$sdev,file=paste0(sdev_out,dataset[x],".csv")))
