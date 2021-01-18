# Download FPKM for TCGA_KIRC 
# Build matrix to do data deconvolution

# Load libraries
rm(list=ls())
load("~/Desktop/Work Documents/ccRCC/TCGA-KIRC/TCGA_KIRC_Target.Rdata")



#retrive file name for RPKM
sample.sheet=read.table("~/Desktop/Work Documents/ccRCC/TCGA-KIRC/gdc_sample_sheet.2020-05-06.tsv", sep="\t", header=T)
head(sample.sheet)
which(sample.sheet$Sample.Type=="Additional - New Primary")
sample.sheet=sample.sheet[-c(528),]
sample.sheet$IDCommon=t(as.data.frame(strsplit(as.character(sample.sheet$File.Name), split='.', fixed = T)))[,1]

target$IDCommon=t(as.data.frame(strsplit(as.character(target$File.Name), split='.', fixed = T)))[,1]
head(target[,c("IDCommon", "Sample.ID")])

MtxBuilder=sample.sheet[,c("IDCommon", "File.ID", "File.Name")]
MtxBuilder=dplyr::right_join(MtxBuilder, target[,c("IDCommon", "Sample.ID")], by='IDCommon')
length(unique(MtxBuilder$Sample.ID))

# Create matrix for RPKM data

file=read.table(paste0("~/Downloads/gdc_download_20200506_170052.213277/",MtxBuilder$File.ID[1],"/",MtxBuilder$File.Name[1]))
colnames(file)=c("EnsemblID",as.character(MtxBuilder$Sample.ID[1]))
data=file
for (i in 2:length(MtxBuilder$Sample.ID)){
  file=read.table(paste0("~/Downloads/gdc_download_20200506_170052.213277/",MtxBuilder$File.ID[i],"/",MtxBuilder$File.Name[i]))
  data=cbind(data,file$V2)
}
colnames(data)=c("EnsemblID",as.character(MtxBuilder$Sample.ID))

matrix<-data
rownames(matrix)=matrix$EnsemblID
matrix=matrix[,-1]
head(matrix[,1:5])

# add gene symbol instead of ENSEMBL
library(biomaRt)
gene.names=data.frame("ensemblID"=rownames(matrix))
gene.names$ensembl <- gsub('\\..+$', '', gene.names$ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#to see what you can retrieve --> listAttributes(ensembl)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = gene.names$ensembl,
                  mart = ensembl )
#write.csv(genemap, "genemap.csv")
#genemap=read.csv("genemap.csv", header = T)
idx <- match( gene.names$ensembl, genemap$ensembl_gene_id )
gene.names$entrez <- genemap$entrezgene_id[ idx ]
gene.names$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
#stock<-gene.names
all(rownames(matrix)==gene.names$ensemblID)
all(!duplicated(gene.names$ensemblID))
all(!duplicated(gene.names$hgnc_symbol)) 
all(!is.na(gene.names$hgnc_symbol))
all(!is.null(gene.names$hgnc_symbol))

which(gene.names$hgnc_symbol=="")
which(is.na(gene.names$hgnc_symbol))

rownames(gene.names)=gene.names$ensemblID
gene.names=gene.names[-c(2,4),]

matrix=matrix[which(rownames(matrix)%in%gene.names$ensemblID),]
all(rownames(matrix)==gene.names$ensemblID)
rownames(matrix)=gene.names$hgnc_symbol

head(matrix)

save(matrix,file="TCGA_KIRC_RPKM_withHGNCsymbol.Rdata")

