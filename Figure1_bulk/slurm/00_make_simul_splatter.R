library(future.apply)
library(pbapply)
plan("multiprocess",workers=30)

data_brca <- read.csv("/shared/projects/mne/elise/bulk/TCGA_BRCA/brca.csv", row.names=1)
data_archs4 <- read.csv("/shared/projects/mne/elise/bulk/ARCHS4/data.csv", row.names=1)
load("/shared/projects/mne/elise/bulk/TCGA_KIRC/TCGA_KIRC_RPKM_withHGNCsymbol.Rdata")
data_kirc <- matrix
rm(matrix)
data <- list(data_brca, data_kirc, data_archs4)
names(data) <- c("TCGA_BRCA","TCGA_KIRC","ARCHS4")

zero_percent <- function(data) { # gene x cell matrix
  return(sum(data==0)/ncol(data)/nrow(data)*100)
}
splatSimDropout <- function(data, dropout_percentage) {
    dropout.mid <- -10
    dropout.shape <- -1
    nCells <- ncol(data)
    nGenes <- nrow(data)
    
    # Generate probabilities based on expression
    drop.prob <- vapply(seq_len(nCells), function(idx) {
        eta <- log(data[,idx])
        return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
    }, as.numeric(seq_len(nGenes)))

    # Decide which counts to keep
    keep <- matrix(rbinom(nCells * nGenes, 1, 1-drop.prob),
                   nrow = nGenes, ncol = nCells)
    counts <- list(data*keep)
    print(paste0("Adding ", round(mean(keep[data!=0]==0)*100), "% of dropout"))
    for (idx in seq(2,length(dropout_percentage))) {
      while (mean(keep[data!=0]==0) < dropout_percentage[idx]/100) {
        dropout.mid <- dropout.mid + 0.5
        drop.prob <- vapply(seq_len(nCells), function(idx) {
        eta <- log(data[,idx])
        return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        }, as.numeric(seq_len(nGenes)))
        keep <- matrix(rbinom(nCells * nGenes, 1, 1-drop.prob),
                   nrow = nGenes, ncol = nCells)
      }
      counts <-  c(counts, list(data*keep))
      print(paste0("Adding ", round(mean(keep[data!=0]==0)*100), "% of dropout"))
    }
    names(counts) <-  dropout_percentage
    return(counts)
}
logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}

dropout_percent <- c(0,10,20,30,40,50,60,70,80,90,95)
data_simul <- future_lapply(data,
                     function(y) splatSimDropout(y,dropout_percent))
sapply(data_simul,
       function(x) sapply(x,zero_percent))
mapply(x=data_simul, y=names(data_simul),
       function(x,y) pbmapply(t=x, z=dropout_percent,
                            function(t,z) write.table(t, file=paste0("/shared/projects/mne/elise/bulk/", y, "/splatter_simu/simul_dropout", z,".csv"))))
print("SUCCESS SIMU SPLATTER")
