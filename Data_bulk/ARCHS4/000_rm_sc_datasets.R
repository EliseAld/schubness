# Compute dropout
setwd("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/ARCHS4/")
data1 <- read.csv("data1.csv", row.names=1)
data2 <- read.csv("data2.csv", row.names=1)
data_archs4 <- list(data1,data2)
data_archs4 <- data.frame(do.call(cbind,data_archs4))
data_brca <- read.csv("/Users/elise/Desktop/GitHub/Hubness_sc/Data_bulk/TCGA_BRCA/brca.csv", row.names=1)

dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
dropout_ref = apply(data_brca,2,function(x) mean(x==0)*100)
hist(dropout_ref)
hist(dropout)

suspicious=which(dropout>=60)
load("human_gsm_meta.rda")
meta=gsmMeta[seq(ncol(data_archs4))]
rm(gsmMeta)
suspicious_meta = meta[suspicious]
suspicious_list = sapply(suspicious_meta, function(x) x$Sample_series_id[1])
suspicious_list = unique(unname(suspicious_list))
write.csv(suspicious_list, file = 'GSE_list_dropout>=60.csv')
#write.csv(names(dropout), file = 'GSM_list.csv')

# Soluce 1
# Go to Python to open the 39 geo urls

# To remove :
data_rm <- list("GSE51254"=sort(names(meta[sapply(meta,
                                                  function(x) c("GSE51254" %in% x$Sample_series_id))]))
                [-c(1,2)],
                "GSE38265"="GSM937713",
                "GSE36552"=names(meta[sapply(meta,
                                             function(x) c("GSE36552" %in% x$Sample_series_id))]),
                "GSE44618"=names(meta[sapply(meta,
                                             function(x) c("GSE44618" %in% x$Sample_series_id))]),
                "GSE49321"=names(meta[sapply(meta,
                                             function(x) c("GSE49321" %in% x$Sample_series_id))]),
                "GSE44183"=names(meta[sapply(meta,
                                             function(x) c("GSE44183" %in% x$Sample_series_id))])
                )
data_rm <- unlist(unname(data_rm))

data_archs4 <- data_archs4[,!c(colnames(data_archs4) %in% data_rm)]
meta <- meta[colnames(data_archs4)]
dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
#ata_archs4 <- data_archs4[,dropout<60]
#dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
hist(dropout)

write.csv(unique(unname(sapply(meta,
                               function(x) x$Sample_series_id[1]))), file = 'GSE_list.csv')


# Soluce 2
# Go to Python to open the 252 geo urls

# To remove :
#data_rm <- list()
#data_rm <- unlist(unname(data_rm))

#data_archs4 <- data_archs4[,!c(colnames(data_archs4) %in% data_rm)]
#dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
#hist(dropout)

# Soluce 3
# Check the meta$Sample_title slot for sRNA, lncRNA, siRNA, CLIP, shRNA, miRNA, etc
dropout[grep("sRNA", sapply(meta,
                  function(x) x$Sample_title))]
dropout[grep("miRNA", sapply(meta,
                             function(x) x$Sample_title))]
dropout[grep("ribosome", sapply(meta,
                                    function(x) x$Sample_title))]
dropout[grep("CLIP", sapply(meta,
                            function(x) x$Sample_title))]
dropout[grep(".rp", sapply(meta,
                            function(x) x$Sample_title))]
dropout[grep("Chr-input", sapply(meta,
                                  function(x) x$Sample_title))]
dropout[grep("small RNA", sapply(meta,
                            function(x) x$Sample_title))]
dropout[grep("^BR", sapply(meta,
                          function(x) x$Sample_title))]
dropout[grep("_3seq", sapply(meta,
                           function(x) x$Sample_title))]
dropout[grep("^CSHL", sapply(meta,
                             function(x) x$Sample_title))]
dropout[grep("METSIM", sapply(meta,
                             function(x) x$Sample_title))]
dropout[grep("-ip", sapply(meta,
                              function(x) x$Sample_title))]
dropout[grep("smallRNA", sapply(meta,
                                function(x) x$Sample_title))]
data_rm <- sort(unique(c(grep("sRNA", sapply(meta,
                                     function(x) x$Sample_title)),
                 grep("miRNA", sapply(meta,
                                      function(x) x$Sample_title)),
                 grep("ribosome", sapply(meta,
                                         function(x) x$Sample_title)),
                 grep("CLIP", sapply(meta,
                                     function(x) x$Sample_title)),
                 grep(".rp", sapply(meta,
                                    function(x) x$Sample_title)),
                 grep("Chr-input", sapply(meta,
                                          function(x) x$Sample_title)),
                 grep("small RNA", sapply(meta,
                                          function(x) x$Sample_title)),
                 grep("^BR", sapply(meta,
                                    function(x) x$Sample_title)),
                 grep("_3seq", sapply(meta,
                                      function(x) x$Sample_title)),
                 grep("^CSHL", sapply(meta,
                                      function(x) x$Sample_title)),
                 grep("METSIM", sapply(meta,
                                       function(x) x$Sample_title)),
                 grep("-ip", sapply(meta,
                                    function(x) x$Sample_title)),
                 grep("smallRNA", sapply(meta,
                                    function(x) x$Sample_title)))))
data_archs4 <- data_archs4[,-data_rm]
meta <- meta[colnames(data_archs4)]
dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
hist(dropout)

data_rm <- list("GSE37686"=names(meta[sapply(meta,
                                             function(x) c("GSE37686" %in% x$Sample_series_id))]),
                "GSE44039"=names(meta[sapply(meta,
                                             function(x) c("GSE44039" %in% x$Sample_series_id))]),
                "GSE40819"=names(meta[sapply(meta,
                                             function(x) c("GSE40819" %in% x$Sample_series_id))]),
                "GSE43550"=names(meta[sapply(meta,
                                             function(x) c("GSE43550" %in% x$Sample_series_id))]),
                "GSE40710"=names(meta[sapply(meta,
                                             function(x) c("GSE40710" %in% x$Sample_series_id))]),
                "GSE39118"=names(meta[sapply(meta,
                                             function(x) c("GSE39118" %in% x$Sample_series_id))]),
                "GSE46579"=names(meta[sapply(meta,
                                             function(x) c("GSE46579" %in% x$Sample_series_id))])
)
data_rm <- unlist(unname(data_rm))
data_archs4 <- data_archs4[,!c(colnames(data_archs4) %in% data_rm)]
meta <- meta[colnames(data_archs4)]
dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
hist(dropout)

write.csv(names(dropout)[dropout>=60], file = 'GSM_list_dropout>=60.csv')

# Last passage manually
data_rm <- c("GSM1207923", "GSM1226164", "GSM1065157", "GSM1226163",
             "GSM862356", "GSM862354", "GSM1067863", "GSM1219398",
             "GSM1071442", "GSM1071443", "GSM862357", "GSM1067861",
             "GSM1219399")
data_archs4 <- data_archs4[,!c(colnames(data_archs4) %in% data_rm)]
meta <- meta[colnames(data_archs4)]
dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
hist(dropout)

data_rm <- sort(unique(grep("single", sapply(meta,
                                             function(x) x$Sample_title))))
data_archs4 <- data_archs4[,-data_rm]
meta <- meta[colnames(data_archs4)]
dropout = apply(data_archs4,2,function(x) mean(x==0)*100)
hist(dropout)

write.csv(data_archs4, file="data.csv")

