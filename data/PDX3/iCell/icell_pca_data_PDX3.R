library(bigSCale)
library(Matrix)

data_sparse <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/Ewing_PDX_normalized/iCell/pdx3_merged.mtx"

# Reduce the size of the data by a factor 2
data_pooled <- iCells(data_sparse, target.cells = 5000)
saveRDS(data_pooled, file="/Users/elise/Desktop/GitHub/Hubness_sc/data/Ewing_PDX_normalized/iCell/data_pooled.rds")

# Reduce the size of the data by a factor 10
data_pooled2 <- iCells(data_sparse, target.cells = 1000)
saveRDS(data_pooled2, file="/Users/elise/Desktop/GitHub/Hubness_sc/data/Ewing_PDX_normalized/iCell/data_pooled2.rds")
