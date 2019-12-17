Hubness_sc
====

Hubness on single-cell RNA sequencing data

Previous steps
-----------
1. Compute the hubness scores for different dimensions of input space (PCA with various number of PCs, initial space), using the data of Guo et al., 2018 (https://www.ncbi.nlm.nih.gov/pubmed/29942094)

2. Demonstrate the progressive appearance of hubs with increasing dimensions

3. Compare the sensitivity of the Lp norms (p comprised between 0.1 and 4) to the dimensionality using the hubs

4. Evaluate the difference in sensitivity as a function of the dimensionality

Current steps
-----------
1. Increase the ranke of the k-parameter (number of neighbors), and the range of PCA-space dimensions

2. Re compute the PCA-transformed matrix

3. Try other methods to compare the distributions (max, nb of hubs, hubs taken as having a nb of k-occurences above 2k [1], etc)

4. Try other type of data (10X, ATAC, etc)

5. Check the transcriptome of "normal cells" vs "hubs"

4. Check the literature for usage of hubness in sc data

Bibliography
-----------
[1] https://www.sciencedirect.com/science/article/pii/S0925231215004336
