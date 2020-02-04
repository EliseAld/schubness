Hubness_sc
====

Hubness on single-cell RNA sequencing data

Goals
-----------
1. Do we have hubs in scRNAseq ?
1a. How to define hubs ?
1b. Are hubs & antihubs different ?
1c. Are hubs sensitive to dropout ?
1d. What happens if we remove hubs ?

2. Can we use hubs to enhance the clustering analysis ?
2a. Choosing a better suited norm
2b. Removing hubs

Previous steps
-----------
1. Compute the hubness scores for different dimensions of input space (PCA with various number of PCs, initial space), using the data of Guo et al., 2018 (https://www.ncbi.nlm.nih.gov/pubmed/29942094)
:white_check_mark:

2. Demonstrate the progressive appearance of hubs with increasing dimensions
:white_check_mark:

3. Compare the sensitivity of the Lp norms (p comprised between 0.1 and 4) to the dimensionality using the hubs
:white_check_mark:

4. Evaluate the difference in sensitivity as a function of the dimensionality
:white_check_mark:

Current steps
-----------
1. Increase the rank of the k-parameter (number of neighbors, from 5 to 200), and the range of PCA-space dimensions (from 2 to 9k)  :white_check_mark:

2. Re compute the PCA-transformed matrix
:white_check_mark:

3. Try other methods to compare the distributions (max, nb of hubs, hubs taken as having a nb of k-occurences >= 2k [1], skewness [1], etc)
:white_check_mark:

4. Try other type of data (10X, ATAC, etc)

5. Check the transcriptome of "normal cells" vs "hubs" and "anti-hubs"

4. Check the literature for usage of hubness in sc data
:white_check_mark:

Bibliography
-----------
[1] https://www.sciencedirect.com/science/article/pii/S0925231215004336

[2] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6284580/
