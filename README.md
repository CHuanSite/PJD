# PJD: R package for Joint and Pairwise Decomposition

## Overview:

Pairwise and Joint Decomposiiton (PJD) is a R package for visualizing biologically structured gene expression matrix environment based on low rank models.

To install this package in R, run the following commands:

```R
library(devtools)
install_github("CHuanSite/PJD")
```

This package implements two categories of algorithms to decompose multiple datasets, concatenate and joint. For each category, there are three available algorithms: Principal Component Analysis (PCA), Independent Component Analysis (ICA) and Nonnegative Matrix Factorization (NMF). For each method, the algorithm takes three arguments, `dataset`, `group` and `comp_num`, specifying which datasets to be used, what is the structure among the datasets and what's the dimension for each component.

In addition to these two categories, two more algorithms is proposed, which are sequential algorithm, `linkedPCA` and `seqPCA`.

Example usage:

```R
library(PJD)
# Simulation the dataset
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
               
## Specify the structure among the datasets
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)

## Concatenated PCA, ICA, NMF
concatPCA_res = concatPCA(dataset, group, comp_num)
concatICA_res = concatICA(dataset, group, comp_num)
concatNMF_res = concatNMF(dataset, group, comp_num)

## Joint PCA, ICA, NMF
jointPCA_res = jointPCA(dataset, group, comp_num)
jointICA_res = jointICA(dataset, group, comp_num)
jointNMF_res = jointNMF(dataset, group, comp_num)

## seqPCA
seqPCA_res = seqPCA(dataset, group, comp_num)
```

To access the component
```R
pairPCA_res$linked_component_list
```

To access the score
```R
pairPCA_res$score_list
```
