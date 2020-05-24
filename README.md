# PJD: R package for Joint and Pairwise Decomposition

## Overview:

Pairwise and Joint Decomposiiton (PJD) is a R package for visualizing biologically structured gene expression matrix environment based on low rank models.

To install this package in R, run the following commands:

```R
library(devtools)
install_github("CHuanSite/PJD")
```

This package implements two categories of algorithms to decompose multiple datasets, pairwise and joint. For each category, there are three available algorithms: Principal Component Analysis (PCA), Independent Component Analysis (ICA) and Nonnegative Matrix Factorization (NMF). For each method, the algorithm takes three arguments, `r dataset`, `r group` and `comp_num`, specifying which datasets to be used, what is the structure among the datasets and what's the dimension for each component.

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

## Pairwise PCA, ICA, NMF
pairPCA_res = pairwisePCA(dataset, group, comp_num)
pairICA_res = pairwiseICA(dataset, group, comp_num)
pairNMF_res = pairwiseNMF(dataset, group, comp_num)

## Joint PCA, ICA, NMF
jointPCA_res = jointPCA(dataset, group, comp_num)
jointICA_res = jointICA(dataset, group, comp_num)
jointNMF_res = jointNMF(dataset, group, comp_num)
```
