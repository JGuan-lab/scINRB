# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data


## 1. Introduction

scINRB, a single-cell gene expression imputation method with network regularization and bulk RNA-seq data, adopts network-regularized non-negative matrix factorization to decompose single-cell gene expression matrix into low-dimensional gene-factor and sample-factor matrices, ensuring that the imputed data maintains the original cell-cell and gene-gene similarities and approaches the gene average expression calculated from bulk RNA-seq data.

The datasets analyzed in the paper are available at: https://zenodo.org/record/8224512

## 2. Installation

### Depends:
    R (>= 4.0.4) 
### Requirements:
    library(MASS)
    source('scINRB.R')
    source('functions.R')
## 3.Quick start


### 3.1 Prepare data
The inputs include scRNA-seq data and bulk RNA-seq data.

    data <- readRDS("data/simulation_data/4_20%.rds")
    data_sc0 <- as.matrix(data$data_dropout)
    data_bulk0 <- as.matrix(data$data_bulk)
    result <- preprocess(data_sc0,data_bulk0)
    data_sc <- result[[1]]
    data_bulk <- result[[2]]

### 3.2 Tuning parameters of scINRB
The default value of the regularization parameter vector is (0.001,0.001,1). The default value of the number of factors in low-dimensional space is 200. Parameters (including three regularization parameters and the number of factors r) can be selected by cross-validation.

    cross_validation(data_sc,data_bulk)
    
### 3.3 Run scINRB
    #run_demo.R
    parameter <- c(0.001,0.001,1) 
    r <- 200
    result <- scINRB(data_sc,data_bulk,parameter,r)
    write.csv(result[[3]], file="scINRB_matrix.csv")

 
