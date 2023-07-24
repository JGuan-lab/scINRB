# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data


## 1.Introduction

scINRB, a single-cell gene expression imputation method with network regularization and bulk RNA-seq data, adopts network-regularized non-negative matrix factorization to decompose single-cell gene expression matrix into low-dimensional gene-factor and sample-factor matrices, ensuring that the imputed data maintains the original cell-cell and gene-gene similarities and approaches the gene average expression calculated from bulk RNA-seq data.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.7501990

## 2.Installation

### Depends:
    R (>= 4.0.4) 
### Requirements:
    library(MASS)
    source('scINRB.R')#“预装“函数
    source('functions.R')#“预装“函数
## 3.Quick start
### 3.1 Prepare data:
The input data includes the input dropout scRNA-seq data and the bulk RNA-seq data
    ...

### 3.2 The parameters used in scINRB:
parameter: the vector of regularization parameters, the default is (0.001,0.001,1).    
r: dimensions of low-dimensional matrix, the default is 200. 
Parameters (including three regularization parameters and the number of factors r) can be selected by cross-validation.
    ...
    
    
### 3.3 Run scINRB:
    #run_demo.R
    setwd(path)
    source('scINRB.R')
    source('functions.R')
    data_sc <- read.csv("data/data_sc.csv",row.names=1)
    data_bulk <- read.csv("data/data_bulk.csv",row.names=1)
    parameter <- c(0.001,0.001,1) 
    r <- 200
    result <- scINRB(data,data_bulk,parameter,r)
    write.csv(result[[3]], file="scINRB_matrix.csv")

 
