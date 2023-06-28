# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data


## Introduction

scINRB, a single-cell gene expression imputation method with network regularization and bulk RNA-seq data, adopts network-regularized non-negative matrix factorization to decompose single-cell gene expression matrix into low-dimensional gene-factor and sample-factor matrices, ensuring that the imputed data maintains the original cell-cell and gene-gene similarities and approaches the gene average expression calculated from bulk RNA-seq data.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.7501990

## Running the tests

### Depends:
    R (>= 4.0.4) 
### Input data:
    data: the input dropout scRNA-seq data.

    data_bulk: the bulk RNA-seq data.
### The parameters used in scINRB:

    Parameters (including three regularization parameters and the number of factors r) can be selected by cross-validation.
    
    parameter: the vector of regularization parameters, the default is (0.001,0.001,1).    
    r: dimensions of low-dimensional matrix, the default is 200. 
    
### Example:
    #run_demo.R
    setwd(path)
    source('scINRB.R')
    source('functions.R')
    data_sc <- read.csv("data_sc.csv",row.names=1)
    data_bulk <- read.csv("data_bulk.csv",row.names=1)
    parameter <- c(0.001,0.001,1) 
    r <- 200
    result <- scINRB(data,data_bulk,parameter,r)
    write.csv(result[[3]], file="scINRB_matrix.csv")

 
