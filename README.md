# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data


## Introduction

For many scRNA-seq data, existing batch data usually exists on the same cell or tissue. Compared with using single-cell data alone, batch RNA-seq data can more accurately estimate the distribution of gene expression between cells. 

Therefore, we use the information of batch RNA-seq data, use network regularized non-negative matrix decomposition to decompose the transcriptional count matrix into two low-dimensional matrices, and propose an interpolation method scINRB to correct the measurement of gene expression in scRNA-seq.

## Running the tests

### Depends:
    R (>= 4.0.4) 
### The parameter used in scINRB:
    data : the input dropout scRNAseq data.
    
    data_bulk : the bulk RNAseq.

    parameter : the vector of parameters , the default is (0.001,0.001,1).
    The first parameter is the value of alpha1 in the mathematical model;
    the second one is the value of alpha2 in the mathematical model;
    the third one is the value of alpha3 in the mathematical model.
    
    r : dimensions of matrix decomposition , the default is 200.  
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

