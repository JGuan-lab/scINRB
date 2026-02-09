# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data

建模：针对基因调控网络构建中基因表达信息缺失及不准确等复杂动态系统建模难题，提出了复杂耦合关联系统的建模方法，建立了融合数据和知识的复杂动态网络系统模型，提升了基因调控网络动态系统模型的准确性和鲁棒性。

## 1. Introduction

scINRB, a single-cell gene expression imputation method with network regularization and bulk RNA-seq data, adopts network-regularized non-negative matrix factorization to decompose single-cell gene expression matrix into low-dimensional gene-factor and sample-factor matrices, ensuring that the imputed data maintains the original cell-cell and gene-gene similarities and approaches the gene average expression calculated from bulk RNA-seq data.

The datasets analyzed in the paper are available at: https://zenodo.org/record/8224512

scINRB corresponds to the following paper:

Yue Kang, Hongyu Zhang, Jinting Guan, scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data, Briefings in Bioinformatics, 25(3): bbae148. https://doi.org/10.1093/bib/bbae148

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

 
