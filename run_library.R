#################################################################################################################
#run simulation data
#################################################################################################################

######## SAVER ##############

run_SAVER_sim <- function(change_rate, scale_num)
{

  # load the raw data
  data <- readRDS(file = paste0('data/simulation_data/',
                                scale_num,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  
  # build the folder saving the imputed data using SAVER
  path <- "imputation_SAVER_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using SAVER
  data_dropout <- as.matrix(data$data_dropout)

  exdata <- saver(data_dropout, ncores = 12, estimates.only = TRUE)
  
  # write the data
  write.table(exdata,
              paste0(path, "SAVER_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
  
}

######## scINRB ##############

run_scINRB_sim <- function(change_rate, scale_num)
{
  # Parameter in the function
  # load the raw data
  data <- readRDS(file = paste0('data/simulation_data/',
                                seed_value,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  
  # build the folder saving the imputed data using scINRB
  path <- "imputation_scINRB_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using scINRB
  data_dropout <- as.matrix(data$data_dropout)
  data_bulk <- as.matrix(data$data_bulk)
  parameter <- c(0.001,0.001,1)
  r <- 200
  exdata <- scINRB(data_dropout,data_bulk,parameter,r)
  
  # write the data
  write.table(exdata,
              paste0(path, "scINRB_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
  
}


######## DrImpute ##############
run_drimpute_sim <- function(change_rate, scale_num){


  # load the raw data
  data <- readRDS(file = paste0('data/simulation_data/',
                                scale_num,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  # build the folder saving the imputed data using DrImpute
  path <- "imputation_drimpute_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data_dropout <- as.matrix(data$data_dropout)
  
  exdata <- DrImpute(data_dropout)
  
  # write the data
  write.table(exdata,
              paste0(path, "drimpute_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## scImpute ##############
run_scimpute_sim <- function(change_rate, scale_num){
  

  # load the data
  data <- readRDS(file = paste0('data/simulation_data/',
                                scale_num,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  
  # build the folder saving the imputed data using scimpute method
  path = "imputation_scimpute_data/"
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using scImpute
  data_dropout = data$data_dropout
  
  write.table(data_dropout,paste0(path, 
                                  "dropout_scimpute_",
                                  change_rate,"_",
                                  scale_num,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  file.remove(paste0(path, "scimpute_", change_rate, "_", scale_num,"_*"))
  
  # run scImpute
  scimpute(
    paste0(path, "dropout_scimpute_",change_rate,"_",scale_num,".csv"),
    infile = "csv",           
    outfile = "csv",         
    out_dir = paste0(path, "scimpute_", change_rate, "_", scale_num,"_"),
    drop_thre = 0.5,         
    Kcluster = 2,
    ncores = 2)       
  # 
  # clean the data
  data_dropout = read.table( file = paste0(path, "scimpute_",
                                           change_rate, "_",
                                           scale_num,
                                           "_scimpute_count.csv") ,
                             header = TRUE, sep=",")
  
  data_dropout$X = NULL
  
  # save the data
  write.table(data_dropout,
              paste0(path, "scimpute_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}


######## VIPER ##############
run_viper_sim <- function(change_rate, scale_num){


  # load the raw data
  data <- readRDS(file = paste0('data/simulation_data/',
                                scale_num,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  # build the folder saving the imputed data using VIPER
  path <- "imputation_viper_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using VIPER
  data_dropout <- as.matrix(data$data_dropout)
  
  exdata <- VIPER(data_dropout, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                  report = FALSE, outdir = NULL, prefix = NULL)
  
  # write the data
  write.table(exdata$imputed,
              paste0(path, "viper_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## SCRABBLE ##############
run_scrabble_sim <- function(change_rate, scale_num){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the data
  data <- readRDS(file = paste0('data/simulation_data/',
                                seed_value,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  
  
  path = "imputation_scrabble_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using SCRABBLE
  data1 = list()
  
  data1[[1]] = data$data_dropout
  
  data1[[2]] = data$data_bulk
  
  # set up the parameters
  parameter = c(1, 1e-06, 1e-04)
  
  # run scrabble
  result = scrabble(data1,
                    parameter = parameter, 
                    nIter = 60,
                    error_out_threshold = 1e-7, 
                    nIter_inner = 100,
                    error_inner_threshold = 1e-5)
  
  # write the data
  write.table(result,
              paste0(path, "scrabble_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F)
  
}



######## knn_smoothing ##############
run_knn_smoothing_sim <- function(change_rate, scale_num){
  
  # load the data

  data <- readRDS(file = paste0('data/simulation_data/',
                                scale_num,
                                '_',
                                change_rate,
                                '%.rds')
  )
  
  
  
  path = "imputation_knn_smoothing_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using knn_smoothing
  data_dropout <- as.matrix(data$data_dropout)
  exdata <- knn_smoothing(data_dropout,k=10)
  
  # write the data
  write.table(exdata,
              paste0(path, "knn_smoothing_",scale_num,'_',change_rate,"%.csv"),
              sep=',',
              row.names = F,
              col.names = F)
  
}


#################################################################################################################
#run true data
#################################################################################################################





######## scINRB ##############

run_scINRB <- function(data_sc,data_bulk,data_name)
{
  
  # build the folder saving the imputed data using scINRB
  path <- "/imputation_true_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using scINRB
  
  parameter <- c(0.001,0.001,1)
  r <- 200
  exdata <- scINRB(data_sc,data_bulk,parameter,r)
  
  # write the data
  # write.table(exdata[1],
  #             paste0(path, "W_",data_name,".csv"),
  #             sep=',',
  #             row.names = F,
  #             col.names = F
  # )
  # write.table(exdata[2],
  #             paste0(path, "H_",data_name,".csv"),
  #             sep=',',
  #             row.names = F,
  #             col.names = F
  # )
  write.table(exdata[3],
              paste0(path, "scINRB_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}
######## bayNorm ##############
run_bayNorm <- function(data_sc,data_name){
  
  
  # build the folder saving the imputed data using DrImpute
  path <- "/imputation_true_data/"
  #dir.create(file.path(path), showWarnings = FALSE)
  data_sc <- as.matrix(data_sc)
  
  data <- SummarizedExperiment::SummarizedExperiment(assays=list(Counts=data_sc))
  exdata <- bayNorm(
    Data=data,
    BETA_vec = NULL,
    mode_version=TRUE,
    mean_version = FALSE,S=20
    ,verbose =FALSE,
    parallel = TRUE)
  
  # write the data
  write.table(exdata,
              paste0(path, "bayNorm_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## knn_smoothing ##############
run_knn_smoothing <- function(data_sc,data_name){
  
  
  # build the folder saving the imputed data using DrImpute
  path <- "/imputation_true_data/"
  #dir.create(file.path(path), showWarnings = FALSE)
  
  data_sc <- as.matrix(data_sc)
  
  exdata <- knn_smoothing(data_sc,k=10)
  
  # write the data
  write.table(exdata,
              paste0(path, "knn_smoothing_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}
######## DrImpute ##############
run_drimpute <- function(data_sc,data_name){
  
  
  # build the folder saving the imputed data using DrImpute
  path <- "/imputation_true_data/"
  
  data_sc <- as.matrix(data_sc)
  
  
  exdata <- DrImpute(data_sc)
  
  # write the data
  write.table(exdata,
              paste0(path, "drimpute_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## scImpute ##############
run_scimpute <- function(data_sc,data_name){
  
  
  # build the folder saving the imputed data using scimpute method
  path = "/imputation_true_data/"
  
  
  # impute the data using scImpute
  
  write.table(data_sc,paste0(path, 
                             "dropout_scimpute_",
                             data_name,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  file.remove(paste0(path, "scimpute_",data_name,"_*"))
  
  # run scImpute
  scimpute(
    paste0(path, "dropout_scimpute_",data_name,".csv"),
    infile = "csv",           
    outfile = "csv",         
    out_dir = paste0(path, "scimpute_", data_name,"_"),
    drop_thre = 0.5,         
    Kcluster = 2,
    ncores = 1)       
  # 
  # clean the data
  data_dropout = read.table( file = paste0(path, "scimpute_",
                                           data_name,
                                           "_scimpute_count.csv") ,
                             header = TRUE, sep=",")
  
  data_dropout$X = NULL
  
  # save the data
  write.table(data_dropout,
              paste0(path, "scimpute_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}



######## VIPER ##############
run_viper <- function(data_sc,data_name){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the raw data
  
  # build the folder saving the imputed data using DrImpute
  path <- "/imputation_true_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using VIPER
  
  
  exdata <- VIPER(data_sc, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                  report = FALSE, outdir = NULL, prefix = NULL)
  
  # write the data
  write.table(exdata$imputed,
              paste0(path, "viper_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## SCRABBLE ##############
run_scrabble <- function(data_sc,data_bulk,data_name){
  
  
  path = "/imputation_true_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data1 = list()
  
  data1[[1]] = data_sc
  
  data1[[2]] = data_bulk
  
  # set up the parameters
  parameter = c(1, 1e-06, 1e-04)
  
  # run scrabble
  result = scrabble(data1,
                    parameter = parameter, 
                    nIter = 60,
                    error_out_threshold = 1e-7, 
                    nIter_inner = 100,
                    error_inner_threshold = 1e-5)
  
  # write the data
  write.table(result,
              paste0(path, "scrabble_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F)
  
}


run_SAVER <- function(data_sc,data_name)
{
  
  
  
  # build the folder saving the imputed data using DrImpute
  path <- "/imputation_true_data/"
  
  #dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  
  
  exdata <- saver(data_sc, ncores = 12, estimates.only = TRUE)
  
  # write the data
  write.table(exdata,
              paste0(path, "SAVER_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
  
}



run_SAVERX <- function(data_sc,data_name)
{
  
  
  path <- "/imputation_true_data/"
  
  
  # impute the data using SAVERX
  #exdata <- saverx(paste0('~/Desktop/真实数据/1/data/processed/cellbench/sc_10x/data_sc.csv'))
  
  exdata <- saverx("/Users/kangkangyue/Desktop/真实数据/1/data/混合数据/data_sc.csv")
  # write the data
  write.table(exdata,
              paste0(path, "SAVERX_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
  
}


run_scRecover <- function(data_sc,data_name)
{
  
  
  path <- "/imputation_true_data/"
  
  
  # impute the data using SAVERX
  exdata <- scRecover(counts = data_sc, Kcluster = 2, outputDir = paste0(commandArgs(trailingOnly = T)[2],'/'), verbose = FALSE)
  
  # write the data
  write.table(exdata,
              paste0(path, "scRecover_",data_name,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
  
}






######## MAGIC ##############
# here we use python script to run magic and we post the python script here.
# -----------------------------------------------------------------------------
# import sys
# import os
# import magic 
# import pandas as pd
# 
# cwd = os.getcwd()
# 
# if not os.path.exists(cwd+"/imputation_magic_data"):
#   os.makedirs(cwd+"/imputation_magic_data")
#
# X =pd.read_csv("simulation_data/simulation_data_dropout_index_"+str(dropout_value)+"_seed_"+str(seed_value)+".txt",sep = ' ', header=None)
#
# magic_operator = magic.MAGIC()
# X_magic = magic_operator.fit_transform(X.T)
#
# out_magic = X_magic.T
# out_magic.to_csv(cwd+"/imputation_magic_data/magic_"+str(dropout_value)+"_"+str(seed_value)+".csv", sep = '\t', header= None)
# -----------------------------------------------------------------------------


