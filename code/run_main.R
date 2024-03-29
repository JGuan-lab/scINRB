# load the libraries library(SAVER)
library(scRecover)
library(SAVERX)
library(SAVER)
library(scImpute)
library(SCRABBLE)
library(DrImpute)
library(viper)
path="../scINRB" 
setwd(path)  
source("run_library.R")
source('scINRB.R')
source('functions.R')


# data in the current folder /simulation_data/
# The following code is used to impute simulation data in batches using different methods.


for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_SAVER_sim(change_rate, scale_num)
    
  }
  
}

for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_scINRB_sim(change_rate, scale_num)
    
  }
  
}

for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_drimpute_sim(change_rate, scale_num)
    
  }
  
}


for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_scimpute_sim(change_rate, scale_num)
    
  }
  
}

for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_viper_sim(change_rate, scale_num)
    
  }
  
}

for(change_rate in c(20,30,40,50,60) ){
  
  for(scale_num in c(1:5)){
    
    run_scrabble_sim(change_rate, scale_num)
    
  }
}

# The following code is used to impute true data using different methods.

data_sc <- read.csv("../scINRB/data/true_data/five_encode/data_sc.csv",row.names=1)
data_bulk <- read.csv("../scINRB/data/true_data/five_encode/data_bulk.csv",row.names=1)
data_label <- read.csv("../scINRB/data/true_data/five_encode/data_label.csv",row.names=1)
data_name <- 'five_encode'


run_scINRB(data_sc,data_bulk,data_name)
run_drimpute(data_sc,data_name)
run_scimpute(data_sc,data_name)
run_viper(data_sc,data_name)
run_scrabble(data_sc,data_bulk,data_name)
run_SAVER(data_sc,data_name)



