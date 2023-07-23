#######Read the data#############
#simulation data
path="../scINRB/data/simulation_data" 
setwd(path)  
data <- readRDS("4_20%.rds")
data_sc <- as.matrix(data$data_dropout)
data_true <- as.matrix(data$data_true)
data_bulk <- as.matrix(data$data_bulk)
#true data
# data_sc <- read.csv("data_sc.csv",row.names=1)
# data_bulk <- read.csv("data_bulk.csv",row.names=1)
###########
source('scINRB.R')
source('functions.R')
library(MASS)
r<-200
parameter <- c(0.001,0.001,1)
result <- scINRB(data_sc,data_bulk,parameter,r)
write.csv(result[[3]], file="scINRB_matrix.csv")
