#install.packages("corrplot")
#install.packages("pheatmap")
rm(list=ls())
library("pheatmap")
library("corrplot")
library("RColorBrewer")
library("reshape2")
library("ggplot2")
library(gridExtra)
library(lattice)
library(aplot)
library(tidyr)
library(Metrics)
path="/Users/kangkangyue/Desktop/scINRB" #声明test2.R所在位置
setwd(path)  #把工作路径设置到path
source('functions.R')#“预装“函数
source('analysis_plot_methods.R')

# gene-gene correlation&cell-cell correlation
pss_c<-list()
pss_g<-list()

for(scale_num in c(3)){

  k=0
  for(change_rate in c(20,30,40,50,60)){
    k=k+1
    data_cor = get_cor_data(change_rate, scale_num)

    data_cell = data_cor[[1]]

    data_gene = data_cor[[2]]

    #相关矩阵转向量
    #as.vector(unlist(x))  矩阵转向量
    c_true<- as.vector(unlist(data_cell[[1]]))
    c_dropout <- as.vector(unlist(data_cell[[2]]))
    c_FPimpute <- as.vector(unlist(data_cell[[3]]))
    c_drimpute <- as.vector(unlist(data_cell[[4]]))
    c_scimpute <- as.vector(unlist(data_cell[[5]]))
    c_scrabble <- as.vector(unlist(data_cell[[6]]))
    c_viper <- as.vector(unlist(data_cell[[7]]))
    c_SAVER <- as.vector(unlist(data_cell[[8]]))
    c_magic <- as.vector(unlist(data_cell[[9]]))
    c_CMF <- as.vector(unlist(data_cell[[10]]))



    g_true<- as.vector(unlist(data_gene[[1]]))
    g_dropout <- as.vector(unlist(data_gene[[2]]))
    g_FPimpute <- as.vector(unlist(data_gene[[3]]))
    g_drimpute <- as.vector(unlist(data_gene[[4]]))
    g_scimpute <- as.vector(unlist(data_gene[[5]]))
    g_scrabble <- as.vector(unlist(data_gene[[6]]))
    g_viper <- as.vector(unlist(data_gene[[7]]))
    g_SAVER <- as.vector(unlist(data_gene[[8]]))
    g_magic <- as.vector(unlist(data_gene[[9]]))
    g_CMF <- as.vector(unlist(data_gene[[10]]))

    
    name<-list('true','dropout','FPimpute','drimpute','scimpute','scrabble','viper','SAVER','magic','CMF')
    #计算相关矩阵的Pearson相似度
    cell_pss <- vector()
    cell_pss[1]<-cor(c_true,c_dropout)
    cell_pss[2]<-cor(c_true,c_FPimpute)
    cell_pss[3]<-cor(c_true,c_drimpute)
    cell_pss[4]<-cor(c_true,c_scimpute)
    cell_pss[5]<-cor(c_true,c_scrabble)
    cell_pss[6]<-cor(c_true,c_viper)
    cell_pss[7]<-cor(c_true,c_SAVER)
    cell_pss[8]<-cor(c_true,c_magic)
    cell_pss[9]<-cor(c_true,c_CMF)

    
    gene_pss <- vector()
    gene_pss[1]<-cor(g_true,g_dropout)
    gene_pss[2]<-cor(g_true,g_FPimpute)
    gene_pss[3]<-cor(g_true,g_drimpute)
    gene_pss[4]<-cor(g_true,g_scimpute)
    gene_pss[5]<-cor(g_true,g_scrabble)
    gene_pss[6]<-cor(g_true,g_viper)
    gene_pss[7]<-cor(g_true,g_SAVER)
    gene_pss[8]<-cor(g_true,g_magic)
    gene_pss[9]<-cor(g_true,g_CMF)
    
    pss_c[[k]]<-cell_pss
    pss_g[[k]]<-gene_pss

  }
  print('数据')
  print(seed_value)
  print('cell pearson')
  print(pss_c)
  print('gene pearson')
  print(pss_g)
  print('___________')
  
}



#Imputded data and true data 矩阵相似性


for(scale_num in c(3)){
  cor_imputed_data<-vector()
  cor_imputed_data_all <-list()
  k=0
  for(change_rate in c(20,30,40,50,60)){
    k=k+1
    cor_imputed_data<- get_cor_imputed_data(change_rate, scale_num)
    cor_imputed_data_all[[k]] <- cor_imputed_data
   
    
    
  }

  print('数据')
  print(seed_value)
  print('impution data pearson')
  print(cor_imputed_data_all)
  print('___________')

}

#Imputded data and true data RMSE
for(scale_num in c(3)){
  k=0
  RMSE_imputed_data_all <-list()
  for(change_rate in c(20,30,40,50,60)){
    k=k+1
    RMSE_imputed_data<- get_RMSE_imputed_data(change_rate, scale_num)
    RMSE_imputed_data_all[[k]] <-RMSE_imputed_data
    
  }
  print('数据')
  print(scale_num)
  print('impution data pearson')
  print(RMSE_imputed_data_all)
  print('___________')
  

  
}


