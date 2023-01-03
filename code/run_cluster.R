Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
library(Rcpp)
Rcpp::sourceCpp("/Users/kangkangyue/Desktop/scINRB//julei.cpp")
source('/Users/kangkangyue/Desktop/scINRB//analysis_methods.R')
library(devtools)
library(fpc)
library(mclust)
library(ngram)

#########Read the simulation data######
data <- readRDS("/Users/kangkangyue/Desktop/scINRB/data/simulation_data/3_60%.rds")
data_dropout <- as.matrix(data$data_dropout)
data_true <- as.matrix(data$data_true)
data_bulk <- as.matrix(data$data_bulk)
data_label <- as.matrix(data[["group"]])
impute_data <- read.csv("/Users/kangkangyue/Desktop/scINRB/imputation_knn_smoothing_data/knn_smoothing_1_60%.csv", header=FALSE)
impute_data <- as.matrix(impute_data)

########Read the true data#######
# data_label <- as.matrix(read.csv("/Users/kangkangyue/Desktop/scINRB/data/five_encode/data_label.csv",row.names=1))
# impute_data <-  read.csv("/Users/kangkangyue/Desktop/论文/kangyue/LPLS9.28/true_data/knn_smoothing/knn_smoothing_混合数据.csv",row.names=1)
# impute_data <- as.matrix(impute_data)
######### 100 averages kmeans clustering#############
H<-t(impute_data)
realcluster <- data_label
ari<-vector(length=100)
sc<-vector(length=100)
nmi<-vector(length=100)
for(j in c(1:100)){
  cluster <- kmeans(H,3)[["cluster"]]
  test=cluster_evalu(H,cluster,realcluster)
  ari[j]<-test[[1]]
  sc[j]<-test[[2]]
  nmi[j]<-test[[3]]
}
cluster_list<-list(ari,sc,nmi)
names(cluster_list)=c('ari','sc','nmi')
ari_ave=mean(ari)
sc_ave=mean(sc)
nmi_ave=mean(nmi)
#test=cluster_evalu(H,cluster,realcluster)
print(ari_ave)
print(sc_ave)
print(nmi_ave)

#write.csv(cluster_list,"/Users/kangkangyue/Desktop/scINRB/imputation_scINDR_data/cluster_3_60%.csv")

######### 100 averages of tsne clustering#############
tsne.coords <- Rtsne(t(impute_data),pca=FALSE,perplexity = 30,theta=0.5,check_duplicates = F)$Y
rownames(tsne.coords) <- rownames(t(impute_data))
colnames(tsne.coords) <- c("tSNE1","tSNE2")
head(tsne.coords)
# Clustering effect after imputation
ari=0
sc=0
nmi=0
for(j in c(1:100)){
  cluster <- kmeans(tsne.coords,2)[["cluster"]]
  test=cluster_evalu(tsne.coords,cluster,realcluster)
  ari=ari+test[[1]]
  sc=sc+test[[2]]
  nmi=nmi+test[[3]]
}
ari_ave=ari/100
sc_ave=sc/100
nmi_ave=nmi/100
#test=cluster_evalu(H,cluster,realcluster)
print(ari_ave)
print(sc_ave)
print(nmi_ave)
# after imputation
#a<-adjustedRandIndex(kmeans(tsne.coords, centers = 3)$cluster, data_label)



######## simulation data kemans clustering
source('/Users/kangkangyue/Desktop/scINRB/analysis_methods.R')
path="/Users/kangkangyue/Desktop/scINRB" 
setwd(path)  
method="deepimpute"

for(scale_num in c(2,3,5,6)){
  print(scale_num)
  for(change_rate in c(60,50,40,30,20)){
    
    result <- get_simulation_data1(method,change_rate,scale_num)
    print(change_rate)
    print(result)
    
  }
}

########simulation data tsne clustering
rm(list=ls()) 
library(Rtsne)
library(ggplot2)
library(devtools)
library(mclust)
library(fpc)
path="/Users/kangkangyue/Desktop/scINRB" 
setwd(path)  
source('analysis_methods.R')
method="deepimpute"

for(scale_num in c(1:5)){
  print(scale_num)
  for(change_rate in c(60,50,40,30,20)){
    
    result <- tsne_simulation_data(method,change_rate, scale_num)
    print(change_rate)
    print(result)
    
  }
  
}

########simulation data : matrix_RMSE\matrix_psss\cell_pss\gene_pss

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
path="/Users/kangkangyue/Desktop/scINRB" 
setwd(path)  
source('methods.R')
source('analysis_methods.R')
method="deepimpute"

for(scale_num in c(1:5)){
  
  print(scale_num)
  
  for(change_rate in c(20,30,40,50,60)){
    
    # result[[1]] = RMSE
    # 
    # result[[2]] = matrix_pss
    # 
    # result[[3]] = cell_pss
    # 
    # result[[4]] = gene_pss
    
    print(change_rate)
    
    print(result)

  }
  
}

