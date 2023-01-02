Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
library(Rcpp)
Rcpp::sourceCpp("/Users/kangkangyue/Desktop/scINRB//julei.cpp")
source('/Users/kangkangyue/Desktop/scINRB//analysis_methods.R')
library(devtools)
library(fpc)
library(mclust)
library(ngram)

#########模拟数据读入######
data <- readRDS("/Users/kangkangyue/Desktop/scINRB/data/simulation_data/3_60%.rds")
data_dropout <- as.matrix(data$data_dropout)
data_true <- as.matrix(data$data_true)
data_bulk <- as.matrix(data$data_bulk)
data_label <- as.matrix(data[["group"]])
impute_data <- read.csv("/Users/kangkangyue/Desktop/scINRB/imputation_knn_smoothing_data/knn_smoothing_1_60%.csv", header=FALSE)
impute_data <- as.matrix(impute_data)

########真实数据读入#######
# data_label <- as.matrix(read.csv("/Users/kangkangyue/Desktop/scINRB/data/five_encode/data_label.csv",row.names=1))
# impute_data <-  read.csv("/Users/kangkangyue/Desktop/论文/kangyue/LPLS9.28/true_data/knn_smoothing/knn_smoothing_混合数据.csv",row.names=1)
# impute_data <- as.matrix(impute_data)
######### kmeans聚类100次均值#############
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

######### tsne聚类100次均值#############
# tSNE降维结果 
# Rtsne输入一个表达量的表格就可以了，每一行为一个细胞，每一列为一个基因
# dims参数设置降维之后的维度，默认值为2，
# perplexity参数的取值必须小于(nrow(data) - 1 )/ 3;
# theta参数取值越大，结果的准确度越低，默认值为0.5，
# max_iter参数设置最大迭代次数。
# pca参数表示是否对输入的原始数据进行PCA分析，然后使用PCA得到的topN主成分进行后续分析
# t-SNE计算量特别大，对于维度较高数据，先采用PCA降维可提高运行效率，默认采用top50的主成分进行后续分析，可通过initial_dims修改这个值。
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



########批量模拟数据kemans聚类 单个方法
source('/Users/kangkangyue/Desktop/scINRB/analysis_methods.R')
path="/Users/kangkangyue/Desktop/scINRB" #声明test2.R所在位置
setwd(path)  #把工作路径设置到path
method="deepimpute"

for(scale_num in c(2,3,5,6)){
  print(scale_num)
  for(change_rate in c(60,50,40,30,20)){
    
    result <- get_simulation_data1(method,change_rate,scale_num)
    print(change_rate)
    print(result)
    
  }
}

########批量模拟数据tsne聚类 单个方法
rm(list=ls()) #清除当前环境中的变量
library(Rtsne)
library(ggplot2)
library(devtools)
library(mclust)
library(fpc)
path="/Users/kangkangyue/Desktop/scINRB" #声明test2.R所在位置
setwd(path)  #把工作路径设置到path
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

########批量模拟数据pss 单个方法

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
source('methods.R')#“预装“函数
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

