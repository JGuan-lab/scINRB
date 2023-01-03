library(splatter)
library(edgeR)

# Parameter in the function
# seed_value: the random seed
# nGenes: the number of genes in the simulation data.
# nCells: the number of cells in the simulation data.
# dropout_mid: control the dropout rate

seed_value=10
nGenes = 1000
nCells = 5000

# Set up the parameters
params = newSplatParams()
params = setParams(params, list(batchCells = nCells ,
                                nGenes = nGenes,
                                group.prob = c(0.20, 0.35, 0.45),
                                de.prob = c(0.045, 0.045, 0.045),
                                de.facLoc = 0.1,
                                de.facScale = 0.4)
)

# Set up the vector of dropout.mid
dropout_mid = 3.2

# Generate the simulation data using Splatter package
sim = splatSimulateGroups(params,
                          dropout.type = "experiment",
                          dropout.shape = -1,
                          dropout.mid = dropout_mid,
                          seed = seed_value)

# genereate the cpm levels of the true simulation data
data_true = cpm(sim@assays@data@listData[["TrueCounts"]])
#data_dropout=cpm(sim@assays@data@listData[["counts"]])

data_dropout = data_true
# generate the dropout data based on the counts in sim
tf<-counts(sim) == 0
tf1<-as.matrix(tf)
data_dropout[tf1] = 0


# calculate the dropout rate
percentage_zeros = round(nnzero(data_dropout == 0, na.counted = NA)/
                           (dim(data_dropout)[1]*dim(data_dropout)[2])*100)


# generate the bulk RNAseq data
data_bulk = data.frame(val = rowMeans(data_true))



data = list()

data$data_bulk = data_bulk

data$data_dropout = data_dropout

data$data_true = data_true

data$percentage_zeros = percentage_zeros

data$group = colData(sim)@listData$Group

saveRDS(data,'3_50%.rds')

########################################################################
#Calculate the zero yield of truet data

nGenes <- nrow(data_true)
nCells <- ncol(data_true)
count_true <- 0
for(i in c(1:nGenes)){
  for(j in c(1:nCells)){
    if(data_true[i,j]==0)
      count_true <- count_true + 1
  }
}

rate_true <- count_true/(nGenes*nCells)
print(rate_true)



########################################################################
#Calculate the zero yield of dropout data

nGenes <- nrow(data_true)
nCells <- ncol(data_true)
count_dropout <- 0
for(i in c(1:nGenes)){
  for(j in c(1:nCells)){
    if(data_dropout[i,j]==0)
      count_dropout <- count_dropout + 1
  }
}

rate_dropout <- count_dropout/(nGenes*nCells)
print(rate_dropout)

########################################################################
#Calculate the zeroing rate of dropout data relative to real data（dropout rate）
data_differ <- data_true-data_dropout
nGenes <- nrow(data_differ)
nCells <- ncol(data_differ)
count <- 0
for(i in c(1:nGenes)){
  for(j in c(1:nCells)){
    if(data_differ[i,j]==0)
      count <- count + 1
  }
}
rate <- count/(nGenes*nCells)
print('changerate:')
print(1-rate)



Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
library(Rcpp)
Rcpp::sourceCpp("julei.cpp")
source('analysis_methods.R')
library(devtools)
library(fpc)
library(mclust)
library(ngram)


########################################################################
#100 kmeans mean clustering of true data
realcluster <- data$group
ari<-vector(length=100)
sc<-vector(length=100)
nmi<-vector(length=100)
H<-t(data_true)
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
print('true cluster:')
print(ari_ave)
print(sc_ave)
print(nmi_ave)

########################################################################
#100 kmeans mean clustering of dropout data
ari1<-vector(length=100)
sc1<-vector(length=100)
nmi1<-vector(length=100)
H1<-t(data_dropout)
for(j in c(1:100)){
  cluster1 <- kmeans(H1,3)[["cluster"]]
  test1=cluster_evalu(H1,cluster1,realcluster)
  ari1[j]<-test1[[1]]
  sc1[j]<-test1[[2]]
  nmi1[j]<-test1[[3]]
}
cluster_list1<-list(ari1,sc1,nmi1)
names(cluster_list1)=c('ari','sc','nmi')
ari_ave1=mean(ari1)
sc_ave1=mean(sc1)
nmi_ave1=mean(nmi1)
print('dropout cluster:')
print(ari_ave1)
print(sc_ave1)
print(nmi_ave1)

