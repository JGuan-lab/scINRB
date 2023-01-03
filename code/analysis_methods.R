get_cor_data <- function(change_rate, scale_num){
  options( warn = -1 )  

  
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  index = rowMeans(data_simulation$data_dropout) > 0 
  
  data_true = data_simulation$data_true
  
  data_true = data_true[index,]
  
  # cell-cell correlation
  
  data_true_cell = cor(as.matrix((data_true)))
  
  data_true_cell[is.na(data_true_cell)] = 0
  
  # gene-gene correlation
  
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  data_true_gene[is.na(data_true_gene)] = 0
  
  data_dropout = data_simulation$data_dropout
  
  data_dropout = data_dropout[index,]
  
  # cell-cell correlation
  
  data_dropout_cell = cor((data_dropout), method = "pearson")
  
  data_dropout_cell[is.na(data_dropout_cell)] = 0
  
  # gene-gene correlation
  
  data_dropout_gene = cor(t((data_dropout)), method = "pearson")
  
  data_dropout_gene[is.na(data_dropout_gene)] = 0
  
  # load imputed data from FPimpute
  data_scINRB = read.table( file = paste0("imputation_scINRB_data/scINRB_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  data_scINRB = data_scINRB[index,]
  
  # cell-cell correlation
  
  data_scINRB_cell = cor((data_scINRB), method = "pearson")
  
  data_scINRB_cell[is.na(data_scINRB_cell)] = 0
  
  # gene-gene correlation
  
  data_scINRB_gene = cor(t((data_scINRB)), method = "pearson")
  
  data_scINRB_gene[is.na(data_scINRB_gene)] = 0
  

  
  # load imputed data from SAVER
  data_SAVER = read.table( file = paste0("imputation_SAVER_data/SAVER_",
                                         scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  data_SAVER = data_SAVER[index,]
  
  # cell-cell correlation
  
  data_SAVER_cell = cor((data_SAVER), method = "pearson")
  
  data_SAVER_cell[is.na(data_SAVER_cell)] = 0
  
  # gene-gene correlation
  
  data_SAVER_gene = cor(t((data_SAVER)), method = "pearson")
  
  data_SAVER_gene[is.na(data_SAVER_gene)] = 0
  
  # load imputed data from CMF

  data_CMF = t(read.table( file = paste0("imputation_CMF_data/CMF_impute_",
                                         scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  ))
  
  data_CMF = data_CMF[index,]
  # cell-cell correlation
  
  data_CMF_cell = cor((data_CMF), method = "pearson")
  
  data_CMF_cell[is.na(data_CMF_cell)] = 0
  
  # gene-gene correlation
  
  data_CMF_gene = cor(t((data_CMF)), method = "pearson")
  
  data_CMF_gene[is.na(data_CMF_gene)] = 0
  
  
  # load the magic results 
  #data = read.delim("~/Desktop/论文/kangyue/LPLS9.28/imputation_magic_data/magic_3_20%.rds", header=FALSE,row.names=1)
  data = read.delim(paste0("imputation_magic_data/magic_",scale_num,"_",change_rate,"%.rds"),
                    header = FALSE,row.names = 1,
                    sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  data_magic = data_magic[index,]
  
  # cell-cell correlation
  
  data_magic_cell = cor((data_magic), method = "pearson")
  
  data_magic_cell[is.na(data_magic_cell)] = 0
  
  # gene-gene correlation
  
  data_magic_gene = cor(t((data_magic)), method = "pearson")
  
  data_magic_gene[is.na(data_magic_gene)] = 0
  
  # load imputed data from scImpute
  
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  data_scimpute = data_scimpute[index,]
  
  # cell-cell correlation
  
  data_scimpute_cell = cor((data_scimpute), method = "pearson")
  
  data_scimpute_cell[is.na(data_scimpute_cell)] = 0
  
  # gene-gene correlation
  
  data_scimpute_gene = cor(t((data_scimpute)), method = "pearson")
  
  data_scimpute_gene[is.na(data_scimpute_gene)] = 0
  
  # load imputed data from Drimpute
  
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  data_drimpute = data_drimpute[index,]
  
  # cell-cell correlation
  
  data_drimpute_cell = cor((data_drimpute), method = "pearson")
  
  data_drimpute_cell[is.na(data_drimpute_cell)] = 0
  
  # gene-gene correlation
  
  data_drimpute_gene = cor(t((data_drimpute)), method = "pearson")
  
  data_drimpute_gene[is.na(data_drimpute_gene)] = 0
  
  
  # load imputed data from Saver
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  )
  
  data_viper = data_viper[index,]
  
  # cell-cell correlation
  
  data_viper_cell = cor((data_viper), method = "pearson")
  
  data_viper_cell[is.na(data_viper_cell)] = 0
  
  # gene-gene correlation
  
  data_viper_gene = cor(t((data_viper)), method = "pearson")
  
  data_viper_gene[is.na(data_viper_gene)] = 0
  
  # load imputed data from scrabble
  
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  data_scrabble = data_scrabble[index,]
  
  
  # cell-cell correlatio
  
  data_scrabble_cell = cor(as.matrix(data_scrabble), method = "pearson")
  
  data_scrabble_cell[is.na(data_scrabble_cell)] = 0
  
  # gene-gene correlation
  
  data_scrabble_gene = cor(t((data_scrabble)), method = "pearson")
  
  data_scrabble_gene[is.na(data_scrabble_gene)] = 0
  

  
  data_cell = list()
  
  data_cell[[1]] = data_true_cell
  
  data_cell[[2]] = data_dropout_cell
  
  data_cell[[3]] = data_FPimpute_cell
  
  data_cell[[4]] = data_drimpute_cell
  
  data_cell[[5]] = data_scimpute_cell
  
  data_cell[[6]] = data_scrabble_cell
  
  data_cell[[7]] = data_viper_cell
  
  data_cell[[8]] = data_SAVER_cell
  
  data_cell[[9]] = data_magic_cell
  
  data_cell[[10]] = data_CMF_cell
  
  data_gene = list()
  
  data_gene[[1]] = data_true_gene
  
  data_gene[[2]] = data_dropout_gene
  
  data_gene[[3]] = data_FPimpute_gene
  
  data_gene[[4]] = data_drimpute_gene
  
  data_gene[[5]] = data_scimpute_gene
  
  data_gene[[6]] = data_scrabble_gene
  
  data_gene[[7]] = data_viper_gene
  
  data_gene[[8]] = data_SAVER_gene
  
  data_gene[[9]] = data_magic_gene
  
  data_gene[[10]] = data_CMF_gene
  

  data_cor = list()
  
  data_cor[[1]] = data_cell
  
  data_cor[[2]] = data_gene
  
  return(data_cor)
  
}




get_cor_imputed_data <- function(change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  data_true = data_simulation$data_true
  
  data_dropout = data_simulation$data_dropout
  
  # load imputed data from scINRB
  data_scINRB = read.table( file = paste0("imputation_scINRB_data/scINRB_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  # load imputed data from SAVER
  data_SAVER = read.table( file = paste0("imputation_SAVER_data/SAVER_",
                                         scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  )
  
  
  
  # load imputed data from CMF
  
  data_CMF = t(read.table( file = paste0("imputation_CMF_data/CMF_impute_",
                                         scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  ))
  
  
 
  
  # load the magic results 
  #data = read.delim("~/Desktop/论文/kangyue/LPLS9.28/imputation_magic_changerate_data/magic_3_20%.rds", header=FALSE,row.names=1)
  data = read.delim(paste0("imputation_magic_data/magic_",scale_num,"_",change_rate,"%.rds"),
                    header = FALSE,row.names = 1,
                    sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  
  
  # load imputed data from scImpute
  
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  
  
  
  # load imputed data from Drimpute
  
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  
  
  # load imputed data from Saver
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  )
  
  
  # load imputed data from scrabble
  
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  )
  
  
  
  
  
  #相关矩阵转向量
  vector_true<- as.vector(unlist(data_true))
  vector_dropout <- as.vector(unlist(data_dropout))
  vector_scINRB <- as.vector(unlist(data_scINRB))
  vector_drimpute <- as.vector(unlist(data_drimpute))
  vector_scimpute <- as.vector(unlist(data_scimpute))
  vector_scrabble <- as.vector(unlist(data_scrabble))
  vector_viper <- as.vector(unlist(data_viper))
  vector_SAVER <- as.vector(unlist(data_SAVER))
  vector_magic <- as.vector(unlist(data_magic))
  vector_CMF <- as.vector(unlist(data_CMF))

  
  imputation_pss <- vector()
  imputation_pss[1]<-cor(vector_true,vector_dropout)
  imputation_pss[2]<-cor(vector_true,vector_scINRB)
  imputation_pss[3]<-cor(vector_true,vector_drimpute)
  imputation_pss[4]<-cor(vector_true,vector_scimpute)
  imputation_pss[5]<-cor(vector_true,vector_scrabble)
  imputation_pss[6]<-cor(vector_true,vector_viper)
  imputation_pss[7]<-cor(vector_true,vector_SAVER)
  imputation_pss[8]<-cor(vector_true,vector_magic)
  imputation_pss[9]<-cor(vector_true,vector_CMF)

  
  return(imputation_pss)
  
}




get_RMSE_imputed_data <- function(change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  data_true = as.matrix(data_simulation$data_true)
  
  data_dropout = as.matrix(data_simulation$data_dropout)
  
  # load imputed data from scINRB
  data_scINRB = as.matrix(read.table( file = paste0("imputation_scINRB_data/scINRB_",
                                                      scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  ))
  
  # load imputed data from SAVER
  data_SAVER = as.matrix(read.table( file = paste0("imputation_SAVER_data/SAVER_",
                                                   scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from CMF
  
  data_CMF = as.matrix(t(read.table( file = paste0("imputation_CMF_data/CMF_impute_",
                                                   scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  )))
  
  
  # load the magic results 
  data = read.delim(paste0("imputation_magic_data/magic_",scale_num,"_",change_rate,"%.rds"),
                    header = FALSE,row.names = 1,
                    sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  data_magic <- as.matrix(data_magic)
  
  # load imputed data from scImpute
  
  data_scimpute = as.matrix(read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                                      scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  ))
  
  
  
  
  # load imputed data from Drimpute
  
  data_drimpute = as.matrix(read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                                      scale_num, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from Saver
  data_viper = as.matrix(read.table( file = paste0("imputation_viper_data/viper_",
                                                   scale_num, "_",
                                         change_rate,
                                         "%.csv") ,
                           header = FALSE, sep=","
  ))
  
  
  # load imputed data from scrabble
  
  data_scrabble = as.matrix(read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            seed_value, "_",
                                            change_rate,
                                            "%.csv") ,
                              header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from netNMFsc_H
  
  data_netNMFsc_H = as.matrix(read.table( file = paste0("imputation_netNMFsc_data/H_",
                                                        scale_num, "_",
                                                        change_rate,
                                                        "%.csv") ,
                                          header = FALSE, sep=","
  ))
  
  data_netNMFsc_W = as.matrix(read.table( file = paste0("imputation_netNMFsc_data/W_",
                                                        scale_num, "_",
                                                        change_rate,
                                                        "%.csv") ,
                                          header = FALSE, sep=","
  ))
  
  data_netNMFsc <- data_netNMFsc_W%*%data_netNMFsc_H
  
  
  imputation_RMSE <- vector()
  imputation_RMSE[1]<-rmse(data_true,data_dropout)
  imputation_RMSE[2]<-rmse(data_true,data_scINRB)
  imputation_RMSE[3]<-rmse(data_true,data_drimpute)
  imputation_RMSE[4]<-rmse(data_true,data_scimpute)
  imputation_RMSE[5]<-rmse(data_true,data_scrabble)
  imputation_RMSE[6]<-rmse(data_true,data_viper)
  imputation_RMSE[7]<-rmse(data_true,data_SAVER)
  imputation_RMSE[8]<-rmse(data_true,data_magic)
  imputation_RMSE[9]<-rmse(data_true,data_CMF)
  imputation_RMSE[10]<-rmse(data_true,data_netNMFsc)
  
  return(imputation_RMSE)
  
}


plot_data <- function(data,name){
  limit <- c(0,5)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  print(dim(data))
  colnames(data) <- NULL
  rownames(data) <- NULL
  longData<-melt(as.matrix(data))
  colnames(longData) <- c("Var1", "Var2","value")
  
  pl <- ggplot(longData, aes_string(x = "Var2", y = "Var1")) +
    geom_raster(aes_string(fill= "value")) +
    scale_colour_gradient2(limits=c(0, 10)) +
    scale_fill_gradientn(colours = c("white", "blue", "red"), limits = c(0,10),values = c(0,0.6,1)) +
    theme_bw()  +
    scale_y_discrete(name ="Genes") +
    ggtitle(name) +
    scale_x_discrete(name ="Cells") +
    theme(panel.grid.major = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          line = element_blank(),
          plot.title = element_text(family = "Helvetica", face = "bold", size = (8)),
          axis.title = element_text(family = "Helvetica", size = (6)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    theme(legend.text=element_text(size=6),legend.title = element_text(size = 6))
  
  return(pl)
}

####cluster_evalu

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::sourceCpp("//Users//kangkangyue//Desktop//论文//kangyue//analysis//julei.cpp")
cluster_evalu<-function(H,cluster,realcluster){
  reallist<-list()
  cluslist<-list()
  for (i in 1:length(unique(realcluster))) {
    reallist[[i]]<-which(realcluster==unique(realcluster)[i])
    reallist[[i]]<-as.vector(reallist[[i]])
  }
  for (i in 1:length(unique(cluster))) {
    cluslist[[i]]<-which(cluster==unique(cluster)[i])
    cluslist[[i]]<-as.vector(cluslist[[i]])
  }
  nmi <- NMI(reallist,cluslist)
  ari <- adjustedRandIndex(realcluster,cluster)
  stats <- cluster.stats(dist(H)^2, cluster)
  sc <- stats$avg.silwidth
  clu_eva<-list(ari,sc,nmi)
  return(clu_eva)
}


get_simulation_data <- function(change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  data_true = as.matrix(data_simulation$data_true)
  
  data_dropout = as.matrix(data_simulation$data_dropout)
  
  data_bulk <- as.matrix(data_simulation$data_bulk)
  
  data_label <- as.matrix(data_simulation[["group"]])
  
  realcluster <- data_label
  
  # load imputed data from scINRB
  
  data_scINRB = as.matrix(read.table( file = paste0("imputation_scINRB_data/scINRB_",
                                                    scale_num, "_",
                                                    change_rate,
                                                    "%.csv") ,
                                      header = FALSE, sep=","
  ))
  
  # load imputed data from SAVER
  data_SAVER = as.matrix(read.table( file = paste0("imputation_SAVER_data/SAVER_",
                                                   scale_num, "_",
                                                   change_rate,
                                                   "%.csv") ,
                                     header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from CMF
  
  data_CMF = as.matrix(t(read.table( file = paste0("imputation_CMF_changerate_data/CMF_impute_",
                                                   scale_num, "_",
                                                   change_rate,
                                                   "%.csv") ,
                                     header = FALSE, sep=","
  )))
  
  
  # load the magic results 
  
  data = read.delim(paste0("imputation_magic_data/magic_",scale_num,"_",change_rate,"%.rds"),
                    header = FALSE,row.names = 1,
                    sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  data_magic <- as.matrix(data_magic)
  
  # load imputed data from scImpute
  
  data_scimpute = as.matrix(read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                                      scale_num, "_",
                                                      change_rate,
                                                      "%.csv") ,
                                        header = FALSE, sep=","
  ))
  
  
  
  
  # load imputed data from Drimpute
  
  data_drimpute = as.matrix(read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                                      scale_num, "_",
                                                      change_rate,
                                                      "%.csv") ,
                                        header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from Saver
  data_viper = as.matrix(read.table( file = paste0("imputation_viper_data/viper_",
                                                   scale_num, "_",
                                                   change_rate,
                                                   "%.csv") ,
                                     header = FALSE, sep=","
  ))
  
  
  # load imputed data from scrabble
  
  data_scrabble = as.matrix(read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                                      scale_num, "_",
                                                      change_rate,
                                                      "%.csv") ,
                                        header = FALSE, sep=","
  ))
  
  
  
  # load imputed data from netNMFsc_H
  
  # data_netNMFsc_H = as.matrix(read.table( file = paste0("imputation_netNMFsc_changerate_data/H_",
  #                                                       scale_num, "_",
  #                                                       change_rate,
  #                                                       "%.csv") ,
  #                                         header = FALSE, sep=","
  # ))
  # 
  # data_netNMFsc_W = as.matrix(read.table( file = paste0("imputation_netNMFsc_changerate_data/W_",
  #                                                       scale_num, "_",
  #                                                       change_rate,
  #                                                       "%.csv") ,
  #                                         header = FALSE, sep=","
  # ))
  # 
  # data_netNMFsc <- data_netNMFsc_W%*%data_netNMFsc_H
  # 
  imputation_data_list <- list()
  
  imputation_data_list[[1]] <- data_true
  
  imputation_data_list[[2]] <- data_dropout
  
  imputation_data_list[[3]] <- data_scINRB
  
  imputation_data_list[[4]] <- data_drimpute
  
  imputation_data_list[[5]] <- data_scimpute
  
  imputation_data_list[[6]] <- data_scrabble
  
  imputation_data_list[[7]] <- data_viper
  
  imputation_data_list[[8]] <- data_SAVER
  
  imputation_data_list[[9]] <- data_magic
  
  imputation_data_list[[10]] <- data_CMF
  
  
  result <- list(imputation_data_list,data_label)
  
  return(result)
  
}




tsen_plot <- function(data,realcluster){
  options( warn = -1 ) 
  name_list <- list()
  name_list[[1]] <- "true data"
  name_list[[2]] <- "dropout data"
  name_list[[3]] <- "scINRB"
  name_list[[4]] <- "drimpute"
  name_list[[5]] <- "scimpute"
  name_list[[6]] <- "scrabble"
  name_list[[7]] <- "viper"
  name_list[[8]] <- "SAVER"
  name_list[[9]] <- "magic"
  name_list[[10]] <- "CMF"
  p <- list()
  asn<-list()
  tsne.coords <- list()
  for(k in 1:10){
    impute_data <- data[[k]]
    tsne.coords[[k]] <- Rtsne(t(impute_data),pca=FALSE,perplexity = 30,theta=0.5,check_duplicates = F)$Y
    #rownames(tsne.coords[[k]]) <- rownames(t(data_true))
    colnames(tsne.coords[[k]]) <- c("tSNE1","tSNE2")
    tsne.coords[[k]] <- as.data.frame(tsne.coords[[k]])
    p[[k]] <- ggplot(tsne.coords[[k]],aes(tSNE1,tSNE2,color=realcluster)) +
      geom_point(size=0.2) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major=element_line(colour=NA),
            panel.grid.minor = element_blank(),
            axis.title =  element_text(size=10,face = "bold"),
      ) +
      labs(title = name_list[[k]],color="Species")
    
    ari=0
    sc=0
    nmi=0
    for(j in c(1:100)){
      cluster <- kmeans(tsne.coords[[k]],3)[["cluster"]]
      test=cluster_evalu(tsne.coords[[k]],cluster,realcluster)
      ari=ari+test[[1]]
      sc=sc+test[[2]]
      nmi=nmi+test[[3]]
    }
    ari_ave=ari/100
    sc_ave=sc/100
    nmi_ave=nmi/100
    asn[[k]] <- c(ari_ave,sc_ave,nmi_ave)
    
    
    
    
  }
  
  result <- list(p,asn)
  
  return(result)
  
}



get_simulation_data1 <- function(methond,change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  data_true = as.matrix(data_simulation$data_true)
  
  data_dropout = as.matrix(data_simulation$data_dropout)
  
  data_bulk <- as.matrix(data_simulation$data_bulk)
  
  data_label <- as.matrix(data_simulation[["group"]])
  
  realcluster <- data_label
  
  # load imputed data from FPimpute
  
  impute_data = as.matrix(read.table( file = paste0("imputation_",method,"_data/",method,"_",
                                                    scale_num, "_",
                                                    change_rate,
                                                    "%.csv") ,
                                      header = FALSE, sep=","
  ))
  H<-t(impute_data)
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
  result <- list(ari_ave,sc_ave,nmi_ave)
  write.csv(cluster_list,paste0("/Users/kangkangyue/Desktop/scINRB/imputation_",method,"_data/cluster_",
                                scale_num,"_",
                                change_rate,
                                "%.csv"))
  return(result)
  
  
  
}


########simulation data tsne clustering 
tsne_simulation_data <- function(methond,change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  data_true = as.matrix(data_simulation$data_true)
  
  data_dropout = as.matrix(data_simulation$data_dropout)
  
  data_bulk <- as.matrix(data_simulation$data_bulk)
  
  data_label <- as.matrix(data_simulation[["group"]])
  
  realcluster <- data_label
  
  # load imputed data from scINRB
  
  impute_data = as.matrix(read.table( file = paste0("imputation_",method,"_data/",method,"_",
                                                    scale_num, "_",
                                                    change_rate,
                                                    "%.csv") ,
                                      header = FALSE, sep=","
  ))
  
  tsne.coords <- Rtsne(t(impute_data),pca=FALSE,perplexity = 30,theta=0.5,check_duplicates = F)$Y
  colnames(tsne.coords) <- c("tSNE1","tSNE2")
  tsne.coords <- as.data.frame(tsne.coords)
  
  ari<-vector(length=100)
  sc<-vector(length=100)
  nmi<-vector(length=100)
  for(j in c(1:100)){
    cluster <- kmeans(tsne.coords,3)[["cluster"]]
    test=cluster_evalu(tsne.coords,cluster,realcluster)
    ari[j]<-test[[1]]
    sc[j]<-test[[2]]
    nmi[j]<-test[[3]]
  }
  cluster_list<-list(ari,sc,nmi)
  names(cluster_list)=c('ari','sc','nmi')
  ari_ave=mean(ari)
  sc_ave=mean(sc)
  nmi_ave=mean(nmi)
  
  
  result <- list(ari_ave,sc_ave,nmi_ave)
  write.csv(cluster_list,paste0("/Users/kangkangyue/Desktop/scINRB/imputation_",method,"_data/cluster_tsne_",
                                scale_num,"_",
                                change_rate,
                                "%.csv"))
  return(result)
  
}

########simulation data matrix pss/RMSE,gene pss/cell pss.

get_cor_result <- function(methond,change_rate, scale_num){
  options( warn = -1 ) 
  
  # load the simulationd data
  data_simulation <- readRDS(file = paste0('data/simulation_data/',
                                           scale_num,
                                           '_',
                                           change_rate,
                                           '%.rds')
  )
  
  
  index = rowMeans(data_simulation$data_dropout) > 0 
  
  data_true = data_simulation$data_true
  
  data_true = data_true[index,]
  
  # cell-cell correlation
  
  data_true_cell = cor(as.matrix((data_true)))
  
  data_true_cell[is.na(data_true_cell)] = 0
  
  # gene-gene correlation
  
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  data_true_gene[is.na(data_true_gene)] = 0
  
  
  
  # load imputed data from imputemethod
  
  impute_data = as.matrix(read.table( file = paste0("imputation_",method,"_data/",method,"_",
                                                    scale_num, "_",
                                                    change_rate,
                                                    "%.csv") ,
                                      header = FALSE, sep=","
  ))
  
  impute_data = impute_data[index,]
  
  # cell-cell correlation
  
  impute_data_cell = cor((impute_data), method = "pearson")
  
  impute_data_cell[is.na(impute_data_cell)] = 0
  
  # gene-gene correlation
  
  impute_data_gene = cor(t((impute_data)), method = "pearson")
  
  impute_data_gene[is.na(impute_data_gene)] = 0
  
  
  vector_true<- as.vector(unlist(data_true))
  vector_impute_data <- as.vector(unlist(impute_data))
  
  c_true<- as.vector(unlist(data_true_cell))
  c_impute_data <- as.vector(unlist(impute_data_cell))

  g_true<- as.vector(unlist(data_true_gene))
  g_impute_data <- as.vector(unlist(impute_data_gene))

  cell_pss<-cor(c_true,c_impute_data)

  gene_pss<-cor(g_true,g_impute_data)
  
  matrix_pss <- cor(vector_true,vector_impute_data)
  
  RMSE <- rmse(data_true,impute_data)
  
  
  
  result <- vector()
  
  result[[1]] = RMSE
  
  result[[2]] = matrix_pss
     
  result[[3]] = cell_pss
    
  result[[4]] = gene_pss
  
  
  return(result)
  
}

