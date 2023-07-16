#####input data#####
path="/../simulation_data_changerate" 
setwd(path)  
data <- readRDS("3_20%.rds")
data_dropout <- as.matrix(data$data_dropout)
data_true <- as.matrix(data$data_true)
data_bulk <- as.matrix(data$data_bulk)
data_label <- as.matrix(data[["group"]])

#####Load the required functions#####
path=".." 
setwd(path) 
source('scINRB.R')
source('fuctions.R')
library(MASS)

#####Generate the zero-one matrix required for cross-verification parameter r#####
#gene 
m = nrow(data_dropout)
#cell
n = ncol(data_dropout)
result <- mm(m,n)

M0<-result[[1]]
m0<-result[[2]]
RMSE <- c(0,0,0,0,0)
min_RMSE_AVE <- Inf  # Initialize with a large value
best_r <- 0  # Initialize with 0

#####Calculate the RMSE mean of cross-validation for each r#####
for(r in c(10,50,100,300,500)){
  
  #####Calculate RMSE under this parameter#####
  for(value in c(1:5)){
    
    M <- M0[[value]]
    m <- m0[[value]]
    parameter <- c(0,0,0)
    newdata <- impute(M,data_dropout,data_bulk,parameter,r)
    x1 <- m*(data_dropout-newdata)
    x2 <- m*(data_true-newdata)
    n <- nrow(data_dropout)*ncol(data_dropout)*0.2
    MSE[value] <- norm(x1,type = "F")/n
    RMSE[value] <- sqrt(MSE[value])
  }
  RMSE_AVE=(RMSE[1]+RMSE[2]+RMSE[3]+RMSE[4]+RMSE[5])/5
  # Update minimum RMSE_AVE and corresponding 'r' if a smaller value is found
  if (RMSE_AVE < min_RMSE_AVE) {
    min_RMSE_AVE <- RMSE_AVE
    best_r <- r
  }
}
# Output the parameter 'r' corresponding to the minimum RMSE_AVE
cat("Best parameter r:", best_r)

###Continue to use cross-validation to select three regularization parameters.###
#####Generate the zero-one matrix required for cross-verification parameters#####
#gene 
m = nrow(data_dropout)
#cell
n = ncol(data_dropout)
result <- mm(m,n)

M0<-result[[1]]
m0<-result[[2]]
RMSE <- c(0,0,0,0,0)
min_RMSE_AVE <- Inf  # Initialize with a large value
best_r <- 0  # Initialize with 0

#####Calculate the RMSE mean of cross-validation for each r#####
for(p1 in c(0.001,0.1,10)){
  p2 <- p1
  
  for(p3 in c(1,10,50,100)){
    #####Calculate RMSE under this parameter#####
    for(value in c(1:5)){
      
      M <- M0[[value]]
      m <- m0[[value]]
      parameter <- c(p1,p2,p3)
      newdata <- impute(M,data_dropout,data_bulk,parameter,best_r)
      x1 <- m*(data_dropout-newdata)
      x2 <- m*(data_true-newdata)
      n <- nrow(data_dropout)*ncol(data_dropout)*0.2
      MSE[value] <- norm(x1,type = "F")/n
      RMSE[value] <- sqrt(MSE[value])
    }
    RMSE_AVE=(RMSE[1]+RMSE[2]+RMSE[3]+RMSE[4]+RMSE[5])/5
    # Update minimum RMSE_AVE and corresponding p1/p2/p3 if a smaller value is found
    if (RMSE_AVE < min_RMSE_AVE) {
      min_RMSE_AVE <- RMSE_AVE
      best_p1 <- p1
      best_p2 <- p2
      best_p3 <- p3
    }

  }
}

# Output the parameter p1\p2\p3 corresponding to the minimum RMSE_AVE
cat("Best parameter p1:", best_p1)
cat("Best parameter p2:", best_p2)
cat("Best parameter p3:", best_p3)
