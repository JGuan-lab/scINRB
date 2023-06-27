#' Runs scINRB
#'
#' @param data the input dropout scRNAseq data.
#'
#' @param data_bulk the bulk RNAseq.
#'
#' @param parameter the vector of parameters. 
#' The first parameter is the value of alpha1 in the mathematical model;
#' the second one is the value of alpha2 in the mathematical model;
#' the third one is the value of alpha3 in the mathematical model.
#'
#' @param nIter the maximum iterations, the default is 2000.
#'
#' @param error_threshold the threshold of the error between the current imputed matrix and the previous one.。
#' Default is 1e-5.
#'
#' @examples
#' # Set up the parameter used in scINRB
#' parameter <- c(0.001,0.001,1)
#' r <- 200
#' 
#' # Run scINRB
#' result <- scINRB(data_dropout,data_bulk,parameter,r)
#'
#' @return A data matrix with the same size of the input scRNAseq data
#'
#' @rdname scINRB
#'
#' @export
#' 
bc = 1e-7
nIter = 2000
error_threshold = 1e-05

scINRB <- function(data,data_bulk,parameter,r,bc = 1e-7,nIter = 2000,error_threshold = 1e-05)
{
  time1<-Sys.time()
  lambda1 <- parameter[1]
  lambda2 <- parameter[2]
  lambda3 <- parameter[3]
  m = nrow(data)#gene 
  n = ncol(data)#cell
  K <- matrix(data=NA, nrow = m, ncol = n)
  row_sum<-rowSums(data)

  #compute sg
  ag=SM(data)
  Ag = ag - diag(m)
  Dg = diag(abs(apply(Ag,2,sum)))
  #det(sqrt(Dg))
  sg = diag(m) - solve(sqrt(Dg)) %*% Ag %*% solve(sqrt(Dg))
  print(paste0('sg is finished'))
  
  #compute sc
  ac=SM(t(data))#similarity
  Ac = ac - diag(n)#adjacency matrix
  Dc = diag(abs(apply(Ac,2,sum)))
  #det(sqrt(Dc))
  sc = diag(n) - solve(sqrt(Dc)) %*% Ac %*% solve(sqrt(Dc))
  print(paste0('sc is finished'))
  
  #initialize H and W
  set.seed(1)
  x1 <- runif(m*r,min=0,max=1)
  W <- matrix(x1,nrow=m,ncol=r)
  set.seed(2)
  x2 <- runif(n*r,min=0,max=1)
  H <- matrix(x2,nrow=n,ncol=r)
  
  #initialize P and Q
  set.seed(3)
  x1 <- runif(m*r,min=0,max=1)
  P <- matrix(x1,nrow=m,ncol=r)
  set.seed(4)
  x2 <- runif(n*r,min=0,max=1)
  Q <- matrix(x2,nrow=n,ncol=r)

  a<-rep(1/n,n)
  a<-matrix(a)
  a<-t(a)

  data_bulk<-matrix(data_bulk)
  data_bulk<-t(data_bulk)
  k <- 0
  error <- 1
  l0 <- 0
  X0 <- matrix(1,m,n)
  print(paste0('begin'))
  #circulate
  while((k < nIter) & (error > error_threshold))
  {
    #print(k)
    data=as.matrix(data)
    t1=2*W%*%t(H)%*%t(a)%*%a%*%H-2*t(data_bulk)%*%a%*%H
    l1=(-2)*data%*%H+2*W%*%t(H)%*%H-2*lambda1*sg%*%W+lambda3*t1-P
    t2=2*t(a)%*%a%*%H%*%t(W)%*%W-2*t(a)%*%data_bulk%*%W
    l2=(-2)*t(data)%*%W+2*H%*%t(W)%*%W+2*lambda2*sc%*%H+lambda3*t2-Q
    
    # update the P、Q
    P <- P+bc*(-W)
    P[P < 0] <- 0
    Q <- Q+bc*(-H)
    Q[Q < 0] <- 0
    
    # update the W、H
    W <- W-bc*l1
    W[W < 0] <- 0
    H <- H-bc*l2
    H[H < 0] <- 0
    #bc<-0.98*bc
    #print(paste0('W H is renewed'))
    list <- list(W,H,W%*%t(H))

    #compute object function
    x1=data-(W%*%t(H))
    b1=norm(x1,type = "F")^2
    x2=sum(diag(t(W)%*%sg%*%W))
    b2=lambda1*x2
    x3=sum(diag(t(H)%*%sc%*%H))
    b3=lambda2*x3
    x4=(a%*%H%*%t(W))-data_bulk
    b4=lambda3*(norm(x4,type = "2")^2)
    l= b1+b2+b3+b4
    #print(paste0('L='))
    #print(l)
    
    #bounce condition——error<err
    if(l0<l)
    {
      flag <- 1   #increase
    }
    else
    {
      flag <- -1   #decrease
    }
    #print(paste0('error is'))
    c1=norm(X0-list[[3]],type = "F")
    c2=norm(X0,type = "F")
    error <- c1/c2
    #print(error)
    #print(flag)
    k<-k+1
    l0=b1+b2+b3+b4
    X0 <- list[[3]]
  }
  print(paste0('Imputation is finished'))

  #replace only non-zero values
  data_Y <- list[[3]]
  mean<-sort(data_Y)[m*n*0.1]
  data_Y[data_Y < mean] <- 0
  data_Y[data_Y < 0] <- 0
  
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(data[i,j]==0)
        data[i,j]<-data_Y[i,j]
    }
  }
  
  for(i in 1:m)
  {
    if(row_sum[i]==0)
      for(j in 1:n)
      {
        K[i,j]<-0
      }
    else
      for(j in 1:n)
      {
        K[i,j]<-1
      }
  }
    
  data_impute <- data*K
    
  list[[3]] <- data_impute
    
  time <- f(time1)
    
  print(time)
    
  return(list)
}
