#####Record Runtime#####
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

#####Calculate sg\sc for scINRB#####
SM <- function(data)
{
  m <- nrow(data)
  n <- ncol(data)
  D <- dist(data,method="euclidean")
  maxD=max(D)
  score <- as.matrix(D)
  score <- score/maxD
  score <- 1/(score+1)
  
}

#####Generate the zero-one matrix required for cross-verification parameters#####
mm <- function(m,n)
{
  x <- gl(5, 0.2*m*n)
  x <- sample(x,length(x))
  M <- matrix(x,nrow=m,ncol=n)
  M1 <- matrix(data=NA, nrow = m, ncol = n)
  M2 <- matrix(data=NA, nrow = m, ncol = n)
  M3 <- matrix(data=NA, nrow = m, ncol = n)
  M4 <- matrix(data=NA, nrow = m, ncol = n)
  M5 <- matrix(data=NA, nrow = m, ncol = n)
  m1 <- matrix(data=NA, nrow = m, ncol = n)
  m2 <- matrix(data=NA, nrow = m, ncol = n)
  m3 <- matrix(data=NA, nrow = m, ncol = n)
  m4 <- matrix(data=NA, nrow = m, ncol = n)
  m5 <- matrix(data=NA, nrow = m, ncol = n)
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(M[i,j]==1)
      {
        M1[i,j]<-0
        m1[i,j]<-1
      }
      else
      {
        M1[i,j]<-1
        m1[i,j]<-0
      }
    }
  }

  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(M[i,j]==2)
      {
        M2[i,j]<-0
        m2[i,j]<-1
      }
      else
      {
        M2[i,j]<-1
        m2[i,j]<-0
      }
    }
  }
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(M[i,j]==3)
      {
        M3[i,j]<-0
        m3[i,j]<-1
      }
      else
      {
        M3[i,j]<-1
        m3[i,j]<-0
      }
    }
  }
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(M[i,j]==4)
      {
        M4[i,j]<-0
        m4[i,j]<-1
      }
      else
      {
        M4[i,j]<-1
        m4[i,j]<-0
      }
    }
  }
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      if(M[i,j]==5)
      {
        M5[i,j]<-0
        m5[i,j]<-1
      }
      else
      {
        M5[i,j]<-1
        m5[i,j]<-0
      }
    }
  }
  list1 <- list(M1,M2,M3,M4,M5)
  list2 <- list(m1,m2,m3,m4,m5)
  return(list(list1,list2))
  
}


