two.mu_unequal <- function(data1,data2){
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)
  if(n1==n2){
    Z <- data1-data2
    mu_z <- apply(Z, 2, mean) 
    sigma1 <- cov(data1)
    sigma2 <- cov(data2)
    #sigma_z <- sigma1+sigma2
    sigma_z <- var(Z)
    TT <- n1*t(mu_z)%*%solve(sigma_z)%*%mu_z
  }
  else{
    if(n1>n2){
      X <- data2
      Y <- data1
      t <- n1
      n1 <- n2
      n2 <- t
    }else{
      X <- data1
      Y <- data2
    }
    sigma1 <- cov(X)
    sigma2 <- cov(Y)
    sum_Y_n <- apply(Y[1:n1,],2,sum)
    sum_Y_m <- apply(Y,2,sum)
    Z <- X-sqrt(n1/n2)*Y[1:n1,]+sum_Y_n/sqrt(n1*n2)-sum_Y_m/n2
    mu_z <- apply(Z,2,mean)
    #sigma_z <- sigma1+n1*sigma2/n2
    sigma_z <- var(Z)
    TT <- n1*t(mu_z)%*%solve(sigma_z)%*%mu_z
  }
  
  Hotelling_dis <- Hotelling$new(n=n1-1,p=p)
  pvalue <- 1-Hotelling_dis$cdf(TT)
  return(list(Hotelling_stat = TT,pvalue = pvalue))
}