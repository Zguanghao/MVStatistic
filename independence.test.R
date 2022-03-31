independence.test <- function(data,split,k){
  n <- nrow(data)
  p <- ncol(data)
  if(missing(split)){split <- rep(1,p)}
  if(missing(k)){k <- p}
  stopifnot(sum(split)==p,length(split)==k)
  
  A <- (n-1)*cov(data)
  index <- cumsum(split)
  index <- c(0,index)
  Ak <- list()
  for(i in 1:k){
    Ak[[i]] <- as.matrix(A[(index[i]+1):index[i+1],(index[i]+1):index[i+1]])
  }
  det_Ak <- sapply(Ak,det)
  Lambda <- (det(A)/cumprod(det_Ak)[k])^((n-1)/2)
  f <- (p^2-sum(split^2))/2
  alpha <- (2*(p^3-sum(split^3))+9*(p^2-sum(split^2)))/(6*(p^2-sum(split^2)))
  rho <- 1-alpha/n
  
  chi_stat <- -2*log(Lambda)
  pvalue <- 1-pchisq(rho*chi_stat,f)
  return(list(chi_stat = chi_stat,pvalue = pvalue))
}
