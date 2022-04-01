#该程序创建列联表类并可进行独立性检验,相合性度量与检验,
#一致性度量与检验.
#author: zheng guanghao

library(R6)
Crosstable <- R6Class("Crosstable",
                public = list(
                  data = NA,
                  nrows = NA,
                  ncols = NA, 
                  row_prop = NA,
                  col_prop = NA,
                  row_order = NA,
                  col_order = NA,
                  crosstable = NA,
                  initialize = function(data,row_prop,col_prop,row_order,col_order){
                    self$nrows <- nrow(data)
                    self$ncols <- ncol(data)
                    if(!missing(data)) self$data <- data
                    if(!missing(row_prop)) 
                      {self$row_prop <- row_prop}
                    else{
                      if(is.null(row.names(data))){
                        self$row_prop <- paste("A",c(1:self$nrows),sep = "")
                      }else{self$row_prop <- row.names(data)}

                    }
                    if(!missing(col_prop)) 
                    {self$col_prop <- col_prop}
                    else{                      
                      if(is.null(colnames(data))){
                      self$col_prop <- paste("B",c(1:self$ncols),sep = "")
                    }else{self$col_prop <- colnames(data)}}
                    
                    
                    if(!missing(row_order))
                    {self$row_order <- row_order}
                    else{self$row_order <- 0}
                    if(!missing(col_order))
                    {self$col_order <- col_order}
                    else{self$col_order <- 0}
                    
                    table <- matrix(rep(0,(self$nrows+1)*(self$ncols+1)),nrow = self$nrows+1)
                    data <- as.matrix(self$data)
                    row_name <- c(self$row_prop,"sum")
                    col_name <- c(self$col_prop,"sum")
                    table[1:self$nrows,1:self$ncols] <- data
                    table[self$nrows+1,1:self$ncols] <- apply(data,2,sum)
                    table[1:self$nrows,self$ncols+1] <- apply(data,1,sum)
                    table[self$nrows+1,self$ncols+1] <- sum(data)
                    crosstable <- as.data.frame(table)
                    row.names(crosstable) <- row_name
                    colnames(crosstable) <- col_name
                    self$crosstable <- crosstable
                  },
                  
                  summary = function(){
                    print(self$crosstable)
                  },
                  
                  chisq.test = function(method = "Pearson",alpha = 0.05){
                    data <- self$data
                    r_sum <- apply(data,1,sum)
                    c_sum <- apply(data,2,sum)
                    all_sum <- sum(data)
                    np_ij <- outer(r_sum,c_sum)/all_sum
                    if(method=="Pearson"){
                    chi_stat <- sum((data-np_ij)^2/np_ij)}
                    else{
                      chi_stat <- -2*sum(data*log(np_ij/data))
                    }
                    df <- (self$nrows-1)*(self$ncols-1)
                    pvalue <- 1-pchisq(chi_stat,df)
                    if(pvalue<alpha){
                      print(paste("alpha =",0.05,", reject the H0."))
                    }else{
                      print(paste("alpha =",0.05,", can't reject the H0."))
                    }
                    return(list(chi_stat = chi_stat,pvalue = pvalue))
                  },
                  
                  corrlation.test = function(alternative = "two.side"){
                    locationG <- list()
                    k <- 1
                    for(i in 1:(self$nrows-1)){
                      for(j in 1:(self$ncols-1)){
                        locationG[[k]] <- c(i,j)
                        k <- k+1
                      }
                    }
                  G_ij <- sapply(locationG,private$calculateG_ij,self$data)
                  G <- sum(G_ij)
                  locationH <- list()
                  k <- 1
                  for(i in 1:(self$nrows-1)){
                    for(j in 2:self$ncols){
                      locationH[[k]] <- c(i,j)
                      k <- k+1
                    }
                  }
                  H_ij <- sapply(locationH,private$calculateH_ij,self$data)
                  H <- sum(H_ij)
                  z <- G-H
                  r_sum <- apply(self$data,1,sum)
                  c_sum <- apply(self$data,2,sum)
                  n <- sum(data)
                  TA <- sum(r_sum*(r_sum-1))/2
                  TB <- sum(c_sum*(c_sum-1))/2
                  TAB <- sum(data*(data-1))/2
                  
                  tau <- z/sqrt((n*(n-1)/2-TA)*(n*(n-1)/2-TB))
                  gamma <- z/(G+H)
                  dBA <- z/(n*(n-1)/2-TA)
                  dAB <- z/(n*(n-1)/2-TB)
                  
                  sigma_z <- (n^3-sum(r_sum^3))*(n^3-sum(c_sum^3))/(9*n^3)
                  U_stat <- z/sqrt(sigma_z)
                  chi_stat <- z^2/sigma_z
                  if(alternative=="two.side"){
                    pvalue <- 1-pchisq(chi_stat,1)
                    stat <- chi_stat
                    test_type <- "Chi_square test"
                  }else if(alternative=="right"){
                    pvalue <- 1-pnorm(U_stat,0,1)
                    stat <- U_stat
                    test_type <- "U test"
                  }else{
                    pvalue <- pnorm(-U_stat,0,1)
                    stat <- U_stat
                    test_type <- "U test"
                  }
                  
                  return(list(Kendall_tau = tau,Gamma = gamma,Somers_d = c(dBA,dAB),
                              test_type = test_type,stat = stat,pvalue = pvalue))
                  },
                  
                  consistency.test = function(){
                    stopifnot(self$nrows==self$ncols)
                    data <- self$data
                    n <- sum(data)
                    q1 <- sum(diag(data))/n
                    r_sum <- apply(data,1,sum)
                    c_sum <- apply(data,2,sum)
                    q2 <- sum((r_sum*c_sum)/n^2)
                    k <- (q1-q2)/(1-q2)
                    sigma_k <- (q2+q2^2-sum(r_sum*c_sum*(r_sum+c_sum)/n^3))/(n-1)/(1-q2)^2
                    U_stat <- k/sqrt(sigma_k)
                    pvalue <- 1-pnorm(U_stat,0,1)
                    return(list(Kappa = k,stat = U_stat,pvalue = pvalue))
                  } 
                  
),
                private = list(
                  calculateG_ij = function(ij,data){
                    return(data[ij[1],ij[2]]*sum(data[(ij[1]+1):self$nrows,(ij[2]+1):self$ncols]))
                  },
                  calculateH_ij = function(ij,data){
                    return(data[ij[1],ij[2]]*sum(data[(ij[1]+1):self$nrows,1:(ij[2]-1)]))
                  }
                )
                      )







