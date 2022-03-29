# This file contains classes of Hotelling's T2 distribution and
# Wilks Lambda distribution, and methods to generate random sa-
# mple(.rvs),calculate the probability(.cdf), calculate the 
# quantile(.quantile).The code below constructs the Hotelling's
# T2 distribution and Wilks Lambda distribution by using F dis-
# tribution,chisquare distribution and the relationship between 
# T2 and Lambda.
#
# author: Zheng guanghao
#
# -------------------------------------------------------------
#                    (n-p+1/np)*T2 ~ F(p,n-p+1)
# -------------------------------------------------------------
#   p  |  l  |      relationship between F and Lambda
# -------------------------------------------------------------
#      |  1  |       (n-p+1)/p * (1-L)/L ~ F(p,n-p+1)        
#      |  2  | (n-p+1)/p * (1-sqrt(L))/sqrt(L) ~ F(2p,2n-2p+2)
#   1  |     |           n/l * (1-L)/L ~ F(l,n)
#   2  |     |  (n-1)/l * (1-sqrt(L))/sqrt(L) ~ F(2l,2(n-1))
#  >2  | >2  |   -rlogL ~ Chisquare(pl), r = n-(p-l+1)/2
# -------------------------------------------------------------
#
# Examples:
#
# >T2 <- Hotelling$new(10,2)
# > T2$cdf(1)
# [1] 0.3487722
# > T2$rvs(2,2)
# [,1]     [,2]
# [1,] 1.298609 2.080543
# [2,] 8.302946 1.000847
# > T2$quantile(0.5)
# [1] 1.66529
# > 
# > L <- Wilks$new(2,10,1)
# > L$cdf(1)
# [1] 1
# > L$cdf(0.8)
# [1] 0.3663574
# > L$rvs(2,2)
# [,1]      [,2]
# [1,] 0.6555047 0.8943807
# [2,] 0.9486193 0.7565949
# > L$quantile(0.5)
# [1] 0.857244
# > 
library(R6)
Hotelling <- R6Class("Hotelling",
                     public = list(
                       n = NA,
                       p = NA,
                       initialize = function(n,p){
                         if(!missing(n)) self$n <- n
                         if(!missing(p)) self$p <- p
                       },
                       cdf = function(x){
                         n <- self$n
                         p <- self$p
                         x_adjust <- ((n-p+1)/(n*p))*x
                         prob <- sapply(x_adjust,pf,df1=p,df2=n-p+1)
                         return(prob)
                       },
                       rvs = function(rows,cols=1){
                         n <- self$n
                         p <- self$p
                         N <- rows*cols
                         r <- rf(N,df1=p,df2=n-p+1)
                         r_adjust <- r*n*p/(n-p+1)
                         return(matrix(r_adjust,nrow = rows,ncol = cols))
                       },
                       quantile = function(prob){
                         n <- self$n
                         p <- self$p
                         q <- sapply(prob,qf,df1=p,df2=n-p+1)
                         q_adjust <- q*n*p/(n-p+1)
                         return(q_adjust)
                       }
                     ),
                  private = list()
)


Wilks <- R6Class("Wilks",
                 public = list(
                   p = NA,
                   n = NA,
                   l = NA,
                   initialize = function(p,n,l){
                     if(!missing(n)) self$n <- n
                     if(!missing(p)) self$p <- p
                     if(!missing(l)) self$l <- l
                     private$set()
                     
                   },
                   para = function(){
                     private$asym_distribution[[1]](1,1,1)
                   },
                   cdf = function(x){
                     x_adjust <- private$adjust_para(x)
                     prob <- sapply(x_adjust,private$asym_distribution[[1]],
                                    private$asym_distribution_para[1],
                                    private$asym_distribution_para[2])
                     return(1-prob)
                   },
                   rvs = function(rows,cols=1){
                     N <- rows*cols
                     r_fun <- private$asym_distribution[[2]]
                     r <- r_fun(N,private$asym_distribution_para[1],
                                  private$asym_distribution_para[2])
                     r_adjust <- private$inv_adjust_para(r)
                     return(matrix(r_adjust,nrow = rows,ncol = cols))
                   },
                   quantile = function(prob){
                     q <- sapply(1-prob,private$asym_distribution[[3]],
                                 private$asym_distribution_para[1],
                                 private$asym_distribution_para[2])
                     q_adjust <- private$inv_adjust_para(q)
                     return(q_adjust)
                   },
                   some_fun = function(){
                     fun_list <- list(private$adjust_para,private$inv_adjust_para)
                     return(fun_list)
                   }
                   
                 ),
                 private = list(
                   adjust_para = NA,
                   inv_adjust_para = NA,
                   asym_distribution = NA,
                   asym_distribution_para = NA,
                   
                  set = function(){
                    if(self$l<3){
                      private$asym_distribution = c(pf,rf,qf)
                      private$adjust_para <- function(lambda){
                        #f
                        A <- (self$n-self$p+1)/self$p
                        A*(1-lambda^(1/self$l))/lambda^(1/self$l)
                      }
                      private$inv_adjust_para <- function(F_stat){
                        #g
                        A <- (self$n-self$p+1)/self$p
                        (F_stat/A+1)^(-self$l)
                      }
                      private$asym_distribution_para <- c(self$l*self$p,self$l*(self$n-self$p-self$l+2))
                    }else if(self$p<3&&self$l>2){
                      private$asym_distribution = c(pf,rf,qf)
                      private$adjust_para <- function(lambda){
                        A <- (self$n-self$p+1)/self$l
                        A*(1-lambda^(1/self$p))/lambda^(1/self$p)
                      }
                      private$inv_adjust_para <- function(F_stat){
                        A <- (self$n-self$p+1)/self$l
                        (F_stat/A+1)^(-self$p)
                      }
                      private$asym_distribution_para <- c(self$p*self$l,self$p*(self$n-self$p+1))
                    }
                    else{
                      private$asym_distribution = c(pchisq,rchisq,qchisq)
                      private$adjust_para <- function(lambda){
                        r <- self$n-(self$p-self$l+1)/2
                        -r*log(lambda)
                      }
                      private$inv_adjust_para <- function(chi_stat){
                        r <- self$n-(self$p-self$l+1)/2
                        exp(-chi_stat/r)
                      }
                      private$asym_distribution_para <- c(self$p*self$l,0)
                    }
                  }
                 )
                 )


















