#' @title A illustration dataset
#' @name mydata
#' @description A dataset used to illustrate the performance of the KCP permutation test.
#' @examples
#' \dontrun{
#' data(mydata)
#' attach(mydata)
#' N <- length(mydata[,1])
#' v <- length(mydata[1,])
#' w <- 15
#' p_VarKCP <- VarKCP(mydata,N,v,w)
#' p_VarKCP
#' p_Vardrop <- Vardrop(mydata,N,v,w,2)
#' p_Vardrop
#' }
#' @useDynLib StatComp22058
NULL

#' @title Running correlations function using R
#' @description Compute the running correlations using R
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param N the number of time points in X
#' @param v the number of variables in X
#' @param w the width of windows. Compute the running correlations in a window.
#' @return a matrix of the running correlations for all variables at these time points.
#' @examples
#' \dontrun{
#' data(mydata)
#' attach(mydata)
#' N <- length(mydata[,1])
#' v <- length(mydata[1,])
#' w <- 15
#' cordata <- CorrelationR(mydata,N,v,w)
#' cordata
#' }
#' @importFrom stats cor
#' @export
CorrelationR <- function(X,N,v,w){
  n <- N-w
  Xcor <- matrix(0,n,v*(v-1)/2)
  for(i in 1:n){
    Xcorpart <- cor(X[i:(i+w),])
    l <- 1
    for(j in 1:(v-1)){
      for(k in (j+1):v){
        Xcor[i,l] = Xcorpart[j,k]
        l <- l+1
      }
    }
  }
  return(Xcor)
}

#' @title The within-phase scatter using R
#' @description Compute the within-phase scatter using R
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param n the number of time points
#' @param a the last observation in the last phase
#' @param b the last observation in the current phase
#' @return the value of within-phase scatter for the current phase
#' @examples 
#' \dontrun{
#' data(mydata)
#' attach(mydata)
#' N <- length(mydata[,1])
#' VR <- VexpressionR(cordata,N,0,N)
#' VR
#' }
#' @import Rcpp
#' @export
VexpressionR <- function(X,n,a,b){
  s <- 0
  h2 <- bandwidthC(X,n)
  for(i in (a+1):(b-1)){
    for(j in (i+1):b){
      s <- s + 2*exp(-sum((X[i,]-X[j,])^2)/(2*h2))
    }
  }
  return((b-a-1)-s/(b-a))
}

#' @title Variance test
#' @description Comparing the average within-phase variance of the original time series and the permuted counterpart when no change points are induced to test the existence of the change point.
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param N the number of time points in X
#' @param v the number of variables in X
#' @param w the width of windows. Compute the running correlations in a window.
#' @return the p-value of test with the null hypothesis that there is no change point in the time series.
#' @examples
#' \dontrun{
#' data(mydata)
#' attach(mydata)
#' N <- length(mydata[,1])
#' v <- length(mydata[1,])
#' w <- 15
#' p_VarKCP <- VarKCP(mydata,N,v,w)
#' p_VarKCP
#' }
#' @import Rcpp
#' @export
VarKCP <- function(X,N,v,w){
  # compute the running correlations for raw data
  Xcor <- CorrelationR(X,N,v,w)
  # the number of time points in Xcor
  n <- N-w
  # the original statistic
  Var0 <- VexpressionC(Xcor,n,0,n)/n 
  
  R <- 999 # the number of permutation
  Var1 <- numeric(R)
  for (i in 1:R){
    n1 <- sample(1:N, size=N, replace = FALSE)
    X1 <- X[n1,]
    X1cor <- CorrelationR(X1,N,v,w)
    # The permutation statistic
    Var1[i] <- VexpressionC(X1cor,n,0,n)/n
  }
  # The significance level of permutation
  p <- mean(c(Var0,Var1)>Var0)
  return(p)
}

#' @title The average within-phase variance
#' @description Given the number of change points, K, compute the average within-phase variance
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param n the number of time points in X
#' @param K the number of change points
#' @param Cpoint contains the locations of the K change points
#' @return The value of average within-phase variance for K+1 phases
#' @import Rcpp
#' @export
Rexpression <- function(Cpoint,X,n,K){
  V <- numeric(K+1) # K+1 is the number of phases
  V[1] <- VexpressionC(X,n,0,Cpoint[1])
  V[K+1] <- VexpressionC(X,n,Cpoint[K],n)
  if(K>1){
    for(i in 2:K){
      V[i] <- VexpressionC(X,n,Cpoint[i-1],Cpoint[i])
    }
  }
  return(sum(V)/n)
}

#' @title Minimize the average within-phase variance
#' @description Compute the minimum of the average within-phase variance with the known number of change points and the unknown locations of them.
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param n the number of time points in X
#' @param K the number of change points
#' @return the minimum of the average within-phase variance
#' @importFrom stats cor runif
#' @importFrom Rsolnp solnp
#' @export
Min_Rexpr <- function(X,n,K){
  Initial <- runif(K,1,n)
  Lower <- numeric(K)+1
  Upper <- numeric(K)+n
  result <- solnp(pars=Initial, fun=Rexpression,
                  LB=Lower, UB=Upper,
                  X=X, n=n, K=K)$values
  return(result[length(result)])
}

#' @title The maximum drop in the minimum of the average within-phase variance
#' @description Given the number of change points, K, compute the maximum drop in the minimum of the average within-phase variance for each change point added.
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param n the number of time points in X
#' @param K the number of change points
#' @return the value of maximum drop variance
#' @import Rcpp
#' @export
Max_Rmindrop <- function(X,n,K){
  Rmins <- Rmindrop <- numeric(K)
  for(i in 1:K){
    Rmins[i] <- Min_Rexpr(X,n,K)
  }
  Rmin0 <- VexpressionC(X,n,0,n)/n
  Rmin <- c(Rmin0,Rmins)
  for(i in 1:K){
    Rmindrop[i] <- Rmin[i+1]-Rmin[i]
  }
  return(max(Rmindrop))
}

#' @title Variance drop test
#' @description Comparing the maximum variance drop of the original time series and the permuted counterpart to test the existence of the change point.
#' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
#' @param N the number of time points in X
#' @param v the number of variables in X
#' @param w the width of windows. Compute the running correlations in a window.
#' @param K the number of change points
#' @return the p-value of test with the null hypothesis that there is no change point in the time series.
#' @examples
#' \dontrun{
#' data(mydata)
#' attach(mydata)
#' N <- length(mydata[,1])
#' v <- length(mydata[1,])
#' w <- 15
#' p_Vardrop <- Vardrop(mydata,N,v,w,2)
#' p_Vardrop
#' }
#' @export
Vardrop <- function(X,N,v,w,K){
  # compute the running correlations for raw data
  Xcor <- CorrelationR(X,N,v,w)
  # the number of time points in Xcor
  n <- N-w
  # the original statistic
  Rmindrop0 <- Max_Rmindrop(Xcor,n,K)
  
  R <- 10 # the number of permutation
  Rmindrops <- numeric(R)
  for (i in 1:R) {
    n1 <- sample(1:N, size=N, replace=FALSE)
    X1 <- X[n1,]
    X1cor <- CorrelationR(X1,N,v,w)
    # The permutation statistic
    Rmindrops[i] <- Max_Rmindrop(X1cor,n,K)
  }
  # The significance level of permutation
  p <- mean(c(Rmindrop0,Rmindrops)>Rmindrop0)
}
