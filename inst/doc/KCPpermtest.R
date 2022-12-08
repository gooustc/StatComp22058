## ----eval=FALSE---------------------------------------------------------------
#  # Sort function
#  NumericVector SortC(NumericVector x){
#    x.sort();
#    return x;
#  }
#  # Median function
#  double medianC(NumericVector A){
#    NumericVector B = SortC(A);
#    int n = B.size();
#    double m;
#    if(n%2 != 0){
#      m = B[n/2];
#    }else{
#      m = (B[n/2]+B[n/2-1])/2.0;
#    }
#    return m;
#  }
#  
#  double bandwidthC(NumericMatrix X,int n){
#    NumericVector d(n*(n-1)/2);
#    int k = 0;
#    for(int i=0; i<(n-1); i++){
#      for(int j=i+1; j<n; j++){
#        d[k] = sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_)));
#        k++;
#      }
#    }
#    return medianC(d);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  CorrelationR <- function(X,N,v,w){
#    n <- N-w
#    Xcor <- matrix(0,n,v*(v-1)/2)
#    for(i in 1:n){
#      Xcorpart <- cor(X[i:(i+w),])
#      l <- 1
#      for(j in 1:(v-1)){
#        for(k in (j+1):v){
#          Xcor[i,l] = Xcorpart[j,k]
#          l <- l+1
#        }
#      }
#    }
#    return(Xcor)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  VexpressionR <- function(X,n,a,b){
#    s <- 0
#    h2 <- bandwidthC(X,n)
#    for(i in (a+1):(b-1)){
#      for(j in (i+1):b){
#        s <- s + 2*exp(-sum((X[i,]-X[j,])^2)/(2*h2))
#      }
#    }
#    return((b-a-1)-s/(b-a))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  double VexpressionC(NumericMatrix X,int n,double a,double b){
#    double s = 0;
#    double h2 = bandwidthC(X,n);
#    for(int i=a; i<(b-1); i++){
#      for(int j=i+1; j<b; j++){
#        s = s + 2*exp(-sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_)))/(2*h2));
#      }
#    }
#    return (b-a-1)-s/(b-a);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  VarKCP <- function(X,N,v,w){
#    # compute the running correlations for raw data
#    Xcor <- CorrelationR(X,N,v,w)
#    # the number of time points in Xcor
#    n <- N-w
#    # the original statistic
#    Var0 <- VexpressionC(Xcor,n,0,n)/n
#  
#    R <- 999 # the number of permutation
#    Var1 <- numeric(R)
#    for (i in 1:R){
#      n1 <- sample(1:N, size=N, replace = FALSE)
#      X1 <- X[n1,]
#      X1cor <- CorrelationR(X1,N,v,w)
#      # The permutation statistic
#      Var1[i] <- VexpressionC(X1cor,n,0,n)/n
#    }
#    # The significance level of permutation
#    p <- mean(c(Var0,Var1)>Var0)
#    return(p)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # The average within-phase variance
#  Rexpression <- function(Cpoint,X,n,K){
#    V <- numeric(K+1) # K+1 is the number of phases
#    V[1] <- VexpressionC(X,n,0,Cpoint[1])
#    V[K+1] <- VexpressionC(X,n,Cpoint[K],n)
#    if(K>1){
#      for(i in 2:K){
#        V[i] <- VexpressionC(X,n,Cpoint[i-1],Cpoint[i])
#      }
#    }
#    return(sum(V)/n)
#  }
#  # Minimize the average within-phase variance
#  Min_Rexpr <- function(X,n,K){
#    Initial <- runif(K,1,n)
#    Lower <- numeric(K)+1
#    Upper <- numeric(K)+n
#    result <- solnp(par=Initial, fun=Rexpression,
#                    LB=Lower, UB=Upper,
#                    X=X, n=n, K=K)$values
#    return(result[length(result)])
#  }
#  # The maximum drop in the minimum of the average within-phase variance
#  Max_Rmindrop <- function(X,n,K){
#    Rmins <- Rmindrop <- numeric(K)
#    for(i in 1:K){
#      Rmins[i] <- Min_Rexpr(X,n,K)
#    }
#    Rmin0 <- VexpressionC(X,n,0,n)/n
#    Rmin <- c(Rmin0,Rmins)
#    for(i in 1:K){
#      Rmindrop[i] <- Rmin[i+1]-Rmin[i]
#    }
#    return(max(Rmindrop))
#  }
#  
#  Vardrop <- function(X,N,v,w,K){
#    # compute the running correlations for raw data
#    Xcor <- CorrelationR(X,N,v,w)
#    # the number of time points in Xcor
#    n <- N-w
#    # the original statistic
#    Rmindrop0 <- Max_Rmindrop(Xcor,n,K)
#  
#    R <- 10 # the number of permutation
#    Rmindrops <- numeric(R)
#    for (i in 1:R) {
#      n1 <- sample(1:N, size=N, replace=FALSE)
#      X1 <- X[n1,]
#      X1cor <- CorrelationR(X1,N,v,w)
#      # The permutation statistic
#      Rmindrops[i] <- Max_Rmindrop(X1cor,n,K)
#    }
#    # The significance level of permutation
#    p <- mean(c(Rmindrop0,Rmindrops)>Rmindrop0)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp22058)
library(microbenchmark)
data(mydata)

N <- length(mydata[,1])# The number of time points
v <- length(mydata[1,])# The number of variables
w <- 15 # Set the width of windows for running correlation computation is 15
cordata <- CorrelationR(mydata,N,v,w)
n <- length(cordata[,1])

tm2 <- microbenchmark(
  VR = VexpressionR(cordata,n,0,n),
  VC = VexpressionC(cordata,n,0,n)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

## ----eval=TRUE----------------------------------------------------------------
p_VarKCP <- VarKCP(mydata,N,v,w)
p_VarKCP

## ----eval=FALSE---------------------------------------------------------------
#  library(Rsolnp)
#  p_Vardrop <- Vardrop(mydata,N,v,w,2)
#  p_Vardrop

