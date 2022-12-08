#include <Rcpp.h>
using namespace Rcpp;

//' @title Sort function using Rcpp
//' @description Sort the values in a NumericVector using Rcpp
//' @param x a NumericVector
//' @return a vector contains the values in the NumericVector x sorted from smallest to largest.
//' @examples
//' \dontrun{
//' A <- c(2,1,4,5)
//' SortC(A)
//' }
//' @export
// [[Rcpp::export]]
NumericVector SortC(NumericVector x){
  x.sort();
  return x;
}

//' @title Median function using Rcpp
//' @description Compute the median of the values in a NumericVector using Rcpp
//' @param A a NumericVector
//' @return the median of the values in the NumericVector A.
//' @examples
//' \dontrun{
//' A <- c(2,1,4,5)
//' medianC(A)
//' }
//' @export
// [[Rcpp::export]]
double medianC(NumericVector A){
  NumericVector B = SortC(A);
  int n = B.size();
  double m;
  if(n%2 != 0){
    m = B[n/2];
  }else{
    m = (B[n/2]+B[n/2-1])/2.0;
  }
  return m;
}

//' @title A Function to compute the bandwidth using Rcpp
//' @description Compute the bandwidth using Rcpp. For a time series dataset, the bandwidth is the median Euclidean distance for the data corresponding to all time points.
//' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
//' @param n the number of time points
//' @return the value of bandwidth for the data in X.
//' @examples
//' \dontrun{
//' data(mydata)
//' attach(mydata)
//' N <- length(mydata[,1])
//' h <- bandwidthC(cordata,N)
//' h
//' }
//' @export
// [[Rcpp::export]]
double bandwidthC(NumericMatrix X,int n){
  NumericVector d(n*(n-1)/2);
  int k = 0;
  for(int i=0; i<(n-1); i++){
    for(int j=i+1; j<n; j++){
      d[k] = sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_)));
      k++;
    }
  }
  return medianC(d);
}

//' @title Running correlations function using Rcpp
//' @description Compute the running correlations using Rcpp
//' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
//' @param N the number of time points in X
//' @param v the number of variables in X
//' @param w the width of windows. Compute the running correlations in a window.
//' @return a matrix of the running correlations for all variables at these time points.
//' @examples
//' \dontrun{
//' data(mydata)
//' attach(mydata)
//' N <- length(mydata[,1])
//' v <- length(mydata[1,])
//' w <- 15
//' cordata <- CorrelationR(mydata,N,v,w)
//' cordata
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix CorrelationC(NumericMatrix X,int N,int v,int w){
  // Obtain the namespace of the stats package.
  // The namespace is assigned to "stats_env".
  Environment stats_env = Environment::namespace_env("stats");
  
  // From the above defined "stats_env", we get function "cor" and 
  // assign it to "cor_in_cpp".
  Function cor_in_cpp = stats_env["cor"];
  
  int n = N-w;
  NumericMatrix Xcor(n,v);
  for(int i=0; i<n; i++){
    NumericMatrix Xcorpart = cor_in_cpp(X(Range(i,(i+w)),_));
    int l = 0;
    for(int j=0; j<v-1; j++){
      for(int k=(j+1); k<v; k++){
        Xcor(i,l) = Xcorpart(j,k);
        l++;
      }
    }
  }
  return Xcor;
}

//' @title The within-phase scatter using Rcpp
//' @description Compute the within-phase scatter using Rcpp
//' @param X a matrix contains multivariate time series data. The rows are time points and the columns are variables.
//' @param n the number of time points
//' @param a the last observation in the last phase
//' @param b the last observation in the current phase
//' @return the value of within-phase scatter for the current phase
//' @examples
//' \dontrun{
//' data(mydata)
//' attach(mydata)
//' N <- length(mydata[,1])
//' VC <- VexpressionC(cordata,N,0,N)
//' VC
//' }
//' @export
// [[Rcpp::export]]
double VexpressionC(NumericMatrix X,int n,double a,double b){
  double s = 0;
  double h2 = bandwidthC(X,n);
  for(int i=a; i<(b-1); i++){
    for(int j=i+1; j<b; j++){
      s = s + 2*exp(-sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_)))/(2*h2));
    }
  }
  return (b-a-1)-s/(b-a);
}
