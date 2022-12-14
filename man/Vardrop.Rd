% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KCPpermR.R
\name{Vardrop}
\alias{Vardrop}
\title{Variance drop test}
\usage{
Vardrop(X, N, v, w, K)
}
\arguments{
\item{X}{a matrix contains multivariate time series data. The rows are time points and the columns are variables.}

\item{N}{the number of time points in X}

\item{v}{the number of variables in X}

\item{w}{the width of windows. Compute the running correlations in a window.}

\item{K}{the number of change points}
}
\value{
the p-value of test with the null hypothesis that there is no change point in the time series.
}
\description{
Comparing the maximum variance drop of the original time series and the permuted counterpart to test the existence of the change point.
}
\examples{
\dontrun{
data(mydata)
attach(mydata)
N <- length(mydata[,1])
v <- length(mydata[1,])
w <- 15
p_Vardrop <- Vardrop(mydata,N,v,w,2)
p_Vardrop
}
}
