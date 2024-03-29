% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xmatrix_function.R
\name{xmatrix_function}
\alias{xmatrix_function}
\title{X-matrix}
\usage{
xmatrix_function(se, sp, nk)
}
\arguments{
\item{se}{sensitivity/PPA}

\item{sp}{specificity/NPA}

\item{nk}{input data vector}
}
\value{
x-matrix
}
\description{
This function calculates x-matrix of the likelihood function
}
\examples{
xmatrix_function(0.99,0.9,nk=c(6535,2816,2128,1929,1699,1399,1288, 41051,9273,8172,7227,6273,5430,4780))
}
