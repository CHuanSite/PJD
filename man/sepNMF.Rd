% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sepNMF.R
\name{sepNMF}
\alias{sepNMF}
\title{Single Data Set Decomposition with Nonnegative Matrix Factorization}
\usage{
sepNMF(dataset, comp_num, weighting = NULL, perturbation = 1e-04)
}
\arguments{
\item{dataset}{A dataframe/matrix to be decomposed}

\item{comp_num}{Number of NMFs to be extracted}

\item{weighting}{Weighting of each dataset, initialized to be NULL}

\item{perturbation}{A small perturbation to ensure nmf works well}
}
\value{
A list of scores and component
}
\description{
Apply NMF (Nonnegative Matrix Factorization) to a single data set
}
\examples{
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
comp_num = 2
res_sepNMF = sepNMF(dataset, comp_num)

}
\keyword{NMF}
\keyword{analysis,}
\keyword{separate}
