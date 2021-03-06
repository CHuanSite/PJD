% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqPCA.rank.R
\name{seqPCA.rank}
\alias{seqPCA.rank}
\title{Sequential Decomposition with PCA and automatic rank selection}
\usage{
seqPCA.rank(dataset, group, total_number = NULL, threshold, backup = 0)
}
\arguments{
\item{dataset}{A list of dataset to be analyzed}

\item{group}{A list of grouping of the datasets, indicating the relationship between datasets}

\item{total_number}{Total number of components will be extracted, if default value is set to NA, then BEMA will be used.}

\item{threshold}{The threshold used to cutoff the eigenvalues}

\item{backup}{A backup variable, which permits the overselection of the components by BEMA}
}
\value{
A list contains the component and the score of each dataset on every component after seqPCA algorithm
}
\description{
Sequential decomposition of several matrices with PCA, the rank selection procedure is automatic based on BEMA
}
\examples{
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
threshold = c(3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5)
res_seqPCA = seqPCA.rank(dataset, group, threshold = threshold)

}
\keyword{PCA}
\keyword{rank,}
\keyword{sequential,}
