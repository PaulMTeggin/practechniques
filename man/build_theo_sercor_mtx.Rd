% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlapping_data.R
\name{build_theo_sercor_mtx}
\alias{build_theo_sercor_mtx}
\title{Build a matrix of theoretical serial correlations}
\usage{
build_theo_sercor_mtx(N, noverlap)
}
\arguments{
\item{N}{The size of the matrix (number of rows and columns)}

\item{noverlap}{The extent of the overlap, e.g. 12 means annual overlaps from monthly data.
1 means no overlap.}
}
\value{
The matrix of theoretical serial correlations.
}
\description{
The rows and columns of the matrix represent indexes to the overlapped data.
All values on the leading diagonal will be 1, since a variable is always perfectly correlated
with itself. Values off the diagonal will depend on the size of the overlap (\code{noverlap}).
}
\examples{
build_theo_sercor_mtx(24, 12)
}
\seealso{
Other Overlapping data functions: 
\code{\link{calc_theo_sercor}()}
}
\concept{Overlapping data functions}
