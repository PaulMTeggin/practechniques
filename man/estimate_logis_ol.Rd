% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_functions.R
\name{estimate_logis_ol}
\alias{estimate_logis_ol}
\title{Estimate location and scale parameters for the logistic distribution by method of moments allowing for potentially overlapping data}
\usage{
estimate_logis_ol(x, noverlap = 1)
}
\arguments{
\item{x}{The data whose mean and standard deviation are required.}

\item{noverlap}{The extent to which the data is overlapped. \code{noverlap == 1} (the default) means no overlap,
and in this case \code{\link{calc_sd_ol_bias_fac}} returns 1, so no adjustment is applied.}
}
\value{
\describe{
\item{location}{the mean of \code{x}}
\item{scale}{the scale parameter derived from standard deviation of \code{x} adjusted for overlapping data bias, if any}
}
}
\description{
This function serves two purposes: (1) collate the location and scale parameters into a single list object,
and (2) allow for bias in the standard deviation arising from the overlapping data
by using \code{\link{calc_sd_ol_bias_fac}}.
}
\details{
The scale parameter is derived from the sample standard deviation by multiplying by \code{sqrt(3) / pi}.
See \href{https://en.wikipedia.org/wiki/Logistic_distribution}{the formula for logistic variance on Wikipedia}.
Note that the adjustment for overlapping data bias is multiplicative to the sample standard deviation
and so applies to the logistic scale parameter in the same way that it does to the Normal sample standard deviation.
}
\examples{
estimate_logis_ol(rlogis(100))
}
\seealso{
Other Estimation functions: 
\code{\link{estimate_exp_ol}()},
\code{\link{estimate_mean_sd_ol}()}
}
\concept{Estimation functions}
