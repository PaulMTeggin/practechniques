% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlapping_data.R
\name{calc_sd_ol_bias_fac}
\alias{calc_sd_ol_bias_fac}
\title{Calculate the overlapped bias factor for sample standard deviation}
\usage{
calc_sd_ol_bias_fac(N, noverlap)
}
\arguments{
\item{N}{The number of overlapped data points.
There will be \code{N + noverlap - 1} in the data prior to taking overlapped samples,
which reduces to \code{N} after overlapping by \code{noverlap}.
E.g. we need 11 more observations than \code{N} in order to have \code{N} annual monthly-overlapped observations.}

\item{noverlap}{The extent of the overlap, e.g. 12 means annual overlaps from monthly data.
1 means no overlap.}
}
\value{
A factor to multiply the sample standard deviation by to give an unbiased estimate.
}
\description{
When overlapped data is used, the sample variance and standard deviation are biased.
This function calculates a factor that corrects for the bias in the sample standard deviation.
It assumes \href{https://en.wikipedia.org/wiki/Bessel\%27s_correction}{Bessel's correction}
has already been applied in the calculation of the sample
standard deviation and removes this to avoid double-counting.
In practice this means this factor can be applied to the results of \code{\link[stats]{sd}}.
}
\details{
The factor is derived in a
\href{https://www.risk.net/risk-management/1509219/error-var-overlapping-intervals}{2009 Risk.net article}
by Sun, Nelken et al. Where \code{noverlap == 1} the factor is 1 and has no numerical effect.
}
\examples{
calc_sd_ol_bias_fac(c(24, 60, 120, 180, 240), 12)
}
\references{
\url{https://www.risk.net/risk-management/1509219/error-var-overlapping-intervals}
}
