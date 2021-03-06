% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getTruth}
\alias{getTruth}
\title{Get the true value of the trend parameter for a given covariance 
weight matrix}
\usage{
getTruth(n = 1e+06, t0 = 6, Upsilon = diag(5), returnL = FALSE)
}
\arguments{
\item{n}{Sample size (should be large to ensure numerical accuracy)}

\item{t0}{Time to compute incidence}

\item{Upsilon}{The covariance weighting matrix}

\item{returnL}{Return the estimated log cumulative incidence ratios}
}
\value{
If \code{returnL = FALSE}, a vector of trend parameter estimates,
the first corresponding to the "intercept" parameter, the second to the "slope"
or trend parameter. If \code{returnL = TRUE}, then a list with named entries
\code{param} (the projection parameters) and \code{L}, a vector of true values 
of L_{j,0}.
}
\description{
Get the true value of the trend parameter for a given covariance 
weight matrix
}
