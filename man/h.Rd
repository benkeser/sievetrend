% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{h}
\alias{h}
\title{Helper function to compute estimate of trend parameter (also returns
estimate of intercept).}
\usage{
h(F, D, j_vec)
}
\arguments{
\item{F}{A vector of cumulative incidence estimates with alternating
treatment = 0 and treatment = 1 for each type.}

\item{D}{A matrix of variance estimates corresponding with the vector 
\code{F}.}

\item{j_vec}{Vector of genetic distances}
}
\description{
Helper function to compute estimate of trend parameter (also returns
estimate of intercept).
}
