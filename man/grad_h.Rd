% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{grad_h}
\alias{grad_h}
\title{Helper function to compute gradient of trend parameter function.}
\usage{
grad_h(F, D, j_vec)
}
\arguments{
\item{F}{A vector of cumulative incidence estimates with alternating
treatment = 0 and treatment = 1 for each type.}

\item{D}{A matrix of variance estimates corresponding with the vector 
\code{F}.}
}
\description{
Helper function to compute gradient of trend parameter function.
}
