% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.apollo.R
\name{summary.apollo}
\alias{summary.apollo}
\title{Prints summary of Apollo model}
\usage{
\method{summary}{apollo}(object, ..., pTwoSided = FALSE)
}
\arguments{
\item{object}{Model object. Estimated model object as returned by function \link{apollo_estimate}.}

\item{...}{further arguments passed to or from other methods.}

\item{pTwoSided}{Logical. Should two-sided p-values be printed instead of one-sided p-values. FALSE by default.
#' @return nothing.}
}
\description{
Receives an estimated model object and prints a summary using the generic summary function.
}
