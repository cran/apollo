% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliaryFunctions.R
\name{aux_validateRows}
\alias{aux_validateRows}
\title{Validates and expands rows if necessary.}
\usage{
aux_validateRows(rows, componentName = NULL, apollo_inputs = NULL)
}
\arguments{
\item{rows}{Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}. Set to \code{"all"} by default if omitted.}

\item{componentName}{Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}
}
\description{
Validates and expands rows if necessary.
}
