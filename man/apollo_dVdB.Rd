% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_dVdB.R
\name{apollo_dVdB}
\alias{apollo_dVdB}
\title{Calculates gradients of utility functions}
\usage{
apollo_dVdB(apollo_beta, apollo_inputs, V)
}
\arguments{
\item{apollo_beta}{Named numeric vector of parameters.}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{V}{List of functions}
}
\value{
Named list. Each element is itself a list of functions: the partial derivatives of the elements of V.
}
\description{
Calculates gradients (derivatives) of utility functions.
}
