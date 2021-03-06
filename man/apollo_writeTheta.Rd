% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_writeTheta.R
\name{apollo_writeTheta}
\alias{apollo_writeTheta}
\title{Writes the vector [beta,ll] to a file called modelname_iterations.csv}
\usage{
apollo_writeTheta(beta, ll, modelName)
}
\arguments{
\item{beta}{vector of parameters to be written.}

\item{ll}{scalar representing the loglikelihood of the whole model.}

\item{modelName}{Character. Name of the model.}
}
\value{
Nothing.
}
\description{
Writes the vector [beta,ll] to a file called modelname_iterations.csv
}
