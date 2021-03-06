% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_insertComponentName.R
\name{apollo_insertComponentName}
\alias{apollo_insertComponentName}
\title{Adds componentName2 to model calls}
\usage{
apollo_insertComponentName(e)
}
\arguments{
\item{e}{An expression or a function. It will usually be apollo_probabilities.}
}
\value{
The original argument 'e' but modified to incorporate a new setting
        called 'componentName2' to every call to apollo_<model> (e.g. 
       apollo_mnl, apollo_nl, etc.).
}
\description{
Adds componentName2 to model calls
}
