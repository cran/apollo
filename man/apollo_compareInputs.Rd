% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_compareInputs.R
\name{apollo_compareInputs}
\alias{apollo_compareInputs}
\title{Compares the content of apollo_inputs to their counterparts in the global environment}
\usage{
apollo_compareInputs(apollo_inputs)
}
\arguments{
\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}
}
\value{
Logical. TRUE if the content of \code{apollo_inputs} is the same than the one in the global environment, FALSE otherwise.
}
\description{
Compares the content of apollo_inputs to their counterparts in the global environment
}
