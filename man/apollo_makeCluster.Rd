% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_makeCluster.R
\name{apollo_makeCluster}
\alias{apollo_makeCluster}
\title{Creates cluster for estimation.}
\usage{
apollo_makeCluster(
  apollo_probabilities,
  apollo_inputs,
  silent = FALSE,
  cleanMemory = FALSE
)
}
\arguments{
\item{apollo_probabilities}{Function. Returns probabilities of the model to be estimated. Must receive three arguments:
\itemize{
  \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
  \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
  \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
}}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{silent}{Boolean. If TRUE, no messages are printed to the terminal. FALSE by default. It overrides \code{apollo_inputs$silent}.}

\item{cleanMemory}{Boolean. If TRUE, it saves apollo_inputs to disc, and removes database and draws from 
the apollo_inputs in .GlobalEnv and the parent environment.}
}
\value{
Cluster (i.e. an object of class cluster from package parallel)
}
\description{
Splits data, creates cluster and loads different pieces of the database on each worker.
}
\details{
Internal use only. Called by \code{apollo_estimate} before estimation. Using multiple cores greatly increases memory consumption.
}
