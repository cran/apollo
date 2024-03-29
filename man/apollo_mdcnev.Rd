% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_mdcnev.R
\name{apollo_mdcnev}
\alias{apollo_mdcnev}
\title{Calculates MDCNEV likelihoods}
\usage{
apollo_mdcnev(mdcnev_settings, functionality)
}
\arguments{
\item{mdcnev_settings}{List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
\itemize{
  \item \strong{\code{alpha}}: Named list. Alpha parameters for each alternative, including for the outside good. As many elements as alternatives.
  \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
  \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'utilities'.
  \item \strong{\code{budget}}: Numeric vector. Budget for each observation.
\item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
  \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
  \item \strong{\code{cost}}: Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
  \item \strong{\code{gamma}}: Named list. Gamma parameters for each alternative, including for the outside good. As many elements as alternatives.
  \item \strong{\code{mdcnevNests}}: Named list. Lambda parameters for each nest. Elements must be named with the nest name. The lambda at the root is fixed to 1, and therefore must be no be defined. The value of the estimated mdcnevNests parameters should be between 0 and 1 to ensure consistency with random utility maximization.
  \item \strong{\code{mdcnevStructure}}: Numeric matrix. One row per nest and one column per alternative. Each element of the matrix is 1 if an alternative belongs to the corresponding nest.
  \item \strong{\code{outside}}: Character. Alternative name for the outside good. Default is "outside"
  \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
  \item \strong{\code{utilities}}: Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
}}

\item{functionality}{Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call \code{apollo_probabilities}, though the user can also call \code{apollo_probabilities} manually with a given functionality for testing/debugging. Possible values are:
\itemize{
  \item \strong{\code{"components"}}: For further processing/debugging, produces likelihood for each model component (if multiple components are present), at the level of individual draws and observations.
  \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
  \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
  \item \strong{\code{"gradient"}}: For model estimation, produces analytical gradients of the likelihood, where possible.
  \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
  \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
  \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
  \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
  \item \strong{\code{"report"}}: Prepares output summarising model and choiceset structure.
  \item \strong{\code{"shares_LL"}}: Produces overall model likelihood with constants only.
  \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
  \item \strong{\code{"zero_LL"}}: Produces overall model likelihood with all parameters at zero.
}}
}
\value{
The returned object depends on the value of argument \code{functionality} as follows.
        \itemize{
          \item \strong{\code{"components"}}: Same as \code{"estimate"}
          \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
          \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the observed consumption for each observation.
          \item \strong{\code{"gradient"}}: Not implemented
          \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
          \item \strong{\code{"prediction"}}: A matrix with one row per observation, and columns indicating means and s.d. of continuous and discrete predicted consumptions.
          \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{mdcnev_settings}.
          \item \strong{\code{"raw"}}: Same as \code{"estimate"}
          \item \strong{\code{"report"}}: Dependent variable overview.
          \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
          \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
          \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
        }
}
\description{
Calculates the likelihoods of a Multiple Discrete Continuous Nested Extreme Value (MDCNEV) model with an outside good and can also perform other operations based on the value of the \code{functionality} argument.
}
