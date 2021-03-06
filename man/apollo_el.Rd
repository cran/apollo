% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_el.R
\name{apollo_el}
\alias{apollo_el}
\title{Calculates Exploded Logit probabilities}
\usage{
apollo_el(el_settings, functionality)
}
\arguments{
\item{el_settings}{List of inputs of the Exploded Logit model. It shoud contain the following.
\itemize{
 \item \strong{\code{"alternatives"}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
 \item \strong{\code{"avail"}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
 \item \strong{\code{"choiceVars"}}: List of numeric vectors. Contain choices for each position of the ranking. The list must be ordered with the best choice first, second best second, etc. It will usually be a list of columns from the database. Use value -1 if a stage does not apply for a given observations (e.g. when some individuals have shorter rankings).
 \item \strong{\code{"V"}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
 \item \strong{\code{"scales"}}: List of vectors. Scale factors of each Logit model. At least one element should be normalized to 1. If omitted, scale=1 for all positions is assumed.
 \item \strong{\code{"rows"}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
 \item \strong{\code{"componentName"}}: Character. Name given to model component.
}}

\item{functionality}{Character. Can take different values depending on desired output.
\itemize{
  \item \code{"estimate"}: Used for model estimation.
  \item \code{"prediction"}: Used for model predictions.
  \item \code{"validate"}: Used for validating input.
  \item \code{"zero_LL"}: Used for calculating null likelihood.
  \item \code{"conditionals"}: Used for calculating conditionals.
  \item \code{"output"}: Used for preparing output after model estimation.
  \item \code{"raw"}: Used for debugging.
}}
}
\value{
The returned object depends on the value of argument \code{functionality} as follows.
        \itemize{
          \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
          \item \strong{\code{"prediction"}}: Not applicable (\code{NA}).
          \item \strong{\code{"validate"}}: Same as \code{"estimate"}
          \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
          \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
          \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
          \item \strong{\code{"raw"}}: Same as \code{"estimate"}
        }
}
\description{
Calculates the probabilities of an Exploded Logit model and can also perform other operations based on the value of the \code{functionality} argument.
The function calculates the probability of a ranking as a product of Multinomial Logit models with gradually reducing availability, where scale differences can be allowed for.
}
