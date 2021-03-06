% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_lc.R
\name{apollo_lc}
\alias{apollo_lc}
\title{Calculates the likelihood of a latent class model}
\usage{
apollo_lc(lc_settings, apollo_inputs, functionality)
}
\arguments{
\item{lc_settings}{List of arguments used by \code{apollo_lc}. It must include the following.
\itemize{
  \item \strong{inClassProb}: List of probabilities. Conditional likelihood for each class. One element per class, in the same order as \code{classProb}.
  \item \strong{classProb}: List of probabilities. Allocation probability for each class. One element per class, in the same order as \code{inClassProb}.
  \item \strong{componentName}: Character. Name given to model component.
}}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{functionality}{Character. Can take different values depending on desired output.
\itemize{
  \item "estimate" Used for model estimation.
  \item "prediction" Used for model predictions.
  \item "validate" Used for validating input.
  \item "zero_LL" Used for calculating null likelihood.
  \item "conditionals" Used for calculating conditionals.
  \item "output" Used for preparing output after model estimation.
  \item "raw" Used for debugging.
  \item "components" Returns \code{P} without changes.
}}
}
\value{
The returned object depends on the value of argument \code{functionality} as follows.
        \itemize{
          \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
          \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all models components, for each class.
          \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but also runs a set of tests on the given arguments.
          \item \strong{\code{"zero_LL"}}: Same as \code{"estimate"}
          \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
          \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
          \item \strong{\code{"raw"}}: Same as \code{"prediction"}
        }
}
\description{
Using the conditional likelihoods of each latent class, as well as teir classification probabilities, calculate the weighted likelihood of the whole model.
}
