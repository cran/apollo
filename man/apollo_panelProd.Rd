% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_panelProd.R
\name{apollo_panelProd}
\alias{apollo_panelProd}
\title{Calculates product across observations from same individual.}
\usage{
apollo_panelProd(P, apollo_inputs, functionality)
}
\arguments{
\item{P}{List of vectors, matrices or 3-dim arrays. Likelihood of the model components.}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

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
Argument \code{P} with (for most functionalities) the original contents after multiplying across observations at the individual level. Shape depends on argument \code{functionality}.
        \itemize{
          \item \strong{\code{"components"}}: Returns \code{P} without changes.
          \item \strong{\code{"conditionals"}}: Returns \code{P} without averaging across draws. Drops all components except \code{"model"}.
          \item \strong{\code{"estimate"}}: Returns \code{P} containing the likelihood of the model after multiplying observations at the individual level. Drops all components except \code{"model"}.
          \item \strong{\code{"gradient"}}: Returns \code{P} containing the gradient of the likelihood after applying the product rule across observations for the same individual.
          \item \strong{\code{"output"}}: Returns \code{P} containing the likelihood of the model after multiplying observations at the individual level.
          \item \strong{\code{"prediction"}}: Returns \code{P} containing the probabilities/likelihoods of all alternatives for all model components averaged across inter-individual draws.
          \item \strong{\code{"preprocess"}}: Returns \code{P} without changes.           
          \item \strong{\code{"raw"}}: Returns \code{P} without changes.
          \item \strong{\code{"report"}}: Returns \code{P} without changes.
          \item \strong{\code{"shares_LL"}}: Returns \code{P} containing the likelihood of the model after multiplying observations at the individual level.
          \item \strong{\code{"validate"}}: Returns \code{P} containing the likelihood of the model averaged across inter-individual draws. Drops all components except \code{"model"}.
          \item \strong{\code{"zero_LL"}}: Returns \code{P} containing the likelihood of the model after multiplying observations at the individual level.
        }
}
\description{
Multiplies likelihood of observations from the same individual, or adds the log of them.
}
\details{
This function should be called inside apollo_probabilities only if the data has a panel structure.
It should be called after apollo_avgIntraDraws if intra-individual draws are used.
}
