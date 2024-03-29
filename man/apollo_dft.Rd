% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_dft.R
\name{apollo_dft}
\alias{apollo_dft}
\title{Calculate DFT probabilities}
\usage{
apollo_dft(dft_settings, functionality)
}
\arguments{
\item{dft_settings}{List of settings for the DFT model. It should contain the following elements.
\itemize{
 \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
 \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
 \item \strong{altStart}: A named list with as many elements as alternatives. Each element can be a scalar or vector containing the starting preference value for the alternative.  
 \item \strong{attrScalings}: A named list with as many elements as attributes, or fewer. Each element is a factor that scale the attribute, and can be a scalar, a vector or a matrix/array. They do not need to add up to one for each observation. \code{attrWeights} and \code{attrScalings} are incompatible, and they should not be both defined for an attribute. Default is 1 for all attributes.
 \item \strong{attrValues}: A named list with as many elements as alternatives. Each element is itself a named list of vectors of the alternative attributes for each observation (usually a column from the database). All alternatives must have the same attributes (can be set to zero if not relevant).
 \item \strong{attrWeights}: A named list with as many elements as attributes, or fewer. Each element is the weight of the attribute, and can be a scalar, a vector with as many elements as observations, or a matrix/array if random. They should add up to one for each observation and draw (if present), and will be re-scaled if they do not. \code{attrWeights} and \code{attrScalings} are incompatible, and they should not be both defined for an attribute. Default is 1 for all attributes.
 \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
 \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
 \item \strong{procPars}: A list containing the four DFT 'process parameters'
    \itemize{
      \item \strong{error_sd}: Numeric scalar or vector. The standard deviation of the the error term in each timestep.
      \item \strong{timesteps}: Numeric scalar or vector. Number of timesteps to consider. Should be an integer bigger than 0.
      \item \strong{phi1}: Numeric scalar or vector. Sensitivity.
      \item \strong{phi2}: Numeric scalar or vector. Process parameter.
    }
 \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
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
          \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
          \item \strong{\code{"gradient"}}: Not implemented.
          \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
          \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
          \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{dft_settings}.
          \item \strong{\code{"raw"}}: Same as \code{"prediction"}
          \item \strong{\code{"report"}}: Choice overview.
          \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
          \item \strong{\code{"validate"}}: Same as \code{"estimate"}
          \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
        }
}
\description{
Calculate probabilities of a Decision Field Theory (DFT) model and can also perform other operations based on the value of the \code{functionality} argument.
}
\section{References}{

Hancock, T.; Hess, S. and Choudhury, C. (2018) Decision field theory: Improvements to current methodology and comparisons with standard choice modelling techniques. Transportation Research 107B, 18 - 40.
Hancock, T.; Hess, S. and Choudhury, C. (Submitted) An accumulation of preference: two alternative dynamic models for understanding transport choices.
Roe, R.; Busemeyer, J. and Townsend, J. (2001) Multialternative decision field theory: A dynamic connectionist model of decision making. Psychological Review 108, 370
}

