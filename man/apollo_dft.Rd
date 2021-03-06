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
  \item \strong{alternatives}: Named numeric vector. Names of alternatives and 
                               their corresponding value in \code{choiceVar}.
  \item \strong{avail}: Named list of numeric vectors or scalars. Availabilities 
                        of alternatives, one element per alternative. Names of 
                        elements must match those in \code{alternatives}. Values 
                        can be 0 or 1.
  \item \strong{choiceVar}: Numeric vector. Contains choices for all observations. 
                            It will usually be a column from the database. Values 
                            are defined in \code{alternatives}.
  \item \strong{attrValues}: A named list with as many elements as alternatives. 
                             Each element is itself a named list of vectors of the 
                             alternative attributes for each observation (usually a 
                             column from the database). All alternatives must have 
                             the same attributes (can be set to zero if not relevant).
  \item \strong{altStart}: A named list with as many elements as alternatives. 
                           Each elment can be a scalar or vector containing the 
                           starting preference value for the alternative.  
  \item \strong{attrWeights}: A named list with as many elements as attributes, 
                              or fewer. Each element is the weight of the attribute, 
                              and can be a scalar, a vector with as many elements as 
                              observations, or a matrix/array if random. They should 
                              add up to one for each observation and draw (if present), 
                              and will be re-scaled if they do not. \code{attrWeights} 
                              and \code{attrScalings} are incompatible, and they should 
                              not be both defined for an attribute. Default is 1 for 
                              all attributes.
  \item \strong{attrScalings}: A named list with as many elements as attributes, 
                               or fewer. Each element is a factor that scale the 
                               attribute, and can be a scalar, a vector or a 
                               matrix/array. They do not need to add up to one 
                               for each observation. \code{attrWeights} and 
                               \code{attrScalings} are incompatible, and they 
                               should not be both defined for an attribute. 
                               Default is 1 for all attributes.
  \item \strong{procPars}: A list containing the four DFT 'process parameters'
    \itemize{
      \item \strong{error_sd}: Numeric scalar or vector. The standard deviation of the the error term in each timestep.
      \item \strong{timesteps}: Numeric scalar or vector. Number of timesteps to consider. Should be an integer bigger than 0.
      \item \strong{phi1}: Numeric scalar or vector. Sensitivity.
      \item \strong{phi2}: Numeric scalar or vector. Process parameter.
    }
 \item \strong{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
 \item \strong{componentName}: Character. Name given to model component.
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
          \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
          \item \strong{\code{"validate"}}: Same as \code{"estimate"}
          \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
          \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
          \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
          \item \strong{\code{"raw"}}: Same as \code{"prediction"}
        }
}
\description{
Calculate probabilities of a Decision Field Theory (DFT) with external thresholds.
}
\section{References}{

Hancock, T.; Hess, S. and Choudhury, C. (2018) Decision field theory: Improvements to current methodology and comparisons with standard choice modelling techniques. Transportation Research 107B, 18 - 40.
Hancock, T.; Hess, S. and Choudhury, C. (Submitted) An accumulation of preference: two alternative dynamic models for understanding transport choices.
Roe, R.; Busemeyer, J. and Townsend, J. (2001) Multialternative decision field theory: A dynamic connectionist model of decision making. Psychological Review 108, 370
}

