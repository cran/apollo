% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_deltaMethod.R
\name{apollo_deltaMethod}
\alias{apollo_deltaMethod}
\title{Delta method for Apollo models}
\usage{
apollo_deltaMethod(model, deltaMethod_settings)
}
\arguments{
\item{model}{Model object. Estimated model object as returned by function \link{apollo_estimate}.}

\item{deltaMethod_settings}{List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
\itemize{
  \item \strong{\code{expression}}: Character vector. A character vector with a single or multiple arbitrary functions of the estimated parameters, as text. 
                    For example: \code{c(VTT="b1/b2*60")}. Each expression can only contain model parameters (estimated or fixed), 
                    numeric values, and operands. At least one of the parameters used needs to not have been fixed in estimation. Variables in the database
                    cannot be included. If the user does not provide a name for an expression, then the expression itself is used in the output. If this setting is provided, then \code{operation}, 
                    \code{parName1}, \code{parName2}, \code{multPar1} and \code{multPar2} are
                    ignored.
  \item \strong{\code{allPairs}}: Logical. If set to TRUE, Delta method calculations are carried out for the ratio and difference for all pairs of parameters and returned as two separate matrices with values and t-ratios. FALSE by default.
  \item \strong{\code{varcov}}: Character. Type of variance-covariance matrix to use in calculations. 
                           It can take values \code{"classical"}, \code{"robust"} and 
                           \code{"bootstrap"}. Default is \code{"robust"}.
  \item \strong{\code{printPVal}}: Logical or Scalar. TRUE or 1 for printing p-values for one-sided test, 2 for printing p-values for two-sided test, FALSE for not printing p-values. FALSE by default.
  \item \strong{\code{operation}}: Character. Function to calculate the delta method for. See details. Not used if \code{expression} is provided.
  \item \strong{\code{parName1}}: Character. Name of the first parameter if \code{operation} is used. See details. Not used if \code{expression} is provided.
  \item \strong{\code{parName2}}: Character. Name of the second parameter if \code{operation} is used. See details. Not used if \code{expression} is provided.. Optional depending on \code{operation}.
  \item \strong{\code{multPar1}}: Numeric scalar. An optional value to scale \code{parName1}. Not used if \code{expression} is provided.
  \item \strong{\code{multPar2}}: Numeric scalar. An optional value to scale \code{parName2}. Not used if \code{expression} is provided.
}}
}
\value{
Matrix containing value, s.e. and t-ratio resulting from the requested expression or operation. This is also printed to screen.
}
\description{
Applies the Delta method to calculate the standard errors of transformations of parameters.
}
\details{
\code{apollo_deltaMethod} can be used in two ways. The first and recommended way is to provide an 
element called \code{expression} inside its argument \code{deltaMethod_settings}. \code{expression} 
should contain the expression or expressions for which the standard error is/are to be calculated, as text. For 
example, to calculate the ratio between parameters b1 and b2, \code{expression=c(vtt="b1/b2")} should be used.

The second method is to provide the name of a specific operation inside \code{deltaMethod_settings}.
The following five operations are supported.
\itemize{
  \item \strong{\code{sum}}: Calculates the s.e. of \code{parName1} + \code{parName2}
  \item \strong{\code{diff}}: Calculates the s.e. of \code{parName1} - \code{parName2} and \code{parName2} - \code{parName1}
  \item \strong{\code{prod}}: Calculates the s.e. of \code{parName1}*\code{parName2}
  \item \strong{\code{ratio}}: Calculates the s.e. of \code{parName1}/\code{parName2} and \code{parName2}/\code{parName1}
  \item \strong{\code{exp}}: Calculates the s.e. of \code{exp(parName1)}
  \item \strong{\code{logistic}}: If only \code{parName1} is provided, it calculates the s.e. of
                  \code{exp(parName1)/(1+exp(parName1))} and \code{1/(1+exp(parName1))}.
                  If \code{parName1} and \code{parName2} are provided, it calculates
                  \code{exp(par_i)/(1+exp(parName1)+exp(parName2))} for i=1, 2, and 3 (par_3 = 1).
  \item \strong{\code{lognormal}}: Calculates the mean and s.d. of a lognormal distribution based on the mean 
                   (\code{parName1}) and s.d. (\code{parName2}) of the underlying normal.
}

By default, \code{apollo_deltaMethod} uses the robust covariance matrix. However, the user can change this through the \code{varcov} setting.
}
