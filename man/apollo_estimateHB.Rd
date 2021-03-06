% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_estimateHB.R
\name{apollo_estimateHB}
\alias{apollo_estimateHB}
\title{Estimates model}
\usage{
apollo_estimateHB(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings = NA
)
}
\arguments{
\item{apollo_beta}{Named numeric vector. Names and values for parameters.}

\item{apollo_fixed}{Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.}

\item{apollo_probabilities}{Function. Returns probabilities of the model to be estimated. Must receive three arguments:
\itemize{
  \item apollo_beta: Named numeric vector. Names and values of model parameters.
  \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
  \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
}}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{estimate_settings}{List. Options controlling the estimation process.
\itemize{
  \item \strong{estimationRoutine}: Character. Estimation method. Can take values "bfgs", "bhhh", or "nr".
                                    Used only if \code{apollo_control$HB} is FALSE. Default is "bfgs".
  \item \strong{maxIterations}: Numeric. Maximum number of iterations of the estimation routine before stopping.
                                Used only if \code{apollo_control$HB} is FALSE. Default is 200.
  \item \strong{writeIter}: Boolean. Writes value of the parameters in each iteration to a csv file. 
                            Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
  \item \strong{hessianRoutine}: Character. Name of routine used to calculate the Hessian of the loglikelihood 
                                 function after estimation. Valid values are \code{"numDeriv"} (default) and 
                                 \code{"maxLik"} to use the routines in those packages, and \code{"none"} to avoid 
                                 estimating the Hessian (and the covariance matrix). Only used if \code{apollo_control$HB=FALSE}.
  \item \strong{printLevel}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. 
                             Ignored if apollo_control$HB is TRUE. Default is 3.
  \item \strong{constraints}: Constraints on parameters to estimate. Should ignore fixed parameters. 
                              See argument \code{constraints} in \link[maxLik]{maxBFGS} for more details.
  \item \strong{scaling}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. 
                          If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but 
                          \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter
                          b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, 
                          specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
 \item \strong{numDeriv_settings}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. 
                                   See argument \code{method.args} in \link[numDeriv]{grad} for more details.
 \item \strong{bootstrapSE}: Numeric. Number of bootstrap samples to calculate standard errors. Default is 0, meaning
                             no bootstrap s.e. will be calculated. Number must zero or a positive integer. Only used
                             if \code{apollo_control$HB} is \code{FALSE}.
 \item \strong{bootstrapSeed}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples.
                               Only used if \code{bootstrapSE>0}. Default is 24.
 \item \strong{silent}: Boolean. If TRUE, no information is printed to the console during estimation. Default is FALSE.
}}
}
\value{
model object
}
\description{
Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
}
\details{
This is the main function of the Apollo package. The estimation process begins by checking the definition of
\code{apollo_probabilities} by estimating it at the starting values. Then it runs the function with argument \code{functionality="validate"}.
If the user requested more than one core for estimation (i.e. \code{apollo_control$nCores>1}), and no bayesian estimation is used
(i.e. \code{apollo_control$HB=FALSE}), then a cluster is created. Using a cluster at least doubles the requires RAM, as the database
must be copied into the cluster.
If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
iterations. If bayesian estimation is used, estimation will finish once the predefined number of iterations are completed.
This functions does not save results into a file nor prints them into the console, so if users want to see and store estimation the results,
they must make sure to call function \code{apollo_modelOutput} and/or \code{apollo_saveOutput} afterwards.
}
