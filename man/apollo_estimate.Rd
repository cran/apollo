% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_estimate.R
\name{apollo_estimate}
\alias{apollo_estimate}
\title{Estimates model}
\usage{
apollo_estimate(
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
  \item \strong{writeIter}: Logical. Writes value of the parameters in each iteration to a csv file. 
                            Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
  \item \strong{hessianRoutine}: Character. Name of routine used to calculate the Hessian of the loglikelihood 
                                 function after estimation. Valid values are \code{"analytic"} (default), 
                                 \code{"numDeriv"} (to use the numeric routine in package numDeric), \code{"maxLik"} 
                                 (to use the numeric routine in packahe maxLik), and \code{"none"} to avoid 
                                 estimating the Hessian and the covariance matrix. Only used if \code{apollo_control$HB=FALSE}.
  \item \strong{printLevel}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. 
                             Ignored if apollo_control$HB is TRUE. Default is 3.
  \item \strong{constraints}: Character vector. Constraints on parameters to estimate. For example \code{c('b1>0', 'b1 + b2>1')}.
                              Only \code{>}, \code{>=} and \code{=} can be used. And they cannot be mixed (e.g. 
                              \code{c(b1=b2, b2>0)} will fail). All parameter names must be on the left side. Fixed 
                              parameters cannot go into constraints. Alternatively, constraints can be defined as in 
                              \link[maxLik]{maxLik}. Constraints can only be used with maximum likelihood estimation, 
                              and the BFGS routine in particular.
  \item \strong{scaling}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. 
                          If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but 
                          \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter
                          b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, 
                          specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
 \item \strong{maxLik_settings}: List. Additional settings for maxLik. See argument \code{control} in \link[maxLik]{maxBFGS}, 
                                 \link[maxLik]{maxBHHH} and \link[maxLik]{maxNM} for more details. Only used for maximum 
                                 likelihood estimation.
 \item \strong{numDeriv_settings}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. 
                                   See argument \code{method.args} in \link[numDeriv]{grad} for more details.
 \item \strong{bootstrapSE}: Numeric. Number of bootstrap samples to calculate standard errors. Default is 0, meaning
                             no bootstrap s.e. will be calculated. Number must zero or a positive integer. Only used
                             if \code{apollo_control$HB} is \code{FALSE}.
 \item \strong{bootstrapSeed}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples.
                               Only used if \code{bootstrapSE>0}. Default is 24.
 \item \strong{silent}: Logical. If TRUE, no information is printed to the console during estimation. Default is FALSE.
 \item \strong{scaleHessian}: Logical. If TRUE, parameters are scaled to 1 for Hessian estimation. Default is TRUE.
}}
}
\value{
model object
}
\description{
Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
}
\details{
This is the main function of the Apollo package. The estimation process begins by running a number of checks on the 
\code{apollo_probabilities} function provided by the user.
If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
iterations. If bayesian estimation is used, estimation will finish once the predefined number of iterations are completed.
By default, this functions writes the estimated parameter values in each iteration to a file in the working directory. Writing 
can be turned off by setting \code{estimate_settings$writeIter} to \code{FALSE}, of by using any estimation algorithm
other than BFGS.
By default, \strong{final results are not written into a file nor printed into the console}, so users must make sure 
to call function \code{apollo_modelOutput} and/or \code{apollo_saveOutput} afterwards.
Users are strongly encouraged to visit www.apolloChoiceModelling.com to download examples on how to use the Apollo package.
The webpage also provides a detailed manual for the package, as well as a user-group to get further help.
}
