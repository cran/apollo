#' Validates the HB_control list of parameters
#'
#' Validates the HB_control list of parameters and sets default values for the omitted ones.
#'
#' This function is only necessary when using bayesian estimation.
#' @param HB_control List. Contains options for bayesian estimation. See \code{?RSGHB::doHB} for details.
#'                   Parameters \code{modelname}, \code{gVarNamesFixed}, \code{gVarNamesNormal},
#'                   \code{gDIST}, \code{svN} and \code{FC} are automatically set based on the
#'                   other arguments of this function.
#' @param hb_dist Numeric vector. Contains a number representing the distribution of each
#'                parameter to be estimated. Possible values are as follows.
#'                \describe{
#'                  \item{-1}{Fixed (as in non-random) parameter.}
#'                  \item{1}{Normal.}
#'                  \item{2}{Positive log-normal.}
#'                  \item{3}{Negative log-normal.}
#'                  \item{4}{Positive censored normal.}
#'                  \item{5}{Negative censored normal.}
#'                  \item{6}{Johnson SB.}
#'                }
#' @param theta_start Named numeric vector. Names and starting values of the parameters.
#' @param theta_fixed Character vector. Names of the fixed parameters (as in those whose
#'                    value is constant throughout estimation).
#' @return Validated HB_control
apollo_validateHBcontrol=function(HB_control, hb_dist, theta_start, theta_fixed){

  if(length(theta_start)!=length(hb_dist)) stop("Argument hb_dist has different length than theta_start.")

  theta_est   = theta_start[!(names(theta_start) %in% theta_fixed)]
  hb_dist_est = hb_dist[!(names(theta_start) %in% theta_fixed)]

  gVarNamesFixed  = names(theta_est[hb_dist_est==0])
  FC              = theta_est[hb_dist_est==0]
  gVarNamesNormal = names(theta_est[hb_dist_est>0])
  gDIST           = hb_dist_est[hb_dist_est>0]
  svN             = theta_est[hb_dist_est>0]

  apollo_control <- get("apollo_control", envir=parent.frame())
  HB_control$modelname=apollo_control$modelName
  HB_control$gVarNamesFixed = gVarNamesFixed
  HB_control$gVarNamesNormal=gVarNamesNormal
  HB_control$gDIST=gDIST
  HB_control$svN=svN
  HB_control$FC=FC

  if(length(HB_control$gVarNamesNormal)==0) {
   HB_control$gVarNamesNormal=NULL
   HB_control$hIW=FALSE
  }
  if(length(HB_control$gVarNamesFixed)==0) HB_control$gVarNamesFixed=NULL
  if(length(HB_control$gDIST)==0) HB_control$gDIST=NULL
  if(length(HB_control$svN)==0) HB_control$svN=NULL
  if(length(HB_control$FC)==0) HB_control$FC=NULL

  cat("All checks on HB_control completed.\n")
  return(HB_control)
}
