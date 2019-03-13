#' Validates the apollo_HB list of parameters
#'
#' Validates the apollo_HB list of parameters and sets default values for the omitted ones.
#'
#' This function is only necessary when using bayesian estimation.
#' @param apollo_HB List. Contains options for bayesian estimation. See \code{?RSGHB::doHB} for details.
#'                   Parameters \code{modelname}, \code{gVarNamesFixed}, \code{gVarNamesNormal},
#'                   \code{gDIST}, \code{svN} and \code{FC} are automatically set based on the
#'                   other arguments of this function.
#'                   It should also include a named character vector called \code{hbDist} identifying 
#'                   the distribution of each parameter to be estimated. Possible values are as follows.
#'                   \itemize{
#'                     \item "DNE": Parameter kept at its starting value (not estimated).
#'                     \item "F": Fixed (as in non-random) parameter.
#'                     \item "N": Normal.
#'                     \item "LN+": Positive log-normal.
#'                     \item "LN-": Negative log-normal.
#'                     \item "CN+": Positive censored normal.
#'                     \item "CN-": Negative censored normal.
#'                     \item "JSB": Johnson SB.
#'                   }
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#'                    value is constant throughout estimation).
#' @param apollo_control List. Options controlling the running of the code. See \link{apollo_validateInputs}.
#' @return Validated apollo_HB
apollo_validateHBControl=function(apollo_HB, apollo_beta, apollo_fixed, apollo_control){
  
  if(!("hbDist" %in% names(apollo_HB))) stop("No 'hbDist' element in apollo_HB.")
  hbDist <- apollo_HB$hbDist
  if(length(apollo_beta)!=length(hbDist)) stop("Argument hbDist has different length than apollo_beta.")
  
  map <- c("F"=0, "N"=1, "LN+"=2, "LN-"=3, "CN+"=4, "CN-"=5, "JSB"=6)
  hbDist_est = hbDist[!(names(apollo_beta) %in% apollo_fixed)]
  hbDist_est = stats::setNames(map[hbDist_est], names(hbDist_est))
  theta_est = apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  
  

  gVarNamesFixed  = names(theta_est[hbDist_est==0])
  FC              = theta_est[hbDist_est==0]
  gVarNamesNormal = names(theta_est[hbDist_est>0])
  gDIST           = hbDist_est[hbDist_est>0]
  svN             = theta_est[hbDist_est>0]

  apollo_HB$modelname=apollo_control$modelName
  apollo_HB$gVarNamesFixed = gVarNamesFixed
  apollo_HB$gVarNamesNormal=gVarNamesNormal
  apollo_HB$gDIST=gDIST
  apollo_HB$svN=svN
  apollo_HB$FC=FC

  if(length(apollo_HB$gVarNamesNormal)==0) {
   apollo_HB$gVarNamesNormal=NULL
   apollo_HB$hIW=FALSE
  }
  if(length(apollo_HB$gVarNamesFixed)==0) apollo_HB$gVarNamesFixed=NULL
  if(length(apollo_HB$gDIST)==0) apollo_HB$gDIST=NULL
  if(length(apollo_HB$svN)==0) apollo_HB$svN=NULL
  if(length(apollo_HB$FC)==0) apollo_HB$FC=NULL

  cat("All checks on apollo_HB completed.\n")
  return(apollo_HB)
}
