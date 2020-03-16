#' Calculates density from a Normal distribution
#'
#' Calculates density from a Normal distribution at a specific value with a specified mean and standard deviation.
#'
#' This function estimates the linear model outcomeNormal = mu + xNormal + epsilon, where epsilon is a random error distributed Normal(0,sigma).
#' If using this function in the context of an Integrated Choice and Latent Variable (ICLV) model with continuous
#' indicators, then \code{outcomeNormal} would be the value of the indicator, \code{xNormal} would be the value of the latent variable (possibly
#' multiplied by a parameter to measure its correlation with the indicator, e.g. xNormal=lambda*LV), and \code{mu} would be
#' an additional parameter to be estimated (the mean of the indicator, which should be fixed to zero if the indicator is
#' centered around its mean beforehand).
#' @param normalDensity_settings List of arguments to the functions. It must contain the following.
#'                               \itemize{
#'                                 \item \code{outcomeNormal}: Numeric vector. Dependant variable.
#'                                 \item \code{xNormal}: Numeric vector. Single explanatory variable.
#'                                 \item \code{mu}: Numeric scalar. Intercept of the linear model.
#'                                 \item \code{sigma}: Numeric scalar. Variance of error component of linear model to be estimated.
#'                                 \item \code{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                 \item \code{componentName}: Character. Name given to model component.
#'                               }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate": Used for model estimation.
#'                        \item "prediction": Used for model predictions.
#'                        \item "validate": Used for validating input.
#'                        \item "zero_LL": Used for calculating null likelihood.
#'                        \item "conditionals": Used for calculating conditionals.
#'                        \item "output": Used for preparing output after model estimation.
#'                        \item "raw": Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the likelihood for each observation.
#'           \item \strong{\code{"prediction"}}: Not implemented. Returns NA.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns NA.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'         }
#' @examples
#' ### Load data
#' xNormal <- runif(100)
#' outcomeNormal <- 1 + 2*xNormal + rnorm(100, mean=0, sd=0.5)
#'
#' ### Parameters
#' b = list(a=1, m=2)
#'
#' ### normalDensity settings
#' normalDensity_settings <- list(
#'    outcomeNormal     = outcomeNormal,
#'    xNormal     = 2*xNormal,
#'    mu    = 1,
#'    sigma = 0.5
#' )
#'
#' ### Compute choice probabilities using MNL model
#' apollo_normalDensity(normalDensity_settings, functionality="estimate")
#' @export
apollo_normalDensity <- function(normalDensity_settings, functionality){
  if(is.null(normalDensity_settings[["componentName"]])) normalDensity_settings[["componentName"]]="NormD"
  componentName = normalDensity_settings[["componentName"]]
  
  if(is.null(normalDensity_settings[["outcomeNormal"]])) stop("The normalDensity_settings list for model component \"",componentName,"\" needs to include an object called \"outcomeNormal\"!")
  if(is.null(normalDensity_settings[["xNormal"]])) stop("The normalDensity_settings list for model component \"",componentName,"\" needs to include an object called \"xNormal\"!")
  if(is.null(normalDensity_settings[["mu"]])) stop("The normalDensity_settings list for model component \"",componentName,"\" needs to include an object called \"mu\"!")
  if(is.null(normalDensity_settings[["sigma"]])) stop("The normalDensity_settings list for model component \"",componentName,"\" needs to include an object called \"sigma\"!")
  if(is.null(normalDensity_settings[["rows"]])) normalDensity_settings[["rows"]]="all"

  outcomeNormal = normalDensity_settings[["outcomeNormal"]]
  xNormal       = normalDensity_settings[["xNormal"]]
  mu            = normalDensity_settings[["mu"]]
  sigma         = normalDensity_settings[["sigma"]]
  rows          = normalDensity_settings[["rows"]]
  
  ### Drop excluded rows
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, length(outcomeNormal))
  if(any(!rows)){
    outcomeNormal <- outcomeNormal[rows]
    xNormal <- apollo_keepRows(xNormal, rows)
    mu <- apollo_keepRows(mu, rows)
    sigma <- apollo_keepRows(sigma, rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){

    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
                                error=function(e) return(list(noValidation=FALSE, noDiagnostics=FALSE)) )

    if(apollo_control$noValidation==FALSE){
      if(is.vector(xNormal)) xlength=length(xNormal)
      if(is.array(xNormal)) xlength=dim(xNormal)[1]
      if(!is.vector(outcomeNormal)) stop("Dependent variable for model component \"",componentName,"\" needs to be one-dimensional!")
      if(xlength!=1 && xlength!=length(outcomeNormal)) stop("Incompatible dimensions for dependent and explanatory variables for model component \"",componentName,"\"!")
      if(length(mu)!=1) stop("Need to use a scalar for mean parameter for model component \"",componentName,"\"!")
      if(length(sigma)!=1) stop("Need to use a scalar for standard deviation parameter for model component \"",componentName,"\"!")
    }
    if(apollo_control$noDiagnostics==FALSE){
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog(title   = paste0("Summary statistics for dependent variable for model component \"",componentName,"\":"), 
                    content = summary(outcomeNormal), apolloLog)
    }
    testL=apollo_normalDensity(normalDensity_settings, functionality="estimate")
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"",componentName,"\"")
    if(any(testL==0)) cat("\nSome observations have zero probability at starting value for model component \"",componentName,"\"")
    return(invisible(testL))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    return(NA)
  }

  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction"){
    return(NA)
  }

  # ############################################### #
  #### functionality="estimate/conditionals/raw" ####
  # ############################################### #

  if(functionality %in% c("estimate", "conditionals", "raw")){
    ans <- stats::dnorm(outcomeNormal-xNormal,mu,sigma)
    if(any(!rows)) ans <- apollo_insertRows(ans, rows, 1)
    return(ans)
  }

  # ################################ #
  #### functionality="output"     ####
  # ################################ #

  if(functionality=="output"){
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog(title   = paste0("Summary statistics for dependent variable for model component \"",componentName,"\":"), 
                  content = summary(outcomeNormal), apolloLog)
    apollo_reportModelTypeLog(modelType="NormD", apolloLog)
    ans <- apollo_normalDensity(normalDensity_settings, functionality="estimate")

    return(ans)
  }

}
