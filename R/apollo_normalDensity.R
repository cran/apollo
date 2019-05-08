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
#'                                 \item outcomeNormal: Numeric vector. Dependant variable.
#'                                 \item xNormal: Numeric vector. Single explanatory variable.
#'                                 \item mu: Numeric scalar. Intercept of the linear model.
#'                                 \item sigma: Numeric scalar. Variance of error component of linear model to be estimated.
#'                       \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
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
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen value for each observation.
#'           \item "prediction": Not applicable.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": Not applicable.
#'           \item "conditionals": Same as "estimate".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "estimate".
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
  if(is.null(normalDensity_settings[["outcomeNormal"]])) stop("The normalDensity_settings list needs to include an object called \"outcomeNormal\"!")
  if(is.null(normalDensity_settings[["xNormal"]])) stop("The normalDensity_settings list needs to include an object called \"xNormal\"!")
  if(is.null(normalDensity_settings[["mu"]])) stop("The normalDensity_settings list needs to include an object called \"mu\"!")
  if(is.null(normalDensity_settings[["sigma"]])) stop("The normalDensity_settings list needs to include an object called \"sigma\"!")
  if(is.null(normalDensity_settings[["rows"]])) normalDensity_settings[["rows"]]="all"

  outcomeNormal=normalDensity_settings[["outcomeNormal"]]
  xNormal=normalDensity_settings[["xNormal"]]
  mu=normalDensity_settings[["mu"]]
  sigma=normalDensity_settings[["sigma"]]
  rows=normalDensity_settings[["rows"]]

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){

    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeNormal))
    outcomeNormal[!rows] = mean(mu)

    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
                                error=function(e) return(list(noValidation=FALSE, noDiagnostics=FALSE)) )

    if(apollo_control$noValidation==FALSE){
      if(is.vector(xNormal)) xlength=length(xNormal)
      if(is.array(xNormal)) xlength=dim(xNormal)[1]
      if(!is.vector(outcomeNormal)) stop("Dependent variable for Normal density model needs to be one-dimensional!")
      if(xlength!=length(outcomeNormal)) stop("Incompatible dimensions for dependent and explanatory variables for Normal density model!")
      if(length(mu)!=1) stop("Need to use a scalar for mean parameter for Normal density model!")
      if(length(sigma)!=1) stop("Need to use a scalar for standard deviation parameter for Normal density model!")
      #cat("\nAll checks passed for Normal density model component.")
    }
    if(apollo_control$noDiagnostics==FALSE){
      #cat("\nSummary statistics for dependent variable for Normal density model component\n")
      #print(summary(outcomeNormal[rows]))
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Summary statistics for dependent variable for Normal density model component:", summary(outcomeNormal[rows]), apolloLog)
    }
    return(invisible(TRUE))
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
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeNormal))
    outcomeNormal[!rows] = mean(mu)
    ans <- stats::dnorm(outcomeNormal-xNormal,mu,sigma)
    ans[!rows] <- 1
    return(ans)
  }

  # ################################ #
  #### functionality="output"     ####
  # ################################ #

  if(functionality=="output"){
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeNormal))

    ## write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
    #apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
    #                            error=function(e){
    #                              cat("apollo_normalDensity could not retrieve apollo_control. No diagnostics in output.\n")
    #                              return(NA)
    #                            } )
    #if(!(length(apollo_control)==1 && is.na(apollo_control))){
    #  fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    #  fileName <- file.path(tempdir(),fileName)
    #  fileConn <- tryCatch( file(fileName, open="at"),
    #                        error=function(e){
    #                          cat('apollo_normalDensity could not write diagnostics to temporary file. No diagnostics in output.\n')
    #                          return(NA)
    #                        })
    #  if(!anyNA(fileConn)){
    #    sink(fileConn)
    #    on.exit({if(sink.number()>0) sink(); close(fileConn)})
    #    if(apollo_control$noDiagnostics==FALSE){
    #      cat("\nSummary statistics for dependent variable for Normal density model component\n")
    #      print(summary(outcomeNormal[rows]))}
    #  }
    #}
    
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Summary statistics for dependent variable for Normal density model component:", summary(outcomeNormal[rows]), apolloLog)
    
    ans <- apollo_normalDensity(normalDensity_settings, functionality="estimate")

    return(ans)
  }

}
