#' Validates apollo_control
#'
#' Validates the options controlling the running of the code \code{apollo_control} and sets default values for the omitted ones.
#'
#' This function should be run before running \code{apollo_validateData}.
#'
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    \itemize{
#'                      \item \code{modelName}: Character. Name of the model. Used when saving the output to files.
#'                      \item \code{modelDescr}: Character. Description of the model. Used in output files.
#'                      \item \code{indivID}: Character. Name of column in the database with each decision maker's ID.
#'                      \item \code{mixing}: Boolean. TRUE for models that include random parameters.
#'                      \item \code{nCores}: Numeric>0. Number of cores to use in calculations of the model likelihood.
#'                      \item \code{seed}: Numeric. Seed for random number generation.
#'                      \item \code{HB}: Boolean. TRUE if using RSGHB for Bayesian estimation of model.
#'                      \item \code{noValidation}: Boolean. TRUE if user does not wish model input to be validated before estimation - FALSE by default.
#'                      \item \code{noDiagnostics}: Boolean. TRUE if user does not wish model diagnostics to be printed - FALSE by default.
#'                      \item \code{weights}: Character. Name of column in database containing weights for estimation.
#'                      \item \code{workInLogs}: Boolean. TRUE for increased numeric precision in models with panel data - FALSE by default.
#'                      \item \code{panelData}: Boolean. TRUE if there are multiple obsrvations (i.e. rows) for each decision maker - Automatically set based on \code{indivID} by default.
#'                    }
#' @param silent Boolean. If TRUE, no messages are printed to screen.
#' @return Validated version of apollo_control, with additional element called panelData set to TRUE for repeated choice data.
#' @export
apollo_validateControl=function(database,apollo_control, silent=FALSE){
  
  if(is.null(apollo_control$debug)) apollo_control$debug <- FALSE
  if(!is.logical(apollo_control$debug) || length(apollo_control$debug)!=1) stop('Setting "debug" in apollo_control should have a logical (boolean) value.')
  debug <- apollo_control$debug
  
  if(is.null(apollo_control$modelName)){
    apollo_control$modelName <- paste("model_", gsub("[: -]", "" , Sys.time(), perl=TRUE), sep="")
    if(!silent) apollo_print(paste0('Model name missing in apollo_control, set to default of "', apollo_control$modelName, '".'))
  }
  
  if(is.null(apollo_control$modelDescr)) apollo_control$modelDescr <- 'No model description provided in apollo_control'
  
  if(is.null(apollo_control$indivID)) stop('Name of column with individual IDs not provided in apollo_control.')
  
  if(is.null(apollo_control$nCores)){
    apollo_control$nCores <- 1
    if(debug) apollo_print("Missing setting for nCores in apollo_control, set to default of 1")
  }
  
  if(is.null(apollo_control$workInLogs)){
    apollo_control$workInLogs <- FALSE
    if(debug) apollo_print("Missing setting for workInLogs in apollo_control, set to default of FALSE")
  }
  
  if(is.null(apollo_control$seed)){
    apollo_control$seed <- 13
    if(debug) apollo_print("Missing setting for seed in apollo_control, set to default of 13")
  }
  
  if(is.null(apollo_control$mixing)&is.null(apollo_control$HB)){
    apollo_control$mixing <- FALSE
    if(debug) apollo_print("Missing setting for mixing in apollo_control, set to default of FALSE")
  }
  
  if(is.null(apollo_control$HB)){
    apollo_control$HB <- FALSE
    if(debug) apollo_print("Missing setting for HB in apollo_control, set to default of FALSE")
  }
  
  if(apollo_control$HB & (!is.null(apollo_control$mixing) && apollo_control$mixing!=FALSE)){
    apollo_control$mixing <- FALSE
    if(!silent) apollo_print("HB set to TRUE in apollo_control, so mixing set to FALSE")
  }
  
  if(apollo_control$HB) apollo_control$mixing <- FALSE
  
  if(apollo_control$HB==TRUE & apollo_control$nCores > 1){
    apollo_control$nCores <- 1
    if(!silent) apollo_print("nCores set to 1 in apollo_control for Bayesian estimation")
  }
  
  if(is.null(apollo_control$noValidation)){
    apollo_control$noValidation <- FALSE
  }
  
  if(is.null(apollo_control$noDiagnostics)){
    apollo_control$noDiagnostics <- FALSE
  }  
  
  if(is.null(apollo_control$panelData)){
  if(length(unique(database[,apollo_control$indivID]))<nrow(database)){
    apollo_control$panelData  = TRUE
    if(!silent) apollo_print(paste0("Several observations per individual detected based on the value of ", apollo_control$indivID, ". Setting panelData in apollo_control set to TRUE.", sep=""))
  } else {
    apollo_control$panelData  = FALSE
    if(debug) apollo_print(paste0("Only one observation per individual detected based on the value of ", apollo_control$indivID, ". Setting panelData in apollo_control set to FALSE.", sep=""))
  }}
  
  # Check that workInLogs is only used with panelData
  if(apollo_control$workInLogs & !apollo_control$panelData){
    apollo_control$workInLogs <- FALSE
    if(!silent) apollo_print("Working in logs is only applicable with panel data; workInLogs set to FALSE in apollo_control.")
  }
  
  if(!is.logical(apollo_control$mixing       )) stop("Setting for mixing in apollo_control should be TRUE or FALSE")
  if(!is.logical(apollo_control$HB           )) stop("Setting for HB in apollo_control should be TRUE or FALSE")
  if(!is.logical(apollo_control$noValidation )) stop("Setting for noValidation in apollo_control should be TRUE or FALSE")
  if(!is.logical(apollo_control$noDiagnostics)) stop("Setting for noDiagnostics in apollo_control should be TRUE or FALSE")
  if(!is.logical(apollo_control$workInLogs   )) stop("Setting for workInLogs in apollo_control should be TRUE or FALSE")

  if(!is.null(apollo_control$weights)){
    w <- apollo_control$weights
    if(length(w)!=1 || !is.character(w) || !(w %in% names(database))) stop("'apollo_control$weights' is not the name of a column in 'database'.")
  }
  
  if((apollo_control$noValidation==TRUE)&!silent) apollo_print("With setting noValidation=TRUE in apollo_control, your model code will not be validated prior to estimation. This may of course be deliberate for large models or models with many components.")
  
  if(is.null(apollo_control$cpp)){
    apollo_control$cpp <- FALSE
    if(debug) apollo_print("Missing setting for cpp in apollo_control, set to default of FALSE.")
  } else {
    if(!is.logical(apollo_control$cpp) || !length(apollo_control$cpp)==1) stop("Setting cpp in apollo_control should be a single logical value!")
  }
  
  if(is.null(apollo_control$analyticGrad)){
    apollo_control$analyticGrad <- TRUE
    if(debug) apollo_print("Missing setting for analyticGrad in apollo_control, set to default of TRUE")
  } else {
    if(!is.logical(apollo_control$analyticGrad) || !length(apollo_control$analyticGrad)==1) stop("Setting analyticGrad in apollo_control should be a single logical value")
  }
  
  if(apollo_control$HB){
    apollo_control$analyticGrad <- FALSE
    if(debug) apollo_print("Analytic gradients cannot be used with setting HB")
  } 
  
  if(apollo_control$workInLogs){
    apollo_control$analyticGrad <- FALSE
    if(debug) apollo_print("Analytic gradients cannot be used with setting workInLogs")
  } 
  
  if(is.null(apollo_control$matrixMult)){
    apollo_control$matrixMult <- FALSE
    if(debug) apollo_print("Missing setting for matrixMult in apollo_control, set to default of FALSE")
  } else {
    if(!is.logical(apollo_control$matrixMult) || !length(apollo_control$matrixMult)==1) stop("Setting matrixMult in apollo_control should be a single logical value")
  }
  if(apollo_control$matrixMult & apollo_control$mixing) stop("Setting matrixMult in apollo_control is only valid for models without mixing!")
  
  
  
  allVars <- c("modelName", "modelDescr", "indivID", "mixing", "nCores", "seed", "HB", 
               "noValidation", "noDiagnostics", "weights", "workInLogs", "panelData", 
               "cpp","subMaxV", "analyticGrad", "matrixMult", "debug")
  unknownVars <- names(apollo_control)[!( names(apollo_control) %in% allVars )]
  if(length(unknownVars)>0){
    apollo_print(paste0("Variable(s) {", paste(unknownVars, collapse=", "), "} were not recognised in apollo_control and will be ignored. Check ?apollo_control for a list of valid control variables."))
  }
  
  if(is.null(apollo_control$subMaxV)) apollo_control$subMaxV <- TRUE
  
  if(!silent) apollo_print("All checks on apollo_control completed.\n")
  return(apollo_control)
  
}
