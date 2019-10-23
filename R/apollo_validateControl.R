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
apollo_validateControl=function(database,apollo_control, silent=FALSE){
  
  if(is.null(apollo_control$modelName)){
    apollo_control$modelName <- "noname_model"
    if(!silent) cat("Model name missing, set to default of 'noname_model'\n")
  }
  
  if(is.null(apollo_control$modelDescr)){
    apollo_control$modelDescr <- "No model description given"
    if(!silent)cat("Model description missing, set to default of 'No model description given'\n")
  }
  
  if(is.null(apollo_control$indivID)){
    stop('Name of column with individual IDs not provided in apollo_control.')
  }
  
  if(is.null(apollo_control$mixing)&is.null(apollo_control$HB)){
    apollo_control$mixing <- FALSE
    if(!silent) cat("Missing setting for mixing, set to default of FALSE\n")
  }
  
  if(is.null(apollo_control$nCores)){
    apollo_control$nCores <- 1
    if(!silent) cat("Missing setting for nCores, set to default of 1\n")
  }
  
  if(is.null(apollo_control$workInLogs)){
    apollo_control$workInLogs <- FALSE
    if(!silent) cat("Missing setting for workInLogs, set to default of FALSE\n")
  }
  
  if(is.null(apollo_control$seed)){
    apollo_control$seed <- 13
    if(!silent) cat("Missing setting for seed, set to default of 13\n")
  }
  
  if(is.null(apollo_control$HB)){
    apollo_control$HB <- FALSE
    if(!silent) cat("Missing setting for HB, set to default of FALSE\n")
  }
  
  if(apollo_control$HB & (!is.null(apollo_control$mixing) && apollo_control$mixing!=FALSE)){
    apollo_control$mixing <- FALSE
    if(!silent) cat("HB set to TRUE, so mixing set to FALSE\n")
  }
  
  if(apollo_control$HB) apollo_control$mixing <- FALSE
  
  if(apollo_control$HB==TRUE & apollo_control$nCores > 1){
    apollo_control$nCores <- 1
    if(!silent) cat("nCores set to 1 for Bayesian estimation\n")
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
    if(!silent) cat("Several observations per individual detected based on the value of ", apollo_control$indivID, ".\n  Setting panelData set to TRUE.\n", sep="")
  } else {
    apollo_control$panelData  = FALSE
    if(!silent) cat("Only one observation per individual detected based on the value of ", apollo_control$indivID, ".\n  Setting panelData set to FALSE.\n", sep="")
  }}
  
  if(!is.logical(apollo_control$mixing       )) stop("Setting for mixing should be TRUE or FALSE")
  if(!is.logical(apollo_control$HB           )) stop("Setting for HB should be TRUE or FALSE")
  if(!is.logical(apollo_control$noValidation )) stop("Setting for noValidation should be TRUE or FALSE")
  if(!is.logical(apollo_control$noDiagnostics)) stop("Setting for noDiagnostics should be TRUE or FALSE")
  if(!is.logical(apollo_control$workInLogs   )) stop("Setting for workInLogs should be TRUE or FALSE")
  
  if(!is.null(apollo_control$weights)){
    w <- apollo_control$weights
    if(length(w)!=1 || !is.character(w) || !(w %in% names(database))) stop("'apollo_control$weights' is not the name of a column in 'database'.")
  }
  
  if((apollo_control$noValidation==TRUE)&!silent) warning("With setting noValidation=TRUE, your model code will not be validated prior to estimation. This may of course be deliberate for large models or models with many components.")
  
  allVars <- c("modelName", "modelDescr", "indivID", "mixing", "nCores", "seed", "HB", "noValidation", "noDiagnostics", "weights", "workInLogs", "panelData")
  unknownVars <- names(apollo_control)[!( names(apollo_control) %in% allVars )]
  if(length(unknownVars)>0){
    warning(" Variable(s) {", paste(unknownVars, collapse=", "), "}\n where not recognised in apollo_control and will be ignored.\n Check ?apollo_control for a list of valid control variables.")
  }
  
  if(!silent) cat("All checks on apollo_control completed.\n")
  return(apollo_control)
  
}
