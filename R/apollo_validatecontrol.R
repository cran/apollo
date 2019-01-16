#' Validates apollo_control
#'
#' Validates the options controlling the running of the code \code{apollo_control} and sets default values for the omitted ones.
#'
#' This function should be run before running \code{apollo_validatedata}.
#'
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    \describe{
#'                      \item{modelName}{Character. Name of the model. Used when saving the output to files.}
#'                      \item{modelDescr}{Character. Description of the model. Used in output files.}
#'                      \item{indivID}{Character. Name of column in the database with each decision maker's ID.}
#'                      \item{mixing}{Boolean. TRUE for models that include random parameters.}
#'                      \item{nCores}{Numeric>0. Number of cores to use in calculations of the model likelihood.}
#'                      \item{seed}{Numeric. Seed for random number generation.}
#'                      \item{HB}{Boolean. TRUE if using RSGHB for Bayesian estimation of model.}
#'                      \item{noValidation}{Boolean. TRUE if user does not wish model input to be validated before estimation - FALSE by default.}
#'                      \item{noDiagnostics}{Boolean. TRUE if user does not wish model diagnostics to be printed - FALSE by default.}
#'                    }
#' @return Validated version of apollo_control, with additional element called panelData set to TRUE for repeated choice data.
apollo_validatecontrol=function(database,apollo_control){

  if(is.null(apollo_control$modelName)){
    apollo_control$modelName <- "noname_model"
    cat("Model name missing, set to default of 'noname_model'\n")
  }

  if(is.null(apollo_control$modelDescr)){
    apollo_control$modelDescr <- "No model description given"
    cat("Model description missing, set to default of 'No model description given'\n")
  }

  if(is.null(apollo_control$indivID)){
    stop('Name of column with individual IDs not provided in apollo_control.')
  }

  if(is.null(apollo_control$mixing)){
    apollo_control$mixing <- FALSE
    cat("Missing setting for mixing, set to default of FALSE\n")
  }

  if(is.null(apollo_control$nCores)){
    apollo_control$nCores <- 1
    cat("Missing setting for nCores, set to default of 1\n")
  }

  if(is.null(apollo_control$seed)){
    apollo_control$seed <- 13
    cat("Missing setting for seed, set to default of 13\n")
  }

  if(is.null(apollo_control$HB)){
    apollo_control$HB <- FALSE
    cat("Missing setting for HB, set to default of FALSE\n")
  }

  if(apollo_control$HB==TRUE & apollo_control$nCores > 1){
    apollo_control$nCores <- 1
    cat("nCores set to 1 for Bayesian estimation\n")
  }

  if(is.null(apollo_control$noValidation)){
    apollo_control$noValidation <- FALSE
    cat("Missing setting for noValidation, set to default of FALSE\n")
  }

  if(is.null(apollo_control$noDiagnostics)){
    apollo_control$noDiagnostics <- FALSE
    cat("Missing setting for noDiagnostics, set to default of FALSE\n")
  }

  if(length(unique(database[,apollo_control$indivID]))<nrow(database)){
    cat("Missing setting for panelData.\n")
    apollo_control$panelData  = TRUE
	cat("Several observations per individual detected based on the value of ", apollo_control$indivID, ".\n  Setting panelData set to TRUE.\n", sep="")
  } else {
    apollo_control$panelData  = FALSE
	cat("Only one observation per individual detected based on the value of ", apollo_control$indivID, ".\n  Setting panelData set to FALSE.\n", sep="")
  }

  if(!(apollo_control$mixing%in%c(TRUE,FALSE))) stop("Setting for mixing should be TRUE or FALSE")
  if(!(apollo_control$HB%in%c(TRUE,FALSE))) stop("Setting for HB should be TRUE or FALSE")
  if(!(apollo_control$noValidation%in%c(TRUE,FALSE))) stop("Setting for noValidation should be TRUE or FALSE")
  if(!(apollo_control$noDiagnostics%in%c(TRUE,FALSE))) stop("Setting for noDiagnostics should be TRUE or FALSE")

  if(apollo_control$noValidation==TRUE) warning("With setting noValidation=TRUE, your model code will not be validated prior to estimation. This may of course be deliberate for large models or models with many components.")

  cat("All checks on apollo_control completed.\n")
  return(apollo_control)

}
