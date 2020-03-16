#' Writes the type of a model inside \code{apolloLog}
#'
#' Writes the type of the model in a list inside apolloLog, which lives inside \code{apollo_inputs}.
#'
#' The variable \code{apolloLog} is an environment created inside \code{apollo_inputs}
#' by \code{apollo_validateInputs}, but re-set by apollo_estimate.
#' As an environment, it can be modified in place, i.e. all changes done within this 
#' function are recorded in apolloLog, even if it belongs to another environment.
#' @param modelType Character. Type of the model, e.g. "MNL" or "MDCEV"
#' @param apolloLog Environment. It contains the character vectors of titles and content.
#' @return TRUE if writing was succesful, FALSE if not.
apollo_reportModelTypeLog <- function(modelType, apolloLog){
  ### Validation
  if(!is.environment(apolloLog)) return(FALSE)
  if(!exists("modelTypeList", envir=apolloLog)) apolloLog$modelTypeList <- list()
  
  ### Update value
  apolloLog$modelTypeList <- c(apolloLog$modelTypeList, modelType)
  return(TRUE)
}