#' Predicts using an estimated model
#' 
#' Calculates apollo_probabilities with functionality="prediction" and extracts one element from the returned list.
#' 
#' Structure of predictions are simplified before returning, e.g. list of vectors are turned into a matrix.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param modelComponent Character. Name of component of apollo_probabilities output to calculate predictions for. Default is "model", i.e. the whole model.
#' @return A vector containing predictions for component \code{modelComponent} of the model described in \code{apollo_probabilities}.
#' @export
apollo_prediction <- function(model, apollo_probabilities, apollo_inputs, modelComponent="model"){
  cat("Updating inputs...")
  apollo_inputs <- apollo_validateInputs(silent=TRUE)
  cat("Done.\n")
  cat("Running predictions from model...")
  apollo_beta=model$estimate
  apollo_fixed=model$apollo_fixed
  
  predictions <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="prediction")
  if(is.null(predictions[[modelComponent]])) stop("A likelihood element with the selected name does not exist in your model!")
  predictions=predictions[[modelComponent]]
  if(anyNA(predictions)){
    cat("\nPredictions do not exist for the selected model component.")
    return(NA)
  }
  if(is.list(predictions)){
    predictions_out           = matrix(unlist(predictions), ncol=length(predictions), byrow=FALSE)
    colnames(predictions_out) = names(predictions)
    cat(" Done.\n")
    return(predictions_out)
  }
  cat(" Done.\n")
  return(predictions)
}