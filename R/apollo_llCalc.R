#' Calculates log-likelihood of all model components
#' 
#' Calculates the log-likelihood of each model component as well as the whole model.
#' 
#' This function calls apollo_probabilities with functionality="output". Then, it reorders the list of likelihoods so that "model" goes first.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If TRUE, no information is printed to the console by the function. Default is FALSE.
#' @return A list of vectors. Each vector corresponds to the log-likelihood of the whole model (first element) or a model component.
#' @export
apollo_llCalc <- function(apollo_beta, apollo_probabilities, apollo_inputs, silent=FALSE){
  apollo_fixed=c()
  
  # Use 'silent' from apollo_inputs if available
  if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent
  
  # Compare apollo_inputs values to those on the global environment
  apollo_compareInputs(apollo_inputs)
  
  apollo_randCoeff  = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars     = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  workInLogs <- apollo_inputs$apollo_control$workInLogs
  
  if(!silent) cat("Calculating LL of each model component...")
  Pout <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="output"),
                   error=function(e) return(NA))
  
  if(!anyNA(Pout) && is.list(Pout)){
    # Give name to unnamed components
    origNames <- names(Pout)
    newNames  <- paste0("component_", 1:length(Pout))
    if(!is.null(origNames)) newNames <- ifelse(origNames!="", origNames, newNames)
    names(Pout) <- newNames
    # Get log of likelihood with "model" first
    tmp <- c("model", newNames[newNames!="model"])
    if(!workInLogs) LLout <- lapply(Pout[tmp], log) else LLout <- Pout[tmp]
    LLout <-lapply(LLout,sum)
    if(!silent) cat("Done.\n")
    return(LLout)
  }
  
  if(!silent) cat("Not applicable to all components.\n")
  return(NA)
}