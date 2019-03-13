#' Applies weights
#' 
#' Applies weights to individual observations in likelihood function.
#' 
#' This function should be called inside \code{apollo_probabilities}, near the end of it, just before \code{return(P)}.
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output of \code{apollo_probabilities}.
#'                      \describe{
#'                        \item{"estimate"}{For model estimation, returns probabilities of chosen alternatives.}
#'                        \item{"prediction"}{For model predictions, returns probabilities of all alternatives.}
#'                        \item{"validate"}{Validates input.}
#'                        \item{"zero_LL"}{Return probabilities with all parameters at zero.}
#'                        \item{"conditionals"}{For conditionals, returns probabilities of chosen alternatives.}
#'                        \item{"output"}{Checks that the model is well defined.}
#'                        \item{"raw"}{For debugging, returns probabilities of all alternatives}
#'                      }
#' @return The likelihood (i.e. probability in the case of choice models) of the model in the appropriate form for the given functionality, 
#' multiplied by individual-specific weights.
#' @export
apollo_weighting=function(P, apollo_inputs, functionality){
  apollo_control=apollo_inputs[["apollo_control"]]
  panelData <- apollo_inputs$apollo_control$panelData
  nIndiv <- length(unique(apollo_inputs$database[, apollo_inputs$apollo_control$indivID]))
  HB <- apollo_inputs$apollo_control$HB
  
  if(!(functionality %in% c("estimate"))) return(P)
  
  if(is.null(P[["model"]])) stop('Element called model is missing in list P!')
  
  if(is.array(P[["model"]])) nPRows <- dim(P[["model"]])[1] else nPRows <- length(P[["model"]])
  if(panelData && nPRows>nIndiv) stop("Need to multiply observations from the same individual! (see ?apollo_panelProd)")
    
    if(is.array(P[["model"]])){
      if(length(dim(P[["model"]]))==3) stop('Need to average over intra-individual draws! (see ?apollo_avgIntraDraws)')
      if(dim(P[["model"]])[2]>1) stop('Need to average over inter-individual draws! (see ?apollo_avgInterDraws)')
      if(dim(P[["model"]])[2]==1) P[["model"]]=as.vector(P[["model"]])
    }

    P_out=P[["model"]]
   
  if(is.null(apollo_inputs$apollo_control$weights)) stop('Call to apollo_weighting performed without a weights variable defined in apollo_control!')

  weights=apollo_inputs$database[,apollo_inputs$apollo_control$weights]
  if(length(weights)==length(P_out)){ 
    weights=weights
  } else {
    indivID <- get(apollo_control$indivID) 
    nObsPerIndiv <- as.vector(table(indivID))
    weights=rowsum(weights, group=indivID)/nObsPerIndiv
  }
  
  if(apollo_inputs$apollo_control$workInLogs==TRUE){
    P_out=weights*P_out
  } else {
    P_out=P_out^weights
  }
  
  return(list("model"=P_out))
}