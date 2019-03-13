#' Checks likelihood
#' 
#' Checks that likelihoods, i.e. probabilities in the case of choice models, are in the appropiate format to be returned.
#' 
#' This function should be called inside \code{apollo_probabilities}, near the end of it, just before \code{return(P)}.
#' This function only performs checks on the shape of P, but does not change its values in any way.
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
#' @return The likelihood (i.e. probability in the case of choice models) of the model in the appropriate form for the given functionality.
#' @export
apollo_prepareProb=function(P, apollo_inputs, functionality){
  panelData <- apollo_inputs$apollo_control$panelData
  nIndiv <- length(unique(apollo_inputs$database[, apollo_inputs$apollo_control$indivID]))
  HB <- apollo_inputs$apollo_control$HB
  
  if(functionality %in% c("prediction", "validate")) return(P)
  
  if(is.null(P[["model"]])) stop('Element called model is missing in list P!')
  
  if(HB) return(P[["model"]])
  
  if(functionality=="estimate"){
    
    if(is.array(P[["model"]])) nPRows <- dim(P[["model"]])[1] else nPRows <- length(P[["model"]])
    if(panelData && nPRows>nIndiv) stop("Need to multiply observations for the same individual! (see ?apollo_panelProd)")
    
    if(is.array(P[["model"]])){
      if(length(dim(P[["model"]]))==3) stop('Need to average over intra-individual draws! (see ?apollo_avgIntraDraws)')
      if(dim(P[["model"]])[2]>1) stop('Need to average over inter-individual draws! (see ?apollo_avgInterDraws)')
      if(dim(P[["model"]])[2]==1) P[["model"]]=as.vector(P[["model"]])
    }
    P_out=P[["model"]]
  } 
  
  if(functionality=="output") {
    P_out=P
    origNames <- names(P_out)
    newNames  <- paste0("component_", 1:length(P_out))
    if(!is.null(origNames)) newNames <- ifelse(origNames!="", origNames, newNames)
    names(P_out) <- newNames
  }
  if(functionality=="raw") P_out=P
  if(functionality=="conditionals") P_out=P[["model"]]
  if(functionality=="zero_LL") P_out=P[["model"]]
  
  
  
  return(P_out)
}