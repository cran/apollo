#' Calculates product of panel observations.
#' 
#' Multiplies likelihood of observations from the same individual, or adds the log of them.
#' 
#' This function should be called inside apollo_probabilities only if the data has a panel structure.
#' It should be called after apollo_avgIntraDraws if intra-individual draws are used.
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \describe{
#'                        \item{"estimate"}{For model estimation, returns probabilities of chosen alternatives.}
#'                        \item{"prediction"}{For model predictions, returns probabilities of all alternatives.}
#'                        \item{"validate"}{Validates input.}
#'                        \item{"zero_LL"}{Return probabilities with all parameters at zero.}
#'                        \item{"conditionals"}{For conditionals, returns probabilities of chosen alternatives.}
#'                        \item{"output"}{Checks that the model is well defined.}
#'                        \item{"raw"}{For debugging, returns probabilities of all alternatives}
#'                      }
#' @return Probabilities at the individual level.
#' @export
apollo_panelProd <- function(P, apollo_inputs, functionality){
  apollo_control = apollo_inputs[["apollo_control"]]
  workInLogs     = apollo_control$workInLogs
  
  if(apollo_control$HB==TRUE) return(P)
  
  if(!apollo_control$panelData) stop('Panel data setting not used, so multiplying over choices not applicable!')
  
  if(functionality %in% c("prediction","raw","validate")) return(P)
  
  inputIsList <- is.list(P)
  indivID <- apollo_inputs$database[, apollo_control$indivID]
  
  if(functionality=="zero_LL"){
    if(inputIsList) P <- P[["model"]]
    Pout = rowsum(log(P), group=indivID)
    if(!workInLogs) Pout <- exp(Pout)
    if(is.matrix(Pout) && ncol(Pout)==1) Pout=as.vector(Pout)
    if(inputIsList) Pout <- list(model=Pout)
    return(Pout)
  }
  
  if(functionality %in% c("estimate", "conditionals")){
    if(inputIsList) P <- P[["model"]]
    if(is.array(P) && length(dim(P))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
    if(is.vector(P) || (is.matrix(P) && !workInLogs) ){
      Pout <- rowsum(log(P), group=indivID)
      if(!workInLogs) Pout <- exp(Pout)
    }
    if(apollo_control$panelData && is.matrix(P) && workInLogs){
      B    <- rowsum(log(P), group=indivID) 
      Bbar <- rowMeans(B) 
      Pout <- Bbar + log( rowMeans(exp(B-Bbar)) ) 
    }
    if(inputIsList) Pout <- list(model=Pout)
    return(Pout)
  }
  
  if(functionality=="output"){
    if(!inputIsList) P <- list(model=P)
    j=1
    Pout=P
    while(j<= length(P)){
      if(is.array(P[[j]]) && length(dim(P[[j]]))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
      if(is.vector(P[[j]]) || (is.matrix(P[[j]]) && !workInLogs) ){
        Pout[[j]] <- rowsum(log(P[[j]]), group=indivID)
        if(!workInLogs) Pout[[j]] <- exp(Pout[[j]])
      }
      if(apollo_control$panelData && is.matrix(P[[j]]) && workInLogs){
        B    <- rowsum(log(P[[j]]), group=indivID) 
        Bbar <- rowMeans(B) 
        Pout[[j]] <- Bbar + log( rowMeans(exp(B-Bbar)) ) 
      }
      j=j+1
    }
    if(!inputIsList) Pout <- Pout[[1]]
    return(Pout)
  }
  
  
}