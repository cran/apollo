#' Calculates product of panel observations.
#'
#' Multiplies likelihood of observations from the same individual, or adds the log of them.
#'
#' This function should be called inside apollo_probabilities only if the data has a panel structure.
#' It should be called after apollo_avgIntraDraws if intra-individual draws are used.
#'
#' @param P List of probabilities. It must contain one element called "model" containing the full likelihood of the model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    See \code{?apollo_validatecontrol} for details.
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
#' @param work_in_logs Boolean. TRUE for higher numeric stability at the expense of computational time. Useful for panel models only. Default is FALSE.
#' @param indivID Numeric vector. Vector with individual's ID. As long as the number of observations.
#' @return Probabilities at the individual level.
apollo_panelProd <- function(P, apollo_control, functionality, work_in_logs, indivID){

  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #

  if(apollo_control$HB==TRUE) return(P)

  # ############################### #
  #### pre-checks                ####
  # ############################### #

  if(!apollo_control$panelData) stop('Panel data setting not used, so multiplying over choices not applicable!')

  # ############################################# #
  #### functionality="prediction/validate/raw" ####
  # ############################################# #

  if(functionality %in% c("prediction","raw","validate")) return(P)

  # ########################################################## #
  #### functionality="estimate/conditionals/zero_LL/output" ####
  # ########################################################## #

  inputIsList <- is.list(P)

  if(functionality=="zero_LL"){
    if(inputIsList) P <- P[["model"]]
    Pout = rowsum(log(P), group=indivID)
    if(!work_in_logs) Pout <- exp(Pout)
    if(is.matrix(Pout) && ncol(Pout)==1) Pout=as.vector(Pout)
    if(inputIsList) Pout <- list(model=Pout)
    return(Pout)
  }

  if(functionality %in% c("estimate", "conditionals")){
    if(inputIsList) P <- P[["model"]]
    if(is.array(P) && length(dim(P))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
    if(is.vector(P) || (is.matrix(P) && !work_in_logs) ){
      Pout <- rowsum(log(P), group=indivID)
      if(!work_in_logs) Pout <- exp(Pout)
    }
    if(apollo_control$panelData && is.matrix(P) && work_in_logs){
      # approach to use if working in logs with mixing
      B    <- rowsum(log(P), group=indivID) # nIndiv x nDraws
      Bbar <- rowMeans(B) # nIndiv x 1
      Pout <- Bbar + log( rowMeans(exp(B-Bbar)) ) # nIndiv x 1
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
      if(is.vector(P[[j]]) || (is.matrix(P[[j]]) && !work_in_logs) ){
        Pout[[j]] <- rowsum(log(P[[j]]), group=indivID)
        if(!work_in_logs) Pout[[j]] <- exp(Pout[[j]])
      }
      if(apollo_control$panelData && is.matrix(P[[j]]) && work_in_logs){
        # approach to use if working in logs with mixing
        B    <- rowsum(log(P[[j]]), group=indivID) # nIndiv x nDraws
        Bbar <- rowMeans(B) # nIndiv x 1
        Pout[[j]] <- Bbar + log( rowMeans(exp(B-Bbar)) ) # nIndiv x 1
      }
      j=j+1
    }
    if(!inputIsList) Pout <- Pout[[1]]
    return(Pout)
  }


}
