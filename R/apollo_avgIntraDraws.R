#' Averages intra-individual draws
#'
#' Averages the likelihood across intra-individual draws.
#'
#' This function performs additional checks on the shape of the probabilities given as arguments.
#' @param P List. Contains the probabilities for each model component, for each inter- and intra-person draw.
#' @param apollo_control List. Contains options for the estimation
#'                    See \link{apollo_validatecontrol} for details.
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
#' @return Average probability over intra-individual draws (shape depends on argument \code{functionality}).
apollo_avgIntraDraws <- function(P, apollo_control, functionality){

  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #

  if(apollo_control$HB==TRUE) return(P)

  # ############################### #
  #### pre-checks                ####
  # ############################### #

  if(!apollo_control$mixing) stop('No mixing used in model!')

  inputIsList <- is.list(P)

  if(inputIsList && functionality!="prediction" && is.null(P[["model"]])) stop('Element called "model" is missing in list P!')

  # ########################################## #
  #### functionality="zero_LL/raw/validate" ####
  # ########################################## #

  if(functionality %in% c("zero_LL","raw","validate")) return(P)

  # ########################################## #
  #### functionality="estimate/conditionals" ####
  # ########################################## #

  if(functionality=="estimate" | functionality=="conditionals"){
    if(!is.list(P)){
      if(is.array(P)){
        if(length(dim(P))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          return( colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3] )
        } else stop('No intra-individual draws present to average over!')
      } else stop('No draws present to average over!')
    } else {
      if(is.array(P[["model"]])){
        if(length(dim(P[["model"]]))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          return( list(model=colSums(aperm(P[["model"]], perm=c(3,1,2)))/dim(P[["model"]])[3] ))
        } else stop('No intra-individual draws present to average over!')
      } else stop('No draws present to average over!')
    }
  }

  # ########################################## #
  #### functionality="prediction"           ####
  # ########################################## #

  if(functionality=="prediction"){
    if(!is.list(P)){
      if(is.array(P)){
        if(length(dim(P))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          output=colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3]
        }}
      return(output)
    } else {
      output_list=P
      j=1
      while(j<= length(P)){
        k=1
        while(k<= length(P[[j]])){
          if(is.array(P[[j]][[k]])){
            if(length(dim(P[[j]][[k]]))==3){
              # Returns a matrix of dimensions nObs x nDrawsInter
              output_list[[j]][[k]]=colSums(aperm(P[[j]][[k]], perm=c(3,1,2)))/dim(P[[j]][[k]])[3]
            }}
          k=k+1}
        j=j+1}
      return(output_list)
    }

  }

  # ########################################## #
  #### functionality="output"               ####
  # ########################################## #

  if(functionality=="output"){
    if(!is.list(P)){
      if(is.array(P)){
        if(length(dim(P))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          output=colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3]
        }}
      return(output)
    } else {
      output_list=P
      j=1
      while(j<= length(P)){
        if(is.array(P[[j]])){
          if(length(dim(P[[j]]))==3){
            # Returns a matrix of dimensions nObs x nDrawsInter
            output_list[[j]]=colSums(aperm(P[[j]], perm=c(3,1,2)))/dim(P[[j]])[3]
          }}
        j=j+1}
      return(output_list)
    }
  }

}
