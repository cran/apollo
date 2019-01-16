#' Averages inter-individual draws
#'
#' Averages the likelihood across inter-individual draws.
#'
#' This function performs additional checks on the shape of the probabilities given as arguments.
#' The shape of \code{P} should be consistent with the value of \code{functionality}, but compliance is
#' assured if \code{P} comes from \link{apollo_mnl}.
#' For \code{functionality} values of "estimate", "zero_LL, "conditionals" and "output", \code{P} should be a list of matrices,
#' with each matrix containing the likelihood of the chosen alternative.
#' For \code{functionality} values of "prediction" and "raw", \code{P} should be a list, with each of
#' its elements a list of matrices, and each matrix containing the likelihood of a different alternative.
#' For the \code{functionality} value of "validate", \code{P} can be anything, as it is not processed.
#' @param P List. Contains the probabilities for each model component, for each inter-person draw.
#'          Intra-person draws must have been averaged over on \code{P}, before giving it as an argument to this function.
#' @param apollo_control List. Contains options for the estimation
#'                    See \link{apollo_validatecontrol} for details.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \describe{
#'                        \item{"estimate"}{Used for model estimation.}
#'                        \item{"prediction"}{Used for model predictions.}
#'                        \item{"validate"}{Used for validating input.}
#'                        \item{"zero_LL"}{Used for getting null likelihood.}
#'                        \item{"conditionals"}{Used for getting conditionals.}
#'                        \item{"output"}{Used for preparing output after model estimation.}
#'                        \item{"raw"}{Used for debugging.}
#'                      }
#' @param indivID Numeric vector. Vector with individual's ID. As long as the number of observations.
#' @return Average probability over inter-individual draws (shape depends on argument \code{functionality}).
#'         \describe{
#'           \item{"estimate"}{Returns P averaged across inter-individual draws for the chosen alternative.}
#'           \item{"prediction"}{Returns P averaged across inter-individual draws for all alternatives.}
#'           \item{"validate"}{Returns P without changes.}
#'           \item{"zero_LL"}{Returns P without changes.}
#'           \item{"conditionals"}{Returns P without changes, but checks its shape.}
#'           \item{"output"}{Returns the same than "estimate", but also prints a summary of estimation data.}
#'           \item{"raw"}{Returns P without changes.}
#'         }
apollo_avgInterDraws <- function(P, apollo_control, functionality, indivID){

  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #

  if(apollo_control$HB==TRUE) return(P)

  # ############################### #
  #### pre-checks                ####
  # ############################### #

  if(!apollo_control$mixing) stop('No mixing used in model!')

  if(apollo_control$panelData){
    nIndiv <- length(unique( indivID ))
    if(is.list(P)){
      if(is.array(P[[1]])) pRows <- dim(P[[1]])[1] else pRows <- length(P[[1]])
    } else {
      if(is.array(P)) pRows <- dim(P)[1] else pRows <- length(P)
    }
    if(nIndiv!=pRows & !(functionality %in% c("zero_LL", "raw", "validate"))) stop("Observations from the same individual must be aggregated (e.g. multiplied) before averaging inter-individual draws.")
  }

  inputIsList <- is.list(P)

  if(inputIsList && functionality!="prediction" && is.null(P[["model"]])) stop('Element called "model" is missing in list P!')

  # ########################################## #
  #### functionality="zero_LL/raw/validate" ####
  # ########################################## #

  if(functionality %in% c("zero_LL","raw","validate")) return(P)

  # ########################################## #
  #### functionality="estimate" ####
  # ########################################## #

  if(functionality=="estimate"){
    if(inputIsList) P <- P[["model"]]
    if(is.vector(P) && !get("work_in_logs", parent.frame()) ) stop('No Inter-individuals draws to average over!')
    if(is.array(P) && length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(is.matrix(P)) P <- rowMeans(P)
    if(inputIsList) P <- list(model=P)
    return(P)
  }

  # ########################################## #
  #### functionality="prediction"           ####
  # ########################################## #

  if(functionality=="prediction"){
    if(!inputIsList){
      if(is.matrix(P)) output=rowMeans(P)
      if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
      return(output)
    }

    if(inputIsList){
      output_list=P
      j=1
      while(j <= length(P)){
        k=1
        while(k<= length(P[[j]])){
          if(is.array(P[[j]][[k]]) && length(dim(P[[j]][[k]]))==3) stop('Intra-individual draws still present to average over!')
          if(is.matrix(P[[j]][[k]])) output_list[[j]][[k]]=rowMeans(P[[j]][[k]])
          k=k+1
        }
        j=j+1
      }
      if(inputIsList) output_list=P
      return(output_list)
    }
  }

  # ########################################## #
  #### functionality="conditionals"         ####
  # ########################################## #

  if(functionality=="conditionals"){
    if(inputIsList) P <- P[["model"]]
    if(!is.array(P)) stop('No draws present to average over!')
    if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(inputIsList) P <- list(model=P)
    return(P)
  }

  # ########################################## #
  #### functionality="output"               ####
  # ########################################## #

  if(functionality=="output"){
    if(!inputIsList) P <- list(model=P)
    output_list=P
    j=1
    while(j<= length(P)){
      if(is.array(P[[j]]) & length(dim(P[[j]]))==3) stop('Intra-individual draws still present to average over!')
      if(is.matrix(P[[j]])) output_list[[j]]=rowMeans(P[[j]])
      j=j+1
    }
    if(!inputIsList) output_list <- output_list[[1]]
    return(output_list)
  }

}
