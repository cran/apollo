# UPDATED
#' Averages across inter-individual draws.
#'
#' Averages individual-specific likelihood across inter-individual draws.
#'
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call apollo_probabilities, though the user can also call apollo_probabilities manually with a given functionality for testing/debugging. Possible values are:
#'                      \itemize{
#'                        \item \strong{\code{"components"}}: For further processing/debugging, produces likelihood for each model component (if multiple components are present), at the level of individual draws and observations.
#'                        \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
#'                        \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"gradient"}}: For model estimation, produces analytical gradients of the likelihood, where possible.
#'                        \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
#'                        \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
#'                        \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
#'                        \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
#'                        \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"zero_LL"}}: Produces overall model likelihood with all parameters at zero.
#'                      }
#' @return Argument \code{P} with (for most functionalities) the original contents averaged over inter-individual draws. Shape depends on argument \code{functionality}.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"conditionals"}}: Returns \code{P} without averaging across draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"estimate"}}: Returns \code{P} containing the likelihood of the model averaged across inter-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"gradient"}}: Returns \code{P} containing the gradient of the likelihood averaged across inter-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"output"}}: Returns \code{P} containing the likelihood of all model components averaged across inter-individual draws.
#'           \item \strong{\code{"prediction"}}: Returns \code{P} containing the probabilities/likelihoods of all alternatives for all model components averaged across inter-individual draws.
#'           \item \strong{\code{"preprocess"}}: Returns \code{P} without changes.           
#'           \item \strong{\code{"raw"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"validate"}}: Returns \code{P} containing the likelihood of the model averaged across inter-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"zero_LL"}}: Returns \code{P} without changes.
#'         }
#' @export
apollo_avgInterDraws <- function(P, apollo_inputs, functionality){
  
  # ###################################################################### #
  #### load and check inputs, prepare variables that are used elsewhere ####
  # ###################################################################### #
  
  apollo_control <- apollo_inputs[["apollo_control"]]
  
  if(apollo_control$HB==TRUE) stop('Function apollo_avgInterDraws should not be used when apollo_control$HB==TRUE!')
  if(!apollo_control$mixing) stop('No mixing used in model!')
  
  inputIsList <- is.list(P)
  
  indivID <- apollo_inputs$database[, apollo_control$indivID]
  nIndiv  <- length(unique(indivID))
  if(inputIsList){
    if(is.array(P[[1]])) pRows <- dim(P[[1]])[1] else pRows <- length(P[[1]])
  } else {
    if(is.array(P)) pRows <- dim(P)[1] else pRows <- length(P)
  }
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality%in%c("components","preprocess","raw", "report")) return(P)
  
  # ############################### #
  #### functionality=="gradient" ####
  # ############################### #
  
  if(functionality=="gradient"){
    if(!is.list(P)) stop("Input P should be a list with at least one component (called model)!")
    if("model" %in% names(P)) P <- P$model
    if(is.null(P$like) || is.null(P$grad)) stop("Missing like and/or grad elements inside components!")
    if(apollo_control$workInLogs && apollo_control$analyticGrad) stop("Setting workInLogs cannot be used in conjunction with analyticGrad!")
    
    if(!is.matrix(P$like)) stop("No inter-draws to average over!")
    P$like <- rowMeans(P$like)
    P$grad <- lapply(P$grad, function(p) if(is.matrix(p)) rowMeans(p) else p)
    return(P)
  }
  
  # ############################## #
  #### functionality=="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    if(inputIsList && is.null(P[["model"]])) stop('Element called model is missing in list P!')
    if(is.list(P)){
      if(any(sapply(P, function(p) is.array(p) && length(dim(p))==3))) stop('Intra-individual draws still present to average over!')
      P <- lapply(P, function(p) if(is.matrix(p)) rowMeans(p) else p)
    }
    if(is.matrix(P)) P <- rowMeans(P) 
    return(P)
  } 
  
  # ########################################### #
  #### functionality=="estimate/validate"    ####
  # ########################################### #
  
  if(functionality %in% c("estimate", "validate")){
    if(nIndiv!=pRows) stop("Observations from the same individual must be combined (i.e. multiplied) before averaging over inter-individual draws.")
    if(inputIsList && is.null(P[["model"]])) stop('Element called model is missing in list P!')
    if(inputIsList) P <- P[["model"]]
    if(is.vector(P) && !apollo_control$workInLogs ) stop('No Inter-individuals draws to average over!')
    if(is.array(P) && length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(is.matrix(P)) P <- rowMeans(P)
    if(inputIsList) P <- list(model=P)
    return(P)
  }
  
  # ############################# #
  #### functionality=="output" ####
  # ############################# #
  
  if(functionality=="output"){
    if(nIndiv!=pRows) stop("Observations from the same individual must be combined (i.e. multiplied) before averaging over inter-individual draws.")
    if(inputIsList && is.null(P[["model"]])) stop('Element called model is missing in list P!')
    if(!inputIsList) P <- list(model=P)
    for(j in 1:length(P)){
      if(is.array(P[[j]]) & length(dim(P[[j]]))==3) stop('Intra-individual draws still present to average over!')
      if(is.matrix(P[[j]])) P[[j]]=rowMeans(P[[j]])
    }
    if(!inputIsList) P <- P[[1]]
    return(P)
  }
  
  # ################################# #
  #### functionality=="prediction" ####
  # ################################# #
  
  if(functionality=="prediction"){
    nInter <- apollo_inputs$apollo_draws$interNDraws
    if(!inputIsList){
      if(is.matrix(P) & ncol(P)==nInter) output=rowMeans(P)
      if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
      return(output)
    } else {
      #output_list=P
      for(j in 1:length(P)){
        if(is.list(P[[j]])){
          for(k in 1:length(P[[j]])){
            if(is.array(P[[j]][[k]]) && length(dim(P[[j]][[k]]))==3) stop('Intra-individual draws still present to average over!')
            if(is.matrix(P[[j]][[k]]) & ncol(P[[j]][[k]])==nInter) P[[j]][[k]]=rowMeans(P[[j]][[k]])
          }}else{
            if(is.array(P[[j]]) && length(dim(P[[j]]))==3) stop('Intra-individual draws still present to average over!')
            if(is.matrix(P[[j]]) && ncol(P[[j]])==nInter) P[[j]]=rowMeans(P[[j]])
          }
      }
      return(P)
    }
  }
  
  # ################################### #
  #### functionality=="conditionals" ####
  # ################################### #
  
  if(functionality=="conditionals"){
    if(nIndiv!=pRows) stop("Observations from the same individual must be combined (i.e. multiplied) before averaging over inter-individual draws.")
    if(inputIsList && is.null(P[["model"]])) stop('Element called model is missing in list P!')
    if(inputIsList) P <- P[["model"]]
    #if(!is.array(P)) stop('No draws present to average over!')
    if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(inputIsList) P <- list(model=P)
    return(P)
  }
}
