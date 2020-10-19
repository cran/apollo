# UPDATED
#' Averages across intra-individual draws.
#'
#' Averages observation-specific likelihood across intra-individual draws.
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
#' @return Argument \code{P} with (for most functionalities) the original contents averaged over intra-individual draws. Shape depends on argument \code{functionality}.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"conditionals"}}: Returns \code{P} containing the likelihood of the model averaged across intra-individual draws. Drops all components except for \code{"model"}.
#'           \item \strong{\code{"estimate"}}: Returns \code{P} containing the likelihood of the model averaged across intra-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"gradient"}}: Returns \code{P} containing the gradient of the likelihood averaged across intra-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"output"}}: Returns \code{P} containing the likelihood of all model components averaged across intra-individual draws.
#'           \item \strong{\code{"prediction"}}: Returns \code{P} containing the probabilities of all alternatives for all model components averaged across intra-individual draws.
#'           \item \strong{\code{"preprocess"}}: Returns \code{P} without changes.           
#'           \item \strong{\code{"raw"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"validate"}}: Returns \code{P} containing the likelihood of the model averaged across intra-individual draws. Drops all components but \code{"model"}.
#'           \item \strong{\code{"zero_LL"}}: Returns \code{P} without changes.
#'         }
#' @export
apollo_avgIntraDraws <- function(P, apollo_inputs, functionality){
  
  # ###################################################################### #
  #### load and check inputs, prepare variables that are used elsewhere ####
  # ###################################################################### #
  
  apollo_control=apollo_inputs[["apollo_control"]]
  
  if(apollo_control$HB==TRUE) stop('Function apollo_avgIntraDraws should not be used when apollo_control$HB==TRUE!')
  if(!apollo_control$mixing) stop('No mixing used in model!')

  isCube <- function(x) is.array(x) && length(dim(x))==3
  inputIsList <- is.list(P)
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality%in%c("components","preprocess","raw", "report")) return(P)
  
  # ########################################### #
  #### functionality=="gradient"             ####
  # ########################################### #

  if(functionality=="gradient"){
    # Checks
    if(!is.list(P)) stop("Input P should be a list with at least one component")
    if(any(sapply(P, function(p) is.null(p$like) || is.null(p$grad)))) stop("Some components are missing the like and/or grad elements")
    if(apollo_control$workInLogs && apollo_control$analyticGrad) stop("workInLogs cannot be used in conjunction with analyticGrad")
    K <- length(P[[1]]$grad) # number of parameters
    if(any(sapply(P, function(p) length(p$grad))!=K)) stop("Dimensions of gradients from different components imply different number of parameters")
    
    # Average intra draws for like and grad of each component
    for(i in 1:length(P)){
      cNam <- apollo_inputs$apolloLog$listOfNames[i]
      test <- !is.null(P[[i]]$like) && !is.null(P[[i]]$grad)
      if(!test) stop("Elements like and/or grad missing for component ", cNam)
      test <- isCube(P[[i]]$like) && is.list(P[[i]]$grad) && any(sapply(P[[i]]$grad,isCube))
      if(!test) stop("Elements like or grad for component ", cNam, " have the wrong dimensions. Maybe there is no need to call apollo_avgIntraDraws")
      P[[i]]$like <- apply(P[[i]]$like, MARGIN=c(1,2), sum)/dim(P[[i]]$like)[3]
      P[[i]]$grad <- lapply(P[[i]]$grad, function(g) if(isCube(g)) apply(g, MARGIN=c(1,2), sum)/dim(g)[3] else g)
    }
    return(P)
  }
  
  # ############################## #
  #### functionality=="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    if(is.list(P)) P <- lapply(P, function(p) if(isCube(p)) apply(p, MARGIN=c(1,2), sum)/dim(p)[3] else p)
    if(isCube(P)) P <- apply(P, MARGIN=c(1,2), sum)/dim(P)[3]
    return(P)
  } 

  # ##################################################### #
  #### functionality=="estimate/conditionals/validate" ####
  # ##################################################### #
  
  if(functionality %in% c("estimate", "conditionals", "validate")){
    if(is.list(P)){
      if(!any(sapply(P, isCube))) stop('No intra-individual draws present to average over on any component!')
      P <- lapply(P, function(p) if(isCube(p)) apply(p, MARGIN=c(1,2), sum)/dim(p)[3] else p)
    } else {
      if(!isCube(P)) stop('No intra-individual draws present to average over!')
      P <- apply(P, MARGIN=c(1,2), sum)/dim(P)[3]
    }
    return(P)
  }
  
  # ########################################### #
  #### functionality=="output"               ####
  # ########################################### #
  
  if(functionality=="output"){
    if(!is.list(P)){
      if(isCube(P)) P <- apply(P, MARGIN=c(1,2), mean)
    } else {
      for(j in 1:length(P)){
        if(isCube(P[[j]])) P[[j]] <- apply(P[[j]], MARGIN=c(1,2), mean)
      }
    }
    return(P)
  }
  
  # ########################################### #
  #### functionality=="prediction"           ####
  # ########################################### #
  
  if(functionality=="prediction"){
    if(!is.list(P)){
      if(isCube(P)) P <- apply(P, MARGIN=c(1,2), mean)
    } else {
      for(j in 1:length(P)){
        for(k in 1:length(P[[j]])){
          if(isCube(P[[j]][[k]])) P[[j]][[k]] <- apply(P[[j]][[k]], MARGIN=c(1,2), mean)
        } 
      }
    }
    return(P)
  }
}
