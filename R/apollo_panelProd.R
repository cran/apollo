# UPDATED
#' Calculates product across observations from same individual.
#' 
#' Multiplies likelihood of observations from the same individual, or adds the log of them.
#' 
#' This function should be called inside apollo_probabilities only if the data has a panel structure.
#' It should be called after apollo_avgIntraDraws if intra-individual draws are used.
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
#'           \item \strong{\code{"gradient"}}: Returns \code{P} containing the gradient of the likelihood after applying the product rule across observations for the same individual.
#'           \item \strong{\code{"output"}}: Returns \code{P} containing the likelihood of all model components averaged across inter-individual draws.
#'           \item \strong{\code{"prediction"}}: Returns \code{P} containing the probabilities/likelihoods of all alternatives for all model components averaged across inter-individual draws.
#'           \item \strong{\code{"preprocess"}}: Returns \code{P} without changes.           
#'           \item \strong{\code{"raw"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"validate"}}: Returns \code{P} containing the likelihood of the model averaged across inter-individual draws. Drops all components except \code{"model"}.
#'           \item \strong{\code{"zero_LL"}}: Returns \code{P} without changes.
#'         }
#' @export
apollo_panelProd <- function(P, apollo_inputs, functionality){

  # ###################################################################### #
  #### load and check inputs, prepare variables that are used elsewhere ####
  # ###################################################################### #

  apollo_control = apollo_inputs[["apollo_control"]]
  if(apollo_control$HB==TRUE) stop('Function apollo_panelProd should not be used when apollo_control$HB==TRUE!')
  if(!apollo_control$panelData) stop('Panel data setting not used, so multiplying over choices using apollo_panelProd not applicable!')
    
  inputIsList <- is.list(P)
  indivID <- apollo_inputs$database[, apollo_control$indivID]

  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality%in%c("components","prediction","preprocess","raw", "report")) return(P)

  # ########################################## #
  #### functionality="gradient"             ####
  # ########################################## #

  if(functionality=="gradient"){
    # Checks
    if(!is.list(P)) stop("Input P should be a list with at least one component (called model)")
    if("model" %in% names(P)) P <- P$model
    if(is.null(P$like) || is.null(P$grad)) stop("Missing like and/or grad elements inside components")
    if(apollo_control$workInLogs && apollo_control$analyticGrad) stop("workInLogs cannot be used in conjunction with analyticGrad")
    
    # Remove zeros from P$like
    P$like[P$like==0] <- 1e-50 #1e-100 #1e-16 #1e-300
    
    # Average gradient respecting product rule
    tmp <- apollo_panelProd(P$like, apollo_inputs, "estimate")
    P$grad <- lapply(P$grad, function(g) tmp*rowsum(g/P$like, group=indivID, reorder=FALSE))
    P$like <- tmp
    return(P)
  }
  
  # ################################################################### #
  #### functionality="estimate/conditionals/validate/output/zero_LL" ####
  # ################################################################### #
  
  if(functionality=='zero_LL'){
    if(!inputIsList) P <- list(model=P)
    if(any(sapply(P, function(p) is.array(p) && length(dim(p)==3)))) stop('Need to average over intra-individual draws first before multiplying over choices!')
    for(j in 1:length(P)){
      test <- is.vector(P[[j]]) && length(P[[j]])==length(indivID)
      test <- test || is.matrix(P[[j]]) && nrow(P[[j]])==length(indivID)
      if(test){
        P[[j]] <- rowsum(log(P[[j]]), group=indivID, reorder=FALSE)
        if(!apollo_control$workInLogs) P[[j]] <- exp(P[[j]])
      }
      if(is.matrix(P[[j]]) && ncol(P[[j]])==1) P[[j]] <- as.vector(P[[j]])
    }
    if(!inputIsList) P <- P[[1]]
    return(P)
  }
  
  if(functionality %in% c("estimate", "conditionals","output","validate")){
    if(inputIsList && is.null(P[["model"]])) stop("The list \"P\" should contain an element called \"model\"!")
    if(functionality %in% c("estimate", "conditionals","validate") && inputIsList) P <- P["model"]  ### only keep the model element, but P is still a list
    if(!inputIsList) P <- list(model=P)
    for(j in 1:length(P)){
      if(is.array(P[[j]]) && length(dim(P[[j]]))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
      if(is.vector(P[[j]]) || (is.matrix(P[[j]]) && !apollo_control$workInLogs) ){
        P[[j]] <- rowsum(log(P[[j]]), group=indivID, reorder=FALSE)
        if(!apollo_control$workInLogs) P[[j]] <- exp(P[[j]])
      }
      if(apollo_control$panelData && is.matrix(P[[j]]) && apollo_control$workInLogs && nrow(P[[j]])==length(indivID)){
        # approach to use if working in logs with mixing
        B    <- rowsum(log(P[[j]]), group=indivID, reorder=FALSE) # nIndiv x nDraws
        Bbar <- apply(B, MARGIN=1, function(r) mean(r[is.finite(r)]) ) # nIndiv x 1 ### FIX: 15/5/2020
        P[[j]] <- ifelse(is.finite(Bbar), Bbar + log( rowMeans(exp(B-Bbar)) ), -Inf) # nIndiv x 1 ### FIX 18/5/2020
      }
      if(is.matrix(P[[j]]) && ncol(P[[j]])==1) P[[j]] <- as.vector(P[[j]])
    }
    if(!inputIsList) P <- P[[1]]
    #if(functionality %in% c("estimate", "conditionals","validate") && inputIsList) Pout <- list(model=Pout) # Removed. Pout is already a list
    return(P)
  }

}