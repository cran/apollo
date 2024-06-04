# UPDATED
#' Combines separate model components.
#' 
#' Combines model components to create likelihood for overall model.
#' 
#' This function should be called inside apollo_probabilities after all model components have been produced.
#' 
#' It should be called before \link{apollo_avgInterDraws}, \link{apollo_avgIntraDraws}, \link{apollo_panelProd} and \link{apollo_prepareProb}, whichever apply, except where these functions are called inside any latent class components of the overall model. 
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call \code{apollo_probabilities}, though the user can also call \code{apollo_probabilities} manually with a given functionality for testing/debugging. Possible values are:
#'                      \itemize{
#'                        \item \strong{\code{"components"}}: For further processing/debugging, produces likelihood for each model component (if multiple components are present), at the level of individual draws and observations.
#'                        \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
#'                        \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"gradient"}}: For model estimation, produces analytical gradients of the likelihood, where possible.
#'                        \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
#'                        \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
#'                        \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
#'                        \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
#'                        \item \strong{\code{"report"}}: Prepares output summarising model and choiceset structure.
#'                        \item \strong{\code{"shares_LL"}}: Produces overall model likelihood with constants only.
#'                        \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"zero_LL"}}: Produces overall model likelihood with all parameters at zero.
#'                      }
#' @param components Character vector. Optional argument. Names of elements in P that should be multiplied to construct the whole model likelihood. If a single element is provided, it is interpreted as a regular expression. Default is to include all components in P.
#' @param asList Logical. Only used if \code{functionality} is \code{"conditionals","estimate","validate","zero_LL"} or \code{"output"}. If \code{TRUE}, it will return a list as described in the 'Value' section. If \code{FALSE}, it will only return a vector/matrix/3-dim array of the product of likelihoods inside P. Default is \code{TRUE}.
#' @return Argument \code{P} with (for most functionalities) an extra element called "model", which is the product of all the other elements. Shape depends on argument \code{functionality}.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"conditionals"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"estimate"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"gradient"}}: Returns \code{P} containing the gradient of the likelihood after applying the product rule across model components.
#'           \item \strong{\code{"output"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"prediction"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"preprocess"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"raw"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"shares_LL"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"validate"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"zero_LL"}}: Returns \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'         }
#' @export
apollo_combineModels=function(P, apollo_inputs, functionality, components=NULL, asList=TRUE){
  
  # ###################################################################### #
  #### load and check inputs, prepare variables that are used elsewhere ####
  # ###################################################################### #
  
  if(length(P)==1) stop("SPECIFICATION ISSUE - No need to call apollo_combineModels for models with only one component!")
  if(!is.null(P[["model"]])) stop("SPECIFICATION ISSUE - A component called model already exists in P before calling apollo_combineModels!")
  
  elements = names(P)
  if(is.null(elements) || length(unique(elements))<length(elements)){
    stop("SPECIFICATION ISSUE - For models using multiple components, all components in P must be named and all names must be unique!")
  }
  
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality %in% c("components","prediction",
                          "preprocess","raw", "report")) return(P)
  
  
  # ############################ #
  #### Drop unused components ####
  # ############################ #
  if(!is.null(components)){
    if(!is.character(components)) stop('SYNTAX ISSUE - Argument "components", if provided, should be a character vector.')
    if(length(components)==1) components <- grep(components, names(P), value=TRUE)
    if(!all(components %in% names(P))) stop('SYNTAX ISSUE - Elements', paste0(components[-which(components %in% names(P))], collapse=', '),
                                            ' do not exist in "P".')
    P <- P[components]
  }
  
  
  # ########################################### #
  #### functionality=="gradient"             ####
  # ########################################### #
  
  if(functionality=="gradient"){
    # Checks
    if(!is.list(P)) stop("INTERNAL ISSUE - Input P should be a list with at least one component")
    if(any(sapply(P, function(p) is.null(p$like) || is.null(p$grad)))) stop("INTERNAL ISSUE - Some components are missing the like and/or grad elements")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("INCORRECT FUNCTION/SETTING USE - workInLogs cannot be used in conjunction with analyticGrad")
    K <- length(P[[1]]$grad) # number of parameters
    if(any(sapply(P, function(p) length(p$grad))!=K)) stop("INTERNAL ISSUE - Dimensions of gradients from different components imply different number of parameters")
    
    # Remove zeros from components' likelihoods
    for(i in 1:length(P)) P[[i]]$like[P[[i]]$like==0] <- 1e-50
    
    # Create "model" component
    model <- list(like=Reduce("*", lapply(P, function(p) p$like)), 
                  grad=list())
    for(k in 1:K) model$grad[[k]] <- model$like*Reduce("+", lapply(P, function(p) p$grad[[k]]/p$like))
    
    return( list(model=model) )
  }
  
  # ############################# #
  #### functionality="hessian" ####
  # ############################# #
  if(functionality=="hessian"){
    
    # Checks
    if(!is.list(P)) stop("INTERNAL ISSUE - Input P should be a list with at least one component")
    if(any(sapply(P, function(p) is.null(p$like) || is.null(p$grad) || is.null(p$hess)))) stop("INTERNAL ISSUE - Some components are missing the like, grad and/or hess elements")
    K <- length(P[[1]]$grad) # number of parameters
    if(any(sapply(P, function(p) length(p$grad))!=K)) stop("INTERNAL ISSUE - Dimensions of gradients from different components imply different number of parameters")
    if(any(sapply(P, function(p) length(p$hess))!=K)) stop("INTERNAL ISSUE - Dimensions of hessian from different components imply different number of parameters")

    # Remove zeros from components' likelihoods
    for(i in 1:length(P)) P[[i]]$like[P[[i]]$like==0] <- 1e-50
    
    # Pre-calculate necessary elements for the hessian
    L <- Reduce("*", lapply(P, function(p) p$like))
    L[L==0] <- 1e-50 # Remove zeros from L
    d1L <- list()
    for(k in 1:K) d1L[[k]] <- L*Reduce("+", lapply(P, function(p) p$grad[[k]]/p$like))
    
    d2L <- vector(mode="list", length=K)
   
    # Loop to calculate hessian
    for(k1 in 1:K){
      d2L[[k1]] <- vector(mode="list", length=K)
      for(k2 in 1:k1){
        part1 <- d1L[[k2]]*d1L[[k1]]/L
        part2 <- L*Reduce("+", lapply(P, function(p) (p[["hess"]][[k1]][[k2]] - p[["grad"]][[k1]]*p[["grad"]][[k2]]/p[["like"]])/p[["like"]]))
        tmp <- part1 + part2
        if(is.matrix(tmp) && dim(tmp)[2]==1) tmp <- as.vector(tmp)
        d2L[[k1]][[k2]] = tmp
        d2L[[k2]][[k1]] = tmp
      }
    }; rm(tmp)
    
    H <- list(
      like = L,    # likelihoods L after combining models
      grad = d1L, # first derivatives of L
      hess = d2L) # second derivatives of L
    return(H)
    
  }
  
  # ############################################################# #
  #### functionality=="conditionals/estimate/validate/zero_LL" ####
  # ############################################################# #
  
  if(functionality %in% c("conditionals","estimate","validate","zero_LL", "shares_LL", "output")){
    if(!apollo_inputs$apollo_control$workInLogs){
      P[["model"]] <- Reduce("*", P)
    } else {
      P[["model"]] <- exp(Reduce("+", lapply(P, log)))
    }
    if(asList) return(P) else return(P[['model']])
  }
}