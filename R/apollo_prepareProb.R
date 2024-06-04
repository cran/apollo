# UPDATED
#' Checks likelihood function
#' 
#' Checks that the likelihood function for the mode is in the appropriate format to be returned.
#' 
#' This function should be called inside \code{apollo_probabilities}, near the end of it, just before \code{return(P)}.
#' This function only performs checks on the shape of P, but does not change its values.
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
#' @return Argument \code{P} with (for most functionalities) the original contents. Output depends on argument \code{functionality}.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"conditionals"}}: Returns only the \code{"model"} component of argument \code{P}.
#'           \item \strong{\code{"estimate"}}: Returns only the \code{"model"} component of argument \code{P}.
#'           \item \strong{\code{"gradient"}}: Returns only the \code{"model"} component of argument \code{P}.
#'           \item \strong{\code{"output"}}: Returns argument \code{P} without any changes to its content, but gives names to unnamed elements.
#'           \item \strong{\code{"prediction"}}: Returns argument \code{P} without any changes.
#'           \item \strong{\code{"preprocess"}}: Returns argument \code{P} without any changes to its content, but gives names to elements corresponding to componentNames.
#'           \item \strong{\code{"raw"}}: Returns argument \code{P} without any changes.
#'           \item \strong{\code{"report"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"shares_LL"}}: Returns argument \code{P} without any changes to its content, but gives names to unnamed elements.
#'           \item \strong{\code{"validate"}}: Returns argument \code{P} without any changes.
#'           \item \strong{\code{"zero_LL"}}: Returns argument \code{P} without any changes to its content, but gives names to unnamed elements.
#'         }
#' @export
apollo_prepareProb=function(P, apollo_inputs, functionality){

  # ###################################################################### #
  #### load and check inputs, prepare variables that are used elsewhere ####
  # ###################################################################### #

  if(!is.list(P)) P=list(model=P)
  ## change 8 July
  #if(is.null(P[["model"]]) && !(functionality %in% c("prediction", "gradient", "preprocess")) ) stop('SYNTAX ISSUE - Element called model is missing in list P!')
  if(is.null(P[["model"]]) && !(functionality %in% c("prediction", "gradient", "hessian", "preprocess", "report", "components","raw")) ) stop('SYNTAX ISSUE - Element called model is missing in list P!')
  ### end change
  panelData <- apollo_inputs$apollo_control$panelData
  nIndiv <- length(unique(apollo_inputs$database[, apollo_inputs$apollo_control$indivID]))
  
  # ############################################### #
  #### functionality="validate" ####
  # ############################################### #
  
  if(functionality=="validate"){
    dimP <- c(1, 1, 1)
    if(is.array(P[["model"]])){
      if(length(dim(P[["model"]]))==2) dimP[1:2] <- dim(P[["model"]])
      if(length(dim(P[["model"]]))==3) dimP      <- dim(P[["model"]])
    } else if(is.vector(P[["model"]])) dimP[1] <- length(P[["model"]]) else {
      stop("SPECIFICATION ISSUE - apollo_probabilities needs to return ",
           "a valid vector of likelihood values.")
    }
    
    if(panelData && dimP[1]>nIndiv && !apollo_inputs$apollo_control$HB) stop("SPECIFICATION ISSUE - Need to ", 
                                         "multiply across observations for the same ", 
                                         "individual! (see ?apollo_panelProd)")
    if(dimP[3]>1) stop("SPECIFICATION ISSUE - Need to average over intra-", 
                       "individual draws! (see ?apollo_avgIntraDraws)")
    if(dimP[2]>1) stop("SPECIFICATION ISSUE - Need to average over inter-", 
                       "individual draws! (see ?apollo_avgInterDraws)")
    return(P)
  }
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality %in% c("components", "prediction", "raw", "report")) return(P)

  # ################### #
  #### HB estimation ####
  # ################### #

  if((apollo_inputs$apollo_control$HB)&&(functionality=="estimate")){
    test = ifelse(is.na(P[["model"]]), TRUE, P[["model"]] < 9.88131291682493e-324)
    if(any(test)){
      tmp <- globalenv()
      if(exists("apollo_HBcensor", envir=tmp)){
        assign("apollo_HBcensor", get("apollo_HBcensor", envir=tmp) + 1, envir=tmp)
      } else assign("apollo_HBcensor", 1, envir=tmp)
    }
    return(P[["model"]])
  } 

  # ########################################## #
  #### functionality="preprocess"           ####
  # ########################################## #
  
  if(functionality=="preprocess"){
    # If there are multiple pre-processing inside a single element of P, promote them to the upper level
    # This can happen in latent class models with multiple components
    isList <- rep(FALSE, length(P))
    for(s in 1:length(P)){
      isList[s] <- is.list(P[[s]]) && all(sapply(P[[s]], is.list)) && all(sapply(P[[s]], function(p) !is.null(p$componentName)))
    }
    if(any(isList)){
      isList <- which(isList)
      for(i in isList) P <- c(P, P[[i]])
      P <- P[-isList]
    }
    
    # Rename pre-processing elements
    for(s in 1:length(P)) if(is.list(P[[s]]) && !is.null(P[[s]]$componentName)){
      names(P)[s]=paste0(P[[s]]$componentName, "_settings")
    }
    
    # Return
    return(P)
  } 

  # ########################################## #
  #### functionality="gradient"             ####
  # ########################################## #
  
  if(functionality=="gradient"){
    if("model" %in% names(P)) P <- P$model
    if(is.null(P$like) || is.null(P$grad)) stop("INTERNAL ISSUE - Missing like and/or grad elements inside components!")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("INTERNAL ISSUE - workInLogs cannot be used in conjunction with analyticGrad!")
    test <- is.array(P$like) && ncol(P$like)==1 && (length(dim(P$like))==2 || (length(dim(P$like))==3 && dim(P$like)[3]==1))
    if(test) P$like <- as.vector(P$like)
    if(is.list(P)) P = do.call(cbind, P$grad)/P$like
    ### New 29 Sept
    #if(!is.null(apollo_inputs$scaling)){
    #  tmp <- apollo_inputs$apollo_beta_names
    #  tmp <- tmp[!(tmp %in% apollo_inputs$apollo_fixed)]
    #  for(i in names(apollo_inputs$scaling)) P[,which(tmp==i)]=P[,which(tmp==i)]/apollo_inputs$scaling[i]
    #}
    ### end new
    return(P)
  }
  
  # ########################################## #
  #### functionality="hessian"              ####
  # ########################################## #
  
  if(functionality=="hessian"){
    if("model" %in% names(P)) P <- P$model
    K <- length(P[["grad"]])
    H <- matrix(0, nrow=K, ncol=K)
    for(k1 in 1:K) for(k2 in 1:k1){
      #H[k1,k2] <- sum((P[["hess"]][,k1,k2] - 
      H[k1,k2] <- sum((P[["hess"]][[k1]][[k2]] - 
                  P[["grad"]][[k1]]*P[["grad"]][[k2]]/P[["like"]])/P[["like"]])
      H[k2,k1] <- H[k1,k2]
    }
    return(H)
    #return(apply(P[["hess"]], MARGIN=c(2,3), sum))
  }
  # ######################################### #
  #### functionality="estimation"          ####
  # ######################################### #
  
  if(functionality=="estimate"){
    
    #### special bit for EM
    if(!is.null(apollo_inputs$EM) && apollo_inputs$EM && !is.null(apollo_inputs$class_specific) && (apollo_inputs$class_specific)>0){
      ### Apply weights
      P = apollo_weighting(P, apollo_inputs, functionality)  
    }

    if(is.array(P[["model"]])) nPRows <- dim(P[["model"]])[1] else nPRows <- length(P[["model"]])
    #if(panelData && nPRows>nIndiv && !apollo_inputs$apollo_control$HB) stop("SPECIFICATION ISSUE - Need to multiply over observations for the same individual! (see ?apollo_panelProd)")
    if(panelData && nPRows>nIndiv && !apollo_inputs$apollo_control$HB && !apollo_inputs$apollo_control$overridePanel) stop("SPECIFICATION ISSUE - Need to multiply over observations for the same individual! (see ?apollo_panelProd)")
    
    if(is.array(P[["model"]])){
      if(length(dim(P[["model"]]))==3) stop('SPECIFICATION ISSUE - Need to average over intra-individual draws! (see ?apollo_avgIntraDraws)')
      if(dim(P[["model"]])[2]>1) stop('SPECIFICATION ISSUE - Need to average over inter-individual draws! (see ?apollo_avgInterDraws)')
      if(dim(P[["model"]])[2]==1) P[["model"]]=as.vector(P[["model"]])
    }
    return(P[["model"]])
  } 
  
  # ############################################## #
  #### functionality="output/zero_LL/shares_LL" ####
  # ############################################## #

  if(functionality%in%c("output","zero_LL","shares_LL")){
    # Give name to unnamed components
    origNames <- names(P)
    newNames  <- paste0("component_", 1:length(P))
    if(!is.null(origNames)) newNames <- ifelse(origNames!="", origNames, newNames)
    names(P) <- newNames
    return(P)
  }

  # ######################################### #
  #### functionality="conditionals"        ####
  # ######################################### #
  
  if(functionality=="conditionals") return(P[["model"]])
  
}