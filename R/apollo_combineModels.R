#' Combines separate model components.
#' 
#' Combines model components to create probability for overall model.
#' 
#' This function should be called inside apollo_probabilities after all model components have been produced.
#' 
#' It should be called before apollo_avgInterDraws, apollo_avgIntraDraws, apollo_panelProd and apollo_prepareProb, whichever apply. 
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"}: For model estimation, returns likelihood of model.
#'                        \item \code{"prediction"}: For model predictions, returns probabilities of all alternatives.
#'                        \item \code{"validate"}: Validates input.
#'                        \item \code{"zero_LL"}: Return probabilities with all parameters at zero.
#'                        \item \code{"conditionals"}: For conditionals, returns likelihood of model.
#'                        \item \code{"output"}: Checks that the model is well defined.
#'                        \item \code{"raw"}: For debugging, returns probabilities of all alternatives.
#'                      }
#' @return Argument \code{P} with an extra element called "model", which is the product of all the other elements. Shape depends on argument \code{functionality}.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: Returns argument \code{P} with an extra component called \code{"model"}, which is the product of all other elements of \code{P}.
#'           \item \strong{\code{"prediction"}}: Returns argument \code{P} without any change.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"zero_LL"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"raw"}}: Returns argument \code{P} without any change.
#'         }
#' @export
apollo_combineModels=function(P, apollo_inputs, functionality){
  
  # ############################### #
  #### pre-checks                ####
  # ############################### #
  
  if(!is.null(P[["model"]])) stop("\nA component called \"model\" already exists in P before calling apollo_combineModels!")
  
  if(length(P)==1) stop("\nNo need to call apollo_combineModels for models with only one component!")
  
  # ########################################## #
  #### functionality="gradient"             ####
  # ########################################## #
  if(functionality=="gradient"){
    if(length(P)==1) return(P[[1]]) else {
      ### If multiple components, implement product rule
      g <- function(b, db){
        #like   <- get("apollo_probabilities", envir=parent.frame())
        #inputs <- get("apollo_inputs", envir=parent.frame(), inherits=TRUE)
        #LL <- tryCatch(get("LL", envir=parent.frame(), inherits=FALSE),
        #               error=function(e) like(b, inputs, functionality="raw")) # this is probably wrong for mixing
        #dLL <- lapply(P, function(p) p(b, db))
        #LL  <- LL[names(dLL)] # remove "model" and make sure both LL and dLL are in the same order
        #LL  <- lapply(LL, log)
        #dLL <- mapply("/", dLL, LL, SIMPLIFY=FALSE)
        #LL  <- Reduce("+", LL)
        #dLL <- Reduce("+", dLL)
        #return(LL*dLL)
        dLL <- Reduce("+", lapply(P, function(p) p(b, db)) )
        return(dLL)
      }
      environment(g) <- new.env(parent=baseenv())
      assign("P", P, envir=environment(g))
      return(g)
    }
  }
  
  
  # ########################################## #
  #### functionality="prediction"           ####
  # ########################################## #
  
  if(functionality=="prediction" | functionality=="raw") return(P)

  # ########################################## #
  #### functionality!="prediction"          ####
  # ########################################## #
  
  elements = names(P)
  if(is.null(elements) || length(unique(elements))<length(elements)){
    stop("For models using multiple components, all components in P must be named and all names must be unique!")
  } 
  
  if(!apollo_inputs$apollo_control$workInLogs){
    P[["model"]] = P[[elements[1]]]
    k = 2
    while(k <= length(elements)){
      P[["model"]] = P[["model"]]*P[[elements[k]]]
      k=k+1
    }
  } else {
    P[["model"]] <- exp(Reduce("+", lapply(P, log)))
  }

  return(P)
}