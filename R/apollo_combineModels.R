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
#'                        \item"estimate": For model estimation, returns likelihood of model.
#'                        \item"prediction": For model predictions, returns probabilities of all alternatives.
#'                        \item"validate": Validates input.
#'                        \item"zero_LL": Return probabilities with all parameters at zero.
#'                        \item"conditionals": For conditionals, returns likelihood of model.
#'                        \item"output": Checks that the model is well defined.
#'                        \item"raw": For debugging, returns probabilities of all alternatives.
#'                      }
#' @return Argument \code{P} with an extra element called "model", which is the product of all the other elements.
#' @export
apollo_combineModels=function(P, apollo_inputs, functionality){
  
  if(!is.null(P[["model"]])) warning("\nA component called \"model\" already exists in P before calling apollo_combineModels!")
  
  if(length(P)==1) warning("\nNo need to call apollo_combineModels for models with only one component!")
  
  if(functionality=="prediction" | functionality=="raw") return(P)

  elements = names(P)
  
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