#' Calculates conditionals of a latent class model.
#' 
#' Calculates posterior expected values (conditionals) of class allocation probabilities for each individual.
#' 
#' This function can only be used with latent class models without continuous heterogeneity.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return A matrix with the posterior class allocation probabilities for each individual.
#' @export
apollo_lcConditionals=function(model, apollo_probabilities, apollo_inputs){
  apollo_beta    = model$estimate
  apollo_fixed   = model$apollo_fixed
  apollo_control = apollo_inputs[["apollo_control"]]
  database       = apollo_inputs[["database"]]
  apollo_lcPars  = apollo_inputs[["apollo_lcPars"]]
  class_prob     = "pi_values" # name of lcpars component with allocation probabilities
  
  ### Validation
  if(apollo_control$mixing) stop("apollo_lcConditionals can only be used for latent class models without continuous random heterogeneity")
  
  cat("Calculating conditionals...")
  ### Get allocation and inClass probs
  lcpars = with(c(apollo_beta, apollo_inputs$database, apollo_inputs$draws), {
    environment(apollo_lcPars) <- environment()
    apollo_lcPars(apollo_beta, apollo_inputs)
  })
  if(is.null(lcpars[[class_prob]])) stop("The lcpars function needs to create an object called \"pi_values\"!")
  L      = apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  classes    = length(lcpars[[class_prob]])
  components = length(L)
  if(components>(classes+1)) stop("apollo_lcConditionals can only be used for latent class models alone (i.e. no hybrid choice)")
  if(components!=(classes+1)) stop("Model should contain one component per class, and an overall model!")
  
  ### Calculate posterior class allocation probs
  post_pi = vector(mode="list", length=classes)
  for(s in 1:classes) post_pi[[s]] = lcpars[[class_prob]][[s]]*L[[s]]/L[["model"]]
  
  ### Prepare output
  conditionals = matrix(unlist(post_pi), ncol = length(post_pi), byrow = FALSE)
  classnames   = paste("Class ",seq(1:classes),sep="")
  rownames(conditionals) = c(unique(database[,apollo_control$indivID]))
  colnames(conditionals) = c(classnames)
  cat("Done.\n")
  
  return(conditionals)
}