#' Calculates conditionals
#' 
#' Calculates posterior expected values (conditionals) of random coefficient models (continuous or discrete mixtures/latent class)
#' 
#' This functions is only meant for use with models using either continuous distributions or latent classes, not both at the same time
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return Depends on whether the model uses continuous mixtures or latent class.
#'         \itemize{
#'           \item If the model contains a continuous mixture, the function returns a list of matrices. 
#'                 Each matrix has dimensions nIndiv x 3. One matrix per random component.
#'                 Each row of each matrix contains the indivID of an individual, and the
#'                 posterior mean and s.d. of this random component for this individual.
#'           \item If the model contains latent classes, the function returns a matrix with 
#'                 the posterior class allocation probabilities for each individual.
#'           \item If the model contains both continuous mixtures and latent classes,
#'                 the function fails.
#'         }
#' @export
apollo_conditionals=function(model, apollo_probabilities, apollo_inputs){
  if(is.null(apollo_inputs$silent)) silent = FALSE else silent = apollo_inputs$silent
  apollo_beta  = model$estimate
  apollo_fixed = model$apollo_fixed
  
  apollo_compareInputs(apollo_inputs)
  
  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_draws     = apollo_inputs[["apollo_draws"]]
  apollo_lcPars    = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  continuous       = apollo_control$mixing
  latentClass      = is.function(apollo_inputs$apollo_lcPars)
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  HB               = apollo_control$HB
  
  if(HB) stop("The function \'apollo_conditionals\' is not applicable for models estimated using HB!") 
  if(continuous&latentClass)stop("The function \'apollo_conditionals\' is not applicable for models combining continuous mixtures with latent class components!")
  if(!(continuous|latentClass)) stop("The function \'apollo_conditionals\' requires either a model using continuous mixtures or a latent class model!")
  if(continuous) return(apollo_mixConditionals(model, apollo_probabilities, apollo_inputs))
  if(latentClass) return(apollo_lcConditionals(model, apollo_probabilities, apollo_inputs))
}