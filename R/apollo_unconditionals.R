#' Returns unconditionals for models with random heterogeneity
#'
#' Returns unconditionals for random parameters in model, both for continuous mixtures and latent class.
#'
#' This functions is only meant for use with models using continuous distributions or latent classes, or both at the same time.
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
#'           \item If the model contains a continuous mixture, it returns a list with one object per
#'                 random coefficient. When using inter-individual draws only, each element will be 
#'                 a matrix with one row per individual, and one column per draw. When using intra-
#'                 individual draws, each element will be a three-dimensional array, with one row 
#'                 per observation, inter-individual draws in the second dimension, and intra-
#'                 individual draws in the third dimension.
#'           \item If the model contains latent classes, it returns a list with as many elements
#'                 as random coefficients in the model, plus one additional element containing 
#'                 the class allocation probabilities.
#'           \item If the model contains both continuous mixing and latent classes, a list with the
#'                 two elements described above will be returned.
#'         }
#' @export
apollo_unconditionals <- function(model, apollo_probabilities, apollo_inputs){
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
  
  if(HB) stop("INCORRECT FUNCTION/SETTING USE - The function \'apollo_unconditionals\' is not applicable for models estimated using HB!") 
  #if(continuous & latentClass) stop("INCORRECT FUNCTION/SETTING USE - The function \'apollo_unconditionals\' is not applicable for models combining continuous mixtures with latent class components!")
  if(!(continuous|latentClass)) stop("INCORRECT FUNCTION/SETTING USE - The function \'apollo_unconditionals\' requires either a model using continuous mixtures and/or a latent class model!")
  if(continuous & latentClass){
    ans <- list()
    ans[["continuous" ]] <- apollo_mixUnconditionals(model, apollo_probabilities, apollo_inputs)
    ans[["latentClass"]] <- apollo_lcUnconditionals(model, apollo_probabilities, apollo_inputs)
    return(ans)
  }
  if(continuous) return(apollo_mixUnconditionals(model, apollo_probabilities, apollo_inputs))
  if(latentClass) return(apollo_lcUnconditionals(model, apollo_probabilities, apollo_inputs))
  
}