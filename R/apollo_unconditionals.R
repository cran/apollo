#' Returns draws for random parameters in a latent class model model
#'
#' Returns draws (unconditionals) for random parameters in model, including interactions with deterministic covariates
#'
#' This functions is only meant for use with continuous distributions
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return List of object, one per random coefficient.
#'         With inter-individual draws only, this will be a matrix, with one row per individual, and one column per draw.
#'         With intra-individual draws, this will be a three-dimensional array, with one row per observation, inter-individual draws in the second dimension, and intra-individual draws in the third dimension.
#' @export
apollo_unconditionals <- function(model, apollo_probabilities, apollo_inputs){
  
  if(is.null(apollo_inputs$silent)) silent = FALSE else silent = apollo_inputs$silent
  apollo_beta  = model$estimate
  apollo_fixed = model$apollo_fixed
  
  #if(!silent) apollo_print("Updating inputs...")
  #apollo_inputs <- apollo_validateInputs(silent=TRUE, recycle=TRUE)
  ### Warn the user in case elements in apollo_inputs are different from those in the global environment
  apollo_compareInputs(apollo_inputs)
  
  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_draws     = apollo_inputs[["apollo_draws"]]
  apollo_lcPars     = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  if(apollo_control$HB) stop("The function \"apollo_unconditionals\" is not applicables for models estimated using HB!") 
  if(is.function(apollo_inputs$apollo_lcPars)) stop("The function \"apollo_unconditionals\" is not applicables for models containing latent class components!")
  
  ### Validate input
  if(!apollo_control$mixing) stop("Sample level random parameters can only be produced for mixture models!")
  if(anyNA(draws)) stop("Random draws have not been specified despite setting mixing=TRUE")
  
  ### Run apollo_randCoeff
  env <- list2env( c(as.list(apollo_beta), apollo_inputs$database, apollo_inputs$draws), 
                   hash=TRUE, parent=parent.frame() )
  environment(apollo_randCoeff) <- env
  randcoeff <- apollo_randCoeff(apollo_beta, apollo_inputs)
  if(any(sapply(randcoeff, is.function))){
    randcoeff = lapply(randcoeff, 
                       function(f) if(is.function(f)){ environment(f) <- env; return(f()) } else { return(f) })
  }
  
  if(apollo_draws$intraNDraws==0){
    nObsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
    nIndiv       <- length(nObsPerIndiv)
    firstRows    <- rep(1, nIndiv)
    for(i in 2:nIndiv) firstRows[i] <- firstRows[i-1] + nObsPerIndiv[i-1]
    j=1
    for(j in 1:length(randcoeff)){
      randcoeff[[j]]=randcoeff[[j]][firstRows,]  
    }
  }
  
  if(!silent) apollo_print("Unconditional distributions computed") 
  return(randcoeff)
}