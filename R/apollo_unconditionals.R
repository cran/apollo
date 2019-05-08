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
#'         With inter-individual draws only, this will be a matrix, with one row per observation, and one column per draw.
#'         With intra-individual draws, this will be a three-dimensional array, with one row per observation, inter-individual draws in the second dimension, and intra-individual draws in the third dimension.
#' @export
apollo_unconditionals <- function(model, apollo_probabilities, apollo_inputs){
  
  apollo_beta=model$estimate
  apollo_fixed=model$apollo_fixed
  
  cat("Updating inputs...")
  apollo_inputs <- apollo_validateInputs(silent=TRUE)
  cat("Done.\n")
  
  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_draws     = apollo_inputs[["apollo_draws"]]
  
  ### Validate input
  if(!apollo_control$mixing) stop("Sample level random parameters can only be produced for mixture models!")
  if(anyNA(draws)) stop("Random draws have not been specified despite setting mixing=TRUE")
  
  toAttach  <- c(as.list(apollo_beta), apollo_inputs$database, apollo_inputs$draws)
  
  randcoeff = with(toAttach, {
    environment(apollo_randCoeff) <- environment()
    apollo_randCoeff(apollo_beta, apollo_inputs)
  } )
  
  if(apollo_draws$intraNDraws==0){
    nObsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
    nIndiv       <- length(nObsPerIndiv)
    firstRows    <- rep(1, nIndiv)
    for(i in 2:nIndiv) firstRows[i] <- firstRows[i-1] + nObsPerIndiv[i-1]
    j=1
    while(j<(length(randcoeff)+1)){
      randcoeff[[j]]=randcoeff[[j]][firstRows,]  
      j=j+1    
    }
  }
  
  cat("Unconditional distributions computed\n") 
  return(randcoeff)
}