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
#' @return List of object, one per random component and one for the class allocation probabilities.
#' @export
apollo_lcUnconditionals <- function(model, apollo_probabilities, apollo_inputs){
  apollo_beta=model$estimate
  apollo_fixed=model$apollo_fixed

  cat("Updating inputs...")
  apollo_inputs <- apollo_validateInputs(silent=TRUE)
  cat("Done.\n")

  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars    = apollo_inputs[["apollo_lcPars"]]
  apollo_draws     = apollo_inputs[["apollo_draws"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  toAttach  <- c(as.list(apollo_beta), apollo_inputs$database)
  if(apollo_control$mixing){
    toAttach  <- c(as.list(apollo_beta), apollo_inputs$database, draws)
    randcoeff = with(toAttach, {
      environment(apollo_randCoeff) <- environment()
      apollo_randCoeff(apollo_beta, apollo_inputs)
    } )
    if(apollo_draws$intraNDraws>0) cat("Your model contains intra-individual draws, so the output will contain draws across individuals and across choices.\n")
    if(apollo_draws$intraNDraws==0) cat("Your model contains only inter-individual draws, so the output will contain draws across individuals only.\n")
  }
  unconditionals = with(toAttach, {
    environment(apollo_lcPars) <- environment()
    apollo_lcPars(apollo_beta, apollo_inputs)
  } )
  cat("Unconditional distributions computed\n")
  return(unconditionals)
}
