#' Returns unconditionals for a latent class model model
#'
#' Returns values for random parameters and class allocation probabilities in a latent class model model.
#'
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return List of object, one per random component and one for the class allocation probabilities.
#' @export
apollo_lcUnconditionals <- function(model, apollo_probabilities, apollo_inputs){
  if(!is.function(apollo_inputs$apollo_lcPars)) stop("SYNTAX ISSUE - This function is for latent class models. For other models use \"apollo_unconditionals\".")
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
  apollo_lcPars    = apollo_inputs[["apollo_lcPars"]]
  apollo_draws     = apollo_inputs[["apollo_draws"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  if(apollo_control$HB) stop("INCORRECT FUNCTION/SETTING USE - The function \'apollo_lcUnconditionals\' is not applicables for models estimated using HB!") 
  
  # Calculate randCoeff if necessary
  toAttach  <- c(as.list(apollo_beta), apollo_inputs$database)
  if(apollo_control$mixing){
    toAttach  <- c(as.list(apollo_beta), apollo_inputs$database, draws)
    randcoeff = with(toAttach, {
      environment(apollo_randCoeff) <- environment()
      apollo_randCoeff(apollo_beta, apollo_inputs)
    } )
    if(apollo_draws$intraNDraws>0) cat("Your model contains intra-individual draws, so the output will contain draws across individuals and across choices.\n")
    if(apollo_draws$intraNDraws==0) cat("Your model contains only inter-individual draws, so the output will contain draws across individuals only.\n")
    toAttach <- c(toAttach, randcoeff)
  }
  
  # Calculate lcPars
  unconditionals = with(toAttach, {
    environment(apollo_lcPars) <- environment()
    apollo_lcPars(apollo_beta, apollo_inputs)
  } )
  
  # If there are no intraDraws, keep only first row of each individual
  nObs <- nrow(database)
  if(any(!is.na(apollo_draws)) && apollo_draws$intraNDraws==0) for(i in 1:length(unconditionals)){
    if(is.list(unconditionals[[i]])){
      for(j in 1:length(unconditionals[[i]])){
        x <- unconditionals[[i]][[j]]
        if(is.array(x)) nRows <- dim(x)[1] else nRows <- length(x)
        if(nRows==nObs) unconditionals[[i]][[j]] <- apollo_firstRow(x, apollo_inputs)
      }
    } else {
      x <- unconditionals[[i]]
      if(is.array(x)) nRows <- dim(x)[1] else nRows <- length(x)
      if(nRows==nObs) unconditionals[[i]] <- apollo_firstRow(x, apollo_inputs)
    }
  }
  
  if(!silent) apollo_print("Unconditional distributions computed")
  return(unconditionals)
}
