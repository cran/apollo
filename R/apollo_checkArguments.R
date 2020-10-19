#' Checks definitions Apollo functions
#' 
#' Checks that the user-defined functions used by Apollo are correctly defined by the user.
#' 
#' It only checks that the functions have the correct definition of inputs. It does not run the functions.
#' 
#' @param apollo_probabilities Function. Likelihood function as defined by the user.
#' @param apollo_randCoeff Function. Defines the random components used inside \code{apollo_probabilities}.
#' @param apollo_lcPars Function. Defines the class allocation probabilities and (optionally) the lists of
#'                      parameters of each class for latent class models.
#' @return Returns (invisibly) TRUE if definitions are correct, and FALSE otherwise.
#' @export
apollo_checkArguments=function(apollo_probabilities=NA,apollo_randCoeff=NA,apollo_lcPars=NA){
  # Check apollo_rpobabilities
  if(is.function(apollo_probabilities)){
    arguments = formals(apollo_probabilities)
    if(!all(names(arguments)==c("apollo_beta", "apollo_inputs", "functionality"))) stop("The arguments for apollo_probabilities need to be apollo_beta, apollo_inputs and functionality")
  } else if(!is.na(apollo_probabilities)) stop("The argument \"apollo_probabilities\" should be a function")
  
  # Check apollo_randCoeff
  if(is.function(apollo_randCoeff)){
    arguments = formals(apollo_randCoeff)
    if(!all(names(arguments)==c("apollo_beta", "apollo_inputs"))) stop("The arguments for apollo_randCoeff need to be apollo_beta and apollo_inputs")
  } else if(!is.na(apollo_randCoeff)) stop("The argument \"apollo_randCoeff\" should be a function")
  
  # Check apollo_lcPars
  if(is.function(apollo_lcPars)){
    arguments = formals(apollo_lcPars)
    if(!all(names(arguments)==c("apollo_beta", "apollo_inputs"))) stop("The arguments for apollo_lcPars need to be apollo_beta and apollo_inputs")
  } else if(!is.na(apollo_lcPars)) stop("The argument \"apollo_lcPars\" should be a function")
  
  return(invisible(TRUE))
}
