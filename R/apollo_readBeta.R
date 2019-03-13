#' Reads parameters from file
#'
#' Reads in parameters from a previously estimated model and copies the values to the given apollo_beta vector, only for those parameters whose name matches.
#'
#' This function will update the values of the parameters in its argument \code{apollo_beta} with the matching values in the file
#' \code{(inputModelName)_estimates.csv}. If there is no match for a given parameter in \code{apollo_beta}, its value will not be updated.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param inputModelName Character. modelName for model from which results are used as starting values.
#' @param overwriteFixed Boolean. TRUE if starting values for fixed parameters should also be updated from input file.
#' @return Named numeric vector. Names and updated starting values for parameters.
#' @export
#' @examples
#' ### Define starting values and fixed parameters
#' apollo_beta <- c(asc1=0, asc2=0, b1=0, b2=0)
#' apollo_fixed <- c("asc1")
#' \dontrun{
#' ## Not run:
#' ### Update starting values
#' apollo_beta <- apollo_readBeta(apollo_beta, apollo_fixed, 
#'                                "oldModelName", overwriteFixed=FALSE)
#' ## End(Not run)
#' }
apollo_readBeta=function(apollo_beta, apollo_fixed, inputModelName, overwriteFixed=FALSE){
  
  tmp <- paste0(inputModelName, "_estimates.csv")
  input_apollo_beta = tryCatch(utils::read.csv(tmp), error=function(e) stop("Could not load ", tmp))

  if(overwriteFixed!=FALSE) overwriteFixed=TRUE

  input_pars=input_apollo_beta[,2]
  names(input_pars)=input_apollo_beta[,1]
  common_pars=names(input_pars)[names(input_pars) %in% names(apollo_beta)]
  common_fixed_pars=names(input_pars)[names(input_pars) %in% apollo_fixed]
  common_var_pars=common_pars[!(common_pars %in% apollo_fixed)]

  if((overwriteFixed==FALSE)&(length(apollo_fixed)>0)){
    if(length(common_var_pars)==0) cat("\nNo parameter names in the input file match those of parameters to be estimated in apollo_beta!")
    if(length(common_var_pars)==1) cat("\nOut of the",length(apollo_beta),"parameters in apollo_beta, 1 parameter was updated with the value from the input file.")
    if(length(common_var_pars)>1) cat("\nOut of the",length(apollo_beta),"parameters in apollo_beta,",length(common_var_pars),"were updated with values from the input file.")
    if(length(common_fixed_pars)==1) cat("\n1 parameter in apollo_beta was kept fixed at its starting value rather than being updated from the input file.")
    if(length(common_fixed_pars)>1) cat("\n",length(common_fixed_pars),"parameters in apollo_beta were kept fixed at their starting values rather than being updated from the input file.")
    apollo_beta[common_var_pars]=input_pars[common_var_pars]
  } else {
    if(length(common_pars)==0) cat("\nNo parameter names in the input file match those in apollo_beta!")
    if(length(common_pars)==1) cat("\nOut of the",length(apollo_beta),"parameters in apollo_beta, 1 parameter was updated with the value from the input file.")
    if(length(common_pars)>1) cat("\nOut of the",length(apollo_beta),"parameters in apollo_beta,",length(common_pars),"were updated with values from the input file.")
    if(length(common_fixed_pars)==1) cat("\nThis includes 1 parameter in apollo_fixed whose value has been set to that from the input file.")
    if(length(common_fixed_pars)>1) cat("\nThis includes",length(common_fixed_pars),"parameters in apollo_fixed whose values have been set to those from the input file.")
    apollo_beta[common_pars]=input_pars[common_pars]
  }
  return(apollo_beta)

}
