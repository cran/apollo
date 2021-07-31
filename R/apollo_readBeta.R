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
apollo_readBeta=function(apollo_beta, apollo_fixed, inputModelName, overwriteFixed=FALSE){
  
  ### Search for output directory
  outputDirectory <- ''
  tmp  <- tryCatch(get('apollo_control', envir=parent.frame(), inherits=FALSE), error=function(e) FALSE)
  test <- is.list(tmp) && !is.null(tmp$outputDirectory) && is.character(tmp$outputDirectory)
  if(test) outputDirectory <- tmp$outputDirectory else {
    tmp <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE), error=function(e) FALSE)
    test <- is.list(tmp) && !is.null(tmp$apollo_control) && !is.null(tmp$apollo_control$outputDirectory)
    test <- test && is.character(tmp$apollo_control$outputDirectory)
    if(test) outputDirectory <- tmp$apollo_control$outputDirectory
  }
  test <- outputDirectory!=''
  test <- test && !(substr(outputDirectory, nchar(outputDirectory), nchar(outputDirectory)) %in% c('/','\\'))
  if(test) outputDirectory <- paste0(outputDirectory, '/')
  
  ### Try to read from file
  filename = paste0(outputDirectory, inputModelName, "_estimates.csv", collapse='')
  if(!file.exists(filename)) filename = paste0(inputModelName,"_estimates.csv", collapse='')
  if(!file.exists(filename)){
    if(outputDirectory=='') stop('File ', filename, ' not found in working directory.') else 
      stop('File ', filename, ' not found in working directory, nor in ', outputDirectory, '.') 
  }
  input_apollo_beta = tryCatch(utils::read.csv(filename), 
                               warning=function(w) x=FALSE,
                               error=function(e) x=FALSE)
  if(is.logical(input_apollo_beta) && input_apollo_beta==FALSE) stop("Could not open file ",filename) 
  
  if(overwriteFixed!=FALSE) overwriteFixed=TRUE

  input_pars=input_apollo_beta[,2]
  names(input_pars)=input_apollo_beta[,1]
  common_pars=names(input_pars)[names(input_pars) %in% names(apollo_beta)]
  common_fixed_pars=names(input_pars)[names(input_pars) %in% apollo_fixed]
  common_var_pars=common_pars[!(common_pars %in% apollo_fixed)]

  if((overwriteFixed==FALSE)&(length(apollo_fixed)>0)){
    if(length(common_var_pars)==0) apollo_print("No parameter names in the input file match those of parameters to be estimated in apollo_beta!")
    if(length(common_var_pars)==1) apollo_print(paste0("Out of the ",length(apollo_beta)," parameters in apollo_beta, 1 parameter was updated with the value from the input file."))
    if(length(common_var_pars)>1) apollo_print(paste0("Out of the ",length(apollo_beta)," parameters in apollo_beta, ",length(common_var_pars)," were updated with values from the input file."))
    if(length(common_fixed_pars)==1) apollo_print("1 parameter in apollo_beta was kept fixed at its starting value rather than being updated from the input file.")
    if(length(common_fixed_pars)>1) apollo_print(paste0(length(common_fixed_pars)," parameters in apollo_beta were kept fixed at their starting values rather than being updated from the input file."))
    apollo_beta[common_var_pars]=input_pars[common_var_pars]
  } else {
    if(length(common_pars)==0) apollo_print("No parameter names in the input file match those in apollo_beta!")
    if(length(common_pars)==1) apollo_print(paste0("Out of the ",length(apollo_beta)," parameters in apollo_beta, 1 parameter was updated with the value from the input file."))
    if(length(common_pars)>1) apollo_print(paste0("Out of the ",length(apollo_beta)," parameters in apollo_beta, ",length(common_pars)," were updated with values from the input file."))
    if(length(common_fixed_pars)==1) apollo_print("This includes 1 parameter in apollo_fixed whose value has been set to that from the input file.")
    if(length(common_fixed_pars)>1) apollo_print(paste0("This includes ",length(common_fixed_pars)," parameters in apollo_fixed whose values have been set to those from the input file."))
    apollo_beta[common_pars]=input_pars[common_pars]
  }
  cat("\n")
  return(apollo_beta)

}
