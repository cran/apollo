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
