apollo_checkArguments=function(apollo_probabilities=NA,apollo_randCoeff=NA,apollo_lcPars=NA){
if(is.function(apollo_probabilities)){
  arguments=formals(apollo_probabilities)
  if(!all(names(arguments)==c("apollo_beta", "apollo_inputs", "functionality"))) stop("The arguments for apollo_probabilities need to be apollo_beta, apollo_inputs and functionality")
}
if(is.function(apollo_randCoeff)){
  arguments=formals(apollo_randCoeff)
  if(!all(names(arguments)==c("apollo_beta", "apollo_inputs"))) stop("The arguments for apollo_randCoeff need to be apollo_beta and apollo_inputs")
}
if(is.function(apollo_lcPars)){
  arguments=formals(apollo_lcPars)
  if(!all(names(arguments)==c("apollo_beta", "apollo_inputs"))) stop("The arguments for apollo_lcPars need to be apollo_beta and apollo_inputs")
}
return(invisible(TRUE))
}
