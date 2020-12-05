#' Calculates conditionals
#' 
#' Calculates posterior expected values (conditionals) of random coefficients, as well as their standard deviations.
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
#' @return List of matrices. Each matrix has dimensions nIndiv x 3. One matrix per random component.
#'         Each row of each matrix contains the indivID of an individual, and the
#'         posterior mean and s.d. of this random component for this individual
#' @export
apollo_conditionals=function(model, apollo_probabilities, apollo_inputs){
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
  
  
  if(is.function(apollo_inputs$apollo_lcPars)) stop("The function \'apollo_conditionals\' is not applicable for models containing latent class components!")
  
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  if(apollo_control$HB) stop("The function \'apollo_conditionals\' is not applicable for models estimated using HB!") 
  
  if(is.null(apollo_control$workInLogs)) apollo_control$workInLogs=FALSE
  if(apollo_control$workInLogs) stop("The function \'apollo_conditionals\' is not applicable for models using the workInLogs setting!") 
  
  if(!apollo_control$mixing) stop("Conditionals can only be estimated for mixture models!")
  if(anyNA(draws)) stop("Random draws have not been specified despite setting mixing=TRUE")
  
  if(apollo_draws$interNDraws==0) stop("This function is only for models that incorporate inter-individual draws!")
  if(apollo_draws$intraNDraws>0) cat("Your model contains intra-individual draws which will be averaged over for conditionals!\n")
  
  
  
  if(!silent) apollo_print("Calculating conditionals...")
  ### Run apollo_randCoeff
  env <- list2env( c(as.list(apollo_beta), apollo_inputs$database, apollo_inputs$draws), 
                   hash=TRUE, parent=parent.frame() )
  environment(apollo_randCoeff) <- env
  randcoeff <- apollo_randCoeff(apollo_beta, apollo_inputs)
  if(any(sapply(randcoeff, is.function))){
    randcoeff = lapply(randcoeff, 
                       function(f) if(is.function(f)){ environment(f) <- env; return(f()) } else { return(f) })
  }
  
  
  ### Get likelihood
  P <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="conditionals")
  
  obsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
  conditionals=list()
  for(j in 1:length(randcoeff)){
    if(length(dim(randcoeff[[j]]))==3) randcoeff[[j]]=colSums(aperm(randcoeff[[j]], perm=c(3,1,2)))/dim(randcoeff[[j]])[3]
    b=randcoeff[[j]]
    b <- rowsum(b, group=database[,apollo_control$indivID], reorder=FALSE)
    b=b/obsPerIndiv
    
    bn=(rowSums(b*P))/rowSums(P)
    bns=sqrt(rowSums(P*(b-bn)^2)/(rowSums(P)))
    conditionals[[names(randcoeff)[j]]]=cbind(unique(database[,apollo_control$indivID]),bn,bns)
    colnames(conditionals[[names(randcoeff)[j]]])=c("ID","post. mean","post. sd")
    rownames(conditionals[[names(randcoeff)[j]]])=c()
  }
  
  if(length(conditionals)==1) conditionals=as.matrix(conditionals[[1]])
  return(conditionals)
}