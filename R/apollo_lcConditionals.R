#' Calculates conditionals for latent class models.
#' 
#' Calculates posterior expected values (conditionals) of class allocation probabilities for each individual.
#' 
#' This function can only be used with latent class models without continuous heterogeneity.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return A matrix with the posterior class allocation probabilities for each individual.
#' @export
apollo_lcConditionals=function(model, apollo_probabilities, apollo_inputs){
  if(!is.function(apollo_inputs$apollo_lcPars)) stop("SYNTAX ISSUE - This function is for latent class models. For other models use \"apollo_conditionals\".")
  
  if(is.null(apollo_inputs$silent)) silent = FALSE else silent = apollo_inputs$silent
  apollo_beta  = model$estimate
  apollo_fixed = model$apollo_fixed
  
  #if(!silent) apollo_print("Updating inputs...")
  #apollo_inputs <- apollo_validateInputs(silent=TRUE, recycle=TRUE)
  ### Warn the user in case elements in apollo_inputs are different from those in the global environment
  apollo_compareInputs(apollo_inputs)
  
  apollo_control = apollo_inputs[["apollo_control"]]
  database       = apollo_inputs[["database"]]
  apollo_lcPars  = apollo_inputs[["apollo_lcPars"]]
  class_prob     = "pi_values" # name of lcpars component with allocation probabilities
  apollo_randCoeff  = apollo_inputs[["apollo_randCoeff"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  if(apollo_control$HB) stop("INCORRECT FUNCTION/SETTING USE - The function \'apollo_lcConditionals\' is not applicables for models estimated using HB!") 

  ### Validation
  if(apollo_control$mixing) stop("SYNTAX ISSUE - apollo_lcConditionals can only be used for latent class models without continuous random heterogeneity")
  
  if(!silent) apollo_print("Calculating conditionals...")
  ### Get allocation and inClass probs
  lcpars = with(c(apollo_beta, apollo_inputs$database, apollo_inputs$draws), {
    environment(apollo_lcPars) <- environment()
    apollo_lcPars(apollo_beta, apollo_inputs)
  })
  if(is.null(lcpars[[class_prob]])) stop("SYNTAX ISSUE - The lcpars function needs to create an object called \"pi_values\"!")
  L      = apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  classes    = length(lcpars[[class_prob]])
  components = length(L)
  if(components>(classes+1)) stop("SYNTAX ISSUE - apollo_lcConditionals can only be used for latent class models alone (i.e. no hybrid choice)")
  if(components!=(classes+1)) stop("SYNTAX ISSUE - Model should contain one component per class, and an overall model!")
  
  ### Calculate posterior class allocation probs
  post_pi = vector(mode="list", length=classes)
  nObsPerIndiv <- apollo_inputs$database[,apollo_inputs$apollo_control$indivID]
  nObsPerIndiv <- as.vector(table(nObsPerIndiv))
  for(s in 1:classes){
    # adjust dimensionality of pi if necessary
    pi <- lcpars[[class_prob]][[s]]
    rL <- ifelse(is.array( L[[s]]), nrow( L[[s]]), length( L[[s]]))
    rP <- ifelse(is.array(pi), nrow(pi), length(pi))
    if(rP!=1 && rL!=1 && rP!=rL){
      if(rP>rL) pi <- apollo_firstRow(pi, apollo_inputs)
      if(rP<rL && is.vector(pi)) pi <- rep(pi, each=nObsPerIndiv)
    }
    if(is.list(L[[s]])){
      if(!is.null(L[[s]]$model)){
        L[[s]]=L[[s]]$model
      }else{
        stop("SPECIFICATION ISSUE: the within class probabilities are lists that do not contain an entry called model!")
      }
    }
      post_pi[[s]] = pi*L[[s]]/L[["model"]]
  }; rm(pi, rL, rP)
  
  ### Prepare output
  conditionals = matrix(unlist(post_pi), ncol = length(post_pi), byrow = FALSE)
  classnames   = paste("Class ",seq(1:classes),sep="")
  conditionals = data.frame(ID=unique(database[,apollo_control$indivID]), conditionals)
  
  return(conditionals)
}