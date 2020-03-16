#' Predicts using an estimated model
#' 
#' Calculates apollo_probabilities with functionality="prediction" and extracts one element from the returned list.
#' 
#' Structure of predictions are simplified before returning, e.g. list of vectors are turned into a matrix.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param prediction_settings List of settings. It can have the following elements.
#'                            \itemize{
#'                              \item \code{modelComponent} Character. Name of component of apollo_probabilities output to calculate predictions for. Default is "model", i.e. the whole model.
#'                              \item \code{silent} Boolean. If TRUE, thsi function won't print any output to screen.
#'                            }
#' @param modelComponent \strong{Deprecated}. Same as \code{modelComponent} inside \code{prediction_settings}.
#' @return A list containing predictions for component \code{modelComponent} of the model described in \code{apollo_probabilities}.
#'         The particular shape of the prediction will depend on the model component.
#' @export
apollo_prediction <- function(model, apollo_probabilities, apollo_inputs, prediction_settings=list(), modelComponent="model"){
  
  if(is.character(prediction_settings) && is.character(modelComponent)){
    modelComponent=prediction_settings
    prediction_settings=list(runs=1)
  }
  
  if(is.character(prediction_settings) && is.list(modelComponent)){
    modelComponent1=prediction_settings
    prediction_settings=modelComponent
    modelComponent=modelComponent1
  }

  if(!exists("prediction_settings")) prediction_settings=list()
  
  if(is.null(prediction_settings$modelComponent)){
    if(exists("modelComponent")){
      prediction_settings$modelComponent=modelComponent
    } else {
      prediction_settings$modelComponent="model"
    }
  }
  if(is.null(prediction_settings$runs)) prediction_settings$runs=1
  if(is.null(prediction_settings$silent)) prediction_settings$silent=FALSE
  
  modelComponent = prediction_settings$modelComponent
  runs           = prediction_settings$runs
  
  if(!is.null(model$apollo_control$HB) && model$apollo_control$HB && runs>1) stop("The calculation of confidence intervals for \'apollo_prediction\' is not applicables for models estimated using HB!") 
  
  #### NEW 16/03
  HB=!is.null(model$apollo_control$HB) && model$apollo_control$HB
  if(!HB && (!is.numeric(model$estimate) | !is.vector(model$estimate) | is.null(names(model$estimate))))  stop("The \'model$estimates\' object should be a named vector for models not estimated using HB!")   
  if(HB){
    nObs  <- nrow(apollo_inputs$database)
    if(any(!(lengths(model$estimate)%in%c(1,nObs)))) stop("\nFor models estimated using HB, the \'model$estimate\' object needs to have one entry per row in the database for each random parameter, and a single value for each fixed parameter. This is not the case here, likely because the prediction data is different from the estimation data. Please restructure \'model$estimate\' accordingly!")
  }
  
  ####
  
  cat("Updating inputs...")
  apollo_inputs <- apollo_validateInputs(silent=TRUE)
  cat("Done.\n")
  
  apollo_randCoeff  = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars     = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  
  if(!prediction_settings$silent) cat("\nRunning predictions from model using parameter estimates")
  apollo_beta  = model$estimate
  apollo_fixed = model$apollo_fixed
  predictions  = apollo_probabilities(apollo_beta, apollo_inputs, functionality="prediction")
  if(is.null(predictions[[modelComponent]])) stop("A likelihood element with the selected name does not exist in your model!")
  predictions=predictions[[modelComponent]]
  if(length(predictions)==1 && is.na(predictions)){
    cat("\nPredictions do not exist for the selected model component.")
    return(NA)
  }
  cat(" Done.\n")
  if(is.list(predictions)){
    predictions = matrix(unlist(predictions), ncol=length(predictions), byrow=FALSE, 
                         dimnames=list(NULL, names(predictions))) 
  }
  predictions = cbind(ID=apollo_inputs$database[,apollo_inputs$apollo_control$indivID],
                      `Choice situation`=apollo_inputs$database$apollo_sequence,
                      predictions)   #### NEW
  
  if(runs>1){
    covar=model$robvarcov
    beta_draws = mvtnorm::rmvnorm(n=runs, mean = apollo_beta[!(names(apollo_beta)%in%apollo_fixed)], sigma = covar)
    fixed_ones = matrix(apollo_beta[(names(apollo_beta)%in%apollo_fixed)],nrow=runs,ncol=length(apollo_fixed),byrow=TRUE)
    colnames(fixed_ones)=apollo_fixed
    beta_draws = cbind(beta_draws,fixed_ones)
    #beta_draws = cbind(beta_draws,apollo_beta[(names(apollo_beta)%in%apollo_fixed)])
    #colnames(beta_draws)=c(names(apollo_beta)[!(names(apollo_beta)%in%apollo_fixed)],names(apollo_beta)[(names(apollo_beta)%in%apollo_fixed)])
    
    cat("\nRunning predictions from vector of model estimates")
    predictions_random=list()
    r=1
    while(r<=runs){
    cat("\n Set of draws ", r, "/", runs, " ", sep="")
      predictions_r  = apollo_probabilities(beta_draws[r,], apollo_inputs, functionality="prediction")
      predictions_r=predictions_r[[modelComponent]]
      cat(" Done.")
      if(is.list(predictions_r)){
        predictions_r = matrix(unlist(predictions_r), ncol=length(predictions_r), byrow=FALSE, 
                             dimnames=list(NULL, names(predictions_r))) 
      }
      predictions_random[[r]] = cbind(ID=apollo_inputs$database[,apollo_inputs$apollo_control$indivID],
                          `Choice situation`=apollo_inputs$database$apollo_sequence,
                          predictions_r)
      r=r+1
    }
    predictions=list(at_estimates=predictions,draws=predictions_random)
  }
  
  if(!prediction_settings$silent){
    cat("\nPredicted demand at model estimates\n")
  if(is.list(predictions)){
   print(colSums(predictions[["at_estimates"]][,-c(1,2,ncol(predictions[["at_estimates"]]))],na.rm=TRUE))
  } else {
   print(colSums(predictions[,-c(1,2,ncol(predictions))],na.rm=TRUE))  
  }
  if(runs>1){
    cat("\nStandard deviation for predicted demand across runs\n")
    #tmp <- sapply(predictions$draws, colSums)
    tmp <- sapply(predictions$draws, function(x) colSums(x,na.rm=TRUE))
    tmp <- apply(tmp, MARGIN=1, sd)
    print(tmp[-c(1,2,ncol(tmp))])
  }}
  return(predictions)
}