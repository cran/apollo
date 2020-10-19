#' Adds covariance matrix to Apollo model
#' 
#' Receives an estimated model object, calculates its Hessian, and classical and robust covariance matrix, and returns the  
#' same model object, but with these additional elements.
#' @param model A model object, as returned by \link{apollo_estimate}
#' @param apollo_inputs List of settings. as returned by \link{apollo_validateInputs}
#' @return model.
#' @export
apollo_addCovariance <- function(model, apollo_inputs){
  #model=apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=list(covarOnly=TRUE))
  L <- apollo_varcov(model$estimate, model$apollo_fixed, list(apollo_probabilities = model$apollo_probabilities,
                                                              apollo_inputs        = apollo_inputs))
  # only use those things that have been updated, keep all the rest from the original model
  for(i in names(L)) model[[i]] <- L[[i]]
  ### we should probably also update the time the model was run, and change the post-estimation time
  
  return(model)
}
