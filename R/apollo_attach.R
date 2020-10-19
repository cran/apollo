#' Attaches predefined variables.
#'
#' Attaches parameters and data to allow users to refer to individual variables by name without reference to the object that contains them. Also applies scaling if in use.
#'
#' This function should be called at the beginning of \code{apollo_probabilities}
#' to make writing the log-likelihood more user-friendly. If used, then \link{apollo_detach}
#' should be called at the end \code{apollo_probabilities}, or more conveniently, 
#' using \link{on.exit}.
#' \code{apollo_attach} attaches \code{apollo_beta}, \code{database}, \code{draws},
#' and the output of \code{apollo_randCoeff} and \code{apollo_lcPars}, if they are
#' defined by the user. The use of \code{apollo_attach} is mandatory in models using scaling.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return Nothing.
#' @export
apollo_attach=function(apollo_beta, apollo_inputs){
  
  # ############################# #
  #### loads and checks inputs ####
  # ############################# #
  
  test <- ((is.vector(apollo_beta) && is.numeric(apollo_beta)) || is.list(apollo_beta)) || !is.null(names(apollo_beta))
  if(!test) stop("The apollo_beta argument needs to be a named numeric vector or list!")
  
  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars    = apollo_inputs[["apollo_lcPars"]]
  
  # ################################## #
  #### Scale and attach apollo_beta ####
  # ################################## #
  
  #if(!is.null(apollo_inputs$scaling) && !is.na(apollo_inputs$scaling)){
  #  r <- names(apollo_beta) %in% names(apollo_inputs$scaling)
  #  r <- names(apollo_beta)[r]
  #  if(is.list(apollo_beta)){
  #    for(j in 1:length(r)){
  #      apollo_beta[[r[j]]] <- apollo_inputs$scaling[r[j]]*apollo_beta[[r[j]]]  
  #    }
  #  }else{
  #    apollo_beta[r] <- apollo_inputs$scaling[r]*apollo_beta[r]
  #  }
  #}
  attach(as.list(apollo_beta))
  attach(database)
  
  # ################################ #
  #### Build and attach randcoeff ####
  # ################################ #

  if(apollo_control$HB==FALSE && apollo_control$mixing){
    if(anyNA(draws)) stop("Random draws have not been specified despite setting apollo_control$mixing==TRUE!")
    if(!is.function(apollo_randCoeff)) stop("apollo_randCoeff function has not been defined despite setting apollo_control$mixing==TRUE!")
    if("draws" %in% search()) detach("draws")
    attach(draws)
    randcoeff = apollo_randCoeff(apollo_beta, apollo_inputs)
    ### FOLLOWING LINE ADDED IN CASE apollo_randCoeff IS A LIST OF FUNCTIONS 8/05/2020
    if(is.list(randcoeff) && any(sapply(randcoeff, is.function)) && (is.null(apollo_inputs$cpp) || !apollo_inputs$cpp) ){
      randcoeff = lapply(randcoeff, function(f) if(is.function(f)){ return(f()) } else { return(f) })
    } 
    if("randcoeff" %in% search()) detach("randcoeff")
    attach(randcoeff)
  }
  
  # ############################# #
  #### Build and attach lcPars ####
  # ############################# #

  if(is.function(apollo_lcPars)){
    lcpars = apollo_lcPars(apollo_beta, apollo_inputs)
    if("lcpars" %in% search()) detach("lcpars")
    ### If class_specific>0, keep only class_specific
    if(!is.null(apollo_inputs[['class_specific']]) && apollo_inputs$class_specific>0){
      if(is.null(lcpars[['pi_values']])) stop('"apollo_lcPars" should return a list with an element called "pi_values" containing the allocation probabilities for each class')
      nClass <- length(lcpars$pi_values)
      if(!all(sapply(lcpars, is.list))) stop('"apollo_lcPars" should return a list, all of whose elements must be lists as well')
      if(!all(sapply(lcpars,length)==nClass)) stop('"apollo_lcPars" should return a list, all of whose elements must be lists with the same length')
      for(i in 1:length(lcpars)) lcpars[[i]] <- lcpars[[i]][apollo_inputs$class_specific]
    }
    attach(lcpars)
  }
  
}
