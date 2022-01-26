#' Calculates gradients of utility functions
#' 
#' Calculates gradients (derivatives) of utility functions.
#' 
#' @param apollo_beta Named numeric vector of parameters.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param V List of functions
#' @return Named list. Each element is itself a list of functions: the partial derivatives of the elements of V.
#' 
#' @export
apollo_dVdB <- function(apollo_beta, apollo_inputs, V){
  # Useful variables
  freeparams <- !(names(apollo_beta) %in% apollo_inputs$apollo_fixed)
  freeparams <- names(apollo_beta[freeparams])
  J <- length(V)
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) ) return(TRUE) else return(FALSE)
  
  # Get apollo_probabilities
  test <- is.list(apollo_inputs) && is.function(apollo_inputs$apollo_probabilities)
  if(test) apollo_probabilities <- apollo_inputs$apollo_probabilities else {
    apollo_probabilities <- tryCatch(get("apollo_probabilities", envir=parent.frame(), inherits=TRUE),
                                     error=function(e) return(NULL))
    # If V contains scaling, but apollo_probability doesn't, then scale apollo_probabilities
    test <- !is.null(apollo_inputs$apollo_scaling)
    test <- test &&  ('apollo_scaling' %in% unlist(sapply(V, function(v) all.vars(body(v)), USE.NAMES=FALSE), use.names=FALSE))
    test <- test && !('apollo_scaling' %in% all.vars(body(apollo_probabilities)))
    if(test) apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
  }
  
  # Extract all variables used in apollo_probabilities
  if(is.null(apollo_probabilities)) return(NULL)
  if(apollo_inputs$apollo_control$debug) apollo_print("dVdB: Extracting variable definitions")
  defs     <- tryCatch(apollo_varList(apollo_probabilities, apollo_inputs), error=function(e) NULL)
  if(is.null(defs)){
    if(apollo_inputs$apollo_control$debug) apollo_print('dVdB: Could not extract variable definitions')
    return(NULL)
  }
  
  # Check there are no indices
  #if(length(grep("[", vars$v[,2], fixed=TRUE))>0){
  #  if(apollo_inputs$apollo_control$debug) apollo_print("dVdB: Analytical gradients cannot be used if utilities contain indices.")
  #  return(NULL)
  #}
  
  # Replace values defined elsewhere
  replaceByDef <- function(e, defs){
    isFunction <- is.function(e)
    if(isFunction){f <- e; e <- body(e)}
    # Case 1: x
    test1 <- is.symbol(e) && (as.character(e) %in% names(defs))
    if(test1) e <- defs[[as.character(e)]]
    # Case 2: L$x or L[['x']] or L[["x"]]
    test2 <- !test1 && is.call(e) && length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('$','[['))
    test2 <- test2 && is.val(e[[3]]) && is.symbol(e[[2]])
    if(test2) tmp  <- paste0(as.character(e[[2]]), '$', as.character(e[[3]]))
    test2 <- test2 && (tmp %in% names(defs))
    if(test2) e <- defs[[tmp]]
    # Case 3: expression
    if( !test1 && !test2 && (is.call(e) || is.expression(e)) ){
      for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- replaceByDef(e[[i]], defs)
    } 
    # Return
    if(isFunction){body(f) <- e; e <- f}
    return(e)
  }
  V <- lapply(V, replaceByDef, defs=defs)
  
  # Symbolic differentiation
  if(apollo_inputs$apollo_control$debug) apollo_print("dVdB: Calculating analytical derivatives")
  dV <- lapply(V, Deriv::Deriv, x=freeparams, combine="list")
  for(k in 1:length(dV)) environment(dV[[k]]) <- new.env(hash=TRUE, parent=baseenv())
  return(dV)
}

#f <- function() b1*x1 + b2*x2
#body(f)[[2]][[2]]
#bf <- body(f)
#bf[[2]][[2]] <- quote(a0 + a1*z1)

#Deriv::Deriv()