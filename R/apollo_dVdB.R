#' Calculates gradients of utility functions
#' 
#' Calculates gradients (derivatives) of utility functions, considering definitions in apollo_randCoeff.
#' 
#' @param apollo_beta Named numeric vector of parameters.
#' @param apollo_inputs List of Apollo main settings.
#' @param V List of functions
#' @return Named list. Each element is itself a list of functions: the partial derivatives of the elements of V.
#' 
#' @export
apollo_dVdB <- function(apollo_beta, apollo_inputs, V){
  freeparams <- !(names(apollo_beta) %in% apollo_inputs$apollo_fixed)
  freeparams <- names(apollo_beta[freeparams])
  J <- length(V)
  
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
  vars     <- apollo_varList(apollo_probabilities, apollo_beta, apollo_inputs, V)
  if(anyNA(vars)) return(NULL)
  
  # Check there are no indices
  #if(length(grep("[", vars$v[,2], fixed=TRUE))>0){
  #  if(apollo_inputs$apollo_control$debug) apollo_print("dVdB: Analytical gradients cannot be used if utilities contain indices.")
  #  return(NULL)
  #}
  
  # Replace values defined elsewhere
  baseVars <- c(vars$b, vars$x, vars$d, 'apollo_inputs', 'apollo_scaling')
  vVars    <- unlist(sapply(V, function(v) all.vars(body(v)), USE.NAMES=FALSE), use.names=FALSE)
  replaceVars <- function(x_old, x_new, expre){
    if(length(expre)==1){
      if(expre==x_old) return(str2lang(x_new))
    } else for(i in 1:length(expre)) expre[[i]] <- replaceVars(x_old, x_new, expre[[i]])
    return(expre)
  }
  while(!all(vVars %in% baseVars)){
    vVars0 <- vVars
    for(j in 1:J){
      varsF <- all.vars(body(V[[j]]))
      if(!is.null(vars$r)) for(k in 1:nrow(vars$r)){
        if(vars$r[k,1] %in% varsF) body(V[[j]]) <- replaceVars(vars$r[k,1], vars$r[k,2], body(V[[j]]))
      }
      if(!is.null(vars$p)) for(k in 1:nrow(vars$p)){
        if(vars$p[k,1] %in% varsF) body(V[[j]]) <- replaceVars(vars$p[k,1], vars$p[k,2], body(V[[j]]))
      }
      vVars <- unlist(sapply(V, function(v) all.vars(body(v)), USE.NAMES=FALSE), use.names=FALSE)
    }
    if(all(vVars0 %in% vVars)){
      if(apollo_inputs$apollo_control$debug) apollo_print("dVdB: Unknown variables inside V")
      return(NULL)
    }
  }
  
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