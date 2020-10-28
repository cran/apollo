#' Scales variables inside a function
#' 
#' It changes the syntax of the function by replacing variable names for their scaled form, 
#' e.g. x --> x*apollo_inputs$apollo_scale[["x"]]. In assignments, it only scales the
#' right side of the assignment.
#' 
#' @param e Function, expression, call or symbol to alter.
#' @param sca Named numeric vector with the scales. The names in these vectors determine which variables should be scaled.
#' @return A function, expression, call or symbol with the corresponding variables scaled.
#' @export
apollo_insertScaling <- function(e, sca){
  # Validate input
  test <- is.function(e) || is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) || is.call(e)
  if(!test) stop('Argument "e" must be a function, a call, a symbol, or a value')
  if(is.function(e)){eOrig <- e; e <- body(e)} else eOrig <- NULL
  if(!is.vector(sca) || !is.numeric(sca) || is.null(names(sca))) stop('Argument "sca" must be a named numeric vector')
  
  # If it is a call, then call recursively to each component
  if(is.call(e)){
    isAssignment <- length(e)==3 && as.character(e[[1]]) %in% c('<-', '=')
    if(isAssignment){
      # If it's an assignment, only modify the right side
      if(!is.null(e[[3]])) e[[3]] <- apollo_insertScaling(e[[3]], sca)
    } else {
      # If it's NOT an assignment, modify everything
      for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- apollo_insertScaling(e[[i]], sca) 
    }
  }
  
  # If it is a symbol
  if(is.symbol(e) && (as.character(e) %in% names(sca)) ){
    e <- str2lang(paste0(as.character(e), '*apollo_inputs$apollo_scaling["', names(sca)[names(sca)==e], '"]'))
  }
  
  # Restore function if needed and return
  if(is.function(eOrig)){
    body(eOrig) <- e
    e <- eOrig
  }
  return(e)
}

#f <- function(x){
#  A <- x + y
#  B <- x^2 + 3
#  return(A + B)
#}
#apollo_insertScaling(f, sca=c(x=2))
