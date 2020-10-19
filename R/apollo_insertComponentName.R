#' Adds componentName2 to model calls
#' 
#' @param e An expression or a function. It will usually be apollo_probabilities.
#' @return The original argument 'e' but modified to incorporate a new setting
#'         called 'componentName2' to every call to apollo_<model> (e.g. 
#'        apollo_mnl, apollo_nl, etc.).
#' @export
apollo_insertComponentName <- function(e){
  # Validate input
  if(is.null(e)) stop('Argument "e" must be a function, a call , or a value')
  if(is.function(e)){eOrig <- e; e <- body(e)} else eOrig <- NULL
  if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e)) return(e)
  if(!is.call(e)) stop('Argument "e" must be a call')
  
  # Figure out if e is a call of the type: x <- apollo_<model>(...)
  test <- FALSE
  # Check if it is an assignment
  test <- length(e)>=3
  test <- test && (e[[1]]=="=" || e[[1]]=="<-")
  # Check that the right hand side of the assignment is a relevant function call
  test <- test && is.call(e[[3]])
  test <- test && (as.character(e[[3]][[1]]) %in% c('apollo_mnl', 'apollo_el', 'apollo_nl', 
                                                    'apollo_cnl',  'apollo_ol', 'apollo_op', 
                                                    'apollo_dft', 'apollo_normalDensity', 
                                                    'apollo_mdcev', 'apollo_mdcnev', 'apollo_lc'))
  
  # If e is NOT of the type: x <- apollo_<model>(...)
  if(!test && length(e)>1) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- apollo_insertComponentName(e[[i]])
  
  # Check that componentName has not been added already
  if(test && !is.null(names(e[[3]])) && names(e[[3]])[2]=='functionality') setPos <- 3 else setPos <- 2
  if(test && !is.null(names(e[[3]][[setPos]])) && ('componentName2' %in% names(e[[3]][[setPos]])) ) test <- FALSE
  
  # If e IS of the type: x <- apollo_<model>(...)
  if(test){
    # get name of variable in the left side
    asignName <- strsplit(as.character(e[[3]][[1]]), split='_')[[1]]
    asignName <- asignName[length(asignName)]
    if(is.symbol(e[[2]])) asignName <- e[[2]]
    if(is.call(e[[2]]) && as.character(e[[2]][[1]]) %in% c('[[','$')) asignName <- e[[2]][[3]]
    # Change expression
    e2 <- str2lang('c(a, componentName2=b)')
    e2[[2]] <- e[[3]][[setPos]]
    e2[[3]] <- asignName
    e[[3]][[setPos]] <- e2
  }
  
  # Return expression (or function)
  if(is.function(eOrig)){
    body(eOrig) <- e
    e <- eOrig
  }
  return(e)
}