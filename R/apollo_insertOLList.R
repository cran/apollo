#' Replaces \code{tau=c(...)} by \code{tau=list(...)} in calls to \code{apollo_ol}
#' 
#' Takes a function, looks for calls to apollo_ol, identifies the corresponding
#' ol_settings, then goes inside the definition of ol_settings and replaces
#' \code{tau=c(...)} for \code{tau=list(...)}.
#' 
#' This only goes one level deep in definitions. For example, it will work
#' correctly in the following cases:
#' \code{
#' ol_settings = list(outcomeOrdered = y1, 
#'                    V   = b1*x1, 
#'                    tau = c(tau11, tau12))
#' P[["OL1"]] = apollo_ol(ol_settings, functionality)
#' P[["OL2"]] = apollo_ol(list(outcomeOrdered=y2, V=b2*x2, tau=c(tau21, tau22)),
#'                        functionality)
#' }
#' But it will not work on the following cases:
#' \code{
#' Tau = c(tau1, tau2, tau3)
#' ol_settings = list(outcomeOrdered = y2,
#'                    V = b2*x2, 
#'                    tau = Tau)
#' P[["OL1"]] = apollo_ol(ol_settings, functionality)
#' P[["OL2"]] = apollo_ol(list(outcomeOrdered=y1, V=b1*x1, tau=Tau), functionality)
#' }
#' 
#' This function is called by apollo_modifyUserDefFunc to allow for analytical
#' gradients when using apollo_ol.
#' 
#' @param f Function. Usually \code{apollo_probabilities}, 
#'          \code{apollo_randCoeff}, or \code{apollo_lcPars}.
#' @return Function \code{f} with \code{tau=c(...)} replaced by 
#'         \code{tau=list(...)}.
#' @export
apollo_insertOLList <- function(f){
  ### Validate inputs
  if(!is.function(f)) stop("SYNTAX ISSUE - Argument 'f' must be a ",
                           "function (usually apollo_probabilities, ", 
                           "apollo_randCoeff, or apollo_lcPars)")
  
  ### Tools
  is.val <- function(e){
    test <- is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) ||
      is.complex(e)
    if(test) return(TRUE) else return(FALSE)
  }
  is.assignment <- function(e) is.call(e) && length(e)==3 && 
    as.character(e[[1]]) %in% c("<-", "=")
  
  ### Search for calls to apollo_ol, and its first argument
  ol_set <- c()
  argName <- function(e){
    # If value, return nothing
    if(is.val(e)) return(NULL)
    # If a call to apollo_ol, return name of the ol_settings variable
    test <- is.call(e) && !is.assignment(e) && length(e)>2 && 
      !is.null(e[[1]]) && as.character(e[[1]])=="apollo_ol"
    if(test){
      firstArg <- e[[2]]
      # If firstArg is c(ol_settings, componentName2=...), keep only ol_settings
      test <- is.call(firstArg) && length(firstArg)>1 && 
        !is.null(e[[1]]) && as.character(firstArg[[1]])=="c"
      if(test) firstArg <- firstArg[[2]]
      # If firstArg is a symbol, add it to ol_set
      if(is.symbol(firstArg)) ol_set <<- c(ol_set, as.character(firstArg))
      return(NULL)
    }
    # If expression or call, look inside each element 
    if(is.expression(e) || is.call(e)){
      for(i in 1:length(e)) if(!is.null(e[[i]])) argName(e[[i]])
      return(NULL)
    }
    stop("INTERNAL ISSUE - Argument 'e' should be a language object")
  }
  e <- body(f)
  argName(e)
  rm(argName)
  
  ### Replace c for list inside ol_settings
  replaceCList <- function(e, insideDef=FALSE){
    # If value, return as is
    if(is.val(e)) return(e)
    # If call to apollo_ol with settings defined in-place
    test <- is.call(e) && length(e)>1 && !is.null(e[[1]]) && 
      as.character(e[[1]])=="apollo_ol" && !is.null(e[[2]]) &&
      !is.symbol(e[[2]]) && length(e[[2]])>1 && !is.null(e[[2]][[1]]) &&
      is.symbol(e[[2]][[1]]) && as.character(e[[2]][[1]]) %in% c("c", "list")
    if(test){
      isList <- as.character(e[[2]][[1]])=="list" && !is.null(names(e[[2]])) &&
        "tau" %in% names(e[[2]])
      isVect <- as.character(e[[2]][[1]])=="c" && !is.null(e[[2]][[2]]) && 
        !is.symbol(e[[2]][[2]]) && length(e[[2]][[2]])>1 && 
        !is.null(e[[2]][[2]][[1]]) && is.symbol(e[[2]][[2]][[1]]) && 
        as.character(e[[2]][[2]][[1]])=="list" && 
        !is.null(names(e[[2]][[2]])) && "tau" %in% names(e[[2]][[2]])
      if(isList || isVect){
        if(isList) ee <- e[[2]] else ee <- e[[2]][[2]]
        iTau <- which(names(ee)=="tau")
        test2 <- is.call(ee[[iTau]]) && length(ee[[iTau]])>1 && 
          !is.null(ee[[iTau]]) && as.character(ee[[iTau]][[1]])=="c"
        if(test2) ee[[iTau]][[1]] <- as.symbol("list")
        if(isList) e[[2]] <- ee else e[[2]][[2]] <- ee
      }
      return(e)
    }
    # If call or function, but not insideDef, look inside each element
    if((is.call(e) || is.expression(e)) && !insideDef){
      for(i in 1:length(e)) if(!is.null(e[[i]])){
        test <- is.assignment(e[[i]]) && is.symbol(e[[i]][[2]]) &&
          (as.character(e[[i]][[2]]) %in% ol_set) &&
          is.call(e[[i]][[3]]) && length(e[[i]][[3]])>1 && 
          !is.null(e[[i]][[3]][[1]]) && as.character(e[[i]][[3]][[1]])=="list" && 
          !is.null(names(e[[i]][[3]])) && ("tau" %in% names(e[[i]][[3]]))
        e[[i]] <- replaceCList(e[[i]], test)
      }
      return(e)
    }
    # If call or function, and insideDef
    if((is.call(e) || is.expression(e)) && insideDef){
      iTau <- which(names(e[[3]])=="tau") # e[[3]] = list(outcomeOrdered=..., )
      test <- is.call(e[[3]][[iTau]]) && length(e[[3]][[iTau]])>1 && 
        !is.null(e[[3]][[iTau]][[1]]) && as.character(e[[3]][[iTau]][[1]])=="c"
      if(test) e[[3]][[iTau]][[1]] <- as.symbol("list")
      return(e)
    }
    stop("INTERNAL ISSUE - Argument 'e' should be a language object")
  }
  e <- replaceCList(e)
  body(f) <- e
  return(f)
}
