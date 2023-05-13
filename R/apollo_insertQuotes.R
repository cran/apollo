#' Introduces quotes into rrm_settings 
#' 
#' Takes a function, looks for the definition of relevant parts of rrm_settings,
#' and introduces quotes on them. This is to facilitate their processing by
#' apollo_rrm under functionality="preprocessing". 
#' 
#' @param f Function. Usually \code{apollo_probabilities}.
#' @return Function \code{f} with relevant expressions turned into character.
#' @export
apollo_insertRRMQuotes <- function(f){
  # Validate inputs
  if(!is.function(f)) stop('INTERNAL ISSUE - Argument "f" should be a function.')
  
  # If there's no call to apollo_rrm, return function as is
  if(!any(grepl("apollo_rrm", deparse(f)))) return(f)
  
  # Check if input is a single value
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) || is.complex(e)) return(TRUE) else return(FALSE)
  
  # Turns any call or value into a character
  # If the argument is a list, it turns its elements into character
  withQuotes <- function(e){
    if(is.character(e)) return(e)
    if(is.val(e)      ) return( deparse(e) )
    isCall <- is.call(e) || is.expression(e)
    if(isCall && length(e)>1 && !is.null(e[[1]]) && e[[1]]=="list"){
      for(i in 2:length(e)) e[[i]] <- withQuotes(e[[i]])
      return(e)
    } else return( paste0(deparse(e), collapse="") )
    #if(isCall) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- withQuotes(e[[i]])
    #return(e)
  }
  
  # Looks for the definition of a list named listName and returns:
  # An expression: If it finds the requested list and it is defined in place,
  #                It modifies the input (e) and returns it with quotes
  # A character  : If it finds the If it couldn't do it, but found that the list is defined 
  #                under a different name
  insertQuotes <- function(e, listName){
    if(is.function(e)){f <- e; e <- body(f)} else f <- NULL
    # Case 1: a value: return immediately
    if(is.val(e)) return(e)
    
    # Check 
    isExp <- is.expression(e) || is.call(e)
    isDef <- isExp && length(e)==3 && is.symbol(e[[1]]) 
    isDef <- isDef && (as.character(e[[1]]) %in% c('<-', '=')) && is.symbol(e[[2]])
    isLis <- isDef && is.call(e[[3]]) && length(e[[3]])>1 && is.symbol(e[[3]][[1]])
    isLis <- isLis && as.character(e[[3]][[1]])=="list"
    
    # Case 2: it's an assignment with listName on the left
    if(isDef && as.character(e[[2]])==listName){
      # Case 2.1: on the right there's a list defined in place
      if(isLis){
        for(i in 2:length(e[[3]])) if(!is.null(e[[3]][[i]])) e[[3]][[i]] <- withQuotes(e[[3]][[i]])
        if(!is.null(f)) body(f) <- e
        return(e)
      }
      # Case 2.2: on the right there's just a symbol
      if( is.symbol(e[[3]]) ) return(as.character(e[[3]]))
    }
    
    # Case 3: an assignment with listName inside a list on the right
    if(isDef && as.character(e[[2]])!=listName && isLis && listName %in% names(e[[3]])){
      i <- which(listName==names(e[[3]]))# - 1
      # Case 3.1: definition is in place
      test <- is.call(e[[3]][[i]]) && length(e[[3]][[i]])>1 && is.symbol(e[[3]][[i]][[1]])
      test <- test && as.character(e[[3]][[i]][[1]])=="list"
      if(test){
        e[[3]][[i]] <- withQuotes(e[[3]][[i]])
        return(e)
      }
      # Case 3.2: definition is just a symbol
      if( is.symbol(e[[3]][[i]]) ) return(as.character(e[[3]][[i]]))
    }
    
    # Case 4: not a definition nor a value: go through it recursively
    if(isExp && !isDef && length(e)>=1){
      for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- insertQuotes(e[[i]], listName)
    }
    if(!is.null(f)){ body(f) <- e; e <- f}
    return(e)
  }
  
  # Insert quotes in all relevant definitions
  for(i in c("rum_inputs", "regret_inputs", "regret_scale")){
    f2 <- insertQuotes(f, i)
    j <- 1
    while(is.character(f2) && j<=10){ f2 <- insertQuotes(f, f2); j <- j + 1}
    if(is.function(f2)) f <- f2
  }
  return(f)
}