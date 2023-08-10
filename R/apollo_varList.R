#' Lists variable names and definitions used inside a function
#' 
#' Returns a list containing the names and definitions of variables in f, \code{apollo_randCoeff} and \code{apollo_lcPars}
#' 
#' It looks for variable definitions inside f, \code{apollo_randCoeff}, and \code{apollo_lcPars}. It returns 
#' them in a list.
#' 
#' @param f A function, usually \code{apollo_probabilities}
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' 
#' @return A list of expressions containing all definitions in f, \code{apollo_randCoeff} and \code{apollo_probabilities}
#' @export
#' @importFrom utils lsf.str
apollo_varList <- function(f, apollo_inputs){
  
  #### Useful functions ####
  ### Checks if argument is a symbol or a constant
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || 
                           is.logical(e) || is.complex(e)) return(TRUE) else return(FALSE)
  
  ### Returns the names of the indices (as a character vector) used inside "for" loops in "e", if any
  getLoopIndices <- function(e){
    if(is.function(e)) e <- body(e)
    ans <- c()
    # Case 1: Is a symbol or a constant
    if(is.val(e)) return(ans)
    # Case 2: is a 'for' loop
    test1 <- is.call(e) || is.expression(e)
    test2 <- test1 && !is.null(e[[1]]) && as.character(e[[1]]) %in% c("for", "while")
    if(test2) ans <- c(ans, as.character(e[[2]]), getLoopIndices(e[[4]]))
    # Case 3: is a call, but not a 'for' loop
    if(test1 && !test2) for(i in 1:length(e)) if(!is.null(e[[i]])) ans <- c(ans, getLoopIndices(e[[i]]))
    # Return
    return(ans)
  }
  
  ### Identifies if an expression is a call to a list element
  isListElem <- function(e){
    test <- (is.call(e) || is.expression(e)) && length(e)>1 && is.symbol(e[[1]])
    test <- test && as.character(e[[1]]) %in% c('[[', '$')
    return(test)
  }
  
  ### If argument is a function definition, returns only the body of the function. e.g. function(x) x^2 --> x^2
  getBody <- function(e){
    isFunDef <- (is.expression(e) || is.call(e)) && length(e)>=1
    isFunDef <- isFunDef && is.symbol(e[[1]]) && as.character(e[[1]])=='function'
    if(isFunDef) return(e[[3]]) else return(e)
  }
  
  ### Replaces expression by their definition
  replaceByDef <- function(e, defs){
    # Case 0: defs is empty, return e
    if(is.null(defs) || (is.list(defs) && length(defs)==0)) return(e)
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
      isDollar <- is.call(e) && length(e)==3 && is.symbol(e[[1]]) && as.character(e[[1]])=='$'
      for(i in 1:length(e)) if(!is.null(e[[i]])){
        emptyIndex <- (i>2 && e[[1]]=="[" && e[[i]]=="")
        if( emptyIndex) e[[i]] <- TRUE
        if(!emptyIndex && !(isDollar && i==3)) e[[i]] <- replaceByDef(e[[i]], defs)
      } 
    } 
    # Return
    return(e)
  }
  
  ### Extract names and definitions inside a function or expression, replacing previous definitions on the fly
  #   It returns the given defs with new elements added or overwritten. So recursive definitions will lead
  #   to a single definition combining all previous ones.
  extractDef <- function(e, defs=NULL){
    # If e is a function, extract its body
    if(is.function(e)) e <- body(e)
    # Case 1: e is a value, return immediately
    if(is.val(e)) return(defs)
    # Variable to store the definitions
    if(is.null(defs)) defs <- list()
    # Case 2: e is an assignment
    if(length(e)==3 && is.symbol(e[[1]]) && (e[[1]]=="=" || e[[1]]=="<-")){
      lTxt <- ""
      le <- e[[2]]
      re <- e[[3]]
      expr   <- (is.expression(re) || is.call(re)) && length(re)>=1
      rList  <- expr && is.symbol(re[[1]]) &&  as.character(re[[1]])=='list'
      rFCall  <- expr && is.symbol(re[[1]]) && (as.character(re[[1]]) %in% as.vector(lsf.str('package:apollo')))
      rFixMat <- expr && is.symbol(re[[1]]) && as.character(re[[1]])=="matrix" && length(re[[2]])==1 && is.numeric(re[[2]])
      # Case 2.1: right side is an empty list, apollo function call, or a fixed matrix; ignore
      if((rList && length(re)==1) || rFCall || rFixMat) return(defs) 
      # Determine name of left side element
      if( length(le)==1 && is.symbol(le) ) lTxt <- as.character(le) # left side is a symbol
      if( isListElem(le) ){ # left side is a list element
        while(isListElem(le)){
          if(!is.val(le[[3]])) return(defs) # give up if element's name is not a value (e.g. is an expression)
          lTxt <- paste0("$", le[[3]], lTxt)
          le <- le[[2]]
        }
        if(is.symbol(le)) lTxt <- paste0(as.character(le), lTxt)
      }; rm(le)
      # If name could not be determined, give up
      if(lTxt=="") return(defs)
      # Case 2.2: right side is a variable, value or expression
      if(!rList) defs[[lTxt]] <- replaceByDef(getBody(re), defs)
      # Case 2.3: right side is a list defined in place, expand into one entry per element of list
      if( rList){
        if(is.null(names(re))) tmp <- 0:(length(re)-1) else tmp <- names(re)
        for(i in 2:length(re)) defs[[paste0(lTxt, "$", tmp[i])]] <- replaceByDef(getBody(re[[i]]), defs)
      }
      # Case 2.4: right side is an Apollo function call, do nothing
      # - - - - - - - - - - - - - - - - - - - - - - - -  do nothing
      # Return defs
      return(defs)
    }
    # Case 3: e is a return statement with a list inside defined in place. e.g: return(list(a=1, b=2))
    test <- length(e)==2 && is.symbol(e[[1]]) && as.character(e[[1]])=='return'
    test <- test && is.call(e[[2]]) && length(e[[2]])>1 && is.symbol(e[[2]][[1]]) && as.character(e[[2]][[1]])=='list'
    test <- test && length(e[[2]])>1
    if(test){
      re <- e[[2]]
      if(is.null(names(re))) tmp <- 0:(length(re)-1) else tmp <- names(re)
      for(i in 2:length(re)) defs[[tmp[i]]] <- replaceByDef(getBody(re[[i]]), defs)
      return(defs)
    }
    # Case 4: e is an expression or call but not an assignment
    if(is.expression(e) || is.call(e)){
      for(i in 1:length(e)) if(!is.null(e[[i]])) defs <- extractDef(e[[i]], defs)
      return(defs)
    }
    # In any other case
    return(defs)
  }
  
  ### Returns the name of the variable returned (will only pick up the first return inside a function)
  identifyReturnVar <- function(e){
    # If e is a function, extract its body
    if(is.function(e)){f <- e; e <- body(e)} else f <- NULL
    # Case 1: e is a value, return immediately
    if(is.val(e)) return(NULL)
    # Case 2: e is a return statement
    isReturn <- is.call(e) && length(e)==2 && is.symbol(e[[1]]) && as.character(e[[1]])=='return'
    if(isReturn){
      # Case 2.1: e returns a single object
      if(is.symbol(e[[2]])) return(as.character(e[[2]]))
      # Case 2.2: e returns an expression (e.g. a list defined inside return)
      if(is.call(e[[2]])) return(NULL)
    }
    # Case 3: e is an expression, but not a return, iterate through e until it finds the first return
    ans <- NULL
    if(is.call(e) || is.expression(e)) for(i in 1:length(e)){
       if(!is.null(e[[i]])) ans <- identifyReturnVar(e[[i]])
       if(!is.null(ans)) return(ans)
    }
    # Case 4: The original e was a function and there is no return in it
    if(is.function(f)){
      e <- e[[length(e)]] # only keep last line of the function
      # Case 4.1: Last line is a symbol, return that symbol as text
      test <- !is.null(e) && is.symbol(e)
      if(test) return(as.character(e))
      # Case 4.2: Last line is a definition, return left side only if its a symbol (i.e. not if its a list element)
      test <- !is.null(e) && length(e)==3 && is.symbol(e[[1]]) && (e[[1]]=="=" || e[[1]]=="<-")
      test <- test && length(e[[2]])==1 && is.symbol(e[[2]])
      if(test) return(as.character(e[[2]]))
    }
    # Case 5: Any other case, return NULL
    return(NULL)
  }
  
  
  
  
  #### Processing of the function ####
  defs <- list()
  
  ### Load definitions inside apollo_randCoeff, except loop indices
  if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
    defs <- c(defs, extractDef(apollo_inputs$apollo_randCoeff)) # extract definitions in apollo_randCoeff
    # Remove the definition of any loop indices
    tmp <- names(defs) %in% getLoopIndices(apollo_inputs$apollo_randCoeff)
    if(any(tmp)) defs <- defs[!tmp]
    # Duplicate the definitions of "return$x", as just "x"
    tmp <- identifyReturnVar(apollo_inputs$apollo_randCoeff)
    if(is.character(tmp) && length(tmp)==1){ # maybe this could work for length(tmp)>1
      tmp2 <- grep(paste0("^",tmp,"\\$"), names(defs), value=TRUE)
      if(length(tmp2)>0){
        defs2 <- defs[tmp2]
        names(defs2) <- gsub(paste0("^",tmp,"\\$"), "", names(defs2))
        defs <- c(defs, defs2)
        rm(defs2)
      }; rm(tmp2)
    }; rm(tmp)
  }
  
  ### Load definitions inside apollo_lcPars, except loop indices
  if(is.function(apollo_inputs$apollo_lcPars)){
    defs <- c(defs, extractDef(apollo_inputs$apollo_lcPars)) # extract definitions in apollo_lcPars
    # Remove the definition of any loop indices
    tmp <- names(defs) %in% getLoopIndices(apollo_inputs$apollo_lcPars)
    if(any(tmp)) defs <- defs[!tmp]
    # Duplicate the definitions of "return$x", as just "x"
    tmp <- identifyReturnVar(apollo_inputs$apollo_lcPars)
    if(is.character(tmp) && length(tmp)==1){ # maybe this could work for length(tmp)>1
      tmp2 <- grep(paste0("^",tmp,"\\$"), names(defs), value=TRUE)
      if(length(tmp2)>0){
        defs2 <- defs[tmp2]
        names(defs2) <- gsub(paste0("^",tmp,"\\$"), "", names(defs2))
        defs <- c(defs, defs2)
        rm(defs2)
      }; rm(tmp2)
    }; rm(tmp)
  }
  
  ### Extract definitions inside f, except loop indices
  defs <- c(defs, extractDef(f)) # extract definitions in f
  # Remove the definition of any loop indices
  tmp <- names(defs) %in% getLoopIndices(f)
  if(any(tmp)) defs <- defs[!tmp]
  # Duplicate the definitions of "return$x", as just "x"
  tmp <- identifyReturnVar(f)
  if(is.character(tmp) && length(tmp)==1){ # maybe this could work for length(tmp)>1
    tmp2 <- grep(paste0("^",tmp,"\\$"), names(defs), value=TRUE)
    if(length(tmp2)>0){
      defs2 <- defs[tmp2]
      names(defs2) <- gsub(paste0("^",tmp,"\\$"), "", names(defs2))
      defs <- c(defs, defs2)
      rm(defs2)
    }; rm(tmp2)
  }; rm(tmp)
  
  return(defs)
}