#' Expands loops in a function or expression
#' 
#' Expands loops replacing the index by its value. It also evaluates \code{paste} and \code{paste0}, and removes \code{get}.
#' 
#' For example, the expression
#' \code{for(j in 1:3) V[[paste0('alt',j)]] = b1*get(paste0('x',j)) + b2*X[,j]}
#' 
#' would be expanded into:
#' 
#' \code{
#' V[[alt1]] = b1*x1 + b2*X[,1]
#' V[[alt2]] = b1*x2 + b2*X[,2]
#' V[[alt3]] = b1*x3 + b2*X[,3]
#' }
#' 
#' @param f function (usually \code{apollo_probabilities}) inside which the name of the components are inserted.
#' @param apollo_inputs List. Main inputs necessary for model estimation. See \link{apollo_validateInputs}.
#' @return A function or an expression (same type as input \code{f})
#' @export
apollo_expandLoop <- function(f, apollo_inputs){
  #### Utilities ####
  
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) ) return(TRUE) else return(FALSE)
  
  replaceByDef <- function(e, defs, rightSide=FALSE){
    isFunction <- is.function(e)
    if(isFunction){f <- e; e <- body(e)}
    # Case 1: x
    test1 <- is.symbol(e) && (as.character(e) %in% names(defs))
    if(test1 && rightSide) e <- defs[[as.character(e)]]
    # Case 2: L$x or L[['x']] or L[["x"]]
    test2 <- !test1 && is.call(e) && length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('$','[['))
    test2 <- test2 && is.val(e[[3]]) && is.symbol(e[[2]])
    if(test2) tmp  <- paste0(as.character(e[[2]]), '$', as.character(e[[3]]))
    test2 <- test2 && (tmp %in% names(defs))
    if(test2 && rightSide) e <- defs[[tmp]]
    # Case 3: expression
    if( !test1 && !test2 && (is.call(e) || is.expression(e)) ){
      isAssign <- length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('<-', '='))
      for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- replaceByDef(e[[i]], defs, rightSide=(rightSide | (isAssign & i==3)))
    } 
    # Return
    if(isFunction){body(f) <- e; e <- f}
    return(e)
  }
  
  replaceIndex <- function(e, jNam, jVal){
    # Case 1: is j
    if(is.symbol(e) && e==jNam) return(jVal)
    # Case 2: is a value but not j
    if(is.val(e)) return(e)
    # Case 3: is an expression
    if(is.expression(e) || is.call(e)) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- replaceIndex(e[[i]], jNam, jVal)
    return(e)
  }
  
  runPaste <- function(e, env=NULL){
    # Case 1: just a value
    if(is.val(e)) return(e)
    # Case 2: a call to paste or paste0
    test1 <- is.expression(e) || is.call(e)
    test1 <- test1 && length(e)>1 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('paste', 'paste0'))
    if(test1){
      e <- tryCatch(eval(e, envir=env), 
                    error=function(em) stop('Could not evaluate ', deparse(e)))
      return(e)
    }
    # Case 3: an expression
    test2 <- (is.call(e) || is.expression(e)) && !test1
    if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- runPaste(e[[i]], env)
    return(e)
  }
  
  rmGet <- function(e, env=NULL){
    # Case 1: just a value
    if(is.val(e)) return(e)
    # Case 2: a call to get
    test1 <- is.expression(e) || is.call(e)
    test1 <- test1 && length(e)>1 && is.symbol(e[[1]]) && as.character(e[[1]])=='get'
    if(test1){
      if(length(all.vars(e))>0) stop('Unknown variables inside "get" (can only use index).')
      test3 <- length(e)==2 && is.character(e[[2]])
      if(!test3) stop('Could not build variable name fetched by "get".')
      e <- as.symbol(e[[2]])
      return(e)
    }
    # Case 3: an expression
    test2 <- (is.call(e) || is.expression(e)) && !test1
    if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- rmGet(e[[i]], env)
    return(e)
  }
  
  evalIndex <- function(e, env=NULL){
    # Case 1: just a value
    if(is.val(e)) return(e)
    # Case 2: a call to '[' or '[['
    test <- (is.expression(e) || is.call(e)) && length(e)==3 && is.symbol(e[[1]]) && as.character(e[[1]]) %in% c('[', '[[')
    if(test){
      # Case 2.1: Index is a value
      if(is.val(e[[3]])) return(e)
      # Case 2.2: Index is an expression without variables
      if(length(all.vars(e[[3]]))==0){
        e[[3]] <- eval(e[[3]])
        return(e)
      }
      # Case 2.3: Index is an expression with variables
      if(length(all.vars(e[[3]]))>0){
        if(!is.null(env)){
          e[[3]] <- tryCatch(eval(e[[3]], envir=env), error=function(e) NULL)
        } else e[[3]] <- NULL
        if(is.null(e[[3]])) stop('expandLoop: Unknown variable inside index')
        return(e)
      }
    }
    # Case 3: an expression
    test2 <- (is.call(e) || is.expression(e)) && !test
    if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- evalIndex(e[[i]], env=env)
    return(e)
  }
  
  expandLoop <- function(e, defs=NULL, env=NULL){
    # Initialise
    isF <- is.function(e)
    if(isF){f <- e; e <- body(e)}
    # Case 1: it is a symbol
    if(is.val(e)) return(e)
    # Case 2: it is a "for" call
    test1 <- is.call(e) || is.expression(e)
    test1 <- test1 && length(e)==4 && is.symbol(e[[1]]) && as.character(e[[1]])=='for'
    if(test1){
      #if(!is.null(defs)) e <- replaceByDef(e, defs)
      jSym <- e[[2]] # index
      if(!is.null(defs)) e[[3]] <- replaceByDef(e[[3]], defs, rightSide=TRUE)
      jValues <- tryCatch(eval(e[[3]], envir=env), # all values that j needs to take
                          error=function(e) stop('Could not evaluate all possible values for index ', jSym, '.'))
      if(!is.null(defs)){ # replace definitions in expression inside loop, except for the index
        tmp <- which(names(defs)==as.character(jSym))
        if(any(tmp)) e[[4]] <- replaceByDef(e[[4]], defs[-tmp]) else e[[4]] <- replaceByDef(e[[4]], defs)
      } 
      ee <- str2lang(paste0('{', paste0(jValues, collapse='; '), '}'))
      for(j in jValues){
        ee[[1+j]] <- replaceIndex(e[[4]], jSym, j)
        ee[[1+j]] <- expandLoop(ee[[1+j]], defs, env) # in case there are nested loops
        ee[[1+j]] <- runPaste(ee[[1+j]], env)
        ee[[1+j]] <- rmGet(ee[[1+j]])
        ee[[1+j]] <- evalIndex(ee[[1+j]], env)
      }
      if(!is.null(defs)) ee <- replaceByDef(ee, defs)
      return(ee)
    }
    # Case 3: It is a call but not a for
    test2 <- (is.call(e) || is.expression(e)) && !test1
    if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- expandLoop(e[[i]], defs, env)
    if(!isF) return(e) else {body(f) <- e; return(f)}
  }
  
  isDef <- function(e, dollar=TRUE){
    # Check if it's a definition
    test1 <- is.call(e) || is.expression(e)
    test1 <- test1 && length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('=', '<-'))
    if(test1){
      # assignment to simple variable (not a list)
      if(is.symbol(e[[2]])) return( as.character(e[[2]]) )
      # assignment to a list element
      test2 <- is.call(e[[2]]) && length(e[[2]])==3 && is.symbol(e[[2]][[1]]) && (as.character(e[[2]][[1]]) %in% c('[[', '$'))
      test2 <- test2 && is.symbol(e[[2]][[2]]) && (is.symbol(e[[2]][[3]]) || is.val(e[[2]][[3]]))
      if(test2 & dollar) return( paste0(as.character(e[[2]][[2]]), '$', as.character(e[[2]][[3]])) )
      if(test2 & !dollar){
        if(is.character(e[[2]][[3]])) tmp <- c('[["', '"]]') else tmp <- c('[[', ']]')
        return( paste0(as.character(e[[2]][[2]]), tmp[1], as.character(e[[2]][[3]]), tmp[2]) )
      } 
      return("")
    } else return("") # return empty character string
  }
  
  simplify <- function(e, defs){
    # Initialise
    isF <- is.function(e)
    if(isF){f <- e; e <- body(e)}
    # Case 1: it is a symbol
    if(is.val(e)) return(e)
    # Case 2: it is a block of curly braces filled with (re)definitions inside
    test1 <- is.call(e) || is.expression(e)
    test1 <- test1 && length(e)>1 && is.symbol(e[[1]]) && as.character(e[[1]])=='{'
    if(test1){
      varNames <- rep("", length(e)-1)
      for(i in 2:length(e)) varNames[i-1] <- isDef(e[[i]])
      test2 <- all(varNames!="" & varNames==varNames[1]) && (varNames[1] %in% names(defs))
      if(test2){
        #ee <- str2lang(paste0(varNames[1], '<-', 0))
        ee <- str2lang(paste0(isDef(e[[2]], dollar=FALSE), '<-', 0))
        ee[[3]] <- defs[[varNames[1]]]
        return(ee)
      }
    }
    # Case 3: It is a call but not as in case 2
    test2 <- (is.call(e) || is.expression(e))
    if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- simplify(e[[i]], defs)
    if(!isF) return(e) else {body(f) <- e; return(f)}
  }
  
  #replaceDollar <- function(e){
  #  # Chek if input is a function
  #  if(is.function(e)){ eOrig <- e; e <- body(e)} else eOrig=NULL
  #  # Case 1: just a value
  #  if(is.val(e)) return(e)
  #  # Case 2: L$x or L[['x']] or L[["x"]]
  #  test2 <- is.call(e) && length(e)==3 && is.symbol(e[[1]]) && as.character(e[[1]])=='$'
  #  test2 <- test2 && is.val(e[[3]]) && is.symbol(e[[2]])
  #  if(test2) tmp  <- paste0(as.character(e[[2]]), '$', as.character(e[[3]]))
  #  test2 <- test2 && (tmp %in% names(defs))
  #  if(test2 && rightSide) e <- defs[[tmp]]
  #  # Case 3: an expression
  #  test2 <- (is.call(e) || is.expression(e)) && !test
  #  if(test2) for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- evalIndex(e[[i]], env=env)
  #  return(e)
  #  # Return
  #  if(!is.null(eOrig)){ body(eOrig) <- q; return(eOrig)} else return(e)
  #}
  
  ### Fetch apollo_beta
  apollo_beta <- tryCatch(get('apollo_beta', envir=parent.frame(1), inherits=FALSE),
                          error=function(e) NULL)
  if(is.null(apollo_beta)) apollo_beta <- tryCatch(get('apollo_beta', envir=globalenv(), inherits=FALSE),
                                                   error=function(e) NULL)
  if(is.null(apollo_beta)) stop('apollo_expandLoop could not fetch apollo_beta.')
  
  #### Process and return ####
  defs  <- apollo_varList(f, apollo_inputs)
  # maybe add apollo_randCoeff and apollo_lcPar?
  test <- anyNA(apollo_inputs$draws)
  if(test) env <- list2env(c(as.list(apollo_beta), apollo_inputs$database), hash=TRUE, parent=parent.frame()) else {
    env <- list2env(c(as.list(apollo_beta), apollo_inputs$database, apollo_inputs$draws), 
                    hash=TRUE, parent=parent.frame())
  }
  env$apollo_inputs <- apollo_inputs
  fNew  <- expandLoop(f, defs, env)
  defs  <- apollo_varList(fNew, apollo_inputs)
  fNew  <- simplify(fNew, defs)
  fNew  <- simplify(fNew, defs) # up to two levels of loops. Maybe I can add a third, but doesn't make much sense
  return(fNew)
}



#rm(list=ls())
#f=function(x){
#  J = 100
#  V = list()
#  for(j in 1:J) V[[paste0("alt",j)]] = b1*get(paste0("x1_",j)) + b2*get(paste0("x2_",j))
#  # 1) Run paste
#  #V[["altj"]] = b1*get("x1_j") + b2*get("x2_j")
#  ## 2) Get rid of get
#  #V[["altj"]] = b1*x1_j + b2*x2_j
#  
#  #X1 = apollo_inputs$database[,paste0('tt', 1:100)]
#  #for(j in 1:J) V[[paste0("alt",j)]] = b1*apollo_inputs$database[,paste0('tt', j)] + b2*X2[,j]
#  ## 1)
#  #V[[paste0("alt",j)]] = b1*apollo_inputs$database[,"ttj"] + b2*X2[,j]
#  return(V)
#}
#defs <- apollo_varList(f, list(apollo_control=list(mixing=FALSE), apollo_randCoeff=NA))
#apollo_expandLoop(f, defs)
#

#e <- expression(for(j in 1:10) V[[paste0("alt",j)]] = b1*get(paste0("x1_",j)) + b2*get(paste0("x2_",j)))[[1]]
#replaceIndex(e, as.symbol('j'), 10)
#runPaste(replaceIndex(e, as.symbol('j'), 5))
#expandLoop(e)
#expandLoop(apollo_probabilities)