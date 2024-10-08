#' Modifies function to make it compatible with analytic gradients
#' 
#' Takes a likelihood function and inserts \code{function ()} before key elements to allow for analytic gradient calculation
#' 
#' It modifies the definition of the following models.
#' \itemize{
#'   \item \strong{\code{apollo_mnl}}: Turns all elements inside \code{mnl_settings$V} into functions.
#'   \item \strong{\code{apollo_ol}}: Turns \code{ol_settings$V} and all elements inside \code{ol_settings$tau} into functions.
#'   \item \strong{\code{apollo_op}}: Turns \code{op_settings$V} and all elements inside \code{op_settings$tau} into functions.
#'   \item \strong{\code{apollo_normalDensity}}: Turns \code{normalDensity_settings$xNormal}, \code{normalDensity_settings$mu} and \code{normalDensity_settings$sigma} into functions.
#' }
#' It can only track a maximum of 3 levels of depth in definitions. For example:
#' \code{
#' V <- list()
#' V[["A"]] <- b1*x1A + b2*x2A
#' V[["B"]] <- b1*x1B + b2*x2B
#' mnl_settings1 <- list(alternatives=c("A", "B"), V = V, choiceVar= Y, avail = 1, componentName="MNL1")
#' P[["MNL1"]] <- apollo_mnl(mnl_settings1, functionality)
#' }
#' But it may not be able to deal with the following:
#' \code{
#' VA <- b1*x1A + b2*x2A
#' V <- list()
#' V[["A"]] <- VA
#' V[["B"]] <- b1*x1B + b2*x2B
#' mnl_settings1 <- list(alternatives=c("A", "B"), V = V, choiceVar= Y, avail = 1, componentName="MNL1")
#' P[["MNL1"]] <- apollo_mnl(mnl_settings1, functionality)
#' }
#' But that might be enough given how apollo_dVdB works.
#' 
#' @param f Function. Expressions inside it will be turned into functions. Usually \code{apollo_probabilities} or 
#'          \code{apollo_randCoeff}.
#' @param like Logical. Must be TRUE if \code{f} is \code{apollo_probabilities}. FALSE otherwise.
#' @param randCoeff Logical. Must be TRUE if \code{f} is \code{apollo_randCoeff}. FALSE otherwise.
#' @param lcPars Logical. Must be TRUE if \code{f} is \code{apollo_lcPars}. FALSE otherwise.
#' @return Function \code{f} but with relevant expressions turned into function definitions.
#' @importFrom utils capture.output
#' @export
apollo_insertFunc <- function(f, like=TRUE, randCoeff=FALSE, lcPars=FALSE){
  # Validate inputs
  if(!is.function(f)) stop('INTERNAL ISSUE - Argument "f" should be a function.')
  
  ## Check if input is a single value
#func  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) || is.complex(e)) return(TRUE) else return(FALSE)
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) || is.complex(e)) return(TRUE) else return(FALSE)
  
  # Vector to store the elements that are yet to be turned into functions
  pendingVars <- c()
  
  ## This functions throws an error if it detects an assignment
  stopIfAssignment <- function(e){
    isCall <- is.call(e)
    isVal  <- is.val(e)
    # Case 0: Wrong input
    if(!isVal && !isCall) stop("INTERNAL ISSUE: Aregument 'e' must be a language object (a call or symbol).")
    # Case 1: a value
    if(isVal) return(FALSE)
    # Case 2: a call
    isAssignment <- !is.null(e) && length(e)==3 && is.symbol(e[[1]])
    isAssignment <- isAssignment && as.character(e[[1]]) %in% c("=", "<-")
    # Case 2.1: an assignment
    if(isAssignment) stop("SYNTAX ISSUE: Utility definitions cannot contain assignments, as in ", deparse(e))
    # Case 2.2: a call with multiple components
    for(i in length(e)) if(!is.null(e[[i]])) stopIfAssignment(e[[i]])
    return(FALSE)
  }
  
  ## Replace definition by function
  # If e is a simple call, it returns 'function() e'
  # If e is a list definition, it returns 'list(function() e[[2]], function() e[[3]], ...)'.
  # If e starts with "function" or "eval" it does not make any changes.
  asFunctionDef <- function(e){
    isCall <- is.call(e)
    isVal  <- is.val(e)
    if(!isCall & !isVal) stop('INTERNAL ISSUE - Argument "e" must be a language object (a call or a symbol).')
    if(isCall && (e[[1]]=="function" || e[[1]]=="eval")) return(e)
    if(isCall && e[[1]]=="list"){
      for(i in 2:length(e)) e[[i]] <- asFunctionDef(e[[i]])
    } else {
      stopIfAssignment(e)
      if(isVal) e <- as.character(e) else e <- paste0(capture.output(print(e)), collapse="")
      e <- str2lang(paste0("function() ", e))
    }
    return(e)
  }
  
  # This function replace expressions on the right side by their definitions
  replaceByDef <- function(e, defs, rightSide=FALSE){
    if(is.function(e)){ f <- e; e <- body(e)} else f <- -1L# f <- NULL
    # Case 1: LV1
    test1 <- rightSide && is.symbol(e) && (as.character(e) %in% names(defs))
    if(test1) e <- defs[[ which(names(defs)==as.character(e))[1] ]]
    # Case 2: L$LV1 or L[['LV1']]
    test2 <- !test1 && rightSide && is.call(e) && length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('$','[['))
    test2 <- test2 && is.val(e[[3]]) && is.symbol(e[[2]])
    if(test2) tmp  <- paste0(as.character(e[[2]]), '$', as.character(e[[3]]))
    test2 <- test2 && (tmp %in% names(defs))
    if(test2) e <- defs[[tmp]]
    # Case 3: expression
    if(!test1 && !test2 && is.call(e)){
      test0 <- length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('<-', '=')) # Is an assignment
      for(i in 1:length(e)) if(!is.null(e[[i]])){
        isFuncArg <- i==2 && is.symbol(e[[i-1]]) && as.character(e[[i-1]])=="function"
        if(!isFuncArg) e[[i]] <- replaceByDef(e[[i]], defs, rightSide=(rightSide | (test0 & i==3)))
      }
    } 
    # Return
    if(is.function(f)){body(f) <- e; return(f)} else return(e)
    #if(is.null(f)) return(e) else {body(f) <- e; return(f)}
  }
  
  ## Searches inside "e" for the definition of listNames, and turns its elements "elemNames" into
  ## functions, as well as all the elements inside elemListNames.
  ## If elemNames="*", then all elements of listNames are turned into functions.
  ## If an elemListNames is defined as another variable, then it adds it to the pendingVars vector
  listElem2Func <- function(e, listNames, elemNames=NULL, elemListNames=NULL){
    # Checks
    if(is.function(e)) e <- body(e)
    if(is.null(elemNames) && is.null(elemListNames)) stop('INTERNAL ISSUE - Arguments "elemNames" and "elemListNames" cannot be both NULL.')
    if(!is.character(listNames)) stop('INTERNAL ISSUE - Argument "listNames" must be a character vector')
    if(!is.null(elemNames) && !is.character(elemNames)) stop('INTERNAL ISSUE - Argument "elemNames", if provided, must be a character vector')
    if(!is.null(elemListNames) && !is.character(elemListNames)) stop('INTERNAL ISSUE - Argument "elemListNames", if provided, must be a character vector')
    # Case 1: value
    if(is.val(e)) return(e)
    # Case 2: L$x <- ... with x a single element
    isAssign <- length(e)>=3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('=', '<-'))
    lList <- isAssign && length(e[[2]])==3 && is.symbol(e[[2]][[1]]) && (as.character(e[[2]][[1]]) %in% c('[[', '$'))
    lList <- lList && is.symbol(e[[2]][[2]]) && (as.character(e[[2]][[2]]) %in% listNames)
    lListElem <- lList && is.symbol(e[[2]][[3]]) && (as.character(e[[2]][[3]]) %in% elemNames)
    lListAll  <- lList && ("*" %in% elemNames) && !(is.symbol(e[[2]][[3]]) && (as.character(e[[2]][[3]]) %in% elemListNames))
    if(lListElem || lListAll) e[[3]] <- asFunctionDef(e[[3]])
    # Case 3: L$x <- ... with x a list
    lListList <- lList && is.symbol(e[[2]][[3]]) && (as.character(e[[2]][[3]]) %in% elemListNames)
    lListList1 <- lListList && is.symbol(e[[3]])
    if(lListList1) pendingVars <<- c(pendingVars, paste0(as.character(e[[3]]), "$*")) # defined somewhere else
    lListList2 <- lListList && is.call(e[[3]]) && length(e[[3]])>=2 && is.symbol(e[[3]][[1]]) && e[[3]][[1]]=="list"
    if(lListList2) e[[3]] <- asFunctionDef(e[[3]]) # defined in place
    # Case 4: L <- list(x=...)
    lListOnly <- isAssign && is.symbol(e[[2]]) && (as.character(e[[2]]) %in% listNames)
    lListOnly <- lListOnly && is.call(e[[3]]) && is.symbol(e[[3]][[1]]) && as.character(e[[3]][[1]])=='list'
    lListOnly <- lListOnly && length(e[[3]])>1  # not assigning empty list
    if(lListOnly){
      tmp <- names(e[[3]]) # names of elements in list definition (including "" at the beginning)
      # L <- list(x=...) when x is a single element
      if("*" %in% elemNames) tmp2 <- 2:length(tmp) else tmp2 <- which(tmp %in% elemNames)
      if(length(tmp2)>0) for(j in tmp2) if(!(tmp[j] %in% elemListNames)) e[[3]][[j]] <- asFunctionDef(e[[3]][[j]])
      # L <- list(x=...) when x is a list 
      tmp2 <- which(tmp %in% elemListNames)
      if(length(tmp2)>0) for(j in tmp2){
        # L <- list(x=y) when x is a list and y is defined somewhere else
        if(is.symbol(e[[3]][[j]])) pendingVars <<- c(pendingVars, paste0(as.character(e[[3]][[j]]), "$*"))
        # L <- list(x=list(...)) when x is a list defined in place
        test <- is.call(e[[3]][[j]]) && length(e[[3]][[j]])>=2 
        test <- test && is.symbol(e[[3]][[j]][[1]]) && as.character(e[[3]][[j]][[1]])=="list"
        if(test) e[[3]][[j]] <- asFunctionDef(e[[3]][[j]])
      }
    }
    # Other cases
    if(!lListElem && !lListAll && !lListList1 && !lListList2 && !lListOnly && is.call(e)){
      for(i in 1:length(e)) if(!is.null(e[[i]])){
        isFuncArg <- i==2 && is.symbol(e[[i-1]]) && as.character(e[[i-1]])=="function"
        if(!isFuncArg) e[[i]] <- listElem2Func(e[[i]], listNames, elemNames, elemListNames)
      }
    }
    # Return
    return(e)
  }
  
  ## Searches for calls to function 'fName'. If the call contains as its arguments the definition 
  ## of a list, searches within that list for elements named 'elemNames' and turn its definition
  ## into functions. Returns the expression 'e' altered in the way just described.
  replaceInModelCall <- function(fName, elemNames=NULL, listNames=NULL, e=NULL){
    if(is.function(e)) e <- body(e)
    if(is.val(e)) return(e)
    if(!is.call(e)) stop('INTERNAL ISSUE - Argument "e" must be a language object')
    if(is.null(elemNames) & is.null(listNames)) return(e)
    if(!is.null(elemNames) && !is.character(elemNames)) stop('INTERNAL ISSUE - Argument "elemNames" must be a character vector')
    if(!is.null(listNames) && !is.character(listNames)) stop('INTERNAL ISSUE - Argument "listNames" must be a character vector')
    # Check if 'e' is a call to fName
    test <- e[[1]]==fName
    if(test) if(!is.null(names(e)) && names(e)[2]=="functionality") setPos <- 3 else setPos <- 2
    # If settings are defined in the same call
    if(test && is.call(e[[setPos]]) && e[[setPos]]!="list()" && e[[setPos]][[1]]=="list"){
      # Turn single elements into functions
      elemPos <- which(names(e[[setPos]]) %in% elemNames)
      if(length(elemPos)>0) for(i in elemPos) e[[setPos]][[i]] <- asFunctionDef(e[[setPos]][[i]])
      # Turn lists into functions
      elemPos <- which(names(e[[setPos]]) %in% listNames)
      if(length(elemPos)>0) for(i in elemPos){
        if(is.symbol(e[[setPos]][[i]])){ # If the list is another variable, change it later
          pendingVars <<- c(pendingVars, paste0(as.character(e[[setPos]][[i]]), "$*"))
        } else e[[setPos]][[i]] <- asFunctionDef(e[[setPos]][[i]])
      }
    }
    # If settings are defined somewhere else
    if(test && is.symbol(e[[setPos]])){
      tmp <- paste0(as.character(e[[setPos]]), "$")
      if(!is.null(elemNames)) pendingVars <<- c(pendingVars, paste0(tmp, elemNames))
      if(!is.null(listNames)) pendingVars <<- c(pendingVars, paste0(tmp, listNames, "$*"))
    } 
    # If the call is as follows fName(c(settingList, componentName2=...), ...)
    test2 <- test && !is.val(e[[setPos]]) && e[[setPos]][[1]]=='c' && !is.null(names(e[[setPos]]))
    test2 <- test2 && names(e[[setPos]])[3]=='componentName2'
    if(test2 && is.symbol(e[[setPos]][[2]])){ # if first element is a symbol (settingList is defined somewhere else)
      tmp <- paste0(as.character(e[[setPos]][[2]]), "$")
      if(!is.null(elemNames)) pendingVars <<- c(pendingVars, paste0(tmp, elemNames))
      if(!is.null(listNames)) pendingVars <<- c(pendingVars, paste0(tmp, listNames, "$*"))
    }
    if(test2 && is.call(e[[setPos]][[2]]) && e[[setPos]][[2]]!='list()' && e[[setPos]][[2]][[1]]=="list"){ # settingList def in place
      # Turn single elements into functions
      elemPos <- which(names(e[[setPos]][[2]]) %in% elemNames)
      if(length(elemPos)>0) for(i in elemPos) e[[setPos]][[2]][[i]] <- asFunctionDef(e[[setPos]][[2]][[i]])
      # Turn lists into functions
      elemPos <- which(names(e[[setPos]][[2]]) %in% listNames)
      if(length(elemPos)>0) for(i in elemPos){
        if(is.symbol(e[[setPos]][[2]][[i]])){ # If the list is another variable, change it later
          pendingVars <<- c(pendingVars, paste0(as.character(e[[setPos]][[2]][[i]]), "$*"))
        } else e[[setPos]][[2]][[i]] <- asFunctionDef(e[[setPos]][[2]][[i]])
      }
    }
    # If 'e' is not a call to 'fName', look one lever deeper
    for(i in 1:length(e)) if(!is.null(e[[i]])){
      isFuncArg <- i==2 && is.symbol(e[[i-1]]) && as.character(e[[i-1]])=="function"
      if(!isFuncArg) e[[i]] <- replaceInModelCall(fName, elemNames, listNames, e[[i]])
    }
    return(e)
  }
  
  processModelDefinition <- function(fName, elemNames=NULL, listNames=NULL, e=NULL){
    # Store original if e is function
    eOrig <- NULL
    if(is.function(e)){
      eOrig <- e
      e <- body(e)
    } 
    if(is.null(elemNames) & is.null(listNames)) stop('INTERNAL ISSUE - Arguments "elemNames" and "listNames" cannot both be NULL')
    
    # Replace elements in function call
    e <- replaceInModelCall(fName, elemNames, listNames, e)
    
    # Process pending variables
    nIter <- 1
    while(length(pendingVars)>0 && nIter<100){
      # Sort
      pend <- sort(pendingVars) # sort by base list name, e.g. c(L1$a, L1$b$*, L2$a, ...)
      pend <- strsplit(pend, "\\$") # this is a list, e.g. list(c(L1, a), c(L1, b, *), c(L2, a), ...)
      if(any(sapply(pend, length)>3)) stop("INTERNAL ISSUE - There are pending elements with more than 2 levels of depth")
      if(length(pend)>1){
        tmp1 <- unique(sapply(pend, function(p) p[1])) # unique list names
        tmp  <- vector(mode="list", length=length(tmp1))
        for(l in 1:length(tmp)){
          tmp[[l]] <- list(tmp1[l],
                           do.call(c, lapply(pend, function(p) if(p[1]==tmp1[l] && length(p)==2) p[2] else NULL)),
                           do.call(c, lapply(pend, function(p) if(p[1]==tmp1[l] && length(p)==3) p[2] else NULL)))
        }
        pend <- tmp
        rm(tmp, tmp1)
      } else {
        pend <- pend[[1]]
        if(length(pend)==2) tmp2 <- pend[2] else tmp2 <- NULL
        if(length(pend)==3) tmp3 <- pend[2] else tmp3 <- NULL
        pend <- list(pend[1], tmp2, tmp3)
        pend <- list(pend)
        rm(tmp2, tmp3)
      } # pend = list( list(L1, elemNamesL1, elemListNamesL2), list(L2, ...), ...)
      
      # Replace in definition
      for(i in 1:length(pend)) e <- listElem2Func(listNames=pend[[i]][[1]], 
                                                  elemNames=pend[[i]][[2]], 
                                                  elemListNames=pend[[i]][[3]], e=e)
      # Remove items dealt with from pendingVars
      done <- c()
      for(l in 1:length(pend)){
        if(!is.null(pend[[l]][[2]])) done <- c(done, paste0(pend[[l]][[1]], "$", pend[[l]][[2]]))
        if(!is.null(pend[[l]][[3]])) done <- c(done, paste0(pend[[l]][[1]], "$", pend[[l]][[3]], "$*"))
      }
      if(length(done)>0) done <- which(pendingVars %in% done)
      if(length(done)>0) pendingVars <<- pendingVars[-done]
      # Increase iteration count
      nIter <- nIter + 1
    }
    
    # If e was function, make it a function again
    if(!is.null(eOrig)){
      body(eOrig) <- e
      e <- eOrig
    }
    
    return(e)
  }
  
  # Check for while loop. If there is one, return without changes.
  containsWhileLoop <- function(e){
    if(is.function(e)) e <- body(e)
    if(is.symbol(e)){ if(as.character(e)=='while') return(TRUE) else return(FALSE) }
    if(is.val(e)) return(FALSE)
    if(!is.call(e)) stop('INTERNAL ISSUE - Argument "e" must be a language object')
    if(is.call(e) && length(e)>0){
      ans <- rep(FALSE, length(e))
      for(i in 1:length(e)) if(!is.null(e[[i]])){
        isFuncArg <- i==2 && is.symbol(e[[i-1]]) && as.character(e[[i-1]])=="function"
        if(!isFuncArg) ans[i] <- containsWhileLoop(e[[i]])
      }
      return( any(ans) )
    }
  }; if(containsWhileLoop(f)) return(f)
  
  if(randCoeff){
    e <- body(f)
    if(is.val(e)) return(f)
    if(!is.call(e)) stop('INTERNAL ISSUE - The body of argument "f" is not a call')
    ## Introduce 'function ()'
    for(i in 1:length(e)){
      test1 <- is.call(e[[i]]) && length(e[[i]])==2 && e[[i]][[1]]=="return"
      # If return is defined in place
      test2 <- test1 && is.call(e[[i]][[2]]) && e[[i]][[2]]!="list()" && e[[i]][[2]][[1]]=="list"
      if(test2) e[[i]][[2]] <- asFunctionDef(e[[i]][[2]])
      # If return is defined somewhere else
      if(test1 && is.symbol(e[[i]][[2]])) e <- listElem2Func(listNames=as.character(e[[i]][[2]]), 
                                                             elemNames="*", elemListNames=c(), e=e)
      # If there are pending elements to convert
      if((test1 | test2) && length(pendingVars)>0) e <- processModelDefinition("31415", "*", c(), e)
      if(test1 | test2) break
    } 
    body(f) <- e
    
    ## Replace LVs by their definitions on other LVs
    # Run function (apollo_randCoeff) with dummy arguments to figure out names of LVs
    rndCoeff <- tryCatch(f(apollo_beta=c(b1=0, b2=0), apollo_inputs=list()), error=function() NULL)
    test <- is.list(rndCoeff) && all(sapply(rndCoeff, is.function)) && length(rndCoeff)>1
    if(!test) return(f)
    # Replace calls in rndCoeff, to make sure I am replacing for definitions without references to other LVs
    for(i in 2:length(rndCoeff)) rndCoeff[[i]] <- replaceByDef(rndCoeff[[i]], rndCoeff[1:(i-1)], rightSide=TRUE)
    # Replace in function and return
    f <- replaceByDef(f, rndCoeff)
    return(f)
  }
  
  if(lcPars){
    f <- processModelDefinition(fName="apollo_classAlloc", elemNames=c(), listNames=c("V", "utilities"), e=f)
    # Replace elements that are functions
    apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame()), 
                              error=function() list(apollo_control=list(mixing=FALSE), apollo_randCoeff=NA, apollo_lcPars=NA))
    defs <- apollo_varList(f, apollo_inputs)
    if(length(defs)>0) f <- replaceByDef(f, defs)
    return(f)
  }
  
  if(like){
    f <- processModelDefinition(fName="apollo_mnl" , elemNames=c()  , listNames=c("V", "utilities") , e=f)
    f <- processModelDefinition(fName="apollo_nl"  , elemNames=c()  , listNames=c("V", "utilities","nlNests") , e=f)
    f <- processModelDefinition(fName="apollo_fmnl", elemNames=c()  , listNames=c("V", "utilities") , e=f)
    f <- processModelDefinition(fName="apollo_ol"  , elemNames=c("V", "utility"), listNames=c("tau"), e=f)
    f <- processModelDefinition(fName="apollo_op"  , elemNames=c("V", "utility"), listNames=c("tau"), e=f)
    f <- processModelDefinition(fName="apollo_normalDensity", elemNames=c("xNormal", "mu", "sigma"), listNames=c(), e=f)
    f <- processModelDefinition(fName="apollo_tobit", elemNames=c("xTobit", "mu", "sigma"), listNames=c(), e=f)
    f <- processModelDefinition(fName="apollo_rrm_2", elemNames=c() , listNames=c("rum_inputs", "regret_inputs", "regret_scale") , e=f)
    f <- processModelDefinition(fName="apollo_rrm_3", elemNames=c() , listNames=c("rum_inputs", "regret_inputs", "regret_scale") , e=f)
    f <- processModelDefinition(fName="apollo_rrm_4", elemNames=c() , listNames=c("rum_inputs", "regret_inputs", "regret_scale") , e=f)
    f <- processModelDefinition(fName="apollo_ownModel", elemNames=c("likelihood") , listNames=c() , e=f)
    # Replace elements that are functions
    apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame()), 
                              error=function() list(apollo_control=list(mixing=FALSE), apollo_randCoeff=NA, apollo_lcPars=NA))
    defs <- apollo_varList(f, apollo_inputs)
    if(length(defs)>0) f <- replaceByDef(f, defs)
    return(f)
  }
  
  
  ## Testing
  #expre <- quote({
  #  x <- 5
  #  L1 <- list(A=1, B=2, C=3)
  #  L2 <- list()
  #  L2$A <- 1
  #  L2[["B"]] <- 2
  #  apollo_model(functionality="estimate", list(a=1, b=2, c=3, d=L1, e=L2, f=list(f1=1, f2=2, f3=3)))
  #  L4 <- list()
  #  L4$A <- 1
  #  L4[["B"]] <- 2
  #  L5 <- list(a=1, b=2, c=3, d=list(A=1, B=2, C=3), e=L4, f=L1)
  #  apollo_model(L5, functionality="gradient")
  #  })
  #expre <- processModelDefinition("apollo_model", c("a", "c"), c("d", "e", "f"), expre)
  #f2 <- processModelDefinition("apollo_mnl", c(), c("V"), apollo_probabilities)
  #f3 <- processModelDefinition("apollo_ol", c("V"), c("tau"), f2)
}