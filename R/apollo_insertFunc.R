#' Modifies function to make it compatible with analytic gradients
#' 
#' Takes a likelihood function and inserts \code{function ()} before key elements to allow for analytic gradient calculation
#' 
#' It modifies the definition of the following models.
#' \itemize{
#'   \item \code{apollo_mnl}: Turns all elements inside \code{mnl_settings$V} into functions.
#'   \item \code{apollo_ol}: Turns \code{ol_settings$V} and all elements inside \code{ol_settings$tau} into functions.
#'   \item \code{apollo_op}: Turns \code{op_settings$V} and all elements inside \code{op_settings$tau} into functions.
#'   \item \code{apollo_normalDensity}: Turns \code{normalDensity_settings$xNormal}, \code{normalDensity_settings$mu} and \code{normalDensity_settings$sigma} into functions.
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
#' @return Function \code{f} but with relevant expressions turned into function definitions.
#' @importFrom utils capture.output
#' @export
apollo_insertFunc <- function(f, like=TRUE, randCoeff=FALSE){
  # Validate inputs
  if(!is.function(f)) stop('Argument "f" should be a function.')
  
  ## Check if input is a single value
  is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) ) return(TRUE) else return(FALSE)
  
  # Vector to store the elements that are yet to be turned into functions
  pendingVars <- c()
  
  ## Replace definition by function
  # If e is a simple call, it returns 'function() e'
  # If e is a list definition, it returns 'list(function() e[[2]], function() e[[3]], ...)'.
  # If e starts with "function" or "eval" it does not make any changes.
  asFunctionDef <- function(e){
    isCall <- is.call(e)
    isVal  <- is.val(e)
    if(!isCall & !isVal) stop('Argument "e" must be a language object (a call or a symbol).')
    if(isCall && (e[[1]]=="function" || e[[1]]=="eval")) return(e)
    if(isCall && e[[1]]=="list"){
      for(i in 2:length(e)) e[[i]] <- asFunctionDef(e[[i]])
    } else {
      if(isVal) e <- as.character(e) else e <- paste0(capture.output(print(e)), collapse="")
      e <- str2lang(paste0("function() ", e))
    }
    return(e)
  }
  
  ## Searches inside "e" for the definition of listNames, and turns its elements "elemNames" into
  ## functions, as well as all the elements inside elemListNames.
  ## If elemNames="*", then all elements of listNames are turned into functions.
  ## If an elemListNames is defined as another variable, then it adds it to the pendingVars vector
  listElem2Func <- function(listNames, elemNames=NULL, elemListNames=NULL, e=NULL, from=1, to=Inf){
    if(is.null(e)) stop('Argument "e" cannot be NULL.')
    if(is.function(e)) e <- body(e)
    if(is.val(e) || from>length(e)) return(e)
    if(!is.call(e)) stop('Argument "e" must be a call')
    if(is.null(elemNames) && is.null(elemListNames)) stop('Arguments "elemNames" and "elemListNames" cannot be both NULL.')
    if(!is.character(listNames)) stop('Argument "listNames" must be a character vector')
    if(!is.null(elemNames) && !is.character(elemNames)) stop('Argument "elemNames", if provided, must be a character vector')
    if(!is.null(elemListNames) && !is.character(elemListNames)) stop('Argument "elemListNames", if provided, must be a character vector')
    to <- min(to, length(e))
    for(i in from:to){
      # Is an assignment
      test0 <- length(e[[i]])>=3
      test0 <- test0 && (e[[i]][[1]]=="=" || e[[i]][[1]]=="<-")
      # Assignment of type L$x <- ...
      test1 <- test0 && length(e[[i]][[2]])==3
      test1 <- test1 && (e[[i]][[2]][[1]]=="[[" || e[[i]][[2]][[1]]=="$")
      test1 <- test1 && as.character(e[[i]][[2]][[2]])[1] %in% listNames
      # L$x <- ... where x is a single elements
      test21 <- test1 && as.character(e[[i]][[2]][[3]])[1] %in% elemNames
      test22 <- test1 && ("*" %in% elemNames) && !(as.character(e[[i]][[2]][[3]])[1] %in% elemListNames)
      if(test21 || test22) e[[i]][[3]] <- asFunctionDef(e[[i]][[3]])
      # L$x <- ... where x is a list defined somewhere else
      test3 <- test1 && as.character(e[[i]][[2]][[3]])[1] %in% elemListNames
      test31<- test3 && is.symbol(e[[i]][[3]])
      if(test31) pendingVars <<- c(pendingVars, paste0(as.character(e[[i]][[3]]), "$*"))
      # L$x <- ... where x is a list defined in place
      test32 <- test3 && is.call(e[[i]][[3]]) && length(e[[i]][[3]])>=2 && e[[i]][[3]][[1]]=="list"
      if(test32) e[[i]][[3]] <- asFunctionDef(e[[i]][[3]])
      # Assignment of type L <- list(x=...)
      test4 <- test0 && as.character(e[[i]][[2]])[1] %in% listNames
      test4 <- test4 && !(is.call(e[[i]][[3]]) && e[[i]][[3]]=="list()") # not assigning empty list
      test4 <- test4 && is.call(e[[i]][[3]]) && e[[i]][[3]][[1]]=="list" # creating a list
      if(test4){
        tmp <- names(e[[i]][[3]]) # names of elements in list definition (including "" at the beginning)
        # L <- list(x=...) when x is a single element
        if("*" %in% elemNames) tmp2 <- 2:length(tmp) else tmp2 <- which(tmp %in% elemNames)
        if(length(tmp2)>0) for(j in tmp2) if(!(tmp[j] %in% elemListNames)) e[[i]][[3]][[j]] <- asFunctionDef(e[[i]][[3]][[j]])
        # L <- list(x=...) when x is a list 
        tmp2 <- which(tmp %in% elemListNames)
        if(length(tmp2)>0) for(j in tmp2){
          # L <- list(x=y) when x is a list and y is defined somewhere else
          if(is.symbol(e[[i]][[3]][[j]])) pendingVars <<- c(pendingVars, paste0(as.character(e[[i]][[3]][[j]]), "$*"))
          # L <- list(x=list(...)) when x is a list defined in place
          test5 <- is.call(e[[i]][[3]][[j]]) && e[[i]][[3]][[j]]!="list()" && length(e[[i]][[3]][[j]])>=2
          test5 <- test5 && e[[i]][[3]][[j]][[1]]=="list"
          if(test5) e[[i]][[3]][[j]] <- asFunctionDef(e[[i]][[3]][[j]])
        }
      }
      # Other situation: look one level deeper
      if(!test1 && !test21 && !test22 && !test31 && !test32 && !test4 && is.call(e[[i]]) && length(e[[i]])>1){
        for(j in 1:length(e[[i]])) if(!is.null(e[[i]][[j]])) e[[i]][[j]] <- listElem2Func(listNames, 
                                                                                          elemNames, 
                                                                                          elemListNames, 
                                                                                          e=e[[i]][[j]])
      }
    }
    return(e)
  }
  
  ## Searches for calls to function 'fName'. If the call contains as its arguments the definition 
  ## of a list, searches within that list for elements named 'elemNames' and turn its definition
  ## into functions. Returns the expression 'e' altered in the way just described.
  replaceInModelCall <- function(fName, elemNames=NULL, listNames=NULL, e=NULL){
    if(is.function(e)) e <- body(e)
    if(is.val(e)) return(e)
    if(!is.call(e)) stop('Argument "e" must be a language object')
    if(is.null(elemNames) & is.null(listNames)) return(e)
    if(!is.null(elemNames) && !is.character(elemNames)) stop('Argument "elemNames" must be a character vector')
    if(!is.null(listNames) && !is.character(listNames)) stop('Argument "listNames" must be a character vector')
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
    for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- replaceInModelCall(fName, elemNames, listNames, e[[i]])
    return(e)
  }
  
  processModelDefinition <- function(fName, elemNames=NULL, listNames=NULL, e=NULL){
    # Store original if e is function
    eOrig <- NULL
    if(is.function(e)){
      eOrig <- e
      e <- body(e)
    } 
    if(is.null(elemNames) & is.null(listNames)) stop('Arguments "elemNames" and "listNames" cannot both be NULL')
    
    # Replace elements in function call
    e <- replaceInModelCall(fName, elemNames, listNames, e)
    
    # Process pending variables
    nIter <- 1
    while(length(pendingVars)>0 && nIter<100){
      # Sort
      pend <- sort(pendingVars) # sort by base list name, e.g. c(L1$a, L1$b$*, L2$a, ...)
      pend <- strsplit(pend, "\\$") # this is a list, e.g. list(c(L1, a), c(L1, b, *), c(L2, a), ...)
      if(any(sapply(pend, length)>3)) stop("There are pending elements with more than 2 levels of depth")
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
    
    # If e was function, make itfunction again
    if(!is.null(eOrig)){
      body(eOrig) <- e
      e <- eOrig
    }
    
    return(e)
  }
  
  if(randCoeff){
    e <- body(f)
    if(is.val(e)) return(f)
    if(!is.call(e)) stop('The body of argument "f" is not a call')
    # Search for name of return variable
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
    
    # Run function (apollo_randCoeff) with dummy arguments to figure out names of LVs
    rndCoeff <- tryCatch(f(apollo_beta=c(b1=0, b2=0), apollo_inputs=list()), error=function() NULL)
    test <- is.list(rndCoeff) && all(sapply(rndCoeff, is.function))
    if(!test) return(f)
    ## replace LVs on the right side by their definitions
    replaceByDef <- function(e, defs, rightSide=FALSE){
      if(is.function(e)){ f <- e; e <- body(e)} else f <- NULL
      # Case 1: LV1
      test1 <- rightSide && is.symbol(e) && (as.character(e) %in% names(defs))
      if(test1) e <- body(defs[[ which(names(defs)==as.character(e))[1] ]])
      # Case 2: L$LV1
      test2 <- !test1 && rightSide && is.call(e) && length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('$','[['))
      test2 <- test2 && is.val(e[[3]]) && (as.character(e[[3]]) %in% names(defs))
      if(test2) e <- body(defs[[ which(names(defs)==as.character(e[[3]]))[1] ]])
      # Case 3: expression
      if(!test1 && !test2 && is.call(e)){
        test0 <- length(e)==3 && is.symbol(e[[1]]) && (as.character(e[[1]]) %in% c('<-', '=')) # Is an assignment
        for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- replaceByDef(e[[i]], defs, rightSide=(rightSide | (test0 & i==3)))
      } 
      # Return
      if(is.null(f)) return(e) else {body(f) <- e; return(f)}
    }
    f <- replaceByDef(f, rndCoeff)
    
    return(f)
  }
  
  if(like){
    f <- processModelDefinition(fName="apollo_mnl", elemNames=c()   , listNames=c("V")  , e=f)
    f <- processModelDefinition(fName="apollo_ol" , elemNames=c("V"), listNames=c("tau"), e=f)
    f <- processModelDefinition(fName="apollo_op" , elemNames=c("V"), listNames=c("tau"), e=f)
    f <- processModelDefinition(fName="apollo_normalDensity", elemNames=c("xNormal", "mu", "sigma"), listNames=c(), e=f)
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