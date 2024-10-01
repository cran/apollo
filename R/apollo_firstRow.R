#' Keeps only the first row for each individual
#' 
#' Given a multi-row input, keeps only the first row for each individual.
#'
#' This a function to keep only the first row of an object per indidividual. It can handle multiple types of components, including scalars, vectors and three-dimensional arrays (cubes).
#' The argument database MUST contain a column called 'apollo_sequence', which is created by \link{apollo_validateData}.
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components (or other object).
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return If \code{P} is a list, then it returns a list where each element has only the first row of each individual.
#'         If \code{P} is a single element, then it returns a single element with only the first row of each individual.
#'         The size of the element is changed only in the first dimension. If input is a scalar, then it returns a vector with the element repeated as many
#'         times as individuals in \code{database}. If the element is a vector, its length will be changed to the number of individuals. If the element is
#'         a matrix, then its first dimension will be changed to the number of individuals, while keeping the size of the second dimension. If the element 
#'         is a cube, then only the first dimension's length is changed, preserving the others.
#' @export
apollo_firstRow=function(P, apollo_inputs){
  ### Functionalities for those where nothing is done
  functionality <- tryCatch(get('functionality', envir=parent.frame(), inherits=FALSE),
                            error=function(e) NULL)
  if(!is.null(functionality) && functionality %in% c('preprocess', 'report')) return(P)
  
  ### Calculate useful values
  firstR <- apollo_inputs$database$apollo_sequence==1
  nObs <- length(firstR)
  K <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                error=function(e) NULL)
  if(!is.null(K)) K <- length(K) - length(apollo_inputs$apollo_fixed)
  
  ### If P is a data.frame
  if(is.data.frame(P)){
    if(nrow(P)==nObs) return(subset(P, firstR)) else {
      stop("INPUT ISSUE - The data.frame passed to apollo_firstRow ",
           "does not have the same number of rows as the database!")
    }
  }
  
  ### Define recursive function
  f <- function(x, validNRow, firstR){
    # If list, apply recursively
    if(is.list(x)) return(lapply(x, f, validNRow=validNRow, firstR=firstR))
    # If not numeric or logical, return as is
    if(!is.numeric(x) && !is.logical(x)) return(x)
    # If not a list and not numeric
    isSca <- length(x)==1
    isVec <- is.vector(x) && !is.array(x)
    isMat <- is.matrix(x)
    isCub <- is.array(x) && !isMat && length(dim(x))==3
    test <- isSca || isVec || isMat || isCub
    if(!test) stop("INPUT ISSUE - Numeric object given to apollo_firstRow ",
                   "is not a scalar, vector, matrix or cube.")
    if(isSca | isVec) nRows <- length(x) else nRows <- nrow(x)
    if(!(nRows %in% validNRow)) stop("INPUT ISSUE - The object passed to ", 
                                     "apollo_firstRow does not have a valid ", 
                                     "number of rows (usually the same as ", 
                                     "the number of rows in the database).")
    if(isSca) x <- rep(x, sum(firstR))
    if(isVec|isMat) x <- subset(x, firstR)
    if(isCub) x = x[firstR,,, drop=FALSE]
    return(x)
  }
  
  ### Keep only
  if(is.null(functionality) || functionality!="gradient") K <- NULL
  validNRow <- c(1,nObs,K)
  P <- f(P, validNRow, firstR)
  return(P)
}