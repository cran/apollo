#' Applies weights
#' 
#' Applies weights to individual observations in likelihood function.
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call \code{apollo_probabilities}, though the user can also call \code{apollo_probabilities} manually with a given functionality for testing/debugging. Possible values are:
#'                      \itemize{
#'                        \item \strong{\code{"components"}}: For further processing/debugging, produces likelihood for each model component (if multiple components are present), at the level of individual draws and observations.
#'                        \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
#'                        \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"gradient"}}: For model estimation, produces analytical gradients of the likelihood, where possible.
#'                        \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
#'                        \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
#'                        \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
#'                        \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
#'                        \item \strong{\code{"report"}}: Prepares output summarising model and choiceset structure.
#'                        \item \strong{\code{"shares_LL"}}: Produces overall model likelihood with constants only.
#'                        \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"zero_LL"}}: Produces overall model likelihood with all parameters at zero.
#'                      }
#' @return The likelihood (i.e. probability in the case of choice models) of the model in the appropriate form for the 
#'         given functionality, multiplied by individual-specific weights.
#' @export
apollo_weighting=function(P, apollo_inputs, functionality){
  
  # ################################## #
  #### preprocess, report           ####
  # ################################## #
  if(functionality %in% c("preprocess", "report")) return(P)
  
  # ####################################### #
  #### Basic checks and useful variables ####
  # ####################################### #
  if(is.null(apollo_inputs$apollo_control$weights)) stop('INCORRECT FUNCTION/SETTING USE - Call to apollo_weighting performed without a weights variable defined in apollo_control!')
  w     <- apollo_inputs$database[, apollo_inputs$apollo_control$weights]
  wInd  <- apollo_firstRow(w, apollo_inputs)
  nObs  <- nrow(apollo_inputs$database)
  indiv <- apollo_inputs$database[, apollo_inputs$apollo_control$indivID]
  nInd  <- length(unique(indiv))
  ### weighting function
  wf <- function(p){
    if(is.list(p)) return( lapply(p, wf) )
    if(all(is.na(p))) return(p) ### added to deal with cases where P fails anyway, e.g. sharesLL for some models
    isVec <- is.vector(p)
    isMat <- is.matrix(p)
    isCub <- is.array(p) && length(dim(p))==3
    if(!isVec && !isMat && !isCub) stop('INTERNAL ISSUE - An element of P is not a numeric vector, matrix or cube (3-dim array)')
    if(isVec) nR <- length(p) else nR <- dim(p)[1]
    #if(nR==nObs){ if(!apollo_inputs$apollo_control$workInLogs) return(p^w) else return(w*p) }
    #if(nR==nObs) return(p^w) 3
    if(nR==nObs){ if(!apollo_inputs$apollo_control$workInLogs) return(p^w) else return(w*p) }# BUG FIX DP 08/03/2023
    if(nR==nInd){ if(!apollo_inputs$apollo_control$workInLogs) return(p^wInd) else return(wInd*p) }
    warning('An element of P did not have as many rows as observations or individuals, so it was not weighted.')
    return(p)
  }
  ### number of rows function
  nRows <- function(p){
    ans <- -1
    if(is.list(p)  ) ans <- sapply(p, nRows)
    if(is.array(p) ) ans <- dim(p)[1]
    if(is.vector(p)) ans <- length(p)
    if(ans==-1) stop('INTERNAL ISSUE - An element of P did not have as many rows as observations or individuals, so it was not weighted.')
    return(ans)
  }
  
  # ############## #
  #### validate #### 
  # ############## #
  # Check dimensionality
  if(functionality=='validate'){
    nRows <- unlist(lapply(P, nRows))
    iInd  <- apollo_firstRow(indiv, apollo_inputs)
    test1 <- all(nRows %in% c(1, nObs, nInd))
    if(!test1) stop('INTERNAL ISSUE - Some elements in "P" have the wrong number of rows (different to 1, nIndiv and nObs).')
    txt <- 'INTERNAL ISSUE - When applying weights at the individual level, weights should be the same for all observations of each individual.'
    if(any(nRows==nInd)) for(i in 1:nInd) if(!all(w[indiv==iInd[i]]==wInd[i])) stop(txt)
    rm(txt)
    P <- lapply(P, wf)
    return(P)
  }
  
  
  # ############################################################## #
  #### estimate, zero_LL, conditionals, output, raw, components #### 
  # ############################################################## #
  if(functionality %in% c('estimate', 'zero_LL','shares_LL', 'conditionals', 'output', 'components')){
    P <- lapply(P, wf)
    return(P)
  }
  
  # ################ #
  #### prediction #### 
  # ################ #
  if(functionality=="prediction"){
    wf2 <- function(p){
      if(is.list(p)) return( lapply(p, wf2) )
      isVec <- is.vector(p)
      isMat <- is.matrix(p)
      isCub <- is.array(p) && length(dim(p))==3
      if(!isVec && !isMat && !isCub) stop('INTERNAL ISSUE - An element of P is not a numeric vector, matrix or cube (3-dim array)')
      if(isVec) nR <- length(p) else nR <- dim(p)[1]
      if(nR==nObs) return(w*p)
      if(nR==nInd) return(wInd*p)
      warning('An element of P did not have as many rows as observations or individuals, so it was not weighted.')
      return(p)
    }
    P <- lapply(P, wf2)
    return(P)
  }
  
  # ############## #
  #### gradient #### 
  # ############## #
  if(functionality=='gradient'){
    # preliminary check
    if(!is.list(P)) stop("INTERNAL ISSUE - Input P should be a list with at least one component (called model)!")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("INTERNAL ISSUE - Setting workInLogs cannot be used in conjunction with analyticGrad!")
    if(!is.null(names(P)) && all(names(P) %in% c('like', 'grad'))){
      # there are no model components
      if(is.null(P$like) || is.null(P$grad)) stop("INTERNAL ISSUE - Missing like and/or grad elements inside components when calculating gradient!")
      nR <- nRows(P$like); W <- 1
      if(nR==nObs) W <- w
      if(nR==nInd) W <- wInd
      tmp <- W*P$like^(W-1)
      P$grad <- lapply(P$grad, '*', tmp)
      P$like <- wf(P$like)
    } else {
      # there are multiple model components
      for(m in 1:length(P)){
        if(is.null(P[[m]]$like) || is.null(P[[m]]$grad)) stop("INTERNAL ISSUE - Missing like and/or grad elements inside components when calculating gradient!")
        nR <- nRows(P$like); W <- 1
        if(nR==nObs) W <- w
        if(nR==nInd) W <- wInd
        tmp <- W*P[[m]]$like^(W-1)
        P[[m]]$grad <- lapply(P[[m]]$grad, '*', tmp)
        P[[m]]$like <- wf(P[[m]]$like)
      }
    }
    return(P)
  }
  
}