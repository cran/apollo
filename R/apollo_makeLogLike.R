#' Creates log-likelihood function.
#'
#' Creates log-likelihood function from the likelihood function apollo_probabilities provided by the user.
#'
#' Internal use only. Called by \code{apollo_estimate} before estimation.
#' The returned function can be single-threaded or multi-threaded based on the model options.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param apollo_estSet List of estimation options. It must contain at least one element called
#'                      estimationRoutine defining the estimation algorithm. See \link{apollo_estimate}.
#' @param cl Cluster as provided by \link[parallel]{makeCluster}. Assumed to be PSock.
#' @return apollo_logLike function.
apollo_makeLogLike <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, apollo_estSet, cl=NA){
  ### Copying from apollo_inputs to current environment
  singleCore <- apollo_inputs$apollo_control$nCores==1
  panelData  <- apollo_inputs$apollo_control$panelData
  modelName  <- apollo_inputs$apollo_control$modelName
  bFixedVal  <- apollo_beta[apollo_fixed]
  workInLogs <- apollo_inputs$apollo_control$workInLogs
  estAlg     <- apollo_estSet$estimationRoutine

  ### Iterations count
  nIter <- 0

  ### SINGLE CORE
  if(singleCore){
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      if(getNIter) return(nIter)
      b <- c(bVar, bFixedVal)
      P <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      if(workInLogs) ll <- P  else ll <- log(P) # condition used to be workInLogs & panelData
      if(sumLL) ll <- sum(ll)
      if(countIter && estAlg=="bfgs" && exists('lastFuncParam', envir=globalenv())){
        lastFuncParam <- get("lastFuncParam", envir=globalenv())
        if(!anyNA(lastFuncParam) && all(lastFuncParam==bVar)) nIter <<- nIter + 1
      }
      if(writeIter) apollo_writeTheta(bVar, sum(ll), modelName)
      return(ll)
    }
  } else {
    ### MULTI CORE

    # Clean up parent environment and check that cluster exists
    rm(apollo_inputs, apollo_probabilities, apollo_beta, apollo_fixed, singleCore)
    if(length(cl)==1 && is.na(cl)) stop("Cluster is missing on 'apollo_makeLogLike' despite nCores>1.")

    # Function to be run on each thread
    apollo_parProb <- function(b){
      P <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      return(P)
    }
    environment(apollo_parProb) <- new.env(parent=parent.env(environment(apollo_parProb)))

    # Loglike
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      if(getNIter) return(nIter)
      b <- c(bVar, bFixedVal)
      P <- parallel::clusterCall(cl=cl, fun=apollo_parProb, b=b)
      P <- unlist(P)
      if(workInLogs) ll <- P  else ll <- log(P) # condition used to be workInLogs & panelData
      if(countIter && estAlg=="bfgs" && exists('lastFuncParam', envir=globalenv())){
        lastFuncParam <- get("lastFuncParam", envir=globalenv())
        if(!anyNA(lastFuncParam) && all(lastFuncParam==bVar)) nIter <<- nIter + 1
      }
      if(writeIter) apollo_writeTheta(bVar, sum(ll), modelName)
      if(sumLL) ll <- sum(ll)
      return(ll)
    }

  }
  
  return(apollo_logLike)
}
