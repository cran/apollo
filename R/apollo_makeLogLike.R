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
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param apollo_estSet List of estimation options. It must contain at least one element called
#'                      estimationRoutine defining the estimation algorithm. See \link{apollo_estimate}.
#' @param cleanMemory Logical. If TRUE, then \code{apollo_inputs$draws} and \code{apollo_inputs$database} are erased
#'                    throughout the calling stack. Used to reduce memory usage in case of multithreading and a large
#'                    database or number o draws.
#' @return apollo_logLike function.
#' @export
apollo_makeLogLike <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                               apollo_estSet=list(estimationRoutine='bgw'), cleanMemory=FALSE){
  
  if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
  if(!is.null(apollo_inputs$apollo_control$debug)) debug <- apollo_inputs$apollo_control$debug else debug <- FALSE
  
  ### Copy (the modified) apollo_probabilities inside apollo_inputs
  apollo_inputs$apollo_probabilities <- apollo_probabilities # Not sure if this is necessary
  
  # # # # # # # # # # # # #
  #### Create workers  ####
  # # # # # # # # # # # # #
  
  if(apollo_inputs$apollo_control$nCores==1) cl <- NA else {
    if(!silent) apollo_print('Creating cluster...')
    # Creates cluster and also deletes database and draws from apollo_inputs in 
    # here and in .GlobalEnv
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent, 
                             cleanMemory=cleanMemory)
    apollo_inputs$apollo_control$nCores <- length(cl)
  }
  
  # # # # # # # # # # # # # # # #
  #### Create apollo_logLike ####
  # # # # # # # # # # # # # # # #
  
  if(debug) apollo_print('Creating apollo_logLike... ')
  
  ### Copying from apollo_inputs to current environment
  singleCore <- apollo_inputs$apollo_control$nCores==1
  panelData  <- apollo_inputs$apollo_control$panelData
  modelName  <- apollo_inputs$apollo_control$modelName
  bFixedVal  <- apollo_beta[apollo_fixed]
  workInLogs <- apollo_inputs$apollo_control$workInLogs
  estAlg     <- apollo_estSet$estimationRoutine
  analyticGrad    <- apollo_inputs$apollo_control$analyticGrad
  apollo_scaling  <- apollo_inputs$apollo_scaling
  outputDirectory <- apollo_inputs$apollo_control$outputDirectory
  
  ### Iterations count
  if(analyticGrad) nIter <- 1 else nIter <- -1
  preLL <- -Inf
  tmp <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
  lastB <- setNames(rnorm(length(tmp)), tmp)
  rm(tmp)

  ### Create apollo_logLike
  if(singleCore){
    ### SINGLE CORE
    
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE, logP=TRUE){
      if(sumLL && !logP) stop("INTERNAL ISSUE - sumLL=TRUE cannot be used alongside logP=FALSE")
      # Return iteration count, if requested
      if(getNIter) return(nIter)
      # Calculate LL
      b  <- c(bVar, bFixedVal)
      ll <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      if(!workInLogs && logP) ll <- log(ll) # condition used to be workInLogs & panelData
      if(workInLogs && !logP) ll <- exp(ll) # 
      sumll <- sum(ll)
      # Keep count of iterations and write them to file, if requested
      if(countIter || writeIter){
        if(analyticGrad){
          # If using analytic gradient, assume a new iteration if the LL improved
          newIter <- is.finite(sumll) && (sumll > preLL)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter){
            if(!logP){
              apollo_writeTheta(b, sum(log(ll)), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)
            }else{
              apollo_writeTheta(b, sum(ll), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)  
            }
          } 
          if(newIter) preLL <<- sumll
        } else {
          # If using numeric gradient, assume new iteration if the parameters haven't changed since last call
          newIter <- all(bVar==lastB)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter){
            if(!logP){
              apollo_writeTheta(b, sum(log(ll)), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)
            }else{
              apollo_writeTheta(b, sum(ll), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)  
            }
          }
          lastB <<- bVar
        }
      }
      # Return LL
      if(sumLL) return(sumll) else return(ll)
    }
  } else {
    ### MULTI CORE
    
    # Clean up parent environment and check that cluster exists
    rm(apollo_inputs, apollo_probabilities, apollo_fixed, apollo_estSet, panelData)
    if(length(cl)==1 && is.na(cl)) stop("INTERNAL ISSUE - Cluster is missing on 'apollo_makeLogLike' despite nCores>1.")

    # Function to be run on each thread
    apollo_parProb <- function(b){
      P <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      return(P)
    }
    environment(apollo_parProb) <- new.env(parent=environment())
    
    # Loglike
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE, logP=TRUE){
      if(sumLL && !logP) stop("INTERNAL ISSUE - sumLL=TRUE cannot be used alongside logP=FALSE")
      # Return iteration count, if requested
      if(getNIter) return(nIter)
      # Calculate LL
      b  <- c(bVar, bFixedVal)
      ll <- parallel::clusterCall(cl=cl, fun=apollo_parProb, b=b)
      ll <- unlist(ll)
      if(!workInLogs && logP) ll <- log(ll) # condition used to be workInLogs & panelData
      if(workInLogs && !logP) ll <- exp(ll) # 
      sumll <- sum(ll)
      # Keep count of iterations and write them to file, if requested
      if(countIter || writeIter){
        if(analyticGrad){
          # If using analytic gradient, assume a new iteration if the LL improved
          newIter <- is.finite(sumll) && (sumll > preLL)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter){
            if(!logP){
              apollo_writeTheta(b, sum(log(ll)), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)
            }else{
              apollo_writeTheta(b, sum(ll), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)  
            }
          }
          if(newIter) preLL <<- sumll
        } else {
          # If using numeric gradient, assume new iteration if the parameters haven't changed since last call
          newIter <- all(bVar==lastB)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter){
            if(!logP){
              apollo_writeTheta(b, sum(log(ll)), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)
            }else{
              apollo_writeTheta(b, sum(ll), modelName, 
                                apollo_scaling, outputDirectory, apollo_beta)  
            }
          }
          lastB <<- bVar
        }
      }
      # Return LL
      if(sumLL) return(sumll) else return(ll)
    }

  }
  
  # # # #  # # # # # # # # # # # # # # # # #
  #### Pre-process apollo_probabilities ####
  # # # #  # # # # # # # # # # # # # # # # #
  
  ### Preprocess model (for optimisation)
  if(debug) apollo_print('Pre-processing model...')
  if(singleCore){
    preprocess_outputs = tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="preprocess"),
                                  error=function(e) NULL)
    if(is.null(preprocess_outputs) & debug) apollo_print('Storing preprocessing failed. It will be repeated in each iteration.')
    apollo_inputs      = c(apollo_inputs, preprocess_outputs)
    assign("apollo_inputs", apollo_inputs, envir=environment(apollo_logLike))
    rm(preprocess_outputs)
  } else if(!anyNA(cl)){
    b <- apollo_beta
    parallel::clusterExport(cl, "b", envir=environment())
    parallel::clusterEvalQ(cl, {
      po = tryCatch(apollo_probabilities(b, apollo_inputs, functionality="preprocess"),
                    error=function(e) NULL)
      if(is.null(po)) apollo_print('Storing preprocessing failed. It will be repeated in each iteration.')
      po = c(apollo_inputs, po)
      tmp <- globalenv()
      assign("apollo_inputs", po, envir=tmp)
      rm(po)
      rm(b, envir=globalenv())
    } )
    rm(b)
  }
  
  #### Pre-processing report
  # fetch reports
  if(debug) apollo_print('Preparing pre-processing report')
  f <- function(){
    name <- grep("_settings$", names(apollo_inputs), value=TRUE)
    grad  <- c(); mType <- c(); countAlt  <- c()
    for(i in name){
      grad  <- c(grad , apollo_inputs[[i]]$gradient)
      mType <- c(mType, apollo_inputs[[i]]$modelType)
      ### 26 April
      countAlt <- c(countAlt, apollo_inputs[[i]]$nAlt) 
      ###
    }
    name <- substr(name, 1, nchar(name)-nchar("_settings"))
    return( list(name=name, grad=grad, mType=mType, countAlt = countAlt) )
  }
  if(!anyNA(cl)) tmp <- parallel::clusterCall(cl, function(f) {environment(f) <- globalenv(); f()}, f)[[1]] else tmp <- f()
  if(length(tmp$name)>0 && length(tmp$name)==length(tmp$mType)) mType <- setNames(tmp$mType, tmp$name) else mType <- NULL
  ### 26 April
  # models with nAlt - exclude LC from this for now, should do as the models are class specific, but will not work if the class specific models are not in P
  if(!is.null(mType)&&!is.null(tmp$countAlt)){
    #countAlt <- setNames(tmp$countAlt, tmp$name[tolower(mType) %in% c("mnl", "nl", "cnl", "el", "dft", "lc", "rrm", "ol", "op")])
    countAlt <- setNames(tmp$countAlt, tmp$name[tolower(mType) %in% c("mnl", "nl", "cnl", "el", "dft", "rrm", "ol", "op")]) 
  }
  ###
  rm(f)
  if(debug){
    if(length(tmp$name)==0) apollo_print("No pre-processing performed") else {
      # Print report
      cat("ComponentName  Type   Gradient Optimisation\n") # cat, as apollo_print messes up the alignment
      for(i in 1:length(tmp$grad)){
        txt <- 14 - nchar(tmp$name[i])
        txt <- ifelse(txt>=0, paste0(tmp$name[i], paste0(rep(" ", txt), collapse="")), substr(tmp$name[i], 1, 14))
        txt <- paste0(txt, " ", substr(paste0(tmp$mType[i], '      '), 1, 6))
        txt <- paste0(txt, " ", ifelse(tmp$grad[i], "analytic\n", "numeric \n"))
        cat(txt) # cat, as apollo_print messes up the alignment
      }; rm(i, txt)
      apollo_print(paste0("Whole model gradient function creation ", ifelse(is.null(tmp$grad), "failed.", "succeeded.")))
    }
    rm(tmp)
  }
  
  
  ### Count number of observations
  if(debug) cat('Counting number of observations...')
  f <- function(){ # this function returns a vector c(component1Name=nObs1, component2Name=nObs2, ...) ignoring LC components
    name    <- grep("_settings$", names(apollo_inputs), value=TRUE)
    nObsTot <- setNames(rep(0, length(name)), name)
    discard <- c()
    for(i in name){
      if(!is.null(apollo_inputs[[i]][['LCNObs']])){ # if it is a latent class component
        nObsTot[i] <- apollo_inputs[[i]]$LCNObs
        discard <- c(discard, apollo_inputs[[i]]$LCCompNames)
      } else nObsTot[i] <- sum(apollo_inputs[[i]]$rows) # if it is a non-LC component
    }
    names(nObsTot) <- substr(name, 1, nchar(name)-nchar("_settings"))
    discard <- which(names(nObsTot) %in% discard)
    if(length(discard)>0) nObsTot <- nObsTot[-discard]
    return(nObsTot)
  }
  if(!anyNA(cl)){
    nObsTot <- parallel::clusterCall(cl, function(f) {environment(f) <- globalenv(); f()}, f)
    if(!all(sapply(nObsTot, length)==length(nObsTot[[1]]))) stop('INTERNAL ISSUE - Model components are inconsistent across workers. Try seting apollo_control$nCores=1')
    nObsTot <- Reduce('+', nObsTot)
  } else nObsTot <- f()
  rm(f)
  if(debug) cat(' Done.\n')
  
  return(apollo_logLike)
}
