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
                               apollo_estSet=list(estimationRoutine='BFGS'), cleanMemory=FALSE){
  
  if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
  if(!is.null(apollo_inputs$apollo_control$debug)) debug <- apollo_inputs$apollo_control$debug else debug <- FALSE
  
  # # # #  # # # # 
  #### Checks ####
  # # # #  # # # # 
  
  ### Check that names of params in apollo_beta, apollo_randCoeff & apollo_lcPars are not re-defined
  tmp <- as.character(body(apollo_probabilities))
  tmp <- gsub("(", "", tmp, fixed=TRUE)
  tmp <- gsub(")", "", tmp, fixed=TRUE)
  tmp <- gsub(" ", "", tmp, fixed=TRUE)
  # check for apollo_beta
  for(i in names(apollo_beta)){
    test <- grep(paste0("^",i,"="), tmp)
    test <- c(test, grep(paste0("^",i,"<-"), tmp))
    if(length(test)>0) stop("Parameter ", i, " from apollo_beta was re-defined ",
                            "inside apollo_probabilities. This is not allowed.")
  }; rm(i, test)
  #check for apollo_randCoeff
  if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
    env <- c(apollo_inputs$database, apollo_inputs$draws, as.list(apollo_beta))
    env <- list2env(env, hash=TRUE, parent=parent.frame())
    rnd <- apollo_inputs$apollo_randCoeff; environment(rnd) <- env
    rnd <- rnd(apollo_beta, apollo_inputs)
    for(i in names(rnd)){
      test <- grep(paste0('^', i, '=|^', i, '<-'), tmp)
      if(length(test)>0) stop("Parameter ", i, " from apollo_randCoeff was re-defined ",
                              "inside apollo_probabilities. This is not allowed.")
    }; rm(env, i, test)
  }
  #check for apollo_lcPars
  if(is.function(apollo_inputs$apollo_lcPars)){
    env <- c(apollo_inputs$database, as.list(apollo_beta))
    if(exists('rnd', inherits=FALSE)) env <- c(env, apollo_inputs$draws, rnd)
    env <- list2env(env, hash=TRUE, parent=parent.frame())
    lcp <- apollo_inputs$apollo_lcPars; environment(lcp) <- env
    lcp <- names(lcp(apollo_beta, apollo_inputs))
    for(i in lcp){
      test <- grep(paste0('^', i, '=|^', i, '<-'), tmp)
      if(length(test)>0) stop("Parameter ", i, " from apollo_lcPars was re-defined ",
                              "inside apollo_probabilities. This is not allowed.")
    }; rm(env, lcp, i, test)
  }; if(exists('rnd')) rm(rnd)
  rm(tmp)
  
  ### Check there are no references to database inside apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- gsub("apollo_inputs$database", " ", tmp, fixed=TRUE)
    tmp <- grep("database", tmp, fixed=TRUE)
    if(length(tmp)>0) stop("The database object is 'attached' and elements should thus be called",
                           " directly in apollo_probabilities without the 'database$' prefix.")
    rm(tmp)
  }
  
  ### Check apollo_weighting is called if apollo_control$weights are defined (unless apollo_inputs$EM is TRUE)
  w <- apollo_inputs$apollo_control[['weights']]
  test <- is.null(apollo_inputs$EM) || (is.logical(apollo_inputs$EM) && !apollo_inputs$EM)
  test <- test && !is.null(w) && !is.null(apollo_inputs$database) && (w %in% names(apollo_inputs$database))
  if(test){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- grep('apollo_weighting', tmp, fixed=TRUE)
    if(length(tmp)==0) stop('When using weights, apollo_weighting should be called inside apollo_probabilities.')
    rm(tmp)
  }; rm(w)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #### Modify apollo_probabilities, apollo_randCoeff & apollo_lcPars ####
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if(debug) apollo_print('Inserting function() in user-defined functions')
  
  ### Insert "function ()" in apollo_probabilities, apollo_randCoeff and apollo_lcPars (if needed)
  if(apollo_inputs$apollo_control$analyticGrad){
    tmp1 <- apollo_insertFunc(apollo_probabilities, like=TRUE)
    test <- !is.null(apollo_inputs$apollo_randCoeff) && is.function(apollo_inputs$apollo_randCoeff)
    if(test) tmp2 <- apollo_insertFunc(apollo_inputs$apollo_randCoeff, randCoeff=TRUE)
    test <- !is.null(apollo_inputs$apollo_lcPars) && is.function(apollo_inputs$apollo_lcPars)
    if(test) tmp3 <- apollo_insertFunc(apollo_inputs$apollo_lcPars, lcPars=TRUE)
    # Test that modified LL actually works
    ll <- tryCatch(tmp1(apollo_beta, apollo_inputs), error=function(e) NULL)
    if(is.null(ll)){
      if(debug) apollo_print('function() could not be inserted correctly in user-defined functions. No analytic gradients possible')
      apollo_inputs$apollo_control$analyticGrad <- FALSE
    } else {
      # Compare result of modified LL against original LL
      tmp4 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs), error=function(e) NULL)
      test <- !is.null(tmp4) && is.numeric(tmp4) && is.numeric(ll) && length(tmp4)==length(ll)
      test <- test && abs(sum(tmp4, na.rm=TRUE) - sum(ll, na.rm=TRUE))/sum(tmp4, na.rm=TRUE) < 0.01
      if(test){
        apollo_probabilities <- tmp1
        if(exists('tmp2', inherits=FALSE)) apollo_inputs$apollo_randCoeff <- tmp2
        if(exists('tmp3', inherits=FALSE)) apollo_inputs$apollo_lcPars    <- tmp3
      }
    }
    rm(ll, test); for(i in paste0('tmp',1:4)) if(exists(i, inherits=FALSE)) rm(i)
  }
  
  ### Copy the modified apollo_probabilities inside apollo_inputs
  apollo_inputs$apollo_probabilities <- apollo_probabilities
  
  # # # # # # # # # # # # #
  #### Create workers  ####
  # # # # # # # # # # # # #
  
  if(apollo_inputs$apollo_control$nCores==1) cl <- NA else {
    if(debug) apollo_print('Creating cluster...')
    # Creates cluster and also deletes database and draws from apollo_inputs in here and in .GlobalEnv
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent, cleanMemory=cleanMemory)
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
  analyticGrad <- apollo_inputs$apollo_control$analyticGrad
  
  ### Iterations count
  if(analyticGrad) nIter <- 1 else nIter <- -1
  preLL <- -Inf
  tmp <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
  lastB <- setNames(rnorm(length(tmp)), tmp)
  rm(tmp)

  ### Create apollo_logLike
  if(singleCore){
    ### SINGLE CORE
    
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      # Return iteration count, if requested
      if(getNIter) return(nIter)
      # Calculate LL
      b  <- c(bVar, bFixedVal)
      ll <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      if(!workInLogs) ll <- log(ll) # condition used to be workInLogs & panelData
      sumll <- sum(ll)
      # Keep count of iterations and write them to file, if requested
      if(countIter || writeIter){
        if(analyticGrad){
          # If using analytic gradient, assume a new iteration if the LL improved
          newIter <- is.finite(sumll) && (sumll > preLL)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter) apollo_writeTheta(b, sumll, modelName)
          if(newIter) preLL <<- sumll
        } else {
          # If using numeric gradient, assume new iteration if the parameters haven't changed since last call
          newIter <- all(bVar==lastB)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter) apollo_writeTheta(b, sumll, modelName)
          lastB <<- bVar
        }
      }
      # Return LL
      if(sumLL) return(sumll) else return(ll)
    }
  } else {
    ### MULTI CORE
    
    # Clean up parent environment and check that cluster exists
    rm(apollo_inputs, apollo_probabilities, apollo_beta, apollo_fixed, apollo_estSet, panelData)
    if(length(cl)==1 && is.na(cl)) stop("Cluster is missing on 'apollo_makeLogLike' despite nCores>1.")

    # Function to be run on each thread
    apollo_parProb <- function(b){
      P <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      return(P)
    }
    environment(apollo_parProb) <- new.env(parent=environment())
    
    # Loglike
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      # Return iteration count, if requested
      if(getNIter) return(nIter)
      # Calculate LL
      b  <- c(bVar, bFixedVal)
      ll <- parallel::clusterCall(cl=cl, fun=apollo_parProb, b=b)
      ll <- unlist(ll)
      if(!workInLogs) ll <- log(ll) # condition used to be workInLogs & panelData
      sumll <- sum(ll)
      # Keep count of iterations and write them to file, if requested
      if(countIter || writeIter){
        if(analyticGrad){
          # If using analytic gradient, assume a new iteration if the LL improved
          newIter <- is.finite(sumll) && (sumll > preLL)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter) apollo_writeTheta(b, sumll, modelName)
          if(newIter) preLL <<- sumll
        } else {
          # If using numeric gradient, assume new iteration if the parameters haven't changed since last call
          newIter <- all(bVar==lastB)
          if(countIter && newIter) nIter <<- nIter + 1
          if(writeIter && newIter) apollo_writeTheta(b, sumll, modelName)
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
    grad  <- c(); mType <- c()
    for(i in name){
      grad  <- c(grad , apollo_inputs[[i]]$gradient)
      mType <- c(mType, apollo_inputs[[i]]$modelType)
    }
    name <- substr(name, 1, nchar(name)-nchar("_settings"))
    return( list(name=name, grad=grad, mType=mType) )
  }
  if(!anyNA(cl)) tmp <- parallel::clusterCall(cl, function(f) {environment(f) <- globalenv(); f()}, f)[[1]] else tmp <- f()
  if(length(tmp$name)>0 && length(tmp$name)==length(tmp$mType)) mType <- setNames(tmp$mType, tmp$name) else mType <- NULL
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
    if(!all(sapply(nObsTot, length)==length(nObsTot[[1]]))) stop('Model components are inconsistent across workers. Try seting apollo_control$nCores=1')
    nObsTot <- Reduce('+', nObsTot)
  } else nObsTot <- f()
  rm(f)
  if(debug) cat(' Done.\n')
  
  return(apollo_logLike)
}
