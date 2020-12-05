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
#' @param cleanMemory Logical. If TRUE, then \code{apollo_inputs$draws} and \code{apollo_inputs$database} are erased
#'                    throughout the calling stack. Used to reduce memory usage in case of multithreading and a large
#'                    database or number o draws.
#' @return apollo_logLike function.
#' @export
apollo_makeLogLike <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                               apollo_estSet, cleanMemory=FALSE){
  
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
  
  ### Checks for scaling and insert them inside apollo_probabilities
  scaling <- setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                      names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  test <- is.null(apollo_inputs$apollo_scaling)
  test <- test || (length(apollo_inputs$apollo_scaling)==1 && is.na(apollo_inputs$apollo_scaling))
  if(!test){
    scaling[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling
    if(any(!(names(scaling) %in% names(apollo_beta)))) stop("Some parameters included in 'scaling' are not included in 'apollo_beta'")
    if(any(names(scaling) %in% apollo_fixed)) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
    if(any(scaling<0)){
      scaling <- abs(scaling)
      txt <- 'Some negative values in "scaling" were replaced by their absolute value'
      if(!silent) apollo_print(paste0('WARNING: ', txt, '.')) else warning(txt)
      rm(txt)
    }
    if(any(scaling<=0)) stop('All terms in "scaling" should be strictly positive!')
  }
  apollo_inputs$apollo_scaling <- scaling
  rm(scaling)
  
  ### Insert scaling if needed
  apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
  if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
    apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
  }
  if(is.function(apollo_inputs$apollo_lcPars)){
    apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
  }
  if(exists('txt', inherits=FALSE)) rm(txt)
  
  ### Insert componentName if missing
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  
  ### Insert "function ()" in apollo_probabilities and apollo_randCoeff (if needed)
  if(apollo_inputs$apollo_control$analyticGrad){
    tmp1 <- apollo_insertFunc(apollo_probabilities, like=TRUE)
    test <- !is.null(apollo_inputs$apollo_randCoeff) && is.function(apollo_inputs$apollo_randCoeff)
    if(test) tmp2 <- apollo_insertFunc(apollo_inputs$apollo_randCoeff, randCoeff=TRUE)
    # Test that modified LL actually works
    ll <- tryCatch(tmp1(apollo_beta, apollo_inputs), error=function(e) NULL)
    if(is.null(ll)){
      if(debug) apollo_print('Modified apollo_probabilities did not work. Defaulting to original one (no analytic gradients possible)')
      apollo_inputs$apollo_control$analyticGrad <- FALSE
    } else {
      apollo_probabilities <- tmp1
      if(exists('tmp2', inherits=FALSE)) apollo_inputs$apollo_randCoeff <- tmp2
    }
    rm(ll, tmp1, test); if(exists('tmp2', inherits=FALSE)) rm(tmp2)
  }
  
  ### Copy the modified apollo_probabilities inside apollo_inputs
  apollo_inputs$apollo_probabilities <- apollo_probabilities
  
  # # # # # # # # # # # # #
  #### Create workers  ####
  # # # # # # # # # # # # #
  
  if(apollo_inputs$apollo_control$nCores==1) cl <- NA else {
    # Creates cluster and also deletes database and draws from apollo_inputs in here and in .GlobalEnv
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent, cleanMemory=cleanMemory)
    apollo_inputs$apollo_control$nCores <- length(cl)
  }
  
  # # # # # # # # # # # # # # # #
  #### Create apollo_logLike ####
  # # # # # # # # # # # # # # # #
  
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

  ### SINGLE CORE
  if(singleCore){
    apollo_logLike <- function(bVar, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      # Return iteration count, if requested
      if(getNIter) return(nIter)
      # Calculate LL
      b  <- c(bVar, bFixedVal)
      ll <- apollo_probabilities(b, apollo_inputs, functionality="estimate")
      if(!workInLogs) ll <- log(ll) # condition used to be workInLogs & panelData
      sumll <- sum(ll)
      # Keep count of iterations and write them to file, if requested
      if(analyticGrad){
        newIter <- is.finite(sumll) && (sumll > preLL)
        if(countIter && newIter) nIter <<- nIter + 1
        if(writeIter && newIter) apollo_writeTheta(bVar, sumll, modelName)
        if(newIter) preLL <<- sumll
      } else if(estAlg=="bfgs" && exists('lastFuncParam', envir=globalenv())) {
        lastFuncParam <- get("lastFuncParam", envir=globalenv())
        newIter       <- !anyNA(lastFuncParam) && length(lastFuncParam)==length(bVar) && all(lastFuncParam==bVar)
        if(countIter && newIter) nIter <<- nIter + 1
        if(writeIter && newIter) apollo_writeTheta(bVar, sumll, modelName)
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
      if(analyticGrad){
        newIter <- is.finite(sumll) && (sumll > preLL)
        if(countIter && newIter) nIter <<- nIter + 1
        if(writeIter && newIter) apollo_writeTheta(bVar, sumll, modelName)
        if(newIter) preLL <<- sumll
      } else if(estAlg=="bfgs" && exists('lastFuncParam', envir=globalenv())) {
        lastFuncParam <- get("lastFuncParam", envir=globalenv())
        newIter       <- !anyNA(lastFuncParam) && length(lastFuncParam)==length(bVar) && all(lastFuncParam==bVar)
        if(countIter && newIter) nIter <<- nIter + 1
        if(writeIter && newIter) apollo_writeTheta(bVar, sumll, modelName)
      }
      # Return LL
      if(sumLL) return(sumll) else return(ll)
    }

  }
  
  # # # #  # # # # # # # # # # # # # # # # #
  #### Pre-process apollo_probabilities ####
  # # # #  # # # # # # # # # # # # # # # # #
  
  ### Preprocess model (for optimisation)
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
  f <- function(){
    name <- grep("_settings$", ls(apollo_inputs), value=TRUE)
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
  rm(f, tmp)
  #if(debug){
  #  if(length(tmp$name)==0) apollo_print("No pre-processing performed") else {
  #    # Print report
  #    cat("ComponentName  Type   Gradient Optimisation\n") # cat, as apollo_print messes up the alignment
  #    for(i in 1:length(tmp$grad)){
  #      txt <- 14 - nchar(tmp$name[i])
  #      txt <- ifelse(txt>=0, paste0(tmp$name[i], paste0(rep(" ", txt), collapse="")), substr(tmp$name[i], 1, 14))
  #      txt <- paste0(txt, " ", substr(paste0(tmp$mType[i], '      '), 1, 6))
  #      txt <- paste0(txt, " ", ifelse(tmp$grad[i], "analytic\n", "numeric \n"))
  #      cat(txt) # cat, as apollo_print messes up the alignment
  #    }; rm(i, txt)
  #    apollo_print(paste0("Whole model gradient function creation ", ifelse(is.null(grad), "failed.", "succeeded.")))
  #  }
  #}; rm(tmp)
  
  ### Count number of observations
  f <- function(){ # this function returns a vector c(component1Name=nObs1, component2Name=nObs2, ...) ignoring LC components
    name    <- grep("_settings$", ls(apollo_inputs), value=TRUE)
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
  
  return(apollo_logLike)
}
