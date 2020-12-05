#' Creates cluster for estimation.
#'
#' Splits data, creates cluster and loads different pieces of the database on each worker.
#' 
#' Internal use only. Called by \code{apollo_estimate} before estimation. Using multiple cores greatly increases memory consumption.
#' 
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If TRUE, no messages are printed to the terminal. FALSE by default. It overrides \code{apollo_inputs$silent}.
#' @param cleanMemory Boolean. If TRUE, it saves apollo_inputs to disc, and removes database and draws from 
#'                    the apollo_inputs in .GlobalEnv and the parent environment.
#' @return Cluster (i.e. an object of class cluster from package parallel)
#' @export
apollo_makeCluster <- function(apollo_probabilities, apollo_inputs, silent=FALSE, cleanMemory=FALSE){
  
  #### Split data ####
  #apollo_control <- apollo_inputs[["apollo_control"]]
  #database       <- apollo_inputs[["database"]]
  ### change 27 July
  ###silent         <- apollo_inputs$silent
  if(silent==FALSE) silent         <- apollo_inputs$silent
  ### end change
  debug          <- apollo_inputs$apollo_control$debug
  
  ### Extract useful elements
  nObs    <- nrow(apollo_inputs$database)
  indivID <- apollo_inputs$database[,apollo_inputs$apollo_control$indivID]
  namesID <- unique(indivID)
  nIndiv  <- length(namesID)
  #nObsID  <- as.vector(table(indivID))
  nObsID  <- rep(0, nIndiv)
  for(n in 1:nIndiv) nObsID[n] <- sum(indivID==namesID[n])
  mixing  <- apollo_inputs$apollo_control$mixing
  nCores  <- apollo_inputs$apollo_control$nCores
  
  ### Figure out how to split the database per individual
  #database <- database[order(indivID),]
  #indivID  <- database[,apollo_control$indivID] # update order
  #namesID  <- unique(indivID)                   # update order
  #apollo_inputs$database <- database            # update order
  #rm(database); gc()
  if(debug) apollo_print(paste0('Attempting to split data into ', nCores, ' pieces.'))
  obj          <- ceiling(nObs/nCores)
  counter      <- 0
  currentCore  <- 1
  assignedCore <- rep(0, nObs)
  assignedCoreIndiv <- rep(0, nIndiv)
  i <- 1
  for(n in 1:nIndiv){
    assignedCore[i:(i+nObsID[n]-1)] <- currentCore
    assignedCoreIndiv[n] <- currentCore
    i <- i + nObsID[n]
    counter <- counter + nObsID[n]
    if(counter>=obj & currentCore<nCores){
      currentCore <- currentCore + 1
      counter <- 0
    }
  }
  nCores <- max(assignedCore)
  if(debug){
    if(nCores!=apollo_inputs$apollo_control$nCores) apollo_print(paste0(nCores, ' workers (threads) will be used for estimation.'))
    coreLoad <- setNames(rep(0, nCores), paste0('worker', 1:nCores))
    for(i in 1:nCores) coreLoad[i] <- sum(assignedCore==i)
    apollo_print(paste0('Obs. per worker (thread): ',paste0(coreLoad, collapse=", ")))
  }
  apollo_inputs$apollo_control$nCores <- nCores
  rm(obj, counter, currentCore, i, n)
  
  
  ### Create names of temp files
  inputPieceFile <- rep("", nCores)
  for(i in 1:nCores) inputPieceFile[i] <- tempfile()
  
  ### Identify elements to split
  # For lists, it looks at its first element
  asLst <- c() # one element per individual
  byObs <- c() # one row per observation
  byInd <- c() # one row per individual
  asIs  <- c() # none of the above
  for(e in ls(apollo_inputs)){ # e is name of element, E is the element
    E <- apollo_inputs[[e]]
    if(is.list(E) && length(E)==nIndiv) asLst <- c(asLst, e) else{
      if(is.list(E)) E <- E[[1]]
      if(is.array(E)) nRows <- dim(E)[1] else nRows <- length(E)
      if(nRows==nObs) byObs <- c(byObs, e)
      if(nRows==nIndiv) byInd <- c(byInd, e)
      if(nRows!=nObs & nRows!=nIndiv) asIs <- c(asIs, e)
    }
  }; rm(e, E, nRows)
  
  ### Create and write to disk each fragment of apollo_inputs
  if(debug) cat("Writing pieces to disk") # do not turn to apollo_print
  for(i in 1:nCores){
    L <- list()
    # Elements to be kept unchanged
    L <- apollo_inputs[ asIs ]
    # Elements to be split by observations
    if(length(byObs)>0) for(e in byObs){
      E <- apollo_inputs[[e]]
      r <- assignedCore==i
      if(is.list(E) & !is.data.frame(E)){
        L[[e]] <- lapply(E, apollo_keepRows, r=r)
      } else L[[e]] <- apollo_keepRows(E, r)
    }
    # Elements to be split by individual
    if(length(byInd)>0) for(e in byInd){
      E <- apollo_inputs[[e]]
      r <- assignedCoreIndiv==i
      if(is.list(E) & !is.data.frame(E)){
        L[[e]] <- lapply(E, apollo_keepRows, r=r)
      } else L[[e]] <- apollo_keepRows(E, r)
    }
    # Elements to be split as list
    if(length(asLst)>0) for(e in asLst){
      E <- apollo_inputs[[e]]
      r <- which(assignedCoreIndiv==i)
      L[[e]] <- E[r]
    }
    # Write to disk
    wroteOK <- tryCatch({
      saveRDS(L, file=inputPieceFile[i])
      TRUE
    }, warning=function(w) FALSE, error=function(e) FALSE)
    if(!wroteOK) stop("Apollo could not write data pieces to disk. This is necessary for multi-core processing to work.")
    if(debug) cat(".")
    if(debug && i==nCores) cat("\n")
  }
  
  if(debug) apollo_print(paste0("Writing completed. ", sum(gc()[,2]), 'MB of RAM in use.'))
  
  
  
  #### Create cluster and load data ####
  
  ### Create cluster
  if(!silent) apollo_print('Preparing workers for multithreading...')
  if(Sys.info()["sysname"] == "Darwin"){
    cl <- parallel::makeCluster(nCores, setup_strategy="sequential")
  } else cl <- parallel::makeCluster(nCores)
  
  ### Delete draws and database from memory in main thread
  if(cleanMemory){
    # No backup of apollo_inputs is created, it is the caller responsibility to do that (e.g. apollo_estimate)
    if(debug) cat('Cleaning memory in main thread...') # do not change to apollo_print
    # Go over calling stack and delete apollo_inputs$database and apollo_inputs$draws wherever found
    for(e in sys.frames()){
      x <- tryCatch(get('apollo_inputs', envir=e, inherits=FALSE), error=function(e) NULL)
      if(!is.null(x)){
        if(!is.null(x$database)) x$database <- NULL
        if(!is.null(x$draws   )) x$draws    <- NULL
        assign('apollo_inputs', x, envir=e)
      }
    }
    tmp <- globalenv()
    x <- tryCatch(get('apollo_inputs', envir=tmp, inherits=FALSE), error=function(e) NULL)
    if(!is.null(x)){
      if(!is.null(x$database)) x$database <- NULL
      if(!is.null(x$draws   )) x$draws    <- NULL
      assign('apollo_inputs', x, envir=tmp)
    }
    gc()
    if(debug) cat(' Done. ',sum(gc()[,2]),'MB of RAM in use.\n',sep='') # do not change to apollo_print
  }
  
  ### Load libraries in the cluster (same as in workspace)
  if(debug) cat('Loading libraries...')  # do not change to apollo_print
  excludePackages<- c('parallel')
  loadedPackages <- search()
  loadedPackages <- loadedPackages[grepl("^(package:)", loadedPackages)]
  loadedPackages <- substr(loadedPackages, start=9, stop=100)
  loadedPackages <- loadedPackages[!(loadedPackages %in% excludePackages)]
  #if(!("apollo" %in% loadedPackages)) loadedPackages <- c(loadedPackages, "apollo")
  if(length(loadedPackages)>0){
    parallel::clusterCall(cl=cl, function(lib, path) {
      .libPaths(path)
      for(i in 1:length(lib)) library(lib[i],character.only=TRUE)
    }, lib=loadedPackages, path=.libPaths())
  }
  ### Copy functions in the workspace to the workers
  # apollo_probabilities is copied separately to ensure it keeps its name
  funcs <- as.vector(utils::lsf.str(envir=.GlobalEnv))
  parallel::clusterExport(cl, funcs, envir=.GlobalEnv)
  parallel::clusterExport(cl, "apollo_probabilities", envir=environment())
  if(debug){
    mbRAM <- sum(gc()[,2])
    if(nCores>1){
      gcClusters <- parallel::clusterCall(cl, gc)
      gcClusters <- Reduce('+',lapply(gcClusters, function(x) sum(x[,2])))
      mbRAM <- mbRAM + gcClusters
    }
    cat(' Done. ',mbRAM,'MB of RAM in use.',sep='') # do not change to apollo_print
  }
  
  ### Load apollo_input pieces in each worker
  if(debug) cat("\nLoading data...") # do not change to apollo_print
  parallel::parLapply(cl, inputPieceFile, fun=function(fileName){
    tmp <- globalenv()
    txt <- " This is necessary for multi-core processing to work."
    if(!file.exists(fileName)) stop(paste0("A piece of data is missing from disk.",txt))
    apollo_inputs_piece <- tryCatch(readRDS(fileName),
                                    warning = function(w) FALSE,
                                    error   = function(e) FALSE)
    if(is.logical(apollo_inputs_piece) && !apollo_inputs_piece) stop(paste0("A piece of data could not be loaded from disk.",txt))
    assign("apollo_inputs", apollo_inputs_piece, envir=tmp)
  })
  unlink(inputPieceFile) # delete tmp files
  
  ### Report memory usage
  if(debug){
    mbRAM1 <- sum(gc()[,2])
    mbRAM2 <- 0
    if(nCores>1){
      mbRAM2 <- parallel::clusterCall(cl, gc)
      mbRAM2 <- Reduce('+',lapply(mbRAM2, function(x) sum(x[,2])))
    }
    cat(' Done. ', mbRAM1+mbRAM2, 'MB of RAM in use\n', sep='') # do not change to apollo_print
  }
  
  # The following two lines enumerate 
  # elements in the workspace of each worker.
  #print(parallel::clusterEvalQ(cl, ls(apollo_inputs$draws)))
  #parallel::clusterCall(cl, ls, envir=.GlobalEnv)
  return(cl)
}