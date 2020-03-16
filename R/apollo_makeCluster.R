#' Creates cluster for estimation.
#'
#' Creates cluster and loads pieces of the database for each worker.
#' 
#' Internal use only. Called by \code{apollo_estimate} before estimation. AT least doubles up memory usage. But during the splitting it uses even more (~250% of the original).
#' 
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If TRUE, it reports progress to the console. Default is FALSE.
#' @param cleanMemory Boolean. If TRUE, it saves apollo_inputs to disc, and removes database and draws from 
#'                    the apollo_inputs in .GlobalEnv and the parent environment.
#' @return Cluster (i.e. an object of class cluster from package parallel)
apollo_makeCluster <- function(apollo_probabilities, apollo_inputs, silent=FALSE, cleanMemory=FALSE){
  
  ### Split data and draws
  inputPieceFile <- apollo_splitDataDraws(apollo_inputs, silent)
  nCores         <- length(inputPieceFile)
  
  ### Create cluster
  if(!silent) cat('\nPreparing workers')
  cl <- parallel::makeCluster(nCores)
  
  ### Delete draws and database from memory in main thread
  if(cleanMemory){
    if(!silent) cat('\n Cleaning memory in main thread...')
    saveRDS(apollo_inputs, file=paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs"))
    apollo_inputs$draws <- NULL
    apollo_inputs$database <- NULL
    assign("apollo_inputs", apollo_inputs, envir=parent.frame())
    tmp <- .GlobalEnv
    if(exists("apollo_inputs", envir=tmp)) assign("apollo_inputs", apollo_inputs, envir=tmp)
    if(!silent) cat(' Done. ',sum(gc()[,2]),'MB of RAM in use.',sep='')
  }
  
  ### Load libraries in the cluster (same as in workspace)
  if(!silent) cat('\n Loading libraries...')
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
  if(!silent){
    mbRAM <- sum(gc()[,2])
    if(nCores>1){
      gcClusters <- parallel::clusterCall(cl, gc)
      gcClusters <- Reduce('+',lapply(gcClusters, function(x) sum(x[,2])))
      mbRAM <- mbRAM + gcClusters
    }
    cat(' Done. ',mbRAM,'MB of RAM in use.',sep='')
  }
  
  ### Load apollo_input pieces in each worker
  if(!silent) cat("\n Loading data...")
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
  if(!silent){
    mbRAM1 <- sum(gc()[,2])
    mbRAM2 <- 0
    if(nCores>1){
      mbRAM2 <- parallel::clusterCall(cl, gc)
      mbRAM2 <- Reduce('+',lapply(mbRAM2, function(x) sum(x[,2])))
    }
    cat(' Done. ', mbRAM1+mbRAM2, 'MB of RAM in use\n', sep='')
  }
  
  # The following two lines enumerate 
  # elements in the workspace of each worker.
  #print(parallel::clusterEvalQ(cl, ls(apollo_inputs$draws)))
  #parallel::clusterCall(cl, ls, envir=.GlobalEnv)
  return(cl)
}