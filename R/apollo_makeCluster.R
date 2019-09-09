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
#' @return Cluster (i.e. an object of class cluster from package parallel)
apollo_makeCluster <- function(apollo_probabilities, apollo_inputs, silent=FALSE){
  
  ### Split data and draws
  LL     <- apollo_splitDataDraws(apollo_inputs, silent)
  nCores <- length(LL)
  
  ### Create cluster
  if(!silent) cat('\nPreparing workers')
  if(!silent) cat('\n Loading libraries...')
  cl <- parallel::makeCluster(nCores)
  
  ### Load libraries in the cluster (same as in workspace)
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
  
  
  ### Report memory usage
  if(!silent){
    mbRAM <- sum(gc()[,2])
    if(nCores>1){
      gcClusters <- parallel::clusterCall(cl, gc)
      gcClusters <- Reduce('+',lapply(gcClusters, function(x) sum(x[,2])))
      mbRAM <- mbRAM + gcClusters
    }
    cat(' Done. ',mbRAM,'MB of RAM in use.',sep='')
  }
  
  ### Copy apollo_inputs (with corresponding database and draws piece) to each workers
  if(!silent) cat("\n Copying data...")
  #parallel::parLapply(cl, LL, fun=function(ll) assign("apollo_inputs", ll, envir=globalenv()) )
  parallel::parLapply(cl, LL, fun=function(ll){
    tmp <- globalenv()
    assign("apollo_inputs", ll, envir=tmp)
  }  )
  
  ### Report memory usage
  if(!silent){
    mbRAM1 <- sum(gc()[,2])
    mbRAM2 <- 0
    if(nCores>1){
      mbRAM2 <- parallel::clusterCall(cl, gc)
      mbRAM2 <- Reduce('+',lapply(mbRAM2, function(x) sum(x[,2])))
    }
    mbRAMmax <- mbRAM1+mbRAM2
    rm(LL)
    mbRAMcurrent <- sum(gc()[,2]) + mbRAM2
    cat(' Done. ', mbRAMcurrent, 'MB of RAM in use (max was ',mbRAMmax,'MB)', sep='')
  }
  
  # The following two lines do the same thing, i.e. enumerate 
  # elements in the workspace of each worker.
  #parallel::clusterEvalQ(cl, ls())
  #parallel::clusterCall(cl, ls, envir=.GlobalEnv)
  return(cl)
}