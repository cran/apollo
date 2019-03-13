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
  
  LL <- apollo_splitDataDraws(apollo_inputs, silent)
  nCores <- length(LL)
  for(i in 1:nCores) LL[[i]]$apollo_control$nCores <- nCores
  
  if(!silent) cat('Creating workers and loading libraries...')
  cl <- parallel::makeCluster(nCores)
  
  excludePackages<- c('parallel')
  loadedPackages <- search()
  loadedPackages <- loadedPackages[grepl("^(package:)", loadedPackages)]
  loadedPackages <- substr(loadedPackages, start=9, stop=100)
  loadedPackages <- loadedPackages[!(loadedPackages %in% excludePackages)]
  if(length(loadedPackages)>0){
    parallel::clusterCall(cl=cl, function(lib) {
      for(i in 1:length(lib)) library(lib[i],character.only=TRUE)
    }, lib=loadedPackages)
  }
  
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
    cat(' Done. ',mbRAM,'Mb of RAM in use.\n',sep='')
  }
  
  if(!silent) cat("Copying data to workers...")
  parallel::parLapply(cl, LL, fun=function(ll){
    tmp <- globalenv()
    assign("apollo_inputs", ll, envir=tmp)
  }  )
  
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
    cat(' Done. ', mbRAMcurrent, 'Mb of RAM in use (max was ',mbRAMmax,'Mb)\n', sep='')
  }
  
  return(cl)
}