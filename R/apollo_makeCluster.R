#' Creates cluster for estimation.
#'
#' Creates cluster and loads pieces of the database in each worker.
#'
#' Internal use only. Called by apollo_estimate before estimation. AT least doubles up memory usage. But during the splitting it uses even more (~250% of the original).
#'
#' @param apollo_control List of model options.
#' @param apollo_probabilities Function. likelihood provided by the user.
#' @param database data.frame Data to be used by model. Gets splitted and copied into the cluster.
#' @param draws List of draws (numeric arrays). Default is NA. Get splitted and copied into the cluster.
#' @param apollo_randcoeff Function. Creates the random coeff or error terms by mixing draws and params. Provided by the user.
#' @param apollo_lcpars Function. Used for latent class models. Not implemented yet.
#' @param silent Boolean. Default is FALSE.
#' @return Cluster (i.e. an object of class cluster from package parallel)
apollo_makeCluster <- function(apollo_control, apollo_probabilities, database,
                            draws=NA, apollo_randcoeff=NA, apollo_lcpars=NA, silent=FALSE){
  LL <- apollo_splitDataDraws(apollo_control, database, draws, silent)
  apollo_control$nCores <- length(LL)

  if(!silent) cat('Creating workers and loading libraries in them.\n')
  cl <- parallel::makeCluster(apollo_control$nCores)

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

  if(!("cmcRcode" %in% loadedPackages)){
    apollo_stuff <- ls(.GlobalEnv)
    apollo_stuff <- grep('apollo_', apollo_stuff, value=TRUE)
    parallel::clusterExport(cl=cl, apollo_stuff)
    rm(apollo_stuff)
  }

  if(!silent){
    mbRAM <- sum(gc()[,2])
    if(apollo_control$nCores>1){
      gcClusters <- parallel::clusterCall(cl, gc)
      gcClusters <- Reduce('+',lapply(gcClusters, function(x) sum(x[,2])))
      mbRAM <- mbRAM + gcClusters
    }
    cat('Workers created and libraries loaded. (',mbRAM,'Mb of RAM in use)\n',sep='')
  }

  if(!silent) cat("Loading data on workers.\n")
  parallel::parLapply(cl, LL, fun=function(ll, prob, ctrl, rndCoeff, lcPars){
    env <- globalenv()
    assign('apollo_probabilities', prob       , envir=env)
    assign('database'  , ll$database, envir=env)
    #assign('apollo_ctrl', ctrl       , envir=env)
    assign('draws', ll$draws   , envir=env)
    #assign('apollo_rndCoeff',rndCoeff, envir=env)
    #assign('apollo_lcpars',    lcPars, envir=env)
  }, prob=apollo_probabilities, ctrl=apollo_control, rndCoeff=apollo_randcoeff, lcPars=apollo_lcpars )

  if(!silent){
    mbRAM1 <- sum(gc()[,2])
    mbRAM2 <- 0
    if(apollo_control$nCores>1){
      mbRAM2 <- parallel::clusterCall(cl, gc)
      mbRAM2 <- Reduce('+',lapply(mbRAM2, function(x) sum(x[,2])))
    }
    mbRAMmax <- mbRAM1+mbRAM2
    rm(LL)
    mbRAMcurrent <- sum(gc()[,2]) + mbRAM2
    cat('Data pieces loaded in workers.\n Max RAM usage was ',mbRAMmax,'Mb, current usage is ',
        mbRAMcurrent,'Mb\n\n', sep='')
  }

  return(cl)
}
