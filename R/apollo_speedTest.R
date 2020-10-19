#' Measures evaluation time of a model
#'
#' Measures the evaluation time of a model for different number of cores and draws.
#'
#' This function evaluates the function \code{apollo_probabilities} several times using different number of threads (a.k.a. processor cores),
#' and draws (if the model uses mixing). Then it plots the estimation time for each combination.
#' Estimation time grows at least linearly with number of draws, while time savings are decreasing with the number of threads.
#' This function can help decide what number of draws and cores to use for estimation, though a high number of draws is always
#' recommended. If the computer will be used for additional activities during estimation, no more than (machine number of cores - 1) should be used.
#' Using more threads than cores available in the machine will lead to reduce performance.
#' The use of additional cores come at the expense of additional memory usage. If R uses more memory than the physical RAM available,
#' then significant slow-downs in processing time can be expected. This function can help avoiding such pitfalls.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \code{apollo_beta}: Named numeric vector. Names and values of model parameters.
#'                            \item \code{apollo_inputs}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \code{functionality}: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param speedTest_settings List containing options for the speed test. The following are valid options.
#'                                   \itemize{
#'                                     \item \code{nDrawsTry}: Numeric vector. Number of inter and intra-person draws to try. Default value is c(50, 100, 200).
#'                                     \item \code{nCoresTry}: Numeric vector. Number of threads to try. Default is from 1 to the detected number of cores.
#'                                     \item \code{nRep}: Numeric scalar. Number of times the likelihood is evaluated for each combination of threads and draws. Default is 10.
#'                                   }
#' @return A matrix with the average time per evaluation for each number of threads and draws combination. A graph is also plotted.
#' @importFrom graphics matplot title legend
#' @export
apollo_speedTest=function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, speedTest_settings=NA){
  
  # # # # # # # # # # #
  #### Initialise  ####
  # # # # # # # # # # #
  
  tmp <- parallel::detectCores(); if(is.na(tmp)) tmp <- 3
  default <- list(nDrawsTry = c(50, 100, 200),
                  nCoresTry = 1:tmp,
                  nRep      = 10)
  if(length(speedTest_settings)==1 && is.na(speedTest_settings)) speedTest_settings <- default
  tmp <- names(default)[ !(names(default) %in% names(speedTest_settings)) ]
  for(i in tmp)  speedTest_settings[[i]] <- default[[i]]
  rm(tmp, default)
  
  apollo_inputs$silent <- TRUE
  
  nDrawsTry        = speedTest_settings[["nDrawsTry"]]
  nCoresTry        = speedTest_settings[["nCoresTry"]]
  nRep             = speedTest_settings[["nRep"]]
  apollo_checkArguments(apollo_probabilities,apollo_inputs$apollo_randCoeff,apollo_inputs$apollo_lcPars)
  
  ### Checks
  if(apollo_inputs$apollo_control$HB) stop("The function 'apollo_speedTest' is not applicable for models estimated using HB.")
  if(apollo_inputs$apollo_control$mixing){
    if(anyNA(apollo_inputs$apollo_draws)) stop("Argument 'apollo_draws' must be provided when estimating mixture models.")
    if(!is.function(apollo_inputs$apollo_randCoeff)) stop("Argument 'apollo_randCoeff' must be provided when estimating mixture models.")
  }
  if(!apollo_inputs$apollo_control$mixing){
    apollo_inputs$apollo_draws <- list(interUnifDraws=c(), interNormDraws=c(), interNDraws=0,
                         intraUnifDraws=c(), intraNormDraws=c(), intraNDraws=0)
    nDrawsTry <- 0
    apollo_inputs$apollo_randCoeff <- NA
  }
  
  usesInter <- (length(apollo_inputs$apollo_draws$interUnifDraws) + length(apollo_inputs$apollo_draws$interNormDraws))>0 &
    apollo_inputs$apollo_draws$interNDraws>0
  usesIntra <- (length(apollo_inputs$apollo_draws$intraUnifDraws) + length(apollo_inputs$apollo_draws$intraNormDraws))>0 &
    apollo_inputs$apollo_draws$intraNDraws>0
  if(!usesInter) apollo_inputs$apollo_draws$interNDraws <- 0
  if(!usesIntra) apollo_inputs$apollo_draws$intraNDraws <- 0
  
  # # # # # # # # # #
  #### Main Loop ####
  # # # # # # # # # #
  
  cat("       ___Draws___              sec/")
  cat("\nnCores inter intra  progress  LLCall RAM(MB)")
  
  if(nRep<=10) stepPrint <- 1 else stepPrint <- floor(nRep/10)
  timeTaken <- matrix(0, nrow=length(nCoresTry), ncol=length(nDrawsTry))
  for(j in 1:length(nDrawsTry)){
    
    for(i in 1:length(nCoresTry)){
      
      apollo_inputs$apollo_control$nCores <- nCoresTry[i]
      if(usesInter) apollo_inputs$apollo_draws$interNDraws <- nDrawsTry[j]
      if(usesIntra) apollo_inputs$apollo_draws$intraNDraws <- nDrawsTry[j]
      
      cat("\n")
      cat(sprintf("%6i",apollo_inputs$apollo_control$nCores), " ", 
          sprintf("%5i",apollo_inputs$apollo_draws$interNDraws), " ", 
          sprintf("%5i",apollo_inputs$apollo_draws$intraNDraws), " ", sep="")
      
      if(apollo_inputs$apollo_control$mixing){
        apollo_inputs$draws <- apollo_makeDraws(apollo_inputs, silent=TRUE)
      } else draws <- NA
      
      #if(apollo_inputs$apollo_control$nCores>1){
      #  cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=TRUE)
      #} else cl <- NA
      
      theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
      theta_fix_val <- apollo_beta[apollo_fixed]
      theta_min <- 0.95*theta_var_val
      theta_max <- 1.05*theta_var_val
      U <- t(apollo_mlhs(nRep, length(theta_var_val), 1))
      deltaTheta <- theta_max - theta_min
      
      apollo_estSet  <- list(estimationRoutine="bfgs", maxIterations=2, writeIter=FALSE, hessianRoutine="none", printLevel=0, silent=FALSE)
      apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, apollo_estSet)
      
      t0 <- Sys.time()
      
      for(k in 1:nRep){
        theta_test <- theta_min + as.vector(U[,k])*deltaTheta
        LL_test    <- apollo_logLike(theta_test, sumLL=TRUE)
        if(k%%stepPrint==0) cat('.')
      }
      
      t1 <- Sys.time()
      
      
      timeTaken[i,j] <- as.numeric(difftime(t1,t0,units='secs'))/nRep
      cat(' ', sprintf("%6.2f", timeTaken[i,j],2),sep='')
      
      mbRAM <- sum(gc()[,2])
      if(apollo_inputs$apollo_control$nCores>1){
        gcClusters <- parallel::clusterCall(environment(apollo_logLike)$cl, gc)
        gcClusters <- Reduce('+',lapply(gcClusters, function(x) sum(x[,2])))
        mbRAM <- mbRAM + gcClusters
      }
      cat(' ', sprintf("%7.1f", mbRAM), sep='')
      
      if(exists('apollo_logLike') && !anyNA(environment(apollo_logLike)$cl)) parallel::stopCluster(environment(apollo_logLike)$cl)
    }
  }
  cat("\n")
  
  colnames(timeTaken) <- paste('draws', nDrawsTry, sep='')
  rownames(timeTaken) <- paste('cores', nCoresTry, sep='')
  
  graphics::matplot(x=nCoresTry, y=timeTaken,
          type='o', lty=1, lwd=2, cex=1.5, fg='gray',
          lab=c(length(nCoresTry), length(nDrawsTry), 1),
          pch=1:length(nDrawsTry), col=1:length(nDrawsTry),
          xlab='Number of cores', ylab='Avg. time per call to apollo_probabilities (sec)') 
  graphics::title(apollo_inputs$apollo_control$modelName)
  if(apollo_inputs$apollo_control$mixing) graphics::legend("topright", legend = colnames(timeTaken),
                                   pch=1:length(nDrawsTry), col=1:length(nDrawsTry), lty=1,
                                   bg='transparent', bty='n', horiz=TRUE) 
  
  cat('\nSummary of results (sec. per call to LL function)','\n', sep="")
  print(round(timeTaken,4))
  invisible(timeTaken)
}

