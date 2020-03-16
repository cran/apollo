#' Bootstrap a model
#'
#' Samples individuals with replacement from the database, and estimates the model in each sample.
#'
#' This function implements a basic block bootstrap. It estimates the model parameters on \code{nRep} number of different samples.
#' Each new sample is constructed by sampling \strong{with replacement} from the original full sample. Each new sample has as many 
#' individuals as the original sample, though some of them may be repeated. Sampling is done at the \strong{individual} level, 
#' therefore if different individuals have different number of observations, each re-sample could have different number of observations.
#' 
#' If the sampling wants to be done at the individual level (not recommended on panel data), then the optional 
#' \code{bootstrap_settings$samples} argument should be provided.
#' 
#' For each sample, only the parameters and loglikelihood are estimated. Standard errors are not calculated (they may be in 
#' future versions). The composition of each re-sample is stored on a file, though it should be consistent across runs.
#' 
#' This function writes three different files to the working directory:
#' \itemize{
#'   \item \code{modelName_bootstrap_params.csv}: Records the estimated parameters, final loglikelihood, and number of observations on each re-sample
#'   \item \code{modelName_bootstrap_samples.csv}: Records the composition of each re-sample.
#'   \item \code{modelName_bootstrap_vcov.csv}: Variance-covariance matrix of the estimated parameters across re-samples.
#' }
#' The first two files are updated throughout the run of this function, while the last one is only written once the function finishes.
#' 
#' When run, this function will look for the first two files above in the working directory. If they are found, the function will
#' attempt to pick up re-sampling from where those files left off. This is useful in cases where the original bootstrapping was 
#' interrupted, or when additional re-sampling wants to be performed.
#' 
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param estimate_settings List. Options controlling the estimation process. See \link{apollo_estimate}.
#'                          \code{hessianRoutine='none'} by default.
#' @param bootstrap_settings List. Options defining the sampling procedure. The following are valid options.
#'                                   \itemize{
#'                                     \item \strong{nRep}: Numeric scalar. Number of times the model must be estimated with different samples. Default is 30.
#'                                     \item \strong{samples}: Numeric matrix or data.frame. Optional argument. Must have as many rows as 
#'                                                             observations in the \code{database}, and as many columns as number of repetitions
#'                                                             wanted. Each column represents a re-sample, and each element the number of times 
#'                                                             that observation must be included in the sample. If this argument is provided, 
#'                                                             then \code{nRep} is ignored. Note that this allows sampling at the observation 
#'                                                             rather than the individual level, which is not recommended for panel data.
#'                                     \item \strong{seed}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples.
#'                                                          Only used if \code{samples} is \code{NA}. Default is 24.
#'                                   }
#' @return Covariance matrix of the \code{nRep} sets of estimated parameters. Also written to file.
#' @export
#' @importFrom maxLik maxLik
#' @importFrom stats cov
apollo_bootstrap <- function(apollo_beta, apollo_fixed,
                               apollo_probabilities, apollo_inputs,
                               estimate_settings=list(estimationRoutine="bfgs",
                                                      maxIterations=200,
                                                      writeIter=FALSE,
                                                      hessianRoutine="none",
                                                      printLevel=2L,
                                                      silent=FALSE),
                               bootstrap_settings=list(nRep=30,
                                                       samples=NA,
                                                       seed=24)){
  
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, hessianRoutine="none", printLevel=2L, silent=FALSE)
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  default <- list(nRep=30, samples=NA, seed=24)
  tmp <- names(default)[!(names(default) %in% names(bootstrap_settings))] # options missing in estimate_settings
  for(i in tmp) bootstrap_settings[[i]] <- default[[i]]
  rm(tmp)
  
  ### Extract values from apollo_inputs
  database         <- apollo_inputs$database
  apollo_inputs$apollo_control$noDiagnostics = TRUE
  apollo_control   <- apollo_inputs$apollo_control
  apollo_draws     <- apollo_inputs$apollo_draws
  apollo_randCoeff <- apollo_inputs$apollo_randCoeff
  apollo_lcPars    <- apollo_inputs$apollo_lcPars
  workInLogs       <- apollo_inputs$apollo_control$workInLogs
  name             <- apollo_control$modelName
  id                     <- database[,apollo_control$indivID]
  if(!is.numeric(id)) id <- as.numeric(as.factor(id))
  database[,apollo_control$indivID] <- id
  rm(id)
  
  ### Extract values from estimate_settings and estimate_settings
  estimationRoutine <- estimate_settings$estimationRoutine
  maxIterations     <- estimate_settings$maxIterations
  nRep              <- bootstrap_settings$nRep
  samples           <- bootstrap_settings$samples
  seed              <- bootstrap_settings$seed
  silent            <- estimate_settings$silent
  
  ### Validate arguments
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  estimationRoutine <- tolower(estimationRoutine)
  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh', or 'nr'.")
  if( (length(apollo_fixed) > 0) & any(!(apollo_fixed %in% names(apollo_beta))) ) stop("Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'")
  maxIterations <- round(maxIterations,0)
  if(maxIterations < 1) stop("Need at least one iteration!")
  if(workInLogs != TRUE) workInLogs=FALSE
  if(is.numeric(nRep) && nRep<0) stop("nRep must be a positive integer.")
  if(apollo_control$mixing){
    if(anyNA(apollo_draws)) stop("Argument 'apollo_draws' must be provided when estimating mixture models.")
    if(!is.function(apollo_randCoeff)) stop("Argument 'apollo_randCoeff' must be provided when estimating mixture models.")
  }
  if(!anyNA(samples)){
    if(is.data.frame(samples)) samples <- as.matrix(samples)
    if(!is.matrix(samples)) stop("The 'samples' argument must be a matrix.")
    if(nrow(samples)!=nrow(database)) stop("The 'samples' matrix must have as many rows as the database.")
    if(any(samples<0)) stop("The 'samples' matrix must only contain non-negative integers.")
    if(ncol(samples)<2) stop("The 'samples' matrix must have at least two columns.")
  }
  
  ### Separate theta into variable and fixed part
  theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  theta_fix_val <- apollo_beta[apollo_fixed]
  
  ### Start clock
  starttime <- Sys.time()
  
  ### Initial report
  if(anyNA(samples)){
    indivs  <- sort(unique(database[,apollo_control$indivID]))
    nIndivs <- length(indivs)
    if(!silent){
      cat("\n ", nRep, " new datasets will be constructed by randomly sampling", sep="")
      cat("\n ", nIndivs, " individuals with replacement from the original dataset.", sep="")
    }
    tmp <- as.vector(table(database[,apollo_control$indivID])) # n obs per indiv
    tmp <- all(tmp==tmp[1]) # TRUE if all indivs have same number of obs
    if(!tmp & !silent){
      cat("\n Not all individuals have the same number of observations,")
      cat("\n therefore not all generated datasets may have the same")
      cat("\n number of observations.")
    }
  } else {
    nRep <- ncol(samples)
    if(!silent){
      cat("\n ", nRep, " new dataset will be constructed by sampling from the original", sep="")
      cat("\n dataset, as described by the 'samples' matrix provided.")
    }
    tmp <- colSums(samples)
    tmp <- all(tmp==tmp[1])
    if(!tmp & !silent) cat("\n Not all samples have the same number of observations.")
  }
  
  ### Get number of LL components in model and create stacks
  if(!silent) cat("\n Preparing bootstrap.")
  llComponents <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  paramStack   <- matrix(0, nrow=nRep, ncol=length(apollo_beta) , dimnames=list(c(), names(apollo_beta)))
  llStack      <- matrix(0, nrow=nRep, ncol=length(llComponents), dimnames=list(c(), paste0("LL_", names(llComponents))) )
  nObsStack    <- rep(0,nRep)
  rm(llComponents)
  
  # Check for previous bootstrap params & samples files
  tmp <- paste0(name, "_bootstrap_params.csv")
  if(file.exists(tmp)) oldParams <- tryCatch(utils::read.csv(tmp), error=function(e) return(NA)) else oldParams <- NA
  if(!anyNA(oldParams)){
    if(!silent) cat("\n File ", tmp, "\n  found in the working directory.", sep="")
    checkCol  <- colnames(oldParams) == c(colnames(paramStack), colnames(llStack), "obs")
    if(!all(checkCol)){
      oldParams <- NA
      if(!silent) cat("\n  But will be ignored as its dimensions are\n  incompatible with current model.")
    }
  }
  tmp <- paste0(name, "_bootstrap_samples.csv")
  if(file.exists(tmp)) oldSamples <- tryCatch(utils::read.csv(tmp), error=function(e) return(NA)) else oldSamples <- NA
  if(!anyNA(oldSamples)){
    if(!silent) cat("\n File ", tmp, "\n  found in the working directory.", sep="")
    checkRow  <- nrow(oldSamples) == nrow(database)
    if(!checkRow){
      oldSamples <- NA
      if(!silent) cat("\n  But will be ignored as its dimensions are\n  incompatible with current model.")
    }
  }
  if(!anyNA(oldParams) & !anyNA(oldSamples)){
    if(!silent) cat("\n New bootstrap repetitions will be added to the existing results.")
  } else {
    oldParams  <- NA
    oldSamples <- NA
    if((exists("checkCol") | exists("checkRow")) & !silent) cat("\n Old bootstrap repetitions will be discarded.")
  }
  
  ### Create samples matrix if necessary
  if(anyNA(samples)){ 
    samples <- matrix(0, nrow=nrow(database), ncol=nRep, dimnames=list(c(), paste0("sample_",1:nRep)))
    if(!anyNA(oldParams) & !anyNA(oldSamples)) seed <- seed + ncol(oldSamples)
    set.seed(seed)
    for(i in 1:nRep){
      newIndivs <- sample(indivs, replace=TRUE)
      for(j in newIndivs){
        r <- which(database[,apollo_control$indivID]==j)
        samples[r,i] <- samples[r,i] + 1
      }
    }
  }
  
  # BOOTSTRAP LOOP
  if(!silent){
    cat("\n Parameters and LL in each repetition will be written to:")
    cat("\n  ", paste0(apollo_control$modelName, "_bootstrap_params.csv"), sep="")
    cat("\n Vectors with the sample in each repetition will be written to:")
    cat("\n  ", paste0(apollo_control$modelName, "_bootstrap_samples.csv"), sep="")
    cat("\n")
  }
  for(i in 1:nRep){
    newObs    <- rep(1:nrow(database), samples[,i])
    newObs1=0*newObs
    j=2
    while(j<=length(newObs)){
     if(newObs[j]==newObs[j-1]) newObs1[j]=newObs1[j-1]+1
    j=j+1}
    database2 <- database[newObs,]
    database2[,apollo_control$indivID] <- database2[,apollo_control$indivID]+(newObs1-1)/1000
    
    apollo_inputs2 <- apollo_validateInputs(database=database2, silent=TRUE)
    apollo_inputs2$apollo_control$noDiagnostics <- TRUE
    if(apollo_control$mixing) draws <- apollo_makeDraws(apollo_inputs2, silent=TRUE)
    
    # Initialize cluster if user asked for it
    if(apollo_control$nCores==1) cl <- NA else {
      cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs2, silent=TRUE)
      apollo_control$nCores <- length(cl)
    }
    
    # Create logLike function
    apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities,
                                         apollo_inputs2, estimate_settings, cl)
    
    # Estimate
    nObsStack[i] <- nrow(database2)
    if(!silent)cat("\nEstimation cycle ", i, " (", nObsStack[i], " obs)\n", sep="")
    model <- maxLik::maxLik(apollo_logLike, start=theta_var_val,
                            method=estimationRoutine, print.level=estimate_settings$printLevel, finalHessian=FALSE,
                            iterlim=maxIterations, countIter=FALSE, writeIter=FALSE, sumLL=FALSE)
    
    # Check convergence
    succesfulEstimation <- FALSE
    if(exists("model")){
      if(estimationRoutine=="bfgs" & model$code==0) succesfulEstimation <- TRUE
      if(estimationRoutine=="bhhh" & (model$code %in% c(2,8)) ) succesfulEstimation <- TRUE
      if(estimationRoutine=="nr" && model$code<=2) succesfulEstimation <- TRUE
    }
    
    # Write results
    if(succesfulEstimation){
      
      # Closes clusters if using multicore
      if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl)
      
      # Store estimated parameters
      temp <- c(model$estimate, apollo_beta[apollo_fixed])
      temp <- temp[names(apollo_beta)]
      paramStack[i,] <- temp
      
      # Store in-sample LL components
      ll <- apollo_probabilities(c(model$estimate, theta_fix_val), apollo_inputs2, functionality="output")
      for(j in 1:ncol(llStack)) llStack[i,j] <- ifelse(workInLogs, sum(ll[[j]]), sum(log(ll[[j]])))
      
      # Save results from bootstrap iteration
      temp <- cbind(paramStack[1:i,,drop=FALSE], llStack[1:i,,drop=FALSE], obs=nObsStack[1:i])
      fileName <- paste(name, "bootstrap_params.csv", sep="_")
      if(!anyNA(oldParams)) temp <- rbind(oldParams, temp)
      tryCatch(utils::write.csv(temp, fileName, row.names=FALSE),
               error=function(e) if(!silent) cat(" Could not write to\n ", fileName,"\n"))
      temp <- samples[,1:i]
      if(!anyNA(oldSamples)) temp <- cbind(oldSamples, temp)
      fileName <- paste(name, "bootstrap_samples.csv", sep="_")
      tryCatch(utils::write.csv(temp, fileName, row.names=FALSE), 
               error=function(e) if(!silent) cat(" Could not write to\n ", fileName))
      if(!silent) cat("Estimation results written to file.\n")
    } else {
      # Report error but continue with next iteration
      if(!silent)cat("\nERROR: Estimation failed in repetition ", i, ".", sep="")
      if(estimationRoutine=="bfgs" & !silent) print(as.matrix(round(get("lastFuncParam", envir=globalenv()),4)))
    }
  }
  
  # Calculate covariance matrix, save it and print it
  if(!anyNA(oldParams)) paramStack <- rbind(oldParams[,1:ncol(paramStack),drop=FALSE], paramStack)
  includeRow <- apply(paramStack, MARGIN=1, function(x) !all(x==0))
  Sigma      <- stats::cov(paramStack[includeRow,])
  fileName   <- paste0(name, "_bootstrap_vcov.csv")
  tryCatch(utils::write.csv(Sigma, fileName, row.names=FALSE),
           error=function(e) cat("\n Could not write to ", fileName))
  
  if(!silent){
    cat("\n Parameters and LL in each repetition written to:")
    cat("\n  ", paste0(name, "_bootstrap_params.csv"), sep="")
    cat("\n Vectors with the sample in each repetition written to:")
    cat("\n  ", paste0(name, "_bootstrap_samples.csv"), sep="")
    cat("\n Covariance matrix of parameters written to\n  ", fileName, sep="")
    cat("\n")
  }
  
  # Stop clock and return results
  endtime   <- Sys.time()
  timeTaken <- difftime(endtime, starttime, units='auto')
  if(!silent) cat("\nBootstrap processing time: ", format(timeTaken), "\n",sep="")
  #output_matrix <- cbind(paramStack, llStack, nObs=nObsStack)
  if(!silent) cat("\nBootstrap covariance matrix produced\n")
  return(Sigma)
}
