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
#'                                     \item \strong{calledByEstimate}: Logical. TRUE if \code{apollo_bootstrap} is called by \code{apollo_estimate}. FALSE by default.
#'                                     \item \strong{recycle}: Logical. If TRUE, the function will look for old output files and append new repetitions to them. If FALSE, output files will be overwritten.
#'                                   }
#' @return List with three elements.
#'         \itemize{
#'           \item \code{estimates}: Matrix containing the parameter estimates for each repetition. As many rows as repetitions and as many columns as parameters.
#'           \item \code{varcov}: Covariance matrix of the estimated parameters across the repetitions.
#'           \item \code{LL}: Vector of final loglikelihoods of each repetition.
#'         }
#'         This function also writes three output files to the working directory, with the following names ('x' represents the model name):
#'         \itemize{
#'           \item \strong{x_bootstrap_params.csv}: Table containing the parameter estimates, loglikelihood, and number of observations for each repetition.
#'           \item \strong{x_bootstrap_samples.csv}: Table containing the description of the sample used in each repetition. Same format than \code{bootstrap_settings$samples}.
#'           \item \strong{x_bootstrap_vcov}: Table containing the covariance matrix of estimated parameters across the repetitions.
#'         }
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
                                                    silent=FALSE,
                                                    maxLik_settings=list()),
                             bootstrap_settings=list(nRep=30,
                                                     samples=NA,
                                                     seed=24,
                                                     calledByEstimate=FALSE,
                                                     recycle=TRUE)){
  
  ### Set missing settings to default values
  if(is.list(estimate_settings) && !is.null(estimate_settings$maxLik_settings)){
    # Copy values in maxLik_settings into estimate_settings${printLevel, iterlim}
    if(!is.null(estimate_settings$maxLik_settings$printLevel)) estimate_settings$printLevel <- estimate_settings$maxLik_settings$printLevel
    if(!is.null(estimate_settings$maxLik_settings$iterlim)) estimate_settings$maxIterations <- estimate_settings$maxLik_settings$iterlim
  }
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=FALSE, 
                  hessianRoutine="none", printLevel=2L, silent=FALSE, maxLik_settings=list())
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  default <- list(nRep=30, samples=NA, seed=24, calledByEstimate=FALSE, recycle=TRUE)
  tmp <- names(default)[!(names(default) %in% names(bootstrap_settings))] # options missing in bootstrap_settings
  for(i in tmp) bootstrap_settings[[i]] <- default[[i]]
  rm(tmp)
  
  ### Write original apollo_inputs to disk before changing it, and make sure to restore it before finishing
  if(!bootstrap_settings$calledByEstimate){
    saveRDS(apollo_inputs, file=paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_bootstrap"))
    on.exit({
      tmp <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_bootstrap")
      apollo_inputs <- tryCatch(readRDS(tmp), error=function(e) NULL)
      tmp2 <- globalenv()
      if(!is.null(apollo_inputs) & exists('apollo_inputs', envir=tmp2)) assign('apollo_inputs', apollo_inputs, envir=tmp2)
      if(file.exists(tmp)) file.remove(tmp)
      rm(tmp)
    })
  }
  
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
  calledByEstimate  <- bootstrap_settings$calledByEstimate
  recycle           <- bootstrap_settings$recycle
  
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
    samples <- samples[, !(colnames(samples) %in% c(apollo_control$indivID, 'apollo_sequence'))]
    if(!is.matrix(samples)) stop("The 'samples' argument must be a matrix.")
    if(nrow(samples)!=nrow(database)) stop("The 'samples' matrix must have as many rows as the database.")
    if(any(samples<0)) stop("The 'samples' matrix must only contain non-negative integers.")
    if(ncol(samples)<2) stop("The 'samples' matrix must have at least two columns.")
  }
  
  ### Start clock
  starttime <- Sys.time()
  
  ### Initial report
  if(anyNA(samples)){
    indivs  <- unique(database[,apollo_control$indivID])
    nIndivs <- length(indivs)
    if(!silent) apollo_print(paste0(nRep, " new datasets will be constructed by randomly sampling ",
                                    nIndivs, " individuals with replacement from the original dataset."))
    tmp <- as.vector(table(database[,apollo_control$indivID])) # n obs per indiv
    tmp <- all(tmp==tmp[1]) # TRUE if all indivs have same number of obs
    if(!tmp & !silent) apollo_print(paste0("Not all individuals have the same number of observations, ",
                                           "therefore not all generated datasets may have the same ",
                                           "number of observations."))
  } else {
    nRep <- ncol(samples)
    if(!silent) apollo_print(paste0(nRep, " new datasets will be constructed by sampling from the original ",
                                    "dataset, as described by the 'samples' matrix provided."))
    tmp <- colSums(samples)
    tmp <- all(tmp==tmp[1])
    if(!tmp & !silent) apollo_print("Not all samples have the same number of observations.")
  }
  
  ### Get number of LL components in model and create stacks
  if(!silent) apollo_print("Preparing bootstrap.")
  llComponents <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  paramStack   <- matrix(NA, nrow=nRep, ncol=length(apollo_beta) , dimnames=list(c(), names(apollo_beta)))
  llStack      <- matrix(0, nrow=nRep, ncol=length(llComponents), dimnames=list(c(), paste0("LL_", names(llComponents))) )
  nObsStack    <- rep(0,nRep)
  rm(llComponents)
  
  # Check for previous bootstrap params & samples files
  fileNameParams <- paste(name, "bootstrap_params.csv", sep="_")
  if(file.exists(fileNameParams)) oldParams <- tryCatch(utils::read.csv(fileNameParams), error=function(e) return(NA)) else oldParams <- NA
  if(!anyNA(oldParams)){
    if(!silent) apollo_print(paste0("File ", fileNameParams, " found in the working directory."))
    checkCol  <- colnames(oldParams) == c(colnames(paramStack), colnames(llStack), "obs")
    if(!all(checkCol)){
      oldParams <- NA
      if(!silent) apollo_print("But will be ignored as its dimensions are incompatible with current model.")
    }
  }
  fileNameSamples <- paste0(name, "_bootstrap_samples.csv")
  if(file.exists(fileNameSamples)) oldSamples <- tryCatch(utils::read.csv(fileNameSamples), error=function(e) return(NA)) else oldSamples <- NA
  if(!anyNA(oldSamples)){
    if(!silent) apollo_print(paste0("File ", fileNameSamples, " found in the working directory."))
    checkRow  <- nrow(oldSamples) == nrow(database)
    if(!checkRow){
      oldSamples <- NA
      if(!silent) apollo_print("But will be ignored as its dimensions are incompatible with current model.")
    }
  }
  if(recycle && !anyNA(oldParams) & !anyNA(oldSamples)){
    if(!silent) apollo_print("New bootstrap repetitions will be added to the existing results.")
  } else {
    oldParams  <- NA
    oldSamples <- NA
    if((exists("checkCol") | exists("checkRow")) & !silent){
      if(recycle) apollo_print("Old bootstrap repetitions will be discarded.")
      if(!recycle) apollo_print("Old bootstrap repetitions will be discarded, as setting 'recycle' is FALSE.")
    } 
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
    apollo_print(paste0("Parameters and LL in each repetition will be written to: ", fileNameParams))
    apollo_print(paste0("Vectors showing sampling rate for each observation in each repetition written to: ", fileNameSamples))
    apollo_print("\n")
  }
  for(i in 1:nRep){
    # Create new database and change indivID for repeated indivs
    newObs  <- rep(1:nrow(database), samples[,i])
    newObs1 <- 0*newObs
    for(j in 2:length(newObs)) if(newObs[j]==newObs[j-1]) newObs1[j]=newObs1[j-1]+1
    database2 <- database[newObs,]
    database2[,apollo_control$indivID] <- database2[,apollo_control$indivID]+(newObs1-1)/1000
    rm(newObs, newObs1)
    database2 <- database2[order(database2[,apollo_control$indivID]),]
    
    # Create apollo_inputs2
    apollo_inputs2 <- apollo_validateInputs(database=database2, silent=TRUE)
    apollo_inputs2$apollo_control$noDiagnostics <- TRUE
    apollo_inputs2$apollo_control$noValidation  <- TRUE
    apollo_inputs2$silent <- TRUE
    if(apollo_control$mixing) draws <- apollo_makeDraws(apollo_inputs2, silent=TRUE)
    
    # Estimate model
    nObsStack[i] <- nrow(database2)
    if(!silent) apollo_print(paste0("Estimation cycle ", i, " (", nObsStack[i], " obs)"))
    estimate_settings$hessianRoutine="none"
    ### change 28 July
    estimate_settings$silent=TRUE
    ### end change
    model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, 
                             apollo_inputs2, estimate_settings)
    
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
      #if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl)
      
      # Store estimated parameters
      temp <- c(model$estimate, apollo_beta[apollo_fixed])
      temp <- temp[names(apollo_beta)]
      paramStack[i,] <- temp
      
      # Store in-sample LL components
      ll <- apollo_probabilities(model$estimate, apollo_inputs2, functionality="output")
      for(j in 1:ncol(llStack)) llStack[i,j] <- ifelse(workInLogs, sum(ll[[j]]), sum(log(ll[[j]])))
      
      # Save results from bootstrap iteration
      temp <- cbind(paramStack[1:i,,drop=FALSE], llStack[1:i,,drop=FALSE], obs=nObsStack[1:i])
      if(!anyNA(oldParams)) temp <- rbind(oldParams, temp)
      tryCatch(utils::write.csv(temp, fileNameParams, row.names=FALSE),
               error=function(e) if(!silent) apollo_print(paste0("Could not write to ", fileNameParams)))
      temp <- samples[,1:i]
      if(!anyNA(oldSamples)) temp <- cbind(oldSamples, temp) else temp <- cbind(database[,c(apollo_control$indivID, 
                                                                                            'apollo_sequence')], temp)
      tryCatch(utils::write.csv(temp, fileNameSamples, row.names=FALSE), 
               error=function(e) if(!silent) apollo_print(paste0("Could not write to ", fileNameSamples)))
      if(!silent && !calledByEstimate) apollo_print("Estimation results written to file.\n")
    } else {
      # Report error but continue with next iteration
      if(!silent) apollo_print(paste0("ERROR: Estimation failed in repetition ", i, "."))
      if(estimationRoutine=="bfgs" & !silent) print(as.matrix(round(get("lastFuncParam", envir=globalenv()),4)))
    }
  }
  
  # Calculate covariance matrix, save it and print it
  if(!anyNA(oldParams)) paramStack <- rbind(oldParams[,1:ncol(paramStack),drop=FALSE], paramStack)
  includeRow <- apply(paramStack, MARGIN=1, function(x) !all(is.na(x)))
  Sigma      <- stats::cov(paramStack[includeRow,])
  # remove fixed params from Sigma (commented out because apollo_estimate includes the fixed ones)
  if(length(apollo_fixed)>0){
  Sigma      <- Sigma[-which(colnames(Sigma) %in% apollo_fixed),-which(colnames(Sigma) %in% apollo_fixed)]
  if(is.null(rownames(Sigma)) && nrow(Sigma)==ncol(Sigma)) rownames(Sigma)=colnames(Sigma)
  }
  fileName   <- paste0(name, "_bootstrap_vcov.csv")
  tryCatch(utils::write.csv(Sigma, fileName, row.names=TRUE),
           error=function(e) apollo_print(paste0("Could not write to ", fileName)))
  
  if(!silent){
    apollo_print(paste0("\n"))
    apollo_print(paste0("Finished bootstrap runs."))
    apollo_print(paste0("Parameters and LL for each repetition written to: ", fileNameParams))
    apollo_print(paste0("Vectors showing sampling rate for each observation in each repetition written to: ", fileNameSamples))
    apollo_print(paste0("Covariance matrix of parameters written to: ", fileName))
    apollo_print(paste0("\n"))
    if(!calledByEstimate){
      apollo_print(paste0("Mean LL across runs: ", round(mean(llStack[,ncol(llStack)]),2)))
      apollo_print("\nMean parameter values across runs: ")
      txt=(data.frame(round(colMeans(paramStack[includeRow,]),4)))
      colnames(txt)="Estimate"
      print(txt)
      apollo_print("\nCovariance matrix across runs:")
      longNames <- colnames(Sigma)
      maxLen    <- max(nchar(longNames))
      for(i in 1:length(longNames)) longNames[i] <- paste0(paste0(rep(" ", maxLen-nchar(longNames[i])), collapse=""), longNames[i])
      tmp <- Sigma #round(Sigma,4)
      colnames(tmp) <- longNames
      if(nrow(tmp)==ncol(tmp)) rownames(tmp) <- longNames
      print(tmp, digits=4)
    }
  }
  
  # Stop clock and return results
  endtime   <- Sys.time()
  timeTaken <- difftime(endtime, starttime, units='auto')
  if(!silent) apollo_print(paste0("Bootstrap processing time: ", format(timeTaken)))
  #output_matrix <- cbind(paramStack, llStack, nObs=nObsStack)
  #if(!silent) apollo_print("Bootstrap covariance matrix produced")
  return(invisible(list(estimates=paramStack[includeRow,],
                        varcov=Sigma,
                        LL=llStack[,ncol(llStack)])))
}
