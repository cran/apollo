#' Out-of-sample fit (LL)
#'
#' Randomly generates estimation and validation samples, estimates the model on the first and 
#' calculates the likelihood for the second, then repeats.
#'
#' A common way to test for overfitting of a model is to measure its fit on a sample not used 
#' during estimation that is, measuring its out-of-sample fit. A simple way to do this is splitting 
#' the complete available dataset in two parts: an estimation sample, and a validation sample. 
#' The model of interest is estimated using only the estimation sample, and then those estimated 
#' parameters are used to measure the fit of the model (e.g. the log-likelihood of the model)
#' on the validation sample. Doing this with only one validation sample, however, may lead to biased 
#' results, as a particular validation sample need not be representative of the population. One way to 
#' minimise this issue is to randomly draw several pairs of estimation and validation samples from the 
#' complete dataset, and apply the procedure to each pair.
#' The splitting of the database into estimation and validaion samples is done at the individual level, 
#' not at the observation level.
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
#' @param outOfSample_settings List. Options defining the sampling procedure. The following are valid options.
#'                                   \describe{
#'                                     \item{nRep}{Numeric scalar. Number of times a different pair of estimation and
#'                                                 validation sets are to be extracted from the full database.
#'                                                 Default is 30.}
#'                                     \item{validationSize}{Numeric scalar. Size of the validation sample. Can be a percentage of the sample (0-1) or the number of individuals in the validation sample (>1). Default is 0.1.}
#'                                   }
#' @return A matrix with the log-likelihood in both the estimation and validation samples. If the 
#'         model has multiple components, the log-likelihood is reported for each of them.
#'         A more complete matrix also containing the estimates is  written to a file 
#'         called <model_name>_outOfSample.csv in the current working directory.
#' @export
#' @importFrom maxLik maxLik
apollo_outOfSample <- function(apollo_beta, apollo_fixed,
                               apollo_probabilities, apollo_inputs,
                               estimate_settings=list(estimationRoutine="bfgs",
                                                      maxIterations=200,
                                                      writeIter=FALSE,
                                                      hessianRoutine="numDeriv",
                                                      printLevel=3L,
                                                      silent=TRUE),
                               outOfSample_settings=list(nRep=10,
                                                         validationSize=0.1)){
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, hessianRoutine="numDeriv", printLevel=3L, silent=FALSE)
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  default <- list(nRep=10, validationSize=0.1)
  tmp <- names(default)[!(names(default) %in% names(outOfSample_settings))] # options missing in estimate_settings
  for(i in tmp) outOfSample_settings[[i]] <- default[[i]]

  # Extract values from apollo_inputs
  database         <- apollo_inputs$database
  apollo_inputs$apollo_control$noDiagnostics = TRUE
  apollo_control   <- apollo_inputs$apollo_control
  apollo_draws     <- apollo_inputs$apollo_draws
  apollo_randCoeff <- apollo_inputs$apollo_randCoeff
  workInLogs       <- apollo_inputs$apollo_control$workInLogs
  
  # Extract values from estimate_settings and estimate_settings
  estimationRoutine <- estimate_settings$estimationRoutine
  maxIterations     <- estimate_settings$maxIterations
  nRep              <- outOfSample_settings$nRep
  validationSize    <- outOfSample_settings$validationSize


  # Validate arguments
  estimationRoutine <- tolower(estimationRoutine)
  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh' or 'nr'.")
  if( (length(apollo_fixed) > 0) & any(!(apollo_fixed %in% names(apollo_beta))) ) stop("Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'")
  maxIterations <- round(maxIterations,0)
  if(maxIterations < 1) stop("Need at least one iteration!")
  if(workInLogs != TRUE) workInLogs=FALSE
  if(validationSize < 0) stop("validationSize must be positive.")
  if(apollo_control$mixing){
    if(anyNA(apollo_draws)) stop("Argument 'apollo_draws' must be provided when estimating mixture models.")
    if(!is.function(apollo_randCoeff)) stop("Argument 'apollo_randCoeff' must be provided when estimating mixture models.")
  }

  # Separate theta into variable and fixed part
  theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  theta_fix_val <- apollo_beta[apollo_fixed]

  # Start clock
  starttime <- Sys.time()

  # Prepare loop
  indivs  <- sort(unique(database[,apollo_control$indivID]))
  nIndivs <- length(indivs)
  if(validationSize < 1) validationSize <- round(outOfSample_settings$validationSize*nIndivs)
  if(!(1<=validationSize & validationSize<nIndivs)) stop("validationSize must be between 1 and (nIndivs-1).")
  cat("\n Number of individuals for estimation : ", nIndivs - validationSize, sep="")
  cat("\n Number of individuals for forecasting: ", validationSize, sep="")
  cat("\n Total number of individuals in sample: ", nIndivs, sep="")
  cat("\n")

  # Get number of LL components in model
  cat("\n Preparing loop.")
  cat("\n")
  llComponents       <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  paramStack         <- matrix(0, nrow=nRep, ncol=length(apollo_beta) , dimnames=list(c(), names(apollo_beta)))
  llInSampleStack    <- matrix(0, nrow=nRep, ncol=length(llComponents), dimnames=list(c(), paste0("inSample_", names(llComponents))) )
  llOutOfSampleStack <- matrix(0, nrow=nRep, ncol=length(llComponents), dimnames=list(c(), paste0("outOfSample_", names(llComponents))) )
  nObsStack <- rep(0,nRep)
  rm(llComponents)

  # BOOTSTRAP LOOP
  cat("\n Result of apollo_outOfSample will be written to:\n  ",paste(apollo_control$modelName, "outOfSample.csv", sep="_"), "\n")
  for(i in 1:nRep){

    # Filter database and create draws
    for(j in 1:i) obsForecast <- database[,apollo_control$indivID] %in% sample(indivs, size=validationSize)
    database2   <- database[!obsForecast,]
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
    cat("\nEstimation cycle ", i, sep="")
    cat("\nUsing ",nrow(database2)," observations\n", sep="")
    nObsStack[i]=nrow(database2)
    model <- maxLik::maxLik(apollo_logLike, start=theta_var_val,
                            method=estimationRoutine, print.level=2, finalHessian=FALSE,
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
      temp           = c(model$estimate, apollo_beta[apollo_fixed])
      temp = temp[names(apollo_beta)]
      paramStack[i,] <- temp

      # Store in-sample LL components
      llin <- apollo_probabilities(c(model$estimate, theta_fix_val), apollo_inputs2, functionality="output")
      for(j in 1:ncol(llInSampleStack)) llInSampleStack[i,j] <- ifelse(workInLogs, sum(llin[[j]]), sum(log(llin[[j]])))

      # Store out-of-sample LL components
      database2 <- database[obsForecast,]
      apollo_inputs2 <- apollo_validateInputs(database=database2, silent=TRUE)
	  apollo_inputs2$apollo_control$noDiagnostics <- TRUE
      llout <- apollo_probabilities(c(model$estimate, theta_fix_val), apollo_inputs2, functionality="output")
      for(j in 1:ncol(llOutOfSampleStack)) llOutOfSampleStack[i,j] <- ifelse(workInLogs, sum(llout[[j]]), sum(log(llout[[j]])))

      # Save results from bootstrap iteration
      utils::write.csv(cbind(paramStack, llInSampleStack, llOutOfSampleStack,inSampleObs=nObsStack,outOfSampleObs=nrow(database)-nObsStack), paste(apollo_control$modelName, "outOfSample.csv", sep="_"), row.names=FALSE)
      cat("Estimation results written to file.\n")
    } else {
      # Report error but continue with next iteration
      cat("\nERROR: Estimation failed in cycle ", i, ".", sep="")
      if(estimationRoutine=="bfgs") print(as.matrix(round(get("lastFuncParam", envir=globalenv()),4)))
    }

  }

  # Stop clock
  endtime   <- Sys.time()
  timeTaken <- difftime(endtime, starttime, units='auto')
  cat("\nProcessing time: ", format(timeTaken), "\n",sep="")

  output_matrix=cbind(llInSampleStack[,"inSample_model"]/nObsStack,llOutOfSampleStack[,"outOfSample_model"]/(nrow(database)-nObsStack))
  colnames(output_matrix)=c("LL per obs in estimation sample","LL per obs in validation sample")

  return(output_matrix)
}
