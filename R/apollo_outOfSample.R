#' Cross-validation of fit (LL)
#'
#' Randomly generates estimation and validation samples, estimates the model on 
#' the first and calculates the likelihood for the second, then repeats.
#'
#' A common way to test for overfitting of a model is to measure its fit on a 
#' sample not used during estimation that is, measuring its out-of-sample fit. 
#' A simple way to do this is splitting the complete available dataset in two 
#' parts: an estimation sample, and a validation sample. 
#' The model of interest is estimated using only the estimation sample, and 
#' then those estimated parameters are used to measure the fit of the model 
#' (e.g. the log-likelihood of the model) on the validation sample. Doing this 
#' with only one validation sample, however, may lead to biased results, as a 
#' particular validation sample need not be representative of the population. 
#' One way to minimise this issue is to randomly draw several pairs of 
#' estimation and validation samples from the complete dataset, and apply the 
#' procedure to each pair.
#' 
#' The splitting of the database into estimation and validation samples is done 
#' at the individual level, not at the observation level. If the sampling wants 
#' to be done at the individual level (not recommended on panel data), then the 
#' optional \code{outOfSample_settings$samples} argument should be provided.
#' 
#' This function writes two different files to the working/output directory:
#' \itemize{
#'   \item \strong{\code{modelName_outOfSample_params.csv}}: Records the 
#'         estimated parameters, final log-likelihood, and number of 
#'         observations on each repetition.
#'   \item \strong{\code{modelName_outOfSample_samples.csv}}: Records the 
#'         sample composition of each repetition.
#' }
#' The first two files are updated throughout the run of this function, while 
#' the last one is only written once the function finishes.
#' 
#' When run, this function will look for the two files above in the 
#' working/output directory. If they are found, the function will attempt to 
#' pick up re-sampling from where those files left off. This is useful in cases 
#' where the original bootstrapping was interrupted, or when additional 
#' re-sampling wants to be performed.
#' 
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in 
#'                     \code{apollo_beta}) of parameters whose value should not 
#'                     change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to 
#'                             be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric 
#'                                  vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List 
#'                                  containing options of the model. See 
#'                                  \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. 
#'                                  Can be either \strong{\code{"components"}}, 
#'                                  \code{"conditionals"}, \code{"estimate"} 
#'                                  (default), \code{"gradient"}, 
#'                                  \code{"output"}, \code{"prediction"}, 
#'                                  \code{"preprocess"}, \code{"raw"}, 
#'                                  \code{"report"}, \code{"shares_LL"}, 
#'                                  \code{"validate"} or \code{"zero_LL"}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function 
#'                      \link{apollo_validateInputs}.
#' @param estimate_settings List. Options controlling the estimation process. 
#'                          See \link{apollo_estimate}.
#' @param outOfSample_settings List. Contains settings for this function. User 
#'                             input is required for all settings except those 
#'                             with a default or marked as optional. 
#'                                   \describe{
#'                                     \item{nRep}{Numeric scalar. Number of 
#'                                                 times a different pair of 
#'                                                 estimation and validation 
#'                                                 sets are to be extracted 
#'                                                 from the full database.
#'                                                 Default is 30.}
#'                                     \item{rmse}{Character matrix with two
#'                                                 columns. Used to calculate 
#'                                                 Root Mean Squared Error (RMSE) 
#'                                                 of prediction. The first 
#'                                                 column must contain the names 
#'                                                 of observed outcomes in the 
#'                                                 database. The second column 
#'                                                 must contain the names of the 
#'                                                 predicted outcomes as 
#'                                                 returned by 
#'                                                 \code{apollo_prediction}. 
#'                                                 If omitted or NULL, no RMSE 
#'                                                 is calculated. This only 
#'                                                 works for models with a 
#'                                                 single component.}
#'                                     \item{samples}{Numeric matrix or 
#'                                                    data.frame. Optional 
#'                                                    argument. Must have as 
#'                                                    many rows as observations 
#'                                                    in the \code{database}, 
#'                                                    and as many columns as 
#'                                                    number of  repetitions 
#'                                                    wanted. Each column 
#'                                                    represents a re-sample, 
#'                                                    and each element must be 
#'                                                    a 0 if the observation 
#'                                                    should be assigned to the 
#'                                                    estimation sample, or 1 
#'                                                    if the observation should 
#'                                                    be assigned to the 
#'                                                    prediction sample. If this 
#'                                                    argument is provided, then 
#'                                                    \code{nRep} and 
#'                                                    \code{validationSize} are 
#'                                                    ignored. Note that this 
#'                                                    allows sampling at the 
#'                                                    observation rather than 
#'                                                    the individual level.}
#'                                     \item{validationSize}{Numeric scalar. 
#'                                                           Size of the 
#'                                                           validation sample. 
#'                                                           Can be a percentage 
#'                                                           of the sample (0-1) 
#'                                                           or the number of 
#'                                                           individuals in the 
#'                                                           validation sample 
#'                                                           (>1). Default is 
#'                                                           0.1.}
#'                                   }
#' @return A matrix with the average log-likelihood per observation for both the 
#'         estimation and validation samples, for each repetition. Two additional 
#'         files with further details are written to the working/output directory.
#' @export
#' @importFrom maxLik maxLik
#' @importFrom utils read.csv
apollo_outOfSample <- function(apollo_beta, apollo_fixed,
                               apollo_probabilities, apollo_inputs,
                               estimate_settings=list(estimationRoutine="bgw",
                                                      maxIterations=200,
                                                      writeIter=FALSE,
                                                      hessianRoutine="none",
                                                      printLevel=3L,
                                                      silent=TRUE),
                               outOfSample_settings=list(nRep=10,
                                                         validationSize=0.1,
                                                         samples=NA, 
                                                         rmse=NULL)){
  ### Set missing settings to default values
  default <- list(estimationRoutine="bgw", maxIterations=200, writeIter=FALSE, 
                  hessianRoutine="none", printLevel=3L, silent=TRUE)
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))]
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  default <- list(nRep=10, validationSize=0.1, samples=NA, rmse=NULL)
  tmp <- names(default)[!(names(default) %in% names(outOfSample_settings))]
  for(i in tmp) outOfSample_settings[[i]] <- default[[i]]
  
  # Start clock
  starttime <- Sys.time()
  
  ### Write original apollo_inputs to disk before changing it, and make sure to 
    # restore it before finishing
  saveRDS(apollo_inputs, 
          file=paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName, 
                      "_outOfSample"))
  on.exit({
    tmp <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,
                  "_outOfSample")
    apollo_inputs <- tryCatch(readRDS(tmp), error=function(e) NULL)
    if(!is.null(apollo_inputs)){
      tmp2 <- globalenv()
      assign('apollo_inputs', apollo_inputs, envir=tmp2)
    } 
    if(file.exists(tmp)) file.remove(tmp)
    rm(tmp)
  })
  
  # Extract values from apollo_inputs
  database         <- apollo_inputs$database
  apollo_inputs$apollo_control$noDiagnostics = TRUE
  apollo_control   <- apollo_inputs$apollo_control
  apollo_draws     <- apollo_inputs$apollo_draws
  apollo_randCoeff <- apollo_inputs$apollo_randCoeff
  apollo_lcPars    <- apollo_inputs$apollo_lcPars
  workInLogs       <- apollo_inputs$apollo_control$workInLogs
  if(!is.null(apollo_inputs$apollo_control$seed)){
    seed <- apollo_inputs$apollo_control$seed + 3
  } else seed <- 13 + 3
  
  # Extract values from estimate_settings and outOfSample_settings
  estimationRoutine <- estimate_settings$estimationRoutine
  maxIterations     <- estimate_settings$maxIterations
  nRep              <- outOfSample_settings$nRep
  validationSize    <- outOfSample_settings$validationSize
  samples           <- outOfSample_settings$samples
  rmse           <- outOfSample_settings$rmse
  
  # Validate arguments
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  estimationRoutine <- tolower(estimationRoutine)
  maxIterations     <- round(maxIterations, 0)
  test <- estimationRoutine %in% c("bfgs","bgw","bhhh", "nr")
  if(!test) stop("SYNTAX ISSUE - Invalid estimationRoutine.",
                 "Use 'bfgs', 'bgw', 'bhhh' or 'nr'.")
  test <- length(apollo_fixed)==0 || all(apollo_fixed %in% names(apollo_beta))
  if(!test) stop("SYNTAX ISSUE - Some parameters included in 'apollo_fixed' ",
                 "are not included in 'apollo_beta'")
  if(maxIterations < 1) stop("SYNTAX ISSUE - Need at least one iteration!")
  if(workInLogs != TRUE) workInLogs <- FALSE
  if(validationSize < 0) stop("SYNTAX ISSUE - validationSize must be positive.")
  test <- !apollo_control$mixing || (!anyNA(apollo_draws) && 
                                       is.function(apollo_randCoeff))
  if(!test) stop("SYNTAX ISSUE - When using mixture models, objects ", 
                 "'apollo_draws' and 'apollo_randCoeff' must be created by ", 
                 "the user")
  if(!anyNA(samples)){
    if(is.data.frame(samples)) samples <- as.matrix(samples)
    test <- is.matrix(samples) && nrow(samples)==nrow(database) && 
      all(samples %in% 0:1) && ncol(samples)>=2
    if(!test) stop("INPUT ISSUE - The 'samples' argument must be a numeric ", 
                   "or logical matrix with the same number of rows as the ", 
                   "database, at least two columns, and only 0/1 elements.")
  }
  indivs  <- unique(database[,apollo_control$indivID])
  nIndivs <- length(indivs)
  if(anyNA(samples)){
    if(validationSize<1) validationSize <- round(validationSize*nIndivs)
    if(validationSize>nIndivs) stop("SYNTAX ISSUE - 'validationSize' must be ", 
                                    "between 1 and (nIndivs-1).")
  }
  apollo_print("Testing likelihood function.")
  llComponents <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
  if(!is.null(rmse)){
    test <- length(llComponents)==1
    if(!test) stop("SYNTAX ISSUE - The model contains more than one component, ", 
                   "but RMSE can only be calculated for single-component models.")
    test <- is.character(rmse)
    test <- test && (is.vector(rmse) || (is.matrix(rmse) && ncol(rmse)==2))
    if(!test) stop("SYNTAX ISSUE - Argument 'rmse' must be a character ", 
                   "matrix with 2 columns, or a character vector.")
    if(is.vector(rmse)) rmse <- matrix(rmse, nrow=length(rmse), ncol=2)
    colnames(rmse) <- c("obs", "for")
    test <- all(rmse[,"obs"] %in% names(apollo_inputs$database)) &&
      all(sapply(apollo_inputs$database[,rmse[,"obs"]], is.numeric))
    if(!test) stop("INPUT ISSUE - All outcomes named in rmse must exist ", 
                   "in the database, and be numeric.")
  }
  
  
  # Initial report
  apollo_print(paste0(nRep, " separate runs will be conducted, each using a ",
                      "random subset of ", round(100*(nIndivs - validationSize)/nIndivs, 0), 
                      "% of individuals for estimation, and the remainder ", 
                      "for validation."))
  if(anyNA(samples)){
    apollo_print("Number of individuals")
    apollo_print(paste0("- for estimation   : ", nIndivs - validationSize))
    apollo_print(paste0("- for forecasting  : ", validationSize))
    apollo_print(paste0("- in sample (total): ", nIndivs))
  } else {
    nRep <- ncol(samples)
    nPre <- colSums(samples)
    if(all(nPre==nPre[1])) nPre <- nPre[1]
    apollo_print("Number of observations")
    apollo_print(paste0("- for estimation : ", ifelse(length(nPre)==1,
                                                      nrow(samples) - nPre, 
                                                      "Changes by sample")))
    apollo_print(paste0("- for forecasting: ", ifelse(length(nPre)==1, nPre, 
                                                      "Changes by sample")))
    apollo_print(paste0("Number of individuals in sample: ", nIndivs))
  }
  
  # Initialise matrices to store results
  cNam               <- names(llComponents)
  paramStack         <- matrix(0, nrow=nRep, ncol=length(apollo_beta) , 
                               dimnames=list(c(), names(apollo_beta)))
  llInSampleStack    <- matrix(0, nrow=nRep, ncol=length(llComponents), 
                               dimnames=list(c(), paste0("inSample_", cNam)))
  llOutOfSampleStack <- matrix(0, nrow=nRep, ncol=length(llComponents), 
                               dimnames=list(c(), paste0("outOfSample_", cNam)))
  nObsStack          <- rep(0,nRep)
  if(anyNA(samples)){
    samples <- matrix(0, nrow=nrow(database), ncol=nRep, 
                      dimnames=list(NULL, paste0("sample_",1:nRep)))
    set.seed(seed)
    for(i in 1:nRep) samples[,i] <- database[,apollo_control$indivID] %in% 
      sample(indivs, size=validationSize)
  }
  if(is.character(rmse)){
    RMSE <- paste0("rmse_", c(rmse[,"obs"], "aggregate"))
    RMSE <- matrix(0, nrow=nRep, ncol=nrow(rmse) + 1, 
                   dimnames=list(NULL, RMSE))
  } else RMSE <- NULL
  rm(llComponents)
  
  # Check if there are previous results. If so, load them
  fileNameParams <- paste0(apollo_control$outputDirectory, 
                           apollo_control$modelName, "_outOfSample_params.csv")
  fileNameSample <- paste0(apollo_control$outputDirectory, 
                           apollo_control$modelName, "_outOfSample_samples.csv")
  nRun <- 0
  if(file.exists(fileNameParams) && file.exists(fileNameSample)){
    apollo_print("Old output files found, they will be recycled.")
    # Read params
    tmp <- as.matrix(read.csv(fileNameParams))
    tmp1 <- tmp[, colnames(tmp) %in% names(apollo_beta), drop=FALSE] # paramStack
    tmp2 <- tmp[, grep('^inSample_'   , colnames(tmp)), drop=FALSE] # llInSampleStack
    tmp3 <- tmp[, grep('^outOfSample_', colnames(tmp)), drop=FALSE] # llOutOfSampleStack
    tmp4 <- as.matrix(read.csv(fileNameSample)) # samples
    tmp5 <- tmp[, grep("^rmse_"       , colnames(tmp)), drop=FALSE] # RMSE
    nRun <- sum(tmp3[,ncol(tmp3)]!=0)
    if(nrow(tmp1)>nRun) tmp1 <- tmp1[1:nRun,] 
    if(nrow(tmp2)>nRun) tmp2 <- tmp2[1:nRun,] 
    if(nrow(tmp3)>nRun) tmp3 <- tmp3[1:nRun,] 
    if(ncol(tmp4)>nRun) tmp4 <- tmp4[, 1:nRun]
    if(ncol(tmp5)>nRun) tmp5 <- tmp5[1:nRun,]
    # Check that result files match
    test <- nRun>0 && ncol(tmp1)==ncol(paramStack)
    test <- test && ncol(tmp2)==ncol(llInSampleStack)
    test <- test && ncol(tmp3)==ncol(llOutOfSampleStack)
    test <- test && nrow(tmp4)==nrow(database)
    if(test){
      apollo_print(paste0(nRun, ' repetitions recovered from old result files. ', 
                          nRep, ' new repetitions will be added.'))
      nRep <- nRun + nRep
      # Expand paramStack, llInSampleStack, llOutOfSampleStack, nObsStack and 
      # samples to fit old results
      paramStack         <- rbind(tmp1, paramStack)
      llInSampleStack    <- rbind(tmp2, llInSampleStack)
      llOutOfSampleStack <- rbind(tmp3, llOutOfSampleStack)
      samples            <- cbind(tmp4, samples)
      colnames(samples)  <- paste0("sample_",1:nRep)
      set.seed(seed)
      for(i in 1:nRep) if(i>nRun) samples[,i] <- database[,apollo_control$indivID] %in% 
        sample(indivs, size=validationSize)
      nObsStack          <- nrow(database) - colSums(samples)
      if(!is.null(RMSE)) RMSE <- rbind(tmp5, RMSE)
    } else { 
      nRun <- 0
      apollo_print("Old result files do not match current model, or contained ", 
                   "no completed runs, so they will be overwritten.")
    }
    rm(tmp, tmp1, tmp2, tmp3, tmp4, test)
  }
  if(nRun>=nRep) stop("INPUT ISSUE - All requested repetitions already exist ", 
                      "in file ", fileNameParams)
  
  # OUT OF SAMPLE LOOP
  apollo_print(paste0("Estimated parameters and log-likelihoods for each ",
                      "sample will be written to: ", fileNameParams))
  apollo_print(paste0("The matrix defining the observations used in each ", 
                      "repetition will be written to: ", fileNameSample))
  apollo_print("\n")
  for(i in (nRun + 1):nRep){
    
    # Filter database and create draws
    database2   <- database[samples[,i]==0,]
    apollo_inputs <- apollo_validateInputs(
      apollo_beta      = apollo_beta, 
      apollo_fixed     = apollo_fixed, 
      database         = database2  , 
      apollo_control   = apollo_control, 
      apollo_HB        = apollo_inputs$apollo_HB, 
      apollo_randCoeff = apollo_inputs$apollo_randCoeff, 
      apollo_lcPars    = apollo_inputs$apollo_lcPars, 
      silent           = TRUE,
      recycle          = TRUE)
    apollo_inputs$apollo_control$noDiagnostics <- TRUE
    apollo_inputs$apollo_control$noValidation  <- TRUE
    
    # Estimate
    apollo_print(paste0('Estimation cycle ', i, ' (', nrow(database2), ' obs.)'))
    nObsStack[i] <- nrow(database2)
    estimate_settings$hessianRoutine="none"
    model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, 
                             apollo_inputs, estimate_settings)
    
    # If apollo_inputs is missing database or draws, restore them
    if(is.null(apollo_inputs$database)) apollo_inputs$database <- database2
    if(is.null(apollo_inputs$draws)){
      if(!apollo_inputs$apollo_control$mixing) apollo_inputs$draws <- NA else {
        apollo_inputs$draws <- apollo_makeDraws(apollo_inputs, silent=TRUE)
      }
    }
    
    # Check convergence
    successfulEstimation <- FALSE
    if(exists("model")){
      test <- estimationRoutine=="bfgs" & model$code==0
      test <- test || ( estimationRoutine=="bgw" && model$code %in% c(3,4,5,6) )
      test <- test || ( estimationRoutine=="bhhh" & (model$code %in% c(2,8)) )
      test <- test || ( estimationRoutine=="nr" && model$code<=2 )
      if(test) successfulEstimation <- TRUE
    }
    
    # Write results
    if(successfulEstimation){
      
      # Store estimated parameters
      temp <- c(model$estimate, apollo_beta[apollo_fixed])
      temp <- temp[names(apollo_beta)]
      paramStack[i,] <- temp
      
      # Store in-sample LL components
      llin <- apollo_probabilities(model$estimate, apollo_inputs, functionality="output")
      for(j in 1:ncol(llInSampleStack)) llInSampleStack[i,j] <- ifelse(workInLogs, 
                                                                       sum(llin[[j]]), 
                                                                       sum(log(llin[[j]])))
      
      # Store out-of-sample LL components
      database2 <- database[samples[,i]>0,]
      apollo_inputs <- apollo_validateInputs(
        apollo_beta      = apollo_beta, 
        apollo_fixed     = apollo_fixed, 
        database         = database2  , 
        apollo_control   = apollo_control, 
        apollo_HB        = apollo_inputs$apollo_HB, 
        apollo_randCoeff = apollo_inputs$apollo_randCoeff, 
        apollo_lcPars    = apollo_inputs$apollo_lcPars, 
        silent           = TRUE,
        recycle          = TRUE)
      apollo_inputs$apollo_control$noDiagnostics <- TRUE
      llout <- apollo_probabilities(model$estimate, apollo_inputs, 
                                    functionality="output")
      for(j in 1:ncol(llOutOfSampleStack)){
        llOutOfSampleStack[i,j] <- ifelse(workInLogs, sum(llout[[j]]), 
                                          sum(log(llout[[j]])))
      } 
      
      # Store RMSE
      if( is.character(rmse) ){
        pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
        if( !is.data.frame(pred) ) stop("INCORRECT FUNCTION/SETTING USE - The ", 
                                        "model uses an incompatible component ", 
                                        "type. RMSE cannot be calculated.")
        test <- all(rmse[,"for"] %in% names(pred))
        if(!test) stop("INPUT ISSUE - The requested forecasted values in ", 
                       "argument 'rmse' are not generated by the model.")
        pred <- pred[, rmse[,"for"]] # keep only relevant forecasts
        obse <- apollo_inputs$database[,rmse[,"obs"]]
        if(nrow(rmse)==1){
          RMSE[i, 1:nrow(rmse)] <- sqrt(mean((pred - obse)^2))
          RMSE[i, 1+nrow(rmse)] <- sum(pred) - sum(obse)
        } else {
          RMSE[i, 1:nrow(rmse)] <- sqrt(colMeans((pred - obse)^2))
          RMSE[i, 1+nrow(rmse)] <- sqrt(mean((colSums(pred) - colSums(obse))^2))
        }
      }
      
      # Save results from cross-validation iteration
      utils::write.csv(cbind(paramStack, 
                             llInSampleStack, 
                             llOutOfSampleStack,
                             inSampleObs=nObsStack, 
                             outOfSampleObs=nrow(database)-nObsStack, 
                             RMSE), 
                       fileNameParams, row.names=FALSE)
      utils::write.csv(samples[,1:i], fileNameSample, row.names=FALSE)
      apollo_print("Estimation results written to file.")
    } else {
      # Report error but continue with next iteration
      apollo_print(paste0("Estimation failed in cycle ", i, "."), type="w")
      if(estimationRoutine=="bfgs") print(as.matrix(round(get("lastFuncParam", 
                                                              envir=globalenv()),4)))
    }
    
  }
  
  # Stop clock
  endtime   <- Sys.time()
  timeTaken <- difftime(endtime, starttime, units='auto')
  apollo_print(paste0("Processing time: ", format(timeTaken)))
  apollo_print("\n")
  
  
  avgObsLL_est = llInSampleStack[,"inSample_model"]/nObsStack
  avgLLObs_val = llOutOfSampleStack[,"outOfSample_model"]/(nrow(database)-nObsStack)
  percentDiff  = 100*(avgLLObs_val-avgObsLL_est)/avgObsLL_est
  M <- cbind(avgObsLL_est, avgLLObs_val, percentDiff)
  rownames(M)=paste0("sample_",1:nrow(M))
  M <- rbind(M, Average = colMeans(M))
  tmp <- colnames(M)
  colnames(M) <- c("LL per obs in estimation sample", 
                   "LL per obs in validation sample", 
                   "% difference")
  apollo_print("Summary of cross-validation:")
  apollo_print(M)
  colnames(M) <- tmp
  
  return(invisible(M))
}
