#' Estimates model
#'
#' Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
#'
#' This is the main function of the Apollo package. The estimation process begins by running a number of checks on the 
#' \code{apollo_probabilities} function provided by the user.
#' If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
#' iterations. If bayesian estimation is used, estimation will finish once the predefined number of iterations are completed.
#' By default, this functions writes the estimated parameter values in each iteration to a file in the working directory. Writing 
#' can be turned off by setting \code{estimate_settings$writeIter} to \code{FALSE}, of by using any estimation algorithm
#' other than BFGS.
#' By default, \strong{final results are not written into a file nor printed into the console}, so users must make sure 
#' to call function \code{apollo_modelOutput} and/or \code{apollo_saveOutput} afterwards.
#' Users are strongly encouraged to visit www.apolloChoiceModelling.com to download examples on how to use the Apollo package.
#' The webpage also provides a detailed manual for the package, as well as a user-group to get further help.
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
#' @param estimate_settings List. Options controlling the estimation process.
#'                                 \itemize{
#'                                   \item \strong{estimationRoutine}: Character. Estimation method. Can take values "bfgs", "bhhh", or "nr".
#'                                                                     Used only if \code{apollo_control$HB} is FALSE. Default is "bfgs".
#'                                   \item \strong{maxIterations}: Numeric. Maximum number of iterations of the estimation routine before stopping.
#'                                                                 Used only if \code{apollo_control$HB} is FALSE. Default is 200.
#'                                   \item \strong{writeIter}: Boolean. Writes value of the parameters in each iteration to a csv file. 
#'                                                             Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
#'                                   \item \strong{hessianRoutine}: Character. Name of routine used to calculate the Hessian of the loglikelihood 
#'                                                                  function after estimation. Valid values are \code{"numDeriv"} (default) and 
#'                                                                  \code{"maxLik"} to use the routines in those packages, and \code{"none"} to avoid 
#'                                                                  estimating the Hessian (and the covariance matrix). Only used if \code{apollo_control$HB=FALSE}.
#'                                   \item \strong{printLevel}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. 
#'                                                              Ignored if apollo_control$HB is TRUE. Default is 3.
#'                                   \item \strong{constraints}: Constraints on parameters to estimate. Should ignore fixed parameters. 
#'                                                               See argument \code{constraints} in \link[maxLik]{maxBFGS} for more details.
#'                                   \item \strong{scaling}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. 
#'                                                           If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but 
#'                                                           \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter
#'                                                           b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, 
#'                                                           specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
#'                                  \item \strong{numDeriv_settings}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. 
#'                                                                    See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                                  \item \strong{bootstrapSE}: Numeric. Number of bootstrap samples to calculate standard errors. Default is 0, meaning
#'                                                              no bootstrap s.e. will be calculated. Number must zero or a positive integer. Only used
#'                                                              if \code{apollo_control$HB} is \code{FALSE}.
#'                                  \item \strong{bootstrapSeed}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples.
#'                                                                Only used if \code{bootstrapSE>0}. Default is 24.
#'                                  \item \strong{silent}: Boolean. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                 }
#' @return model object
#' @export
#' @importFrom numDeriv hessian grad
#' @importFrom sandwich sandwich
#' @importFrom maxLik maxLik
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov runif
apollo_estimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
  
  ### Detach things if necessary
  apollo_detach()
    
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="numDeriv", printLevel=3L, constraints=NULL, 
                  numDeriv_settings=list(), scaling=NA, bootstrapSE=0, bootstrapSeed=24, silent=FALSE)
  if(length(estimate_settings)==1 && is.na(estimate_settings)) estimate_settings <- default
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  
  ### Extract variables from pollo_input
  #database          = apollo_inputs[["database"]]
  apollo_control    = apollo_inputs[["apollo_control"]]
  #draws             = apollo_inputs[["draws"]]
  apollo_randCoeff  = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars     = apollo_inputs[["apollo_lcPars"]]
  apollo_HB         = apollo_inputs[["apollo_HB"]]
  workInLogs        = apollo_control$workInLogs
  estimationRoutine = tolower( estimate_settings[["estimationRoutine"]] )
  maxIterations     = estimate_settings[["maxIterations"]]
  writeIter         = estimate_settings[["writeIter"]]
  hessianRoutine    = estimate_settings[["hessianRoutine"]]
  printLevel        = estimate_settings[["printLevel"]]
  silent            = estimate_settings[["silent"]]
  numDeriv_settings = estimate_settings[["numDeriv_settings"]]
  constraints       = estimate_settings[["constraints"]]
  scaling           = estimate_settings[["scaling"]]
  bootstrapSE       = estimate_settings[["bootstrapSE"]]
  bootstrapSeed     = estimate_settings[["bootstrapSeed"]]
  
  # ################################## #
  #### initial processing & testing ####
  # ################################## #
  
  ### Validation of input
  apollo_checkArguments(apollo_probabilities, apollo_randCoeff, apollo_lcPars)
  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh' or 'nr'.")
  if( !(hessianRoutine %in% c("numDeriv", "maxLik", "none")) ) stop("Invalid hessianRoutine. Use 'numDeriv', 'maxLik' or 'none'.")
  
  # Check apollo_beta and apollo_fixed
  if(!is.numeric(apollo_beta) | !is.vector(apollo_beta) | is.null(names(apollo_beta))) stop("The \"apollo_beta\" argument needs to be a named vector")
  if(length(apollo_fixed)>0 && !is.character(apollo_fixed)) stop("'apollo_fixed' is not an empty vector nor a vector of names.")
  if(length(unique(names(apollo_beta)))<length(apollo_beta)) stop("The \"apollo_beta\" argument contains duplicate elements")
  if(length(unique(apollo_fixed))<length(apollo_fixed)) stop("The \"apollo_fixed\" argument contains duplicate elements")
  if(!all(apollo_fixed %in% names(apollo_beta))) stop("Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'.")
  
  # Check bootstrap settings
  if(!is.numeric(bootstrapSE) || length(bootstrapSE)!=1 || bootstrapSE<0) stop("'bootstrapSE' is not zero or a positive integer.")
  bootstrapSE <- as.integer(bootstrapSE); estimate_settings$bootstrapSE <- bootstrapSE
  if(!is.numeric(bootstrapSeed) || length(bootstrapSeed)!=1 || bootstrapSeed<=0) stop("'bootstrapSeed' is not a positive integer.")
  bootstrapSeed <- as.integer(bootstrapSeed); estimate_settings$bootstrapSeed <- bootstrapSeed
  
  # Check printLevel
  if(!is.integer(printLevel)) printLevel <- as.integer(round(printLevel,0))
  if(printLevel<0L){ printLevel <- 0L; estimate_settings$printLevel <- 0L }
  if(3L<printLevel){ printLevel <- 3L; estimate_settings$printLevel <- 3L }
  
  # Check maxIterations
  if(maxIterations<0) stop("Need positive number of iterations!")
  maxIterations     = round(maxIterations,0)
  estimate_settings$maxIterations = maxIterations
  
  # Check constraints
  if(!is.null(constraints) && apollo_control$HB) stop("Constraints cannot be used with Bayesian estimation.")
  if(!is.null(constraints) && estimationRoutine!="bfgs"){
    estimationRoutine="bfgs"
    apollo_inputs$estimationRoutine = estimationRoutine
    if(!silent) cat("WARNING: Estimation routine changed to 'BFGS'.\n         Only 'BFGS' supports constrained optimization.\n")
  }
  
  # Check writeIter
  if(estimationRoutine!="bfgs" & writeIter==TRUE){
    writeIter = FALSE
    estimate_settings$writeIter = writeIter
    txt <- "witeIter set to FALSE. Writing parameters values at each iteration is only available for BFGS estimation method."
    if(!silent) cat(txt, "\n", sep="") else warning(txt)
    rm(txt)
  }
  
  # create temporary copy of starting values for use later
  temp_start=apollo_beta
  
  # Checks for scaling
  if(length(scaling)>0 && !is.na(scaling)){
    if(any(!(names(scaling) %in% names(apollo_beta)))) stop("Some parameters included in 'scaling' are not included in 'apollo_beta'")
    #scaling <- scaling[!(names(scaling) %in% apollo_fixed)] # remove fixed parameters from scaling
    if(any((names(scaling) %in% apollo_fixed))) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
    if(any(scaling<0)){
      scaling <- abs(scaling)
      txt <- "WARNING: Some values in 'scaling' were negative,\n         they were replaced by their absolute value.\n"
      if(!silent) cat(txt) else warning("Some negative values in 'scaling' were replaced by their absotulte value")
    }
    if(any(!(scaling>0))) stop("All terms in in 'scaling' should be strictly positive!")
    txt <- "During estimation, parameters will be scaled\n using the values in estimate_settings$scaling"
    if(!silent) cat(txt,"\n", sep="") else warning(txt)
    # Copy scaling into apollo_inputs
    apollo_inputs$scaling = scaling
  }
  
  ### Call apollo_estimateHB
  if(apollo_control$HB){
    model <- apollo_estimateHB(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings)
    return(model)
  }
  
  # Scale starting parameters
  if(length(scaling)>0 && !is.na(scaling)){
    r <- names(apollo_beta) %in% names(apollo_inputs$scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- 1/apollo_inputs$scaling[r]*apollo_beta[r]
  }
  
  ### Start clock & apolloLog
  starttime <- Sys.time()
  apollo_inputs$apolloLog <- new.env(parent=emptyenv())
  
  ### Create useful variables
  nObs  <- nrow(apollo_inputs$database)
  indiv <- unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  nObsPerIndiv <- table(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  nObsPerIndiv <- setNames(as.vector(nObsPerIndiv), names(nObsPerIndiv))
  nObsPerIndiv <- nObsPerIndiv[nObsPerIndiv>0]
  
  ### Calculate analytical gradient if required
  if(!is.null(apollo_control$analyticGrad) && apollo_control$analyticGrad){
    g    <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="gradient")
    grad <- function(b, countIter, writeIter, sumLL, getNIter){
      b <- c(b, temp_start[apollo_fixed])
      b <- b[names(temp_start)]
      db <- get("apollo_inputs", envir=parent.frame(), inherits=TRUE)$database # maybe from .GlobalEnv to make it work with multicore?
      return( g(b, db) )
    }
  } else grad = NULL
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  
  if(!apollo_control$noValidation && !apollo_control$HB){
    if(!silent) cat("Validating inputs of likelihood function (apollo_probabilities)\n")
    apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")
    if(!silent & !apollo_control$noDiagnostics) cat("\n", apollo_printLog(apollo_inputs$apolloLog), sep="")
  }
  
  ### Check there are no references to database inside apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- gsub("apollo_inputs$database", " ", tmp, fixed=TRUE)
    tmp <- grep("database", tmp, fixed=TRUE)
    if(length(tmp)>0) stop("The database object is 'attached' and elements should thus be called directly in apollo_probabilities without the 'database$' prefix.")
  }
  
  ### Create multi-core version of likelihood (if needed)
  if(!apollo_control$HB){
    if(apollo_control$nCores==1) cl <- NA else {
      # Creates cluster and also deletes database and draws from apollo_inputs in here and in .GlobalEnv
      cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent, cleanMemory=TRUE)
      apollo_control$nCores <- length(cl)
      apollo_inputs$apollo_control$nCores <- length(cl)
    }
    on.exit(if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl))
    apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                         apollo_probabilities, apollo_inputs, 
                                         estimate_settings, cl=cl)
  }
  
  # Get starting LL
  beta0 <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  testLL = apollo_logLike(beta0)
  
  if(!apollo_control$noValidation){
    if(!silent) cat("\nTesting likelihood function ")
    # Check there are no invalid values
    if(any(!is.finite(testLL))){
      cat("\nLog-likelihood calculation fails at starting values!\n")
      cat("Affected individuals:\n")
      LLtemp <- data.frame(ID=indiv, LL=testLL)
      LLtemp = subset(LLtemp,!is.finite(LLtemp[,2]))
      colnames(LLtemp) <- c("ID","LL")
      print(LLtemp, row.names=FALSE)
      stop('Not proceeding with estimation!')
    } else {if(!silent) cat(".")}
    ### Test for unused parameters
    beta1   = beta0 + 0.001*runif(length(beta0))
    base_LL = apollo_logLike(beta1)
    #if(!workInLogs) base_LL = log(base_LL)
    if(any(!is.finite(base_LL))){
      cat("\nDuring testing, Apollo added disturbances smaller than 0.001 to all starting values\n")
      cat("This led to a log-likelihood calculation failure!\n")
      cat("Affected individuals:\n")
      LLtemp <- data.frame(ID=indiv, LL=base_LL)
      LLtemp = subset(LLtemp,!is.finite(LLtemp[,2]))
      colnames(LLtemp) <- c("ID","LL")
      print(LLtemp, row.names=FALSE)
      stop('Not proceeding with estimation!')
    } else {if(!silent) cat(".")}
    base_LL = sum(base_LL)
    for(p in names(beta0)){
      beta1p <- beta1 - (names(beta1)==p)*0.001
      beta1m <- beta1 + (names(beta1)==p)*0.001
      test1_LL = apollo_logLike(beta1p)
      test2_LL = apollo_logLike(beta1m)
      #if(!workInLogs){test1_LL <- log(test1_LL); test2_LL <- log(test2_LL)}
      test1_LL = sum(test1_LL)
      test2_LL = sum(test2_LL)
      if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      if(base_LL==test1_LL & base_LL==test2_LL) stop("Parameter ",p," does not influence the log-likelihood of your model!")
      if(!silent) cat(".")
    }
    if(!silent) cat("Done.\n")
    rm(beta1, beta1p, beta1m, base_LL, test1_LL, test2_LL)
  }
  rm(beta0)
  
  # ################################## #
  #### classical estimation         ####
  # ################################## #
  
  ### Split parameters between variable and fixed
  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  beta_fix_val <- apollo_beta[apollo_fixed]
  
  ### Main (classical) estimation
  if(exists("lastFuncParam")){
    tmp <- globalenv()
    assign("lastFuncParam", rep(0, length(beta_var_val)), envir=tmp)
    rm(tmp)
  }
  if(!silent) cat("\n\nStarting main estimation\n") else printLevel=0
  model <- maxLik::maxLik(apollo_logLike, grad=grad, start=beta_var_val,
                          method=estimationRoutine, finalHessian=FALSE,
                          control=list(printLevel=printLevel, iterlim=maxIterations),
                          constraints=constraints,
                          countIter=TRUE, writeIter=writeIter, sumLL=FALSE, getNIter=FALSE)
  
  # Checks if main estimation was successful
  succesfulEstimation <- FALSE
  if(exists("model")){
    if(estimationRoutine=="bfgs" & model$code==0) succesfulEstimation <- TRUE
    if(estimationRoutine=="bhhh" & (model$code %in% c(2,8)) ) succesfulEstimation <- TRUE
    if(estimationRoutine=="nr" && model$code<=2) succesfulEstimation <- TRUE
  }
  
  # Print estimated parameters
  if(exists("model")& !silent){
    cat("\nEstimated parameters:\n")
    tmp <- c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
    if(!is.null(apollo_inputs$scaling)) for(i in names(apollo_inputs$scaling)){
      tmp[i] <- apollo_inputs$scaling[i]*tmp[i]
    }
    print(as.matrix(round(tmp,4)))
    rm(tmp)
    cat("\n")
  }
  
  # If estimation failed, return whatever can be salvaged
  if(!succesfulEstimation){
    txt <- "ERROR: Estimation failed. No covariance matrix to compute."
    if(!silent) cat(txt,"\n", sep="") else warning(txt)
    if(exists("model")){
      if(length(scaling)>0 && !is.na(scaling)) model$scaling <- scaling
      return(model)
    } else stop("Sorry, no estimated model to return.\n")
  }
  
  ### Calculation of the Hessian
  if(hessianRoutine!="none"){
    success_nd <- FALSE
    nNA_nd <- -1
    
    # If Hessian is to be calculated with numDeriv
    if(hessianRoutine=="numDeriv"){
      if(!silent) cat("Computing covariance matrix using numDeriv package.\n (this may take a while)\n")
      # Create closure to keep track of progress and estimate hessian using numDeriv
      sumLogLike <- function(k, silent){
        i <- 0
        I <- 2+8*( k*(k+1)/2 )
        step <- ceiling(I/20)
        
        function(theta){
          if(i==0 & !silent) cat(' 0%')
          tmp <- apollo_logLike(theta, countIter=FALSE, writeIter=FALSE, sumLL=TRUE)
          i <<- i+1
          if(i%%step==0 & !silent){
            if(i%%(5*step)==0) cat(i/(5*step)*25,'%',sep='') else cat('.')
          }
          if(i==I & !silent) cat('100%\n')
          return(tmp)
        }
      }
      # Estimate hessian using numDeriv
      H <- tryCatch(numDeriv::hessian(func=sumLogLike(length(model$estimate), silent),
                                      x=model$estimate, method.args=numDeriv_settings),
                    error = function(e) return(NA))
      # Check if estimation was succesful
      if(length(H)==1 && anyNA(H)){
        if(!silent) cat(" ERROR: Hessian calculation using numDeriv failed.\n")
        hessianRoutine <- "maxLik"
      } else {
        success_nd <- TRUE
        nNA_nd <- sum(is.na(H))
        if(nNA_nd>0 & !silent) cat(" Some (",nNA_nd,") NA values found in numDeriv Hessian.\n", sep="")
        if(success_nd && nNA_nd==0 && !silent) cat(" Hessian calculated with numDeriv will be used.\n")
      }
    }
    
    # If Hessian is to be calculated with maxLik
    if(hessianRoutine=="maxLik" | nNA_nd>0){
      success_ml <- FALSE
      nNA_ml <- -1
      if(!silent) cat("Computing covariance matrix using maxLik package.\n (this may take a while, no progress bar displayed)\n")
      # Estimate Hessian using maxLik
      model2 <- tryCatch(maxLik::maxLik(apollo_logLike, start=model$estimate, method=estimationRoutine, print.level=0,
                                        finalHessian=TRUE, iterlim=2, countIter=FALSE, writeIter=FALSE, sumLL=FALSE),
                         error=function(e) return(NA))
      # Check the hessian was correctly estimated with maxLik
      if(length(model2)==1 && anyNA(model2)){ if(!silent) cat(" ERROR: Hessian calculation using maxLik failed.\n") } else {
        success_ml <- TRUE
        nNA_ml <- sum(is.na(model2$hessian))
        if(nNA_ml>0 & !silent) cat(" Some (",nNA_ml,") NA values found in maxLik Hessian.\n", sep="")
        if(hessianRoutine=="maxLik" | nNA_ml<nNA_nd){
          H <- model2$hessian
          if(!silent) cat(" Hessian calculated with maxLik will be used.\n")
        }
      }
    }
    
    if(success_nd || success_ml){
      rownames(H) <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
      colnames(H) <- rownames(H)
    } else H <- NULL
  } else H <- NULL
  
  ### Copy Hessian to model, and checks if s.e. cannot be calculated
  model$hessian <- H
  if(is.null(model$hessian) & hessianRoutine!="none"){
    if(!silent) cat(" ERROR: Hessian could not be calculated.\n Postprocessing aborted.\n")
  }
  if(!is.matrix(try(solve(model$hessian),silent=T))){
    if(!silent) cat(' ERROR: Singular Hessian, cannot calculate s.e.\n')
    tryCatch({
      colnames(model$hessian) <- names(beta_var_val)
      rownames(model$hessian) <- names(beta_var_val)
      utils::write.csv(model$hessian, paste(apollo_control$modelName, "hessian.csv", sep="_"))
      if(!silent) cat(" Hessian written to", paste(apollo_control$modelName, "hessian.csv", sep="_"), "\n")
    }, error=function(e) if(!silent) cat(" Could not write hessian to a file.\n"))
  }
  eigValue <- tryCatch(eigen(H)$values, error=function(e) return(NA))
  if(!anyNA(eigValue)){
    eigValue <- eigValue[which.min(abs(eigValue))]
    model$eigenpos=any(eigValue>0)
    if(model$eigenpos) cat(" Warning: Some eigenvalues of hessian are positive,\n          indicating potential problems!\n") 
    model$eigen=round(eigValue, 6)
    cat(" Min abs eigenvalue of hessian    : ",abs(model$eigen),"\n") 
  }
  
  
  ### Calculate s.e.
  dummyVCM           <- matrix(NA, nrow=length(model$estimate), ncol=length(model$estimate))
  rownames(dummyVCM) <- names(model$estimate)
  colnames(dummyVCM) <- names(model$estimate)
  model$varcov     <- tryCatch(stats::vcov(model), error=function(e) return(dummyVCM))
  if(is.null(model$varcov)) model$varcov <- dummyVCM
  model$robvarcov  <- tryCatch(sandwich::sandwich(model), error=function(e) return(dummyVCM))
  model$se         <- sqrt(diag(model$varcov))
  model$robse      <- sqrt(diag(model$robvarcov))
  model$corrmat    <- tryCatch(model$varcov/(model$se%*%t(model$se)), error=function(e) return(dummyVCM))
  model$robcorrmat <- tryCatch(model$robvarcov/(model$robse%*%t(model$robse)) , error=function(e) return(dummyVCM))
  
  ### Calculate probabilities to identify outliers
  P <- exp(apollo_logLike(model$estimate))
  if(apollo_control$panelData){
    model$avgCP <- setNames(P^(1/nObsPerIndiv), names(nObsPerIndiv))
  } else model$avgCP <- setNames(P, 1:length(P))
  
  #### Calculate gradient
  #if(!silent) cat("Calculating gradient norm... ")
  #model$grad <- tryCatch(numDeriv::grad(apollo_logLike, model$estimate, countIter=FALSE, writeIter=FALSE, sumLL=TRUE),
  #                       error=function(e) return(NA))
  #if(!silent) if(length(model$grad)==1 && is.na(model$grad)) cat("ERROR\n") else cat(round( sqrt(sum(model$grad^2)) ,4), "\n", sep="")
  
  ### Stores number of iterations
  model$nIter <- ifelse(estimationRoutine=="bfgs", apollo_logLike(NA, getNIter=TRUE), model$iterations)
  
  ### Delete apollo_logLike and close cluster to free resources
  rm(apollo_logLike)
  if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl)
  rm(cl)
  
  ### Restore apollo_inputs including database and draws
  if(apollo_control$nCores>1){
    if(!silent) cat('Restoring data in main thread...')
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs")
    apollo_inputs <- readRDS(file=fileName)
    tmp <- .GlobalEnv
    if(exists("apollo_inputs", envir=tmp)) assign("apollo_inputs", apollo_inputs, envir=tmp)
    unlink(fileName); rm(fileName)
    if(!silent) cat(' Done. ',sum(gc()[,2]),'MB of RAM in use\n',sep='')
  }
  
  ### Calculate bootstrap s.e.
  if(bootstrapSE>0){
    if(!silent) cat("Starting bootstrap calculation of standard errors.")
    tmp <- list(estimationRoutine=estimationRoutine, maxIterations=maxIterations,
                writeIter=FALSE, hessianRoutine="none", printLevel=printLevel,
                silent=silent)
    model$bootvarcov <- apollo_bootstrap(apollo_beta, apollo_fixed,
                                         apollo_probabilities, apollo_inputs,
                                         estimate_settings=tmp,
                                         bootstrap_settings=list(nRep=bootstrapSE,
                                                                  seed=bootstrapSeed))
    model$bootse <- sqrt(diag(model$bootvarcov))
    model$bootse[apollo_fixed] <- NA
    model$bootcorrmat <- tryCatch(model$bootvarcov/(model$bootse%*%t(model$bootse)), error=function(e) return(dummyVCM))
  }
  
  ### Calculate Zero LL
  if(!silent) cat("Calculating LL(0)... ")
  model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
                        error=function(e) return(NA))
  if(is.list(model$LL0)) model$LL0=model$LL0[c("model",names(model$LL0)[names(model$LL0)!="model"])]
  if(is.list(model$LL0)){
    if(!workInLogs) model$LL0=sapply(model$LL0,function(x) sum(log(x)))
    if(workInLogs) model$LL0=sapply(model$LL0,sum)
  } else {
    model$LL0 <- ifelse( workInLogs, sum(model$LL0), sum(log(model$LL0)) )
  }
  if(!silent) cat("Done.\n")
  
  
  ### Restore fixed parameters to model$estimate and harmonise with model$se, model$robse
  temp           = c(model$estimate, apollo_beta[apollo_fixed])
  model$estimate = temp[names(apollo_beta)]
  temp           = c(model$se, apollo_beta[apollo_fixed])
  model$se       = temp[names(apollo_beta)]
  temp           = c(model$robse, apollo_beta[apollo_fixed])
  model$robse    = temp[names(apollo_beta)]
  model$se[apollo_fixed]    = NA
  model$robse[apollo_fixed] = NA
  
  ### Get LL at optimum for each component (functionality="output" is used inside apollo_llCalc)
  apollo_inputs$apolloLog <- new.env(parent=emptyenv()) # Log is emptied
  if(!silent) cat("Calculating LL of each model component... ")
  LLout <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="output"),
                   error=function(e) return(NA))
  
  if(!anyNA(LLout) && is.list(LLout)){
    # Give name to unnamed components (already done in apollo_prepareProb)
    #origNames <- names(Pout)
    #newNames  <- paste0("component_", 1:length(Pout))
    #if(!is.null(origNames)) newNames <- ifelse(origNames!="", origNames, newNames)
    #names(Pout) <- newNames
    # Get log of likelihood with "model" first
    #tmp <- c("model", newNames[newNames!="model"])
    LLout <- LLout[c("model",names(LLout)[names(LLout)!="model"])]
    if(!workInLogs) LLout <- lapply(LLout, log)
    LLout <-lapply(LLout,sum)
    if(!silent) cat("Done.\n")
    model$Pout  <- LLout
    #model$LLout <- unlist(lapply(model$Pout, sum))
    model$LLout <- unlist(LLout)
  } else{
    model$Pout  <- LLout
    model$LLout <- list(NA)
    if(!silent) cat("Not applicable to all components.\n")
  }
  
  if(!silent) cat("\n", apollo_printLog(apollo_inputs$apolloLog, book=2), sep="")
  
  ### use pre-scaling values as starting values in output
  model$bootstrapSE       <- bootstrapSE
  model$apollo_beta       <- temp_start
  model$LLStart           <- sum(testLL)
  model$startTime         <- starttime
  model$apollo_control    <- apollo_control
  model$nObs              <- nObs
  model$nIndivs           <- length(indiv)
  model$apollo_draws      <- apollo_inputs$apollo_draws
  model$apollo_randCoeff  <- apollo_randCoeff
  model$apollo_lcPars     <- apollo_lcPars
  model$apollo_fixed      <- apollo_fixed
  model$estimationRoutine <- estimationRoutine
  model$apolloLog         <- apollo_inputs$apolloLog
  
  ### Apply scaling to estimates and covariance matrix
  if(length(scaling)>0 && !is.na(scaling)) model$scaling <- scaling
  if(!is.null(apollo_inputs$scaling)) for(i in names(apollo_inputs$scaling)){
    s <- apollo_inputs$scaling[i]
    model$estimate[i] <- s*model$estimate[i]
    model$apollo_beta[i] <- s*model$apollo_beta[i]
    model$varcov[i,] <- s*model$varcov[i,]
    model$varcov[,i] <- s*model$varcov[,i]
    #if(i %in% names(model$gradient)) model$gradient[i] <- s*model$gradient[i]
    model$robvarcov[i,] <- s*model$robvarcov[i,]
    model$robvarcov[,i] <- s*model$robvarcov[,i]
    model$se[i] <- s*model$se[i]
    model$robse[i] <- s*model$robse[i]
  }
  
  ### Stop clock
  endtime                 <- Sys.time()
  model$timeTaken         <- as.numeric(difftime(endtime,starttime,units='secs'))
  
  return(model)
}
