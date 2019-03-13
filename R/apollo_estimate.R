#' Estimates model
#'
#' Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
#'
#' This is the main function of the Apollo package. The estimation process begins by checking the definition of
#' \code{apollo_probabilities} by estimating it at the starting values. Then it runs the function with argument \code{functionality="validate"}.
#' If the user requested more than one core for estimation (i.e. \code{apollo_control$nCores>1}), and no bayesian estimation is used
#' (i.e. \code{apollo_control$HB=FALSE}), then a cluster is created. Using a cluster at least doubles the requires RAM, as the database
#' must be copied into the cluster.
#' If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
#' iterations. If bayesian estimation is used, estimation will finish once the predefined number of iterations are completed.
#' This functions does not save results into a file nor prints them into the console, so if users want to see and store estimation the results,
#' they must make sure to call function \code{apollo_modelOutput} and/or \code{apollo_saveOutput} afterwards.
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
#'                                   \item estimationRoutine: Character. Estimation method. Can take values "bfgs", "bhhh", or "nr".
#'                                                            Used only if \code{apollo_control$HB} is FALSE. Default is "bfgs".
#'                                  \item maxIterations: Numeric. Maximum number of iterations of the estimation routine before stopping.
#'                                                       Used only if \code{apollo_control$HB} is FALSE. Default is 200.
#'                                  \item writeIter: Boolean. Writes value of the parameters in each iteration to a csv file. Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
#'                                  \item hessianRoutine: Character. Name of routine used to calculate the Hessian of the loglikelihood function after estimation. Valid values are \code{"numDeriv"} (default) and \code{"maxLik"} to use the routines in those packages, and \code{"none"} to avoid estimating the Hessian (and the covariance matrix). Only used if \code{apollo_control$HB=FALSE}.
#'                                  \item printLevel: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. Ignored if apollo_control$HB is TRUE. Default is 3.
#'                                  \item numDeriv_settings: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                                  \item silent: Boolean. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                 }
#' @return model object
#' @export
#' @importFrom numDeriv hessian grad
#' @importFrom sandwich sandwich
#' @importFrom RSGHB doHB
#' @importFrom maxLik maxLik
apollo_estimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){

  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, hessianRoutine="numDeriv", printLevel=3L, numDeriv_settings=list(), silent=FALSE)
  if(length(estimate_settings)==1 && is.na(estimate_settings)) estimate_settings <- default
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] 
  for(i in tmp) estimate_settings[[i]] <- default[[i]]

  database          = apollo_inputs[["database"]]
  apollo_control    = apollo_inputs[["apollo_control"]]
  draws             = apollo_inputs[["draws"]]
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


  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh' or 'nr'.")
  if( !(hessianRoutine %in% c("numDeriv", "maxLik", "none")) ) stop("Invalid hessianRoutine. Use 'numDeriv', 'maxLik' or 'none'.")
  if((length(apollo_fixed)>0) & any(!(apollo_fixed %in% names(apollo_beta)))) stop("Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'")
  if(!is.integer(printLevel)) printLevel <- as.integer(round(printLevel,0))
  if(printLevel<0L) printLevel <- 0L
  if(3L<printLevel) printLevel <- 3L
  if(maxIterations<0) stop("Need positive number of iterations!")
  maxIterations=round(maxIterations,0)
  if(estimationRoutine!="bfgs" & writeIter==TRUE){
    writeIter = FALSE
    estimate_settings[["writeIter"]] = FALSE
    apollo_inputs$apollo_estimate$writeIter <- FALSE
    cat("witeIter set to FALSE. Writing parameters values at each iteration is only available for BFGS estimation method.\n")
  }

  tempOutputFile <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
  tempOutputFile <- file.path(tempdir(),tempOutputFile)
  if(file.exists(tempOutputFile)){
    tryCatch( file.remove(tempOutputFile),
              error=function(e) cat("Could not delete old file ",tempOutputFile,".\n", sep=""))
  }

  starttime <- Sys.time()

  if(!silent & !apollo_control$noValidation) cat("Testing probability function (apollo_probabilities)\n")
  apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")
  testLL = apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate")
  if(anyNA(testLL)) stop('Log-likelihood calculation fails at starting values!')

  if(apollo_control$HB){

    tmp <- tryCatch( get("apollo_fixed", envir=globalenv()), error=function(e) 1 )
    if( length(tmp)>0 && any(tmp %in% c(apollo_HB$gVarNamesFixed, apollo_HB$gVarNamesFixed)) ) stop("apollo_fixed seems to have changed since calling apollo_inputs.")

    gFix <- apollo_HB$gVarNamesFixed
    gNor <- apollo_HB$gVarNamesNormal
    apollo_HB_likelihood=function(fc,b){
      if(is.null(gFix)) fc1 <- NULL else fc1 <- stats::setNames(as.list(fc)     , gFix)
      if(is.null(gNor)) b1  <- NULL else b1  <- stats::setNames(as.data.frame(b), gNor)
      if(length(apollo_fixed)==0) fp <- NULL else fp  <-  stats::setNames( as.list(apollo_beta[apollo_fixed]), apollo_fixed )
      P <- apollo_probabilities(apollo_beta=c(fc1,b1,fp), apollo_inputs, functionality="estimate")
      return(P)
    }

    model <- RSGHB::doHB(apollo_HB_likelihood, database, apollo_HB)

    model$apollo_HB   <- apollo_HB
    model$apollo_beta <- apollo_beta
    if(workInLogs) model$LLStart <- sum(testLL) else model$LLStart <- sum(log(testLL))
    model$LL0         <- sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL")))
    model$startTime   <- starttime
    model$apollo_control <- apollo_control
    model$nObs        <- nrow(database)
    model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
    endtime           <- Sys.time()
    model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
    model$apollo_fixed <- apollo_fixed
    model$estimationRoutine <- "Hierarchical Bayes"

    return(model)
  }

  if(apollo_control$nCores==1) cl <- NA else {
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent)
    apollo_control$nCores <- length(cl)
    apollo_inputs$apollo_control$nCores <- length(cl)
  }
  on.exit(if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl))
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities,
                                       apollo_inputs, estimate_settings, cl=cl)

  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  beta_fix_val <- apollo_beta[apollo_fixed]


  if(exists("lastFuncParam")){
    tmp <- globalenv()
    assign("lastFuncParam", rep(0, length(beta_var_val)), envir=tmp)
	rm(tmp)
  }
  if(!silent) cat("\n\nStarting main estimation\n") else printLevel=0
  model <- maxLik::maxLik(apollo_logLike, start=beta_var_val,
                          method=estimationRoutine, finalHessian=FALSE,
                          control=list(printLevel=printLevel, iterlim=maxIterations),
                          countIter=TRUE, writeIter=writeIter, sumLL=FALSE, getNIter=FALSE)

  succesfulEstimation <- FALSE
  if(exists("model")){
    if(estimationRoutine=="bfgs" & model$code==0) succesfulEstimation <- TRUE
    if(estimationRoutine=="bhhh" & (model$code %in% c(2,8)) ) succesfulEstimation <- TRUE
    if(estimationRoutine=="nr" && model$code<=2) succesfulEstimation <- TRUE
  }

  if(exists("model")& !silent){
    cat("\nEstimated values:\n")
    tmp <- c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
    print(as.matrix(round(tmp,4)))
    rm(tmp)
    cat("\n")
  }

  if(!succesfulEstimation){
    cat("ERROR: Estimation failed. No covariance matrix to compute.\n")
    if(exists("model")){
      print(as.matrix(model$estimate, ncol=1))
      return(model)
    } else stop("Sorry, no estimated model to return.\n")
  }

  if(hessianRoutine!="none"){
    success_nd <- FALSE
    nNA_nd <- -1

    if(hessianRoutine=="numDeriv"){
      if(!silent) cat("Computing covariance matrix using numDeriv package.\n (this may take a while)\n")
      sumLogLike <- function(k, silent){
        i <- 0
        I <- 2+8*( k*(k+1)/2 )
        step <- ceiling(I/20)

        function(theta){
          if(i==0 & !silent) cat('0%')
          tmp <- apollo_logLike(theta, countIter=FALSE, writeIter=FALSE, sumLL=TRUE)
          i <<- i+1
          if(i%%step==0 & !silent){
            if(i%%(5*step)==0) cat(i/(5*step)*25,'%',sep='') else cat('.')
          }
          if(i==I & !silent) cat('100%\n')
          return(tmp)
        }
      }
      H <- tryCatch(numDeriv::hessian(func=sumLogLike(length(model$estimate), silent),
                                      x=model$estimate, method.args=numDeriv_settings),
                    error = function(e) return(NA))
      if(length(H)==1 && anyNA(H)){
        if(!silent) cat("ERROR: Hessian calculation using numDeriv failed.\n")
        hessianRoutine <- "maxLik"
      } else {
        success_nd <- TRUE
        nNA_nd <- sum(is.na(H))
        if(nNA_nd>0 & !silent) cat("Some (",nNA_nd,") NA values found in numDeriv Hessian.\n", sep="")
        if(success_nd && nNA_nd==0 && !silent) cat("Hessian calculated with numDeriv will be used.\n")
      }
    }

    if(hessianRoutine=="maxLik" | nNA_nd>0){
      success_ml <- FALSE
      nNA_ml <- -1
      if(!silent) cat("Computing covariance matrix using maxLik package.\n (this may take a while, no progress bar displayed)\n")
      model2 <- tryCatch(maxLik::maxLik(apollo_logLike, start=model$estimate, method=estimationRoutine, print.level=0,
                                        finalHessian=TRUE, iterlim=2, countIter=FALSE, writeIter=FALSE, sumLL=FALSE),
                         error=function(e) return(NA))
      if(length(model2)==1 && anyNA(model2)){ if(!silent) cat("ERROR: Hessian calculation using maxLik failed.\n") } else {
        success_ml <- TRUE
        nNA_ml <- sum(is.na(model2$hessian))
        if(nNA_ml>0 & !silent) cat("Some (",nNA_ml,") NA values found in maxLik Hessian.\n", sep="")
        if(hessianRoutine=="maxLik" | nNA_ml<nNA_nd){
          H <- model2$hessian
          if(!silent) cat("Hessian calculated with maxLik will be used.\n")
        }
      }
    }

    if(success_nd || success_ml){
      rownames(H) <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
      colnames(H) <- rownames(H)
    } else H <- NULL
  } else H <- NULL

  model$hessian <- H
  if(is.null(model$hessian) & hessianRoutine!="none"){
    if(!silent) cat("ERROR: Hessian could not be calculated. Postprocessing aborted.\n")
  }
  if(!is.matrix(try(solve(model$hessian),silent=T))){
    if(!silent) cat('ERROR: Singular Hessian, cannot calculate s.e.\n')
    tryCatch({
      colnames(model$hessian) <- names(beta_var_val)
      rownames(model$hessian) <- names(beta_var_val)
      utils::write.csv(model$hessian, paste(apollo_control$modelName, "hessian.csv", sep="_"))
      if(!silent) cat("Hessian written to", paste(apollo_control$modelName, "hessian.csv", sep="_"), "\n")
    }, error=function(e) if(!silent) cat("Could not write hessian to a file.\n"))
  }

  dummyVCM           <- matrix(NA, nrow=length(model$estimate), ncol=length(model$estimate))
  rownames(dummyVCM) <- names(model$estimate)
  colnames(dummyVCM) <- names(model$estimate)
  model$varcov    <- tryCatch(stats::vcov(model), error=function(e) return(dummyVCM))
  model$robvarcov <- tryCatch(sandwich::sandwich(model), error=function(e) return(dummyVCM))
  model$se        <- sqrt(diag(model$varcov))
  model$robse     <- sqrt(diag(model$robvarcov))
  model$corrmat    <- tryCatch(model$varcov/(model$se%*%t(model$se)), error=function(e) return(dummyVCM))
  model$robcorrmat <- tryCatch(model$robvarcov/(model$robse%*%t(model$robse)) , error=function(e) return(dummyVCM))


  P <- exp(apollo_logLike(model$estimate))
  if(apollo_control$panelData){
    nObsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
  } else nObsPerIndiv <- 1
  model$avgCP <- P^(1/nObsPerIndiv)
  names(model$avgCP) <- unique(database[,apollo_control$indivID])


  if(!silent) cat("Calculating LL(0)... ")
  model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
                        error=function(e) return(NA))
  if(!anyNA(model$LL0)){
    model$LL0 <- sum(log(model$LL0))
    if(!silent) cat(round(model$LL0,2), "\n", sep="")
  } else {
  if(!silent) cat("No LL0 for some components.\n")
  }


  temp           = c(model$estimate, apollo_beta[apollo_fixed])
  model$estimate = temp[names(apollo_beta)]
  temp           = c(model$se, apollo_beta[apollo_fixed])
  model$se       = temp[names(apollo_beta)]
  temp           = c(model$robse, apollo_beta[apollo_fixed])
  model$robse    = temp[names(apollo_beta)]
  model$se[apollo_fixed]    = NA
  model$robse[apollo_fixed] = NA

  model$Pout  <- apollo_llCalc(model$estimate, apollo_probabilities, apollo_inputs,silent)
  model$LLout <- unlist(lapply(model$Pout, sum))



  model$apollo_beta <- c(beta_var_val, beta_fix_val)[names(apollo_beta)]
  if(workInLogs) model$LLStart <- sum(testLL) else model$LLStart <- sum(log(testLL))
  model$startTime   <- starttime
  model$nIter       <- ifelse(estimationRoutine=="bfgs", apollo_logLike(NA, getNIter=TRUE), model$iterations)
  model$apollo_control <- apollo_control
  model$nObs        <- nrow(database)
  model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
  model$apollo_draws <- apollo_inputs$apollo_draws
  model$apollo_randCoeff<-apollo_randCoeff
  model$apollo_lcPars   <- apollo_lcPars
  endtime           <- Sys.time()
  model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
  model$apollo_fixed <- apollo_fixed
  model$estimationRoutine <- estimationRoutine

  return(model)
}
