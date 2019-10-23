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
#' @importFrom RSGHB doHB
#' @importFrom maxLik maxLik
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov
#' 
apollo_estimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
  
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="numDeriv", printLevel=3L, constraints=NULL, 
                  numDeriv_settings=list(), scaling=NA, bootstrapSE=0, bootstrapSeed=24, silent=FALSE)
  if(length(estimate_settings)==1 && is.na(estimate_settings)) estimate_settings <- default
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  
  ### Extract variables from pollo_input
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
  constraints       = estimate_settings[["constraints"]]
  scaling           = estimate_settings[["scaling"]]
  bootstrapSE       = estimate_settings[["bootstrapSE"]]
  bootstrapSeed     = estimate_settings[["bootstrapSeed"]]
  
  # ################################## #
  #### initial processing & testing ####
  # ################################## #
  
  ### Validation of input
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh' or 'nr'.")
  if( !(hessianRoutine %in% c("numDeriv", "maxLik", "none")) ) stop("Invalid hessianRoutine. Use 'numDeriv', 'maxLik' or 'none'.")
  if(!is.numeric(apollo_beta) | !is.vector(apollo_beta) | is.null(names(apollo_beta))) stop("The \"apollo_beta\" argument needs to be a named vector")
  if(length(apollo_fixed)>0 && !is.character(apollo_fixed)) stop("'apollo_fixed' is not an empty vector nor a vector of names.")
  if(length(unique(names(apollo_beta)))<length(apollo_beta)) stop("The \"apollo_beta\" argument contains duplicate elements")
  if(length(unique(apollo_fixed))<length(apollo_fixed)) stop("The \"apollo_fixed\" argument contains duplicate elements")
  if(!all(apollo_fixed %in% names(apollo_beta))) stop("Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'.")
  if(!is.numeric(bootstrapSE) || length(bootstrapSE)!=1 || bootstrapSE<0) stop("'bootstrapSE' is not zero or a positive integer.")
  bootstrapSE <- as.integer(bootstrapSE)
  if(!is.numeric(bootstrapSeed) || length(bootstrapSeed)!=1 || bootstrapSeed<=0) stop("'bootstrapSeed' is not a positive integer.")
  bootstrapSeed <- as.integer(bootstrapSeed)
  ### create temporary copy of starting values for use later
  temp_start=apollo_beta
  if(length(scaling)>0 && !is.na(scaling)){
    if(any(!(names(scaling) %in% names(apollo_beta)))) stop("Some parameters included in 'scaling' are not included in 'apollo_beta'")
    if(any((names(scaling) %in% apollo_fixed))) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
    if(any(!(scaling>0))) stop("All terms in in 'scaling' should be strictly positive!")
    txt <- "During estimation, parameters will be scaled using the\n values in estimate_settings$scaling"
    if(!silent) cat(txt,"\n", sep="") else warning(txt)
    apollo_inputs$scaling = scaling
    r <- names(apollo_beta) %in% names(apollo_inputs$scaling)
    r <- names(apollo_beta)[r]
    if(!apollo_control$HB){
      # classical version
      apollo_beta[r] <- 1/apollo_inputs$scaling[r]*apollo_beta[r]
    }else{
      # bayesian version
      if(!is.null(apollo_HB$gVarNamesFixed)){
        r <- ( names(apollo_beta) %in% names(apollo_inputs$scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
        r <- names(apollo_beta)[r]
        apollo_HB$FC[r] <- 1/apollo_inputs$scaling[r]*apollo_HB$FC[r]
      }
      if(!is.null(apollo_HB$gVarNamesNormal)){
        r <- ( names(apollo_beta) %in% names(apollo_inputs$scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
        r <- names(apollo_beta)[r]
        dists_normal=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==1])
        dists_lnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==2])
        dists_lnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==3])
        dists_cnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==4])
        dists_cnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==5])
        dists_sb=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==6])
        if(length(dists_normal)>0) apollo_HB$svN[dists_normal] <- 1/apollo_inputs$scaling[dists_normal]*apollo_HB$svN[dists_normal]
        if(length(dists_lnp)>0) apollo_HB$svN[dists_lnp] <- -log(apollo_inputs$scaling[dists_lnp])+apollo_HB$svN[dists_lnp]
        if(length(dists_lnn)>0) apollo_HB$svN[dists_lnn] <- -log(apollo_inputs$scaling[dists_lnn])+apollo_HB$svN[dists_lnn]
        if(length(dists_cnp)>0) apollo_HB$svN[dists_cnp] <- 1/apollo_inputs$scaling[dists_cnp]*apollo_HB$svN[dists_cnp]
        if(length(dists_cnn)>0) apollo_HB$svN[dists_cnn] <- 1/apollo_inputs$scaling[dists_cnn]*apollo_HB$svN[dists_cnn]
        if(length(dists_sb)>0){
          names(apollo_HB$gMINCOEF)=names(apollo_HB$svN)
          names(apollo_HB$gMAXCOEF)=names(apollo_HB$svN)
          apollo_HB$gMINCOEF[dists_sb] <- 1/apollo_inputs$scaling[dists_sb]*apollo_HB$gMINCOEF[dists_sb]
          apollo_HB$gMAXCOEF[dists_sb] <- 1/apollo_inputs$scaling[dists_sb]*apollo_HB$gMAXCOEF[dists_sb]
        }
      }
    }
  }
  if(!is.integer(printLevel)) printLevel <- as.integer(round(printLevel,0))
  if(printLevel<0L) printLevel <- 0L
  if(3L<printLevel) printLevel <- 3L
  if(maxIterations<0) stop("Need positive number of iterations!")
  maxIterations=round(maxIterations,0)
  if(!is.null(constraints) && apollo_control$HB) stop("Constraints cannot be used with Bayesian estimation.")
  if(!is.null(constraints) && estimationRoutine!="bfgs"){
    estimationRoutine="bfgs"
    warning("Estimation routine changed to 'BFGS'. Only 'BFGS' supports constrained optimization.")
  }
  if(estimationRoutine!="bfgs" & writeIter==TRUE){
    writeIter = FALSE
    estimate_settings[["writeIter"]] = FALSE
    apollo_inputs$apollo_estimate$writeIter <- FALSE
    txt <- "witeIter set to FALSE. Writing parameters values at each iteration is only available for BFGS estimation method."
    if(!silent) cat(txt, "\n", sep="") else warning(txt)
    rm(txt)
  }
  
  ### Start clock & apolloLog
  starttime <- Sys.time()
  apollo_inputs$apolloLog <- new.env(parent=emptyenv())
  
  ### Validate and test probability
  if(!silent & !apollo_control$noValidation) cat("Testing probability function (apollo_probabilities)\n")
  if(!apollo_control$HB){
    ### Validation for classical estimation
    apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")
    if(!silent & !apollo_control$noDiagnostics) cat("\n", apollo_printLog(apollo_inputs$apolloLog), sep="")
    testLL = apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate")
    if(!workInLogs) testLL=log(testLL)
    # Maybe here we could return the value of the likelihood and print and error wuth cat, instead of simply stopping
    if(anyNA(testLL)) stop('Log-likelihood calculation fails at starting values!')
  } else {
    ### Validation using HB estimation
    apollo_test_beta=apollo_beta
    if(!is.null(apollo_HB$gVarNamesFixed)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
      r <- names(apollo_beta)[r]
      apollo_test_beta[r] <- apollo_HB$FC[r]
    }
    if(!is.null(apollo_HB$gVarNamesNormal)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
      r <- names(apollo_beta)[r]
      dists_normal=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==1])
      dists_lnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==2])
      dists_lnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==3])
      dists_cnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==4])
      dists_cnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==5])
      dists_sb=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==6])
      if(length(dists_normal)>0) apollo_test_beta[dists_normal] <- apollo_HB$svN[dists_normal]
      if(length(dists_lnp)>0) apollo_test_beta[dists_lnp] <- exp(apollo_HB$svN[dists_lnp])
      if(length(dists_lnn)>0) apollo_test_beta[dists_lnn] <- -exp(apollo_HB$svN[dists_lnn])
      if(length(dists_cnp)>0) apollo_test_beta[dists_cnp] <- apollo_HB$svN[dists_cnp]*(apollo_HB$svN[dists_cnp]>0)
      if(length(dists_cnn)>0) apollo_test_beta[dists_cnn] <- apollo_HB$svN[dists_cnn]*(apollo_HB$svN[dists_cnn]<0)
      if(length(dists_sb)>0){
        names(apollo_HB$gMINCOEF)=names(apollo_HB$svN)
        names(apollo_HB$gMAXCOEF)=names(apollo_HB$svN)
        apollo_test_beta[dists_sb] <- apollo_HB$gMINCOEF[dists_sb]+(apollo_HB$gMAXCOEF[dists_sb]-apollo_HB$gMINCOEF[dists_sb])/(1+exp(-apollo_HB$svN[dists_sb]))
      }
    }
    apollo_probabilities(apollo_test_beta, apollo_inputs, functionality="validate")
    if(!silent & !apollo_control$noDiagnostics ) cat(apollo_printLog(apollo_inputs$apolloLog))
    testLL = apollo_probabilities(apollo_test_beta, apollo_inputs, functionality="estimate")
    if(!workInLogs) testLL=log(testLL)
    # Maybe here we could return the value of the likelihood and print and error with cat, instead of simply stopping
    if(anyNA(testLL)) stop('Log-likelihood calculation fails at starting values!')
  }
  
  ### Test for unused parameters
  #if(!silent & !apollo_control$noValidation) cat("Testing probability function (apollo_probabilities)\n")
  if(!apollo_control$HB){
    apollo_beta_base=apollo_beta+0.001
    base_LL=apollo_probabilities(apollo_beta_base, apollo_inputs, functionality="estimate")
    if(workInLogs) base_LL=sum(base_LL) else base_LL=sum(log(base_LL))
    freeparams=apollo_beta_base[!names(apollo_beta_base)%in%apollo_fixed]
    for(p in names(freeparams)){
      apollo_beta_test1=apollo_beta_base
      apollo_beta_test2=apollo_beta_base
      apollo_beta_test1[p]=apollo_beta_test1[p]-0.001
      apollo_beta_test2[p]=apollo_beta_test2[p]+0.001
      test1_LL=apollo_probabilities(apollo_beta_test1, apollo_inputs, functionality="estimate")
      test2_LL=apollo_probabilities(apollo_beta_test2, apollo_inputs, functionality="estimate")
      if(workInLogs){
        test1_LL=sum(test1_LL)
        test2_LL=sum(test2_LL)
      } else{
        test1_LL=sum(log(test1_LL))
        test2_LL=sum(log(test2_LL))
      }
      if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      if(base_LL==test1_LL & base_LL==test2_LL) stop("Parameter ",p," does not influence the log-likelihood of your model!")
    }
  } else {
    ### Validation using HB estimation
    apollo_test_beta=apollo_beta
    if(!is.null(apollo_HB$gVarNamesFixed)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
      r <- names(apollo_beta)[r]
      apollo_test_beta[r] <- apollo_HB$FC[r]
    }
    if(!is.null(apollo_HB$gVarNamesNormal)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
      r <- names(apollo_beta)[r]
      dists_normal=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==1])
      dists_lnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==2])
      dists_lnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==3])
      dists_cnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==4])
      dists_cnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==5])
      dists_sb=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==6])
      if(length(dists_normal)>0) apollo_test_beta[dists_normal] <- apollo_HB$svN[dists_normal]
      if(length(dists_lnp)>0) apollo_test_beta[dists_lnp] <- exp(apollo_HB$svN[dists_lnp])
      if(length(dists_lnn)>0) apollo_test_beta[dists_lnn] <- -exp(apollo_HB$svN[dists_lnn])
      if(length(dists_cnp)>0) apollo_test_beta[dists_cnp] <- apollo_HB$svN[dists_cnp]*(apollo_HB$svN[dists_cnp]>0)
      if(length(dists_cnn)>0) apollo_test_beta[dists_cnn] <- apollo_HB$svN[dists_cnn]*(apollo_HB$svN[dists_cnn]<0)
      if(length(dists_sb)>0){
        names(apollo_HB$gMINCOEF)=names(apollo_HB$svN)
        names(apollo_HB$gMAXCOEF)=names(apollo_HB$svN)
        apollo_test_beta[dists_sb] <- apollo_HB$gMINCOEF[dists_sb]+(apollo_HB$gMAXCOEF[dists_sb]-apollo_HB$gMINCOEF[dists_sb])/(1+exp(-apollo_HB$svN[dists_sb]))
      }
    }
    apollo_beta_base=apollo_test_beta+0.001
    base_LL=apollo_probabilities(apollo_beta_base, apollo_inputs, functionality="estimate")
    if(workInLogs) base_LL=sum(base_LL) else base_LL=sum(log(base_LL))
    freeparams=apollo_beta_base[!names(apollo_beta_base)%in%apollo_fixed]
    for(p in names(freeparams)){
      apollo_beta_test1=apollo_beta_base
      apollo_beta_test2=apollo_beta_base
      apollo_beta_test1[p]=apollo_beta_test1[p]-0.001
      apollo_beta_test2[p]=apollo_beta_test2[p]+0.001
      test1_LL=apollo_probabilities(apollo_beta_test1, apollo_inputs, functionality="estimate")
      test2_LL=apollo_probabilities(apollo_beta_test2, apollo_inputs, functionality="estimate")
      if(workInLogs){
        test1_LL=sum(test1_LL)
        test2_LL=sum(test2_LL)
      } else{
        test1_LL=sum(log(test1_LL))
        test2_LL=sum(log(test2_LL))
      }
      if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      if(base_LL==test1_LL & base_LL==test2_LL) stop("Parameter ",p," does not influence the log-likelihood of your model!")    }
  }  
  
  # ################################## #
  #### HB estimation                ####
  # ################################## #
  
  if(apollo_control$HB){
    
    ### Check that apollo_fixed hasn't changed since calling apollo_validateInputs
    tmp <- tryCatch( get("apollo_fixed", envir=globalenv()), error=function(e) 1 )
    if( length(tmp)>0 && any(tmp %in% c(apollo_HB$gVarNamesFixed, apollo_HB$gVarNamesFixed)) ) stop("apollo_fixed seems to have changed since calling apollo_inputs.")
    
    ### Function masking apollo_probabilities and compatible with RSGHB
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
    ### use pre-scaling values as starting values in output 
    model$apollo_beta <- apollo_test_beta
    model$LLStart <- sum(testLL)
    ### model$LL0         <- sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL")))
    if(workInLogs) model$LL0 <- sum((apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"))) else model$LL0 <- sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL")))
    
    model$startTime   <- starttime
    model$apollo_control <- apollo_control
    model$nObs        <- nrow(database)
    model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
    endtime           <- Sys.time()
    model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
    model$apollo_fixed <- apollo_fixed
    model$estimationRoutine <- "Hierarchical Bayes"
    
    if(!is.null(model$F)){
      tmp <- coda::geweke.diag(model$F[,2:(ncol(model$F))], frac1=0.1, frac2=0.5)[[1]]
      names(tmp) <- model$params.fixed
      model$F_convergence=tmp
    }
    if(!is.null(model$A)){
      tmp <- coda::geweke.diag(model$A[,2:(ncol(model$A))], frac1=0.1, frac2=0.5)[[1]]
      model$A_convergence=tmp
    }
    if(!is.null(model$D)){
      # This assumes the matrix is square
      tmp <- c()
      for(i in 1:dim(model$D)[1]) for(j in 1:i){
        if(i==1 & j==1) Dmatrix <- as.matrix(model$D[i,j,]) else Dmatrix <- cbind(Dmatrix, as.vector(model$D[i,j,]))
        tmp <- c(tmp, paste(colnames(model$A)[i+1],colnames(model$A)[j+1], sep="_"))
      }
      colnames(Dmatrix) <- tmp
      tmp <- coda::geweke.diag(Dmatrix, frac1=0.1, frac2=0.5)[[1]]
      model$D_convergence=tmp
    }
    
    if(length(apollo_HB$gVarNamesFixed)>0 | length(model$apollo_fixed)>0){
      if(length(apollo_HB$gVarNamesFixed)>0){
        non_random=matrix(0,nrow=length(apollo_HB$gVarNamesFixed),2)
        non_random[,1]=colMeans(model$F)[2:ncol(model$F)]
        non_random[,2]=apply(model$F,FUN=stats::sd,2)[2:ncol(model$F)]
        rownames(non_random)=apollo_HB$gVarNamesFixed}
      if(length(model$apollo_fixed)>0){
        if(length(apollo_HB$gVarNamesFixed)>0){
          non_random=rbind(non_random,cbind(matrix(model$apollo_beta[model$apollo_fixed]),NA))
          rownames(non_random)[(length(apollo_HB$gVarNamesFixed)+1):nrow(non_random)]=model$apollo_fixed
        } else{
          non_random=cbind(matrix(model$apollo_beta[model$apollo_fixed]),NA)
          rownames(non_random)=model$apollo_fixed
        }
      }
      colnames(non_random)=c("Mean","SD")
      originalOrder <- names(model$apollo_beta)[names(model$apollo_beta) %in% rownames(non_random)]
      model$chain_non_random=non_random[originalOrder,,drop=FALSE]
    }
    
    apollo_HB$gVarNamesFixed <- model$params.fixed
    apollo_HB$gVarNamesNormal <- model$params.vary
    if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
      random_mean     = matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
      random_mean[,1] = colMeans(model$A)[2:ncol(model$A)]
      random_mean[,2] = apply(model$A,FUN=stats::sd,2)[2:ncol(model$A)]
      rownames(random_mean)=apollo_HB$gVarNamesNormal
      colnames(random_mean)=c("Mean","SD")
      model$chain_random_mean=random_mean
      
      random_cov_mean           = apply(model$D,FUN=mean,c(1,2))
      random_cov_sd             = apply(model$D,FUN=stats::sd,c(1,2))
      rownames(random_cov_mean) = apollo_HB$gVarNamesNormal
      colnames(random_cov_mean) = apollo_HB$gVarNamesNormal
      model$chain_random_cov_mean=random_cov_mean
      
      rownames(random_cov_sd) = apollo_HB$gVarNamesNormal
      colnames(random_cov_sd) = apollo_HB$gVarNamesNormal
      model$chain_random_cov_sd=random_cov_sd
      
      posterior=matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
      posterior[,1]=colMeans(model$C)[3:ncol(model$C)]
      posterior[,2]=apply(model$C,FUN=stats::sd,2)[3:ncol(model$C)]
      rownames(posterior)=apollo_HB$gVarNamesNormal
      model$posterior_mean=posterior
      
      ### create matrix of draws from distributions
      
      draws=10000
      covMat=random_cov_mean
      meanA=random_mean[,1]
      pars = length(meanA)
      covMat=as.matrix(covMat)
      Ndraws=mvtnorm::rmvnorm(draws,meanA,covMat,method="chol")
      i=1
      while(i<(pars+1)){
        if(apollo_HB$gDIST[i]==6){
          Ndraws[,i]=apollo_HB$gMINCOEF+(apollo_HB$gMAXCOEF[i]-apollo_HB$gMINCOEF[i])*1/(1+exp(-Ndraws[,i]))
        }
        if(apollo_HB$gDIST[i]==5){
          Ndraws[,i]=(Ndraws[,i]<0)*Ndraws[,i]
        }
        if(apollo_HB$gDIST[i]==4){
          Ndraws[,i]=(Ndraws[,i]>0)*Ndraws[,i]
        }
        if(apollo_HB$gDIST[i]==3){
          Ndraws[,i]=-exp(Ndraws[,i])
        }
        if(apollo_HB$gDIST[i]==2){
          Ndraws[,i]=exp(Ndraws[,i])
        }
        i=i+1
      }
      
    }
    
    if(length(scaling)>0 && !is.na(scaling)){
      s=1
      while(s<=length(scaling)){
        ss=names(scaling)[s]
        if(ss%in%colnames(model$C)) model$C[,ss]=scaling[s]*model$C[,ss]
        if(ss%in%colnames(model$Csd)) model$Csd[,ss]=scaling[s]*model$Csd[,ss]
        if(ss%in%colnames(model$F)) model$F[,ss]=scaling[s]*model$F[,ss]
        if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){if(ss%in%colnames(Ndraws)) Ndraws[,ss]=scaling[s]*Ndraws[,ss]}
        if(ss%in%rownames(model$chain_non_random)) model$chain_non_random[ss,]=scaling[s]*model$chain_non_random[ss,]
        if(ss%in%rownames(model$posterior_mean)) model$posterior_mean[ss,]=scaling[s]*model$posterior_mean[ss,]
        s=s+1
      }
      model$scaling <- scaling
    }
    
    if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
      model$random_coeff_summary=cbind(colMeans(Ndraws),apply(Ndraws,2,sd))
      colnames(model$random_coeff_summary)=c("Mean","SD")
      if(length(apollo_HB$gVarNamesNormal)>1){
        model$random_coeff_covar=cov(Ndraws)
        model$random_coeff_corr=cor(Ndraws)
      }
    }
    
    ### produce model$estimate
    
    panelData <- apollo_control$panelData
    indivID   <- database[,apollo_control$indivID]
    nObs <- length(indivID)
    if(!panelData) indivID <- 1:nObs
    nIndiv <- length(unique(indivID))
    obsPerIndiv <- as.vector(table(indivID))
    
    if(is.null(model$chain_non_random)){
      fc1 <- NULL
    }else{
      fc1 <- stats::setNames(as.list(model$chain_non_random[,1]),names(model$chain_non_random[,1]))
    }
    
    if(is.null(model$C)){
      b1 <- NULL
    }else{
      M=model$C[,-c(1,2),drop=FALSE]
      M1 <- matrix(0, nrow=nObs, ncol=ncol(M))
      r1 <- 1
      for(i in 1:nIndiv){
        r2 <- r1 + obsPerIndiv[i] - 1
        M1[r1:r2,] <- matrix(as.vector(M[i,]), nrow=r2-r1+1, ncol=ncol(M), byrow=TRUE)
        r1 <- r2 + 1
      }
      b1  <- stats::setNames(as.data.frame(M1), colnames(M))
    } 
    model$estimate=c(fc1,b1)
    
    # Report number of times the probs have been censored
    if(exists("HBcensor", envir=apollo_inputs$apolloLog)){
      txt <- paste0(" Please note that in at least some iterations RSGHB has\n",
                    " avoided numerical issues by left censoring the\n",
                    " probabilities. This has the side effect of zero or\n",
                    " negative probabilities not leading to failures!", collapse="")
      apollo_addLog(title="WARNING: RSGHB has censored the probabilities", content=txt, apollo_inputs$apolloLog)
    }
    if(!silent) cat("\n", apollo_printLog(apollo_inputs$apolloLog), sep="")
    model$apolloLog <- apollo_inputs$apolloLog
    
    return(model)
  }
  
  # ################################## #
  #### classical estimation         ####
  # ################################## #
  
  ### Create cluster (if needed)
  if(apollo_control$nCores==1) cl <- NA else {
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=silent)
    apollo_control$nCores <- length(cl)
    apollo_inputs$apollo_control$nCores <- length(cl)
  }
  on.exit(if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl))
  
  ### Create loglike function
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                       apollo_probabilities, apollo_inputs, 
                                       estimate_settings, cl=cl)
  
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
  model <- maxLik::maxLik(apollo_logLike, start=beta_var_val,
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
    cat("\nEstimated values:\n")
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
        #if(!is.null(numDeriv_settings) && !is.null(numDeriv_settings$d)) I <- 2+(2*numDeriv_settings$d)*( k*(k+1)/2 )
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
  
  ### Calculate s.e.
  dummyVCM           <- matrix(NA, nrow=length(model$estimate), ncol=length(model$estimate))
  rownames(dummyVCM) <- names(model$estimate)
  colnames(dummyVCM) <- names(model$estimate)
  model$varcov     <- tryCatch(stats::vcov(model), error=function(e) return(dummyVCM))
  model$robvarcov  <- tryCatch(sandwich::sandwich(model), error=function(e) return(dummyVCM))
  model$se         <- sqrt(diag(model$varcov))
  model$robse      <- sqrt(diag(model$robvarcov))
  model$corrmat    <- tryCatch(model$varcov/(model$se%*%t(model$se)), error=function(e) return(dummyVCM))
  model$robcorrmat <- tryCatch(model$robvarcov/(model$robse%*%t(model$robse)) , error=function(e) return(dummyVCM))
  
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
  
  #### Calculate gradient
  #if(!silent) cat("Calculating gradient norm... ")
  #model$grad <- tryCatch(numDeriv::grad(apollo_logLike, model$estimate, countIter=FALSE, writeIter=FALSE, sumLL=TRUE),
  #                       error=function(e) return(NA))
  #if(!silent) if(length(model$grad)==1 && is.na(model$grad)) cat("ERROR\n") else cat(round( sqrt(sum(model$grad^2)) ,4), "\n", sep="")
  
  
  ### Calculate probabilities to identify outliers
  P <- exp(apollo_logLike(model$estimate))
  if(apollo_control$panelData){
    nObsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
  } else nObsPerIndiv <- 1
  model$avgCP <- P^(1/nObsPerIndiv)
  names(model$avgCP) <- unique(database[,apollo_control$indivID])
  
  
  ### Calculate Zero LL
  if(!silent) cat("Calculating LL(0)... ")
  model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
                        error=function(e) return(NA))
  if(!anyNA(model$LL0)){
    model$LL0 <- ifelse( workInLogs, sum(model$LL0), sum(log(model$LL0)) )
    if(!silent) cat(round(model$LL0,2), "\n", sep="")
  } else if(!silent) cat("No LL0 for some components.\n")
  
  
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
  if(!silent) cat("Calculating LL of each model component...")
  Pout <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="output"),
                   error=function(e) return(NA))
  if(!anyNA(Pout) && is.list(Pout)){
    # Give name to unnamed components
    origNames <- names(Pout)
    newNames  <- paste0("component_", 1:length(Pout))
    if(!is.null(origNames)) newNames <- ifelse(origNames!="", origNames, newNames)
    names(Pout) <- newNames
    # Get log of likelihood with "model" first
    tmp <- c("model", newNames[newNames!="model"])
    if(!workInLogs) LLout <- lapply(Pout[tmp], log) else LLout <- Pout[tmp]
    LLout <-lapply(LLout,sum)
    if(!silent) cat("Done.\n")
    model$Pout  <- LLout
    model$LLout <- unlist(lapply(model$Pout, sum))
  } else{
    model$Pout  <- Pout
    model$LLout <- list(NA)
    if(!silent) cat("Not applicable to all components.\n")
  }
  
  if(!silent) cat("\n", apollo_printLog(apollo_inputs$apolloLog), sep="")
  
  ### use pre-scaling values as starting values in output
  ### OLD VERSION model$apollo_beta <- c(beta_var_val, beta_fix_val)[names(apollo_beta)]
  model$bootstrapSE <- bootstrapSE
  model$apollo_beta <- temp_start
  model$LLStart     <- sum(testLL)
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
  model$apolloLog <- apollo_inputs$apolloLog
  
  ### Apply scaling to estimates and covariance matrix
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
  
  
  
  if(length(scaling)>0 && !is.na(scaling)) model$scaling <- scaling
  return(model)
}
