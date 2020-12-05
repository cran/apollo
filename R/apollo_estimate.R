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
#'                                   \item \strong{writeIter}: Logical. Writes value of the parameters in each iteration to a csv file. 
#'                                                             Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
#'                                   \item \strong{hessianRoutine}: Character. Name of routine used to calculate the Hessian of the loglikelihood 
#'                                                                  function after estimation. Valid values are \code{"analytic"} (default), 
#'                                                                  \code{"numDeriv"} (to use the numeric routine in package numDeric), \code{"maxLik"} 
#'                                                                  (to use the numeric routine in packahe maxLik), and \code{"none"} to avoid 
#'                                                                  estimating the Hessian and the covariance matrix. Only used if \code{apollo_control$HB=FALSE}.
#'                                   \item \strong{printLevel}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. 
#'                                                              Ignored if apollo_control$HB is TRUE. Default is 3.
#'                                   \item \strong{constraints}: Character vector. Constraints on parameters to estimate. For example \code{c('b1>0', 'b1 + b2>1')}.
#'                                                               Only \code{>}, \code{>=} and \code{=} can be used. And they cannot be mixed (e.g. 
#'                                                               \code{c(b1=b2, b2>0)} will fail). All parameter names must be on the left side. Fixed 
#'                                                               parameters cannot go into constraints. Alternatively, constraints can be defined as in 
#'                                                               \link[maxLik]{maxLik}. Constraints can only be used with maximum likelihood estimation, 
#'                                                               and the BFGS routine in particular.
#'                                   \item \strong{scaling}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. 
#'                                                           If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but 
#'                                                           \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter
#'                                                           b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, 
#'                                                           specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
#'                                  \item \strong{maxLik_settings}: List. Additional settings for maxLik. See argument \code{control} in \link[maxLik]{maxBFGS}, 
#'                                                                  \link[maxLik]{maxBHHH} and \link[maxLik]{maxNM} for more details. Only used for maximum 
#'                                                                  likelihood estimation.
#'                                  \item \strong{numDeriv_settings}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. 
#'                                                                    See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                                  \item \strong{bootstrapSE}: Numeric. Number of bootstrap samples to calculate standard errors. Default is 0, meaning
#'                                                              no bootstrap s.e. will be calculated. Number must zero or a positive integer. Only used
#'                                                              if \code{apollo_control$HB} is \code{FALSE}.
#'                                  \item \strong{bootstrapSeed}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples.
#'                                                                Only used if \code{bootstrapSE>0}. Default is 24.
#'                                  \item \strong{silent}: Logical. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                  \item \strong{scaleHessian}: Logical. If TRUE, parameters are scaled to 1 for Hessian estimation. Default is TRUE.
#'                                 }
#' @return model object
#' @export
#' @importFrom numDeriv hessian grad
#' @importFrom maxLik maxLik
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov runif
apollo_estimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
  # #################### #
  #### Loading Inputs ####
  # #################### #
  
  ### First checkpoint
  time1 <- Sys.time()
  
  ### Detach things if necessary
  apollo_detach()
  
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="analytic", printLevel=3L, constraints=NULL, maxLik_settings=NULL, 
                  numDeriv_settings=list(), scaling=NA, bootstrapSE=0, bootstrapSeed=24, silent=FALSE, 
                  scaleHessian=TRUE)
  if(length(estimate_settings)==1 && is.na(estimate_settings)) estimate_settings <- default
  if(is.null(estimate_settings$maxLik_settings)){
    estimate_settings$maxLik_settings <- list(printLevel=3, iterlim=200)
    if(!is.null(estimate_settings$printLevel)) estimate_settings$maxLik_settings$printLevel <- estimate_settings$printLevel
    if(!is.null(estimate_settings$maxIterations)) estimate_settings$maxLik_settings$iterlim <- estimate_settings$maxIterations
  } else {
    test <- !is.null(estimate_settings$maxLik_settings$iterlim) && !is.null(estimate_settings$maxIterations)
    if(test) estimate_settings$maxLik_settings$iterlim <- estimate_settings$maxIterations
  }
  tmp <- names(default)[!(names(default) %in% names(estimate_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings[[i]] <- default[[i]]
  rm(default, tmp)
  if(exists("i")) rm(i)
  apollo_inputs$apollo_scaling <- estimate_settings$scaling
  apollo_inputs$scaling <- NULL
  
  
  
  ### Warn the user in case elements in apollo_inputs are different from those in the global environment
  apollo_compareInputs(apollo_inputs)
  
  ### Extract variables from apollo_input
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
  printLevel        = estimate_settings[["printLevel"]]
  silent            = estimate_settings[["silent"]]
  constraints       = estimate_settings[["constraints"]]
  #scaling           = estimate_settings[["scaling"]]
  bootstrapSE       = estimate_settings[["bootstrapSE"]]
  bootstrapSeed     = estimate_settings[["bootstrapSeed"]]
  maxLik_settings   = estimate_settings[["maxLik_settings"]]
  
  if(maxIterations==0){
    apollo_control$noValidation=TRUE
    apollo_inputs$apollo_control$noValidation=TRUE
    estimate_settings$writeIter = FALSE
    writeIter = FALSE ## we're using both the one in estimate_settings and a local one, we should fix that
  } 
  
  apollo_inputs$silent <- estimate_settings$silent
  silent <- apollo_inputs$silent
  debug  <- apollo_inputs$apollo_control$debug
  
  
  # ########################## #
  #### Validation of inputs ####
  # ########################## #
  
  ### Validation of input
  apollo_checkArguments(apollo_probabilities, apollo_randCoeff, apollo_lcPars)
  if( !(estimationRoutine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimationRoutine. Use 'bfgs', 'bhhh' or 'nr'.")
  if( !(estimate_settings$hessianRoutine %in% c('analytic', 'numDeriv', 'maxLik', 'none')) ) stop("Invalid hessianRoutine. Use 'analytic', 'numDeriv', 'maxLik' or 'none'.")
  
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
    estimationRoutine               = "bfgs"
    apollo_inputs$estimationRoutine = estimationRoutine
    if(!silent) apollo_print("WARNING: Estimation routine changed to 'BFGS'. Only 'BFGS' supports constrained optimization.\n")
  }
  if(is.vector(constraints) && is.character(constraints)){
    nCon <- length(constraints)
    bVar <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
    nVar <- length(bVar)
    bVar <- list2env(setNames(split(diag(nVar), rep(1:nVar,each=nVar)), bVar))
    bVal <- list2env(as.list(apollo_beta))
    A <- matrix(0, nrow=nCon, ncol=nVar); b <- rep(0, nCon); mid0 <- ''
    for(i in 1:nCon){
      # turn constraint into expression, and split into left, right and middle
      e <- tryCatch(str2lang(constraints[i]), error=function() NULL)
      test <- is.null(e) || !is.call(e) || length(e)!=3
      if(test) stop('Constraint "', constraints[i], '" is not a valid linear constraint expression.')
      mid <- e[[1]]; lef <- e[[2]]; rig <- e[[3]]
      # Checks
      test <- is.symbol(mid) && (as.character(mid) %in% c('>=', '>', '='))
      if(!test) stop('Constraint "', constraints[i], '" does not contain one (and only one) of the following: >=, >, or =.')
      test <- mid0=='' | as.character(mid)==mid0; mid0 <- as.character(mid)
      if(!test) stop('All constraints must involve the same operator >=, >, or =. Mixing is not allowed.')
      test <- length(all.vars(rig))==0
      if(!test) stop('The right side of constraint "', constraints[i],'" should only contain numeric values.')
      test <- all(all.vars(lef) %in% ls(bVar))
      if(!test) stop('All the variables in the left side of constraint "', constraints[i],
                     '" should be in apollo_beta. Fixed parameters cannot go into constraints.')
      test <- eval(e, envir=bVal)
      if(!test) stop('Starting values of parameters do not satisfy constraint "', constraints[i],'".')
      # Fill A & b
      A[i,] <- eval(lef, envir=bVar)
      b[i]  <- eval(rig)
    }
    if(mid0 %in% c('>', '>=')) constraints <- list(ineqA=A, ineqB=b) else constraints <- list(eqA=A, eqB=b)
    rm(nCon, bVar, nVar, A, b, mid0, i, e, test, mid, lef, rig, bVal)
  }
  
  # Check writeIter
  if(estimationRoutine!="bfgs" & writeIter==TRUE){
    writeIter                   = FALSE
    estimate_settings$writeIter = writeIter
    txt <- "witeIter set to FALSE. Writing parameters values at each iteration is only available for BFGS estimation method."
    if(!silent) apollo_print(txt) else warning(txt)
    rm(txt)
  }
  
  # create temporary copy of starting values for use later
  temp_start = apollo_beta
  
  # Validate scaling
  scaling <- setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                      names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  test <- is.null(apollo_inputs$apollo_scaling)
  test <- test || (length(apollo_inputs$apollo_scaling)==1 && is.na(apollo_inputs$apollo_scaling))
  if(test) apollo_inputs$apollo_scaling <- scaling else {
    scaling[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling
    #if(any(!(names(scaling) %in% names(apollo_beta)))) stop("Some parameters included in 'scaling' are not included in 'apollo_beta'")
      if(anyDuplicated(names(scaling))){
        txt <- paste0(names(scaling)[duplicated(names(scaling))], collapse=", ")
        txt <- paste0("The \"scaling\" setting contains duplicate elements (", txt, ").")
        stop(txt)
      }
      if(!all(names(scaling) %in% names(apollo_beta))){
        txt <- names(scaling)[!(names(scaling) %in% names(apollo_beta))]
        txt <- paste0(txt, collapse=", ")
        txt <- paste0("Some parameters included in 'scaling' (", txt, ") are not included in 'apollo_beta'.")
        stop(txt)
      }

    if(any(names(scaling) %in% apollo_fixed)) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
    if(any(scaling<0)){
      scaling <- abs(scaling)
      txt <- 'Some negative values in "scaling" were replaced by their absolute value'
      if(!silent) apollo_print(paste0('WARNING: ', txt, '.')) else warning(txt)
    }
    if(any(scaling<=0)) stop('All terms in "scaling" should be strictly positive!')
    txt <- "During estimation, parameters will be scaled using the values in estimate_settings$scaling"
    if(!all(scaling==1)){ if(!silent) apollo_print(txt) else warning(txt)}
    rm(txt)
    apollo_inputs$apollo_scaling <- scaling
  }
  
  ### Recommend using multi-core if appropriate
  if(!apollo_control$HB && apollo_inputs$apollo_control$mixing && apollo_inputs$apollo_control$nCores==1 && !silent){
    n   <- nrow(apollo_inputs$database)
    tmp <- apollo_inputs$apollo_draws$interNDraws
    if(!is.null(tmp) && is.numeric(tmp) && tmp>0) n <- n*tmp
    tmp <- apollo_inputs$apollo_draws$interNDraws
    if(!is.null(tmp) && is.numeric(tmp) && tmp>0) n <- n*tmp
    tmp <- paste0('WARNING: You can use multiple processor cores to speed up estimation. ',
                  'To do so, specify the desired number of cores using the setting nCores ',
                  'inside "apollo_control", for example "nCores=2". This computer has ',
                  parallel::detectCores(), ' available cores.')
    if(parallel::detectCores()>2) tmp <- paste0(tmp, ' We recommend using no more ',
                                                'than ', parallel::detectCores()-1, ' cores.')
    if(n>1e5) apollo_print(tmp, highlight=TRUE)
    rm(n, tmp)
  }
  
  ### Create useful variables
  nObs  <- nrow(apollo_inputs$database)
  indiv <- unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  nObsPerIndiv <- rep(0, length(indiv))
  for(n in 1:length(indiv)) nObsPerIndiv[n] <- sum(apollo_inputs$database[,apollo_inputs$apollo_control$indivID]==indiv[n])
  names(nObsPerIndiv) <- indiv
  #nObsPerIndiv <- table(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  #nObsPerIndiv <- setNames(as.vector(nObsPerIndiv), names(nObsPerIndiv))
  #nObsPerIndiv <- nObsPerIndiv[nObsPerIndiv>0]
  
  # ############# #
  #### Back-up ####
  # ############# #
  
  ### Store unaltered version of apollo_probabilities, apollo_randCoeff and apollo_lcPars
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_lcPars
  
  ### If multicore is in use, save original apollo_inputs to disk before changing it, and make sure to close cluster
  if(apollo_inputs$apollo_control$nCores>1){
    saveRDS(apollo_inputs, file=paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs"))
    on.exit({
      tmp <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs")
      if(file.exists(tmp)) file.remove(tmp)
      rm(tmp)
    })
  }
  
  # ############# #
  #### Call HB ####
  # ############# #
  
  ### Call apollo_estimateHB
  if(apollo_control$HB){
    ### Insert scaling if needed
    apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
    test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
    if(test) apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
    test <- is.function(apollo_inputs$apollo_lcPars)
    if(test) apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
    ### Call apollo_estimateHB
    preHBtime=as.numeric(difftime(Sys.time(),time1,units='secs'))
    model <- apollo_estimateHB(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings)
    model$timeTaken <- as.numeric(model$timeTaken+preHBtime)
    model$timePre   <- as.numeric(model$timePre+preHBtime)
    return(model)
  }
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  
  # Scale starting parameters
  if(length(apollo_inputs$apollo_scaling)>0 && !anyNA(apollo_inputs$apollo_scaling)){
    r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- apollo_beta[r]/apollo_inputs$apollo_scaling[r]
  }
  
  ### Scale apollo_probabilities, apollo_randCoeff and apollo_lcPars, and names components in apollo_probabilities
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
  test <- is.function(apollo_inputs$apollo_lcPars)
  if(test) apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
  
  ### Validate likelihood function
  test <- apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$mixing
  if(test && !silent) apollo_print(paste0('CAUTION: In models using mixing and workInLogs=TRUE, Apollo assumes (and does not check) that ',
                                          'the call to apollo_panelProd is immediately followed by a call to apollo_avgInterDraws.'))
  # Validate LL at starting values
  if(!silent) apollo_print(c("\n", "Testing likelihood function..."))
  #if(covarOnly){ silentOriginal = silent; apollo_inputs$silent = TRUE }
  testLL <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")
  testLL <- testLL[["model"]]
  if(!apollo_inputs$apollo_control$workInLogs) testLL <- log(testLL)
  if(!silent & any(!is.finite(testLL))){
    apollo_print("Log-likelihood calculation fails at starting values!")
    apollo_print("Affected individuals:")
    LLtemp <- data.frame(ID=indiv, LL=testLL)
    LLtemp <- subset(LLtemp,!is.finite(LLtemp[,2]))
    colnames(LLtemp) <- c("ID","LL")
    print(LLtemp, row.names=FALSE)
  }
  #if(covarOnly) apollo_inputs$silent = silentOriginal
  if(any(!is.finite(testLL))) stop('Log-likelihood calculation fails at values close to the starting values!')
  
  ### Restore unaltered functions
  apollo_probabilities <- apollo_probabilities_ORIG
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_inputs$apollo_randCoeff <- apollo_randCoeff_ORIG
  if(is.function(apollo_inputs$apollo_lcPars)) apollo_inputs$apollo_lcPars <- apollo_lcPars_ORIG
  rm(test)
  
  
  # ################################## #
  #### Pre-process                  ####
  # ################################## #
  if(!silent) apollo_print(c("\n","Pre-processing likelihood function..."))
  
  ### Create multi-core version of likelihood (if needed)
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                       apollo_probabilities, apollo_inputs, 
                                       estimate_settings, cleanMemory=TRUE)
  on.exit({
    test <- apollo_control$nCores>1 && exists('apollo_logLike', inherits=FALSE)
    test <- test && !anyNA(environment(apollo_logLike)$cl)
    if(test) parallel::stopCluster(environment(apollo_logLike)$cl)
  }, add=TRUE)
  
  ### Create gradient function if required (and possible)
  if(!is.null(apollo_inputs$apollo_control$analyticGrad) && apollo_inputs$apollo_control$analyticGrad){
    # apollo_makeGrad will create gradient function ONLY IF all components have analytical gradient.
    grad <- apollo_makeGrad(apollo_beta, apollo_fixed, apollo_logLike, validateGrad=TRUE)
    if(is.null(grad) && debug) apollo_print("Analytical gradients could not be calculated for all components, numerical gradients will be used.")
    if(is.null(grad)){
      apollo_inputs$apollo_control$analyticGrad <- FALSE # Not sure this is necessary, but it can't hurt
      environment(apollo_logLike)$analyticGrad  <- FALSE # fix for iteration counting
      environment(apollo_logLike)$nIter         <- 1     # fix for iteration counting
    }
  } else grad <- NULL
  
  # ############################################## #
  #### Test all parameters influence likelihood ####
  # ############################################## #
  
  ### Check that likelihood is sensitive to changes in all parameters
  if(!apollo_control$noValidation){
    if(!silent) cat("\nTesting influence of parameters")
    beta0   = apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    beta1   = beta0 + 0.001*runif(length(beta0))
    base_LL = apollo_logLike(beta1)
    if(any(!is.finite(base_LL))){
      if(!silent){
        cat("\n")
        apollo_print(paste0("During testing, Apollo added disturbances smaller than 0.001 ",
                            "to all starting values. This led to a log-likelihood calculation failure!"))
        apollo_print("Affected individuals:")
        LLtemp <- subset( data.frame(ID=indiv, LL=base_LL), !is.finite(base_LL))
        colnames(LLtemp) <- c("ID","LL")
        if(!silent) print(LLtemp, row.names=FALSE)
      }; stop('Log-likelihood calculation fails at values close to the starting values!')
    }
    base_LL = sum(base_LL)
    for(p in names(beta0)){
      beta1p <- beta1 - (names(beta1)==p)*0.001
      beta1m <- beta1 + (names(beta1)==p)*0.001
      test1_LL = sum( apollo_logLike(beta1p) )
      test2_LL = sum( apollo_logLike(beta1m) )
      if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      if(base_LL==test1_LL & base_LL==test2_LL) stop("Parameter ",p," does not influence the log-likelihood of your model!")
      if(!silent) cat(".")
    }
    rm(beta0, beta1, beta1p, beta1m, base_LL, test1_LL, test2_LL, p)
  }
  
  
  # ################################## #
  #### Classical estimation         ####
  # ################################## #
  
  ### Second checkpoint
  time2 <- Sys.time()
  
  ### Split parameters between variable and fixed
  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  beta_fix_val <- apollo_beta[apollo_fixed]
  
  ### Preparations for iteration writing and counting
  if(estimate_settings$writeIter){
    # Remove modelName_iterations if it exists
    tmp <- paste0(apollo_inputs$apollo_control$modelName, "_iterations.csv")
    txt <- paste0('Could not delete old ', tmp, ' file. New iterations will be written after old ones.')
    if(file.exists(tmp)) tryCatch(file.remove(tmp), error=function(e) apollo_print(txt))
    # Create modelName_iterations with starting values if using analytic gradient
    if(is.function(grad)) apollo_writeTheta(beta_var_val, sum(testLL), apollo_inputs$apollo_control$modelName)
    rm(tmp, txt)
  }
  # Reset lastFuncParam if it exists
  if(exists("lastFuncParam")){
    tmp <- globalenv()
    assign("lastFuncParam", rep(0, length(beta_var_val)), envir=tmp)
    rm(tmp)
  }
  
  ### Main (classical) estimation
  if(!silent) apollo_print(c("\n", "Starting main estimation")) else maxLik_settings$printLevel=0
  model <- maxLik::maxLik(apollo_logLike, grad=grad, start=beta_var_val,
                          method=estimationRoutine, finalHessian=FALSE,
                          control=maxLik_settings,
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
    apollo_print("Estimated parameters:")
    tmp <- c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
    tmp[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling*tmp[names(apollo_inputs$apollo_scaling)]
    tmp <- matrix(tmp, nrow=length(tmp), ncol=1, dimnames=list(names(tmp), 'Estimate'))
    apollo_print(tmp); rm(tmp)
    apollo_print("\n")
  }
  
  # If estimation failed, return whatever can be salvaged
  if(!succesfulEstimation){
    if(exists("model")){
      if(!silent) apollo_print("Estimation failed. No covariance matrix to compute.", highlight=TRUE)
      #if(length(scaling)>0 && !is.na(scaling)) model$scaling <- scaling
      #return(model)
    } else stop("ERROR: Estimation failed, no estimated model to return.")
  }
  ### Stores number of iterations
  model$nIter <- apollo_logLike(NA, getNIter=TRUE)
  
  ### Calculate probabilities to identify outliers
  indLL <- apollo_logLike(model$estimate)
  P     <- exp(indLL)
  if(apollo_control$panelData){
    model$avgLL <- setNames(indLL/nObsPerIndiv, names(nObsPerIndiv))
    model$avgCP <- setNames(P^(1/nObsPerIndiv), names(nObsPerIndiv))
  } else {
    model$avgLL <- setNames(indLL, 1:length(P))
    model$avgCP <- setNames(P, 1:length(P))
  }
  
  ### Store functions and other relevant things
  model$apollo_randCoeff     <- apollo_randCoeff_ORIG
  model$apollo_lcPars        <- apollo_lcPars_ORIG
  model$apollo_probabilities <- apollo_probabilities_ORIG
  model$apollo_fixed         <- apollo_fixed
  model$modelTypeList        <- environment(apollo_logLike)$mType
  model$nObsTot              <- environment(apollo_logLike)$nObsTot
  
  ### Third checkpoint
  time3 <- Sys.time()
  
  # ############################ #
  #### Hessian & covar matrix ####
  # ############################ #
  
  if(!succesfulEstimation) estimate_settings$hessianRoutine <- 'none'
  b <- c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
  b[names(apollo_inputs$apollo_scaling)] <- b[names(apollo_inputs$apollo_scaling)]*apollo_inputs$apollo_scaling
  varcov <- apollo_varcov(apollo_beta=b, apollo_fixed, 
                          varcov_settings=list(hessianRoutine    = estimate_settings$hessianRoutine,
                                               scaleBeta         = estimate_settings$scaleHessian, 
                                               numDeriv_settings = estimate_settings$numDeriv_settings,
                                               apollo_logLike    = apollo_logLike,
                                               apollo_grad       = grad))
  if(is.list(varcov)) for(i in names(varcov)) model[[i]] <- varcov[[i]]
  rm(varcov)
  
  
  ### Close cluster and delete apollo_logLike
  if(!environment(apollo_logLike)$singleCore){
    parallel::stopCluster(environment(apollo_logLike)$cl)
  }; rm(apollo_logLike)
  
  ### Restore apollo_inputs including database and draws
  if(apollo_control$nCores>1 & is.null(apollo_inputs$database)){
    if(debug) cat('Restoring data in main thread...')
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs")
    apollo_inputs <- tryCatch(readRDS(file=fileName), error=function(e) NULL)
    if(!is.null(apollo_inputs)){
      tmp <- .GlobalEnv
      assign("apollo_inputs", apollo_inputs, envir=tmp)
      unlink(fileName); rm(fileName)
      if(debug) cat(' Done. ',sum(gc()[,2]),'MB of RAM in use\n',sep='') # do not change to apollo_print
    } else if(debug) cat(' Failed.\n')
  }
  
  ### Calculate bootstrap s.e. if required
  if(bootstrapSE>0){
    if(!silent) apollo_print("\nStarting bootstrap calculation of standard errors.")
    # If using multiple cores, save copy of apollo_inputs
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName, "_inputs_extra")
    if(apollo_inputs$apollo_control$nCores>1) saveRDS(apollo_inputs, file=fileName)
    # Call apollo_bootstrap
    tmp <- list(estimationRoutine=estimationRoutine, maxIterations=maxIterations,
                writeIter=FALSE, hessianRoutine="none", printLevel=printLevel,
                maxLik_settings=maxLik_settings, silent=silent)
    tmpPar=c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
    model$bootvarcov <- apollo_bootstrap(tmpPar, apollo_fixed,
                                         apollo_probabilities, apollo_inputs,
                                         estimate_settings=tmp,
                                         bootstrap_settings=list(nRep=bootstrapSE,
                                                                  seed=bootstrapSeed,
                                                                 calledByEstimate=TRUE))$varcov
    model$bootse <- sqrt(diag(model$bootvarcov))
    bVar         <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    dummyVCM     <- matrix(NA, nrow=length(bVar), ncol=length(bVar), dimnames=list(names(bVar), names(bVar)))
    model$bootcorrmat <- tryCatch(model$bootvarcov/(model$bootse%*%t(model$bootse)), error=function(e) return(dummyVCM))
    # update number of bootstrap repetitions
    bootstrapSE <- tryCatch(nrow(read.csv(paste0(apollo_inputs$apollo_control$modelName, '_bootstrap_params.csv'))),
                            error=function(e) bootstrapSE)
    if(length(apollo_fixed)>0){
      model$bootse <- c(model$bootse, rep(NA, length(apollo_fixed)))
      names(model$bootse) <- c(colnames(model$bootvarcov), apollo_fixed)
      model$bootse <- model$bootse[names(apollo_beta)]
    }
    # If using multiple cores, restore apollo_inputs
    if(apollo_inputs$apollo_control$nCores>1) apollo_inputs <- tryCatch(readRDS(file=fileName), error=function(e) NULL)
    if(is.null(apollo_inputs)) stop('apollo_inputs could not be restored from disk after bootstrap')
    rm(bVar, dummyVCM)
  }
  
  # #################### #
  #### Prepare output ####
  # #################### #
  
  ### Scale apollo_probabilities, apollo_randCoeff and apollo_lcPars, and names components in apollo_probabilities
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
  test <- is.function(apollo_inputs$apollo_lcPars)
  if(test) apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
  
  ### Restore fixed parameters to model$estimate
  temp           = c(model$estimate, apollo_beta[apollo_fixed])
  model$estimate = temp[names(apollo_beta)]
  
  ### Save and print overview
  model$componentReport <- apollo_probabilities(model$estimate, apollo_inputs, functionality="report")
  if(!silent) for(r in model$componentReport) if(!is.null(r$param) && length(r$param)>0) for(j in r$param) cat(j, '\n', sep='')
  if(exists('r')) rm(r); if(exists('j')){cat('\n'); rm(j)}
  
  ### Calculate Zero LL
  if(!silent) apollo_print("Calculating LL(0)... ")
  model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
                        error=function(e) return(NA))
  if(is.list(model$LL0)) model$LL0=model$LL0[c("model",names(model$LL0)[names(model$LL0)!="model"])]
  if(is.list(model$LL0)){
    if(!workInLogs) model$LL0=sapply(model$LL0,function(x) sum(log(x)))
    if(workInLogs) model$LL0=sapply(model$LL0,sum)
  } else model$LL0 <- ifelse( workInLogs, sum(model$LL0), sum(log(model$LL0)) )
  
  ### Get LL at optimum for each component
  if(!silent) apollo_print("Calculating LL of each model component... ")
  LLout <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="output"),
                    error=function(e){ apollo_print("Could not complete validation using estimated parameters."); return(NA) }
  )
  if(!anyNA(LLout) && is.list(LLout)){
    LLout <- LLout[c("model",names(LLout)[names(LLout)!="model"])]
    if(!workInLogs) LLout <- lapply(LLout, log)
    LLout <-lapply(LLout,sum)
    model$Pout  <- LLout
    model$LLout <- unlist(LLout)
  } else{
    model$Pout  <- LLout
    model$LLout <- list(NA)
    if(!silent) apollo_print("LL could not be calculated for all components.")
  }
  
  ### use pre-scaling values as starting values in output
  model$bootstrapSE          <- bootstrapSE
  model$apollo_beta          <- temp_start
  model$LLStart              <- sum(testLL)
  model$startTime            <- time1
  model$apollo_control       <- apollo_control
  model$nObs                 <- nObs
  model$nIndivs              <- length(indiv)
  model$apollo_draws         <- apollo_inputs$apollo_draws
  model$estimationRoutine    <- estimationRoutine
  model$scaling              <- apollo_inputs$apollo_scaling
  
  ### De-scale parameters
  model$estimate[names(model$scaling)] <- model$estimate[names(model$scaling)]*model$scaling
  
  ### Fourth checkpoint
  time4 <- Sys.time()
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))
  model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))
  model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))
  
  return(model)
}
