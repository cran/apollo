#' Estimates model
#'
#' Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
#'
#' This is the main function of the Apollo package. The estimation process begins by running a number of checks on the 
#' \code{apollo_probabilities} function provided by the user.
#' If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
#' iterations. If Bayesian estimation is used, estimation will finish once the predefined number of iterations are completed.
#' By default, this functions writes the estimated parameter values in each iteration to a file in the working/output directory. Writing 
#' can be turned off by setting \code{estimate_settings$writeIter} to \code{FALSE}.
#' By default, \strong{final results are not written into a file nor printed to the console}, so users must make sure 
#' to call function \link{apollo_modelOutput} and/or \link{apollo_saveOutput} afterwards.
#' Users are strongly encouraged to visit \url{http://www.apollochoicemodelling.com/} to download examples on how to use the Apollo package.
#' The webpage also provides a detailed manual for the package, as well as a user-group to get further help.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param estimate_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                                 \itemize{
#'                                  \item \strong{\code{bgw_settings}}: List. Additional arguments to the BGW optimisation method. See \link[bgw]{bgw_mle} for more details.
#'                                  \item \strong{\code{bootstrapSE}}: Numeric. Number of bootstrap samples to calculate standard errors. Default is 0, meaning no bootstrap s.e. will be calculated. Number must zero or a positive integer. Only used if \code{apollo_control$estMethod!="HB"}.
#'                                  \item \strong{\code{bootstrapSeed}}: Numeric scalar (integer). Random number generator seed to generate the bootstrap samples. Only used if \code{bootstrapSE>0}. Default is 24.
#'                                  \item \strong{\code{constraints}}: Character vector. Linear constraints on parameters to estimate. For example \code{c('b1>0', 'b1 + 2*b2>1')}. Only \code{>}, \code{<} and \code{=} can be used. Inequalities cannot be mixed with equality constraints, e.g. \code{c(b1-b2=0, b2>0)} will fail. All parameter names must be on the left side. Fixed parameters cannot go into constraints. Alternatively, constraints can be defined as in \link[maxLik]{maxLik}. Constraints can only be used with maximum likelihood estimation and the BFGS routine in particular.
#'                                  \item \strong{\code{estimationRoutine}}: Character. Estimation method. Can take values "bfgs", "bgw", "bhhh", or "nr". Used only if \code{apollo_control$HB} is FALSE. Default is "bfgs".
#'                                  \item \strong{\code{hessianRoutine}}: Character. Name of routine used to calculate the Hessian of the log-likelihood function after estimation. Valid values are \code{"analytic"} (default), \code{"numDeriv"} (to use the numeric routine in package numDeric), \code{"maxLik"} (to use the numeric routine in packahe maxLik), and \code{"none"} to avoid calculating the Hessian and the covariance matrix. Only used if \code{apollo_control$HB=FALSE}.
#'                                  \item \strong{\code{maxIterations}}: Numeric. Maximum number of iterations of the estimation routine before stopping. Used only if \code{apollo_control$HB} is FALSE. Default is 200.
#'                                  \item \strong{\code{maxLik_settings}}: List. Additional settings for maxLik. See argument \code{control} in \link[maxLik]{maxBFGS}, \link[maxLik]{maxBHHH} and \link[maxLik]{maxNM} for more details. Only used for maximum likelihood estimation.
#'                                  \item \strong{\code{numDeriv_settings}}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                                  \item \strong{\code{printLevel}}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. Ignored if apollo_control$HB is TRUE. Default is 3.
#'                                  \item \strong{\code{scaleAfterConvergence}}: Logical. Used to increase numerical precision of convergence. If TRUE, parameters are scaled to 1 after convergence, and the estimation is repeated from this new starting values. Results are reported scaled back, so it is a transparent process for the user. Default is TRUE.
#'                                  \item \strong{\code{scaleHessian}}: Logical. If TRUE, parameters are scaled to 1 for Hessian estimation. Default is TRUE.
#'                                  \item \strong{\code{scaling}}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
#'                                  \item \strong{\code{silent}}: Logical. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                  \item \strong{\code{validateGrad}}: Logical. If TRUE, the analytical gradient (if used) is compared to the numerical one. Default is TRUE.
#'                                  \item \strong{\code{writeIter}}: Logical. Writes value of the parameters in each iteration to a csv file. Works only if \code{estimation_routine="bfgs"}. Default is TRUE.
#'                                 }
#' @return model object
#' @export
#' @importFrom numDeriv hessian grad
#' @importFrom maxLik maxLik
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov runif
#' @importFrom bgw bgw_mle
apollo_estimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
  # #################### #
  #### Loading Inputs ####
  # #################### #
  
  test <- is.vector(apollo_beta) && is.function(apollo_probabilities) && is.list(apollo_inputs)
  if(!test) stop('SYNTAX ISSUE - Arguments apollo_beta, apollo_fixed, apollo_probabilities ", 
                 "and apollo_inputs must be provided.')
  
  ### First checkpoint
  time1 <- Sys.time()
  
  ### Detach things if necessary
  apollo_detach()
  
  ### Set missing settings to default values
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="analytic", printLevel=3L, constraints=NULL, 
                  maxLik_settings=NULL, numDeriv_settings=list(), scaling=NA, 
                  bootstrapSE=0, bootstrapSeed=24, silent=FALSE, 
                  scaleHessian=TRUE, scaleAfterConvergence=TRUE,
                  validateGrad=TRUE,
                  bgw_settings=list())
  prtLvlMan <- is.list(estimate_settings) && !is.null(estimate_settings$printLevel)
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
  if(!prtLvlMan && tolower(estimate_settings$estimationRoutine)=="bgw") estimate_settings$printLevel <- 2
  rm(default, tmp, prtLvlMan)
  
  ### Set missing BGW settings to default values
  default <- list(maxIterations    = estimate_settings$maxIterations,
                  maxFunctionEvals = 3*estimate_settings$maxIterations,
                  silent           = estimate_settings$silent,
                  outputDirectory  = apollo_inputs$apollo_control$outputDirectory,
                  printLevel       = estimate_settings$printLevel, 
                  printNonDefaultSettings = FALSE, 
                  printStartingValues     = FALSE,
                  printFinalResults       = FALSE)
  if(length(estimate_settings$bgw_settings)==1 && is.na(estimate_settings$bgw_settings)) estimate_settings$bgw_settings <- default
  tmp <- names(default)[!(names(default) %in% names(estimate_settings$bgw_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings$bgw_settings[[i]] <- default[[i]]
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
  bootstrapSE       = estimate_settings[["bootstrapSE"]]
  bootstrapSeed     = estimate_settings[["bootstrapSeed"]]
  maxLik_settings   = estimate_settings[["maxLik_settings"]]
  scaleAfterConvergence = estimate_settings[['scaleAfterConvergence']]
  validateGrad = estimate_settings[['validateGrad']]
  bgw_settings      = estimate_settings[["bgw_settings"]]
  
  # Extract silent and debug
  apollo_inputs$silent <- estimate_settings$silent
  silent <- apollo_inputs$silent
  debug  <- apollo_inputs$apollo_control$debug
  
  
  # ########################## #
  #### Validation of inputs ####
  # ########################## #
  
  ### If maxIterations are zero, then do not do any validation
  if(maxIterations==0){
    apollo_control$noValidation=TRUE
    apollo_inputs$apollo_control$noValidation=TRUE
    estimate_settings$writeIter = FALSE
    writeIter = FALSE ## we're using both the one in estimate_settings and a local one, we should fix that
  }
  
  ### Validation of input
  apollo_checkArguments(apollo_probabilities, apollo_randCoeff, apollo_lcPars)
  if( !(estimationRoutine %in% c("bfgs", "bgw", "bhhh", "nr")) ) stop("SYNTAX ISSUE - Invalid estimationRoutine. Use 'bfgs', 'bgw', 'bhhh' or 'nr'.")
  if( !(estimate_settings$hessianRoutine %in% c('analytic', 'numDeriv', 'maxLik', 'none')) ) stop("SYNTAX ISSUE - Invalid hessianRoutine. Use 'analytic', 'numDeriv', 'maxLik' or 'none'.")
  
  # Check apollo_beta and apollo_fixed
  if(!is.numeric(apollo_beta) | !is.vector(apollo_beta) | is.null(names(apollo_beta))) stop("INPUT ISSUE - The \"apollo_beta\" argument needs to be a named vector")
  if(length(apollo_fixed)>0 && !is.character(apollo_fixed)) stop("INPUT ISSUE - 'apollo_fixed' is not an empty vector nor a vector of names.")
  if(length(unique(names(apollo_beta)))<length(apollo_beta)) stop("INPUT ISSUE - The \"apollo_beta\" argument contains duplicate elements")
  if(length(unique(apollo_fixed))<length(apollo_fixed)) stop("INPUT ISSUE - The \"apollo_fixed\" argument contains duplicate elements")
  if(!all(apollo_fixed %in% names(apollo_beta))) stop("INPUT ISSUE - Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'.")
  
  # Check bootstrap settings
  if(!is.numeric(bootstrapSE) || length(bootstrapSE)!=1 || bootstrapSE<0) stop("SYNTAX ISSUE - 'bootstrapSE' is not zero or a positive integer.")
  bootstrapSE <- as.integer(bootstrapSE); estimate_settings$bootstrapSE <- bootstrapSE
  if(!is.numeric(bootstrapSeed) || length(bootstrapSeed)!=1 || bootstrapSeed<=0) stop("SYNTAX ISSUE - 'bootstrapSeed' is not a positive integer.")
  bootstrapSeed <- as.integer(bootstrapSeed); estimate_settings$bootstrapSeed <- bootstrapSeed
  
  # Check printLevel
  if(!is.integer(printLevel)) printLevel <- as.integer(round(printLevel,0))
  if(printLevel<0L){ printLevel <- 0L; estimate_settings$printLevel <- 0L }
  if(3L<printLevel){ printLevel <- 3L; estimate_settings$printLevel <- 3L }
  
  # Check maxIterations
  if(maxIterations<0) stop("SYNTAX ISSUE - Need positive number of iterations!")
  maxIterations     = round(maxIterations,0)
  estimate_settings$maxIterations = maxIterations
  
  # create temporary copy of starting values for use later
  temp_start = apollo_beta
  
  # Check constraints
  if(!is.null(constraints) && apollo_control$HB) stop("INCORRECT FUNCTION/SETTING USE - Constraints cannot be used with Bayesian estimation.")
  if(!is.null(constraints) && estimationRoutine!="bfgs"){
    estimationRoutine               = "bfgs"
    apollo_inputs$estimationRoutine = estimationRoutine
    if(!silent) apollo_print("Estimation routine changed to 'BFGS'. Only 'BFGS' supports constrained optimization.", type="w")
  }
  if(is.vector(constraints) && is.character(constraints)){
    apollo_constraints <- constraints # copy them so that they can be stored in the model object later
    nCon <- length(constraints)
    bVar <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
    nVar <- length(bVar)
    bVar <- list2env(setNames(split(diag(nVar), rep(1:nVar,each=nVar)), bVar))
    bVal <- list2env(as.list(apollo_beta))
    A <- matrix(0, nrow=nCon, ncol=nVar, dimnames=list(NULL, names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]))
    b <- rep(0, nCon)
    mid0 <- ''
    for(i in 1:nCon){
      # turn constraint into expression, and split into left, right and middle
      e <- tryCatch(str2lang(constraints[i]), error=function(e) NULL)
      test <- is.null(e) || !is.call(e) || length(e)!=3
      if(test) stop('SYNTAX ISSUE - Constraint "', constraints[i], '" is not a valid linear constraint expression.')
      mid <- e[[1]]; lef <- e[[2]]; rig <- e[[3]]
      # Checks
      test <- is.symbol(mid) && (as.character(mid) %in% c(">", "=", "<"))
      if(!test) stop('SYNTAX ISSUE - Constraint "', constraints[i], '" does not contain one (and only one) of the following: >, <, or =.')
      test <- c(mid0, as.character(mid)); test <- mid0=="" | all(test %in% c("<", ">")) | all(test=="=")
      #test <- mid0=='' | as.character(mid)==mid0; mid0 <- as.character(mid)
      if(!test) stop('SYNTAX ISSUE - All constraints must be either equalities or inequealities, but not a mix of them.')
      mid0 <- as.character(mid)
      test <- length(all.vars(rig))==0
      if(!test) stop('SYNTAX ISSUE - The right side of constraint "', constraints[i],'" should only contain numeric values.')
      test <- all(all.vars(lef) %in% ls(bVar))
      if(!test) stop('INCORRECT FUNCTION/SETTING USE - All the variables in the left side of constraint "', constraints[i],
                     '" should be in apollo_beta. Fixed parameters cannot go into constraints.')
      if(as.character(mid)=='=') e[[1]] <- as.symbol('==')
      test <- eval(e, envir=bVal)
      if(!test) stop('SYNTAX ISSUE - Starting values of parameters do not satisfy constraint "', constraints[i],'".')
      # Fill A & b
      A[i,] <- eval(lef, envir=bVar)*ifelse(mid0=="<",-1, 1)
      b[i]  <- -eval(rig)*ifelse(mid0=="<",-1, 1)
    }
    #if(!all(apollo_inputs$apollo_scaling==1)) for(a in ls(bVar)) A[,a] <- A[,a]*apollo_inputs$apollo_scaling[a]
    if(mid0 %in% c(">", "<")) constraints <- list(ineqA=A, ineqB=b) else constraints <- list(eqA=A, eqB=b)
    rm(nCon, bVar, nVar, A, b, mid0, i, e, test, mid, lef, rig, bVal)
  }
  hasEqConst   <- length(constraints)==2 && all(names(constraints) %in% c('eqA', 'eqB'))
  hasIneqConst <- length(constraints)==2 && all(names(constraints) %in% c('ineqA', 'ineqB'))
  
  ### Recommend using multi-core if appropriate
  if(!apollo_control$HB && apollo_inputs$apollo_control$mixing && apollo_inputs$apollo_control$nCores==1 && !silent){
    n   <- nrow(apollo_inputs$database)
    tmp <- apollo_inputs$apollo_draws$interNDraws
    if(!is.null(tmp) && is.numeric(tmp) && tmp>0) n <- n*tmp
    tmp <- apollo_inputs$apollo_draws$interNDraws
    if(!is.null(tmp) && is.numeric(tmp) && tmp>0) n <- n*tmp
    tmp <- paste0('You can use multiple processor cores to speed up estimation. ',
                  'To do so, specify the desired number of cores using the setting nCores ',
                  'inside "apollo_control", for example "nCores=2". This computer has ',
                  parallel::detectCores(), ' available cores.')
    if(parallel::detectCores()>2) tmp <- paste0(tmp, ' We recommend using no more ',
                                                'than ', parallel::detectCores()-1, ' cores.')
    if(n>1e5) apollo_print(tmp, pause=5, type="i")
    rm(n, tmp)
  }
  
  ### Check if any of the fixed parameters are not zero
  if(length(apollo_fixed)>0){
    fixed_params=apollo_beta[apollo_fixed]
    if(any(!(fixed_params%in%c(0,1)))){
      tmp=names(fixed_params[!(fixed_params%in%c(0,1))])
      one=(length(tmp)==1)
      txt <- paste0('Element', ifelse(one,' ', 's '), paste0(tmp, collapse=', '), 
                    ' in \'apollo_fixed\' ', ifelse(one,'is ', 'are '),' constrained to ', ifelse(one,'a value ', 'values '),' other than zero or one. This may be intentional. ',
                    'If not, stop this function by pressing the "Escape" key and adjust the starting values accordingly.')
      apollo_print(txt, pause=5, type="w")
    }
  }

  ### Create useful variables
  nObs  <- nrow(apollo_inputs$database)
  indiv <- unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  nObsPerIndiv <- rep(0, length(indiv))
  for(n in 1:length(indiv)) nObsPerIndiv[n] <- sum(apollo_inputs$database[,apollo_inputs$apollo_control$indivID]==indiv[n])
  names(nObsPerIndiv) <- indiv
  
  
  # ########################################################## #
  #### Enhance user-defined functions and back-up originals ####
  # ########################################################## #
  
  ### Store unaltered version of apollo_probabilities, apollo_randCoeff and 
    # apollo_lcPars
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_lcPars
  
  ### Check and modify user-defined functions
  L <- apollo_modifyUserDefFunc(apollo_beta, apollo_fixed, apollo_probabilities,
                                apollo_inputs, validate=TRUE, 
                                noModification=apollo_control$noModification)
  
  ### Update functions if modification was successful, or fall back if not
  #if(L$success){
    apollo_probabilities           <- L$apollo_probabilities
    apollo_inputs$apollo_randCoeff <- L$apollo_randCoeff
    apollo_inputs$apollo_lcPars    <- L$apollo_lcPars
    apollo_inputs$apollo_scaling   <- L$apollo_scaling
    apollo_inputs$manualScaling    <- L$manualScaling
  #} else {
  # apollo_inputs$apollo_scaling   <- L$apollo_scaling
  # apollo_inputs$manualScaling    <- all(L$apollo_scaling==1)
  # estimate_settings$scaleHessian <- FALSE
  # estimate_settings$scaleAfterConvergence <- FALSE
  # scaleAfterConvergence          <- estimate_settings[['scaleAfterConvergence']]
  #}; 
  rm(L)
  
  ### If multicore is in use, save apollo_inputs to disk before changing it, 
    # and make sure to delete later
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
    # Scale starting parameters (this should always happen)
    if(length(apollo_inputs$apollo_scaling)>0 && !anyNA(apollo_inputs$apollo_scaling)){
      r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
      r <- names(apollo_beta)[r]
      apollo_beta[r] <- apollo_beta[r]/apollo_inputs$apollo_scaling[r]
    }
    
    preHBtime = as.numeric(difftime(Sys.time(),time1,units='secs'))
    model <- apollo_estimateHB(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings)
    model$timeTaken <- as.numeric(model$timeTaken+preHBtime)
    model$timePre   <- as.numeric(model$timePre+preHBtime)
    model$apollo_probabilities <- apollo_probabilities_ORIG
    return(model)
  }
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  
  # Scale starting parameters (this should always happen)
  if(length(apollo_inputs$apollo_scaling)>0 && !anyNA(apollo_inputs$apollo_scaling)){
    r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- apollo_beta[r]/apollo_inputs$apollo_scaling[r]
  }
  
  ### Validate likelihood function
  test <- apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$mixing
  if(test && !silent) apollo_print(paste0('CAUTION: In models using mixing and workInLogs=TRUE, Apollo assumes (and does not check) that ',
                                          'the call to apollo_panelProd is immediately followed by a call to apollo_avgInterDraws.'), type="w")
  # Validate LL at starting values
  if(!silent) apollo_print(c("\n", "Testing likelihood function..."))
  testLL <- apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")
  test <- is.list(testLL) && !is.null(names(testLL)) && 'model' %in% names(testLL) && is.numeric(testLL[['model']])
  if(!test) stop('CALCULATION ISSUE - Log-likelihood calculation fails!')
  testLL <- testLL[["model"]]
  if(!apollo_inputs$apollo_control$workInLogs) testLL <- log(testLL)
  if(!silent & any(!is.finite(testLL))){
    cat("\n")
    apollo_print("Log-likelihood calculation fails at starting values!")
    apollo_print("Affected individuals:")
    LLtemp <- data.frame(ID=indiv, LL=testLL)
    LLtemp <- subset(LLtemp,!is.finite(LLtemp[,2]))
    colnames(LLtemp) <- c("ID","LL")
    print(LLtemp, row.names=FALSE)
  }
  if(any(!is.finite(testLL))) stop('CALCULATION ISSUE - Log-likelihood calculation fails at values close to the starting values!')
  
  
  # ################################## #
  #### Pre-process                  ####
  # ################################## #
  if(!silent) apollo_print(c("\n","Pre-processing likelihood function..."))
  
  ### Create multi-core version of likelihood (if needed)
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                       apollo_probabilities, apollo_inputs, 
                                       apollo_estSet=estimate_settings, cleanMemory=TRUE)
  on.exit({
    test <- apollo_control$nCores>1 && exists('apollo_logLike', inherits=FALSE)
    test <- test && !anyNA(environment(apollo_logLike)$cl)
    if(test) parallel::stopCluster(environment(apollo_logLike)$cl)
  }, add=TRUE)
  
  ### Create gradient function if required (and possible)
  if(!is.null(apollo_inputs$apollo_control$analyticGrad) && apollo_inputs$apollo_control$analyticGrad){
    # apollo_makeGrad will create gradient function ONLY IF all components have analytical gradient.
    grad <- apollo_makeGrad(apollo_beta, apollo_fixed, apollo_logLike, validateGrad)
    if(is.null(grad)){
      if(!silent) apollo_print(paste0("Analytical gradients could not be calculated for all components, ", 
                                    "numerical gradients will be used."))
      apollo_inputs$apollo_control$analyticGrad <- FALSE # Not sure this is necessary, but it can't hurt
      environment(apollo_logLike)$analyticGrad  <- FALSE # fix for iteration counting
      environment(apollo_logLike)$nIter         <- 1     # fix for iteration counting
    }
  } else grad <- NULL
  
  ### Extract useful values from apollo_logLike
  nObsTot       <- environment(apollo_logLike)$nObsTot
  modelTypeList <- environment(apollo_logLike)$mType
  countAlt      <- environment(apollo_logLike)$countAlt
  
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
      }; stop('CALCULATION ISSUE - Log-likelihood calculation fails at values close to the starting values!')
    }
    base_LL = sum(base_LL)
    if(!is.null(grad)){
      test_gradient=colSums(grad(beta1))
      if(any(test_gradient==0)){
        tmp=names(beta1)[which(test_gradient==0)]
        if(length(tmp)==1){
          stop("SPECIFICATION ISSUE - Parameter ",tmp," does not influence the log-likelihood of your model!")
        }else{
          stop("SPECIFICATION ISSUE - Parameters ",paste0(tmp,collapse=", ")," do not influence the log-likelihood of your model!")
        }
      }
    }else{
    for(p in names(beta0)){
      #beta1p <- beta1 - (names(beta1)==p)*0.001
      beta1m <- beta1 + (names(beta1)==p)*0.001
      #test1_LL = sum( apollo_logLike(beta1p) )
      test2_LL = sum( apollo_logLike(beta1m) )
      #if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      #if(base_LL==test1_LL & base_LL==test2_LL) stop("SPECIFICATION ISSUE - Parameter ",p," does not influence the log-likelihood of your model!")
      if(base_LL==test2_LL) stop("SPECIFICATION ISSUE - Parameter ",p," does not influence the log-likelihood of your model!")
      if(!silent) cat(".")
    }
    rm(beta0, beta1, beta1m, base_LL, test2_LL, p)#, beta1p, test1_LL)
    }
  }
  
  
  # ################################## #
  #### Classical main estimation    ####
  # ################################## #
  
  ### Second checkpoint
  time2 <- Sys.time()
  
  ### Split parameters between variable and fixed
  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  beta_fix_val <- apollo_beta[apollo_fixed]
  
  ### Preparations for iteration writing and counting
  if(estimate_settings$writeIter){
    # Remove modelName_iterations if it exists
    tmp <- paste0(apollo_inputs$apollo_control$outputDirectory,apollo_inputs$apollo_control$modelName, "_iterations.csv")
    txt <- paste0('Could not delete old ', tmp, ' file. New iterations will be written after old ones.')
    if(file.exists(tmp)) tryCatch(file.remove(tmp), error=function(e) apollo_print(txt))
    if(is.function(grad)) apollo_writeTheta(apollo_beta, sum(testLL), apollo_inputs$apollo_control$modelName)
    rm(tmp, txt)
  }
  
  ### Scale constraints if necessary
  if(hasEqConst) for(a in colnames(constraints$eqA)) constraints$eqA[,a] <- constraints$eqA[,a]*apollo_inputs$apollo_scaling[a]
  if(hasIneqConst) for(a in colnames(constraints$ineqA)) constraints$ineqA[,a] <- constraints$ineqA[,a]*apollo_inputs$apollo_scaling[a]
  
  ### Main (classical) estimation
  if(!silent) apollo_print(c("\n", "Starting main estimation")) else maxLik_settings$printLevel=0
  hasEqConst   <- length(constraints)==2 && all(names(constraints) %in% c('eqA', 'eqB'))
  if(hasEqConst & !silent) apollo_print("Your model uses equality constraints, no estimation progress will be printed to screen.", type="i")
  ## Check to see if panelProd is to be overridden
  #test <- is.na(apollo_inputs$apollo_lcPars) & 
  #        is.na(apollo_inputs$apollo_draws) & 
  #        is.na(apollo_inputs$apollo_randCoeff) & 
  #        !(apollo_inputs$apollo_control$preventOverridePanel)
  test <- !is.function(apollo_inputs$apollo_lcPars) && 
    !apollo_inputs$apollo_control$mixing && 
    !apollo_inputs$apollo_control$preventOverridePanel # BUG FIX (DP 10/11/2022)
  if(test){
    if(environment(apollo_logLike)$singleCore){
      environment(apollo_logLike)$apollo_inputs$apollo_control$overridePanel <- TRUE
    } else parallel::clusterEvalQ(cl=environment(apollo_logLike)$cl,
                                  apollo_inputs$apollo_control$overridePanel <- TRUE)
  }

    
  ##
  if(estimationRoutine=="bgw"){
    f <- function(b) apollo_logLike(b, logP=FALSE, countIter=FALSE, writeIter=FALSE)
    model <- bgw::bgw_mle(calcR=f, betaStart=beta_var_val, calcJ=grad, bgw_settings=bgw_settings)
    model$nIter <- model$iterations
    rm(f)
  }else{
    model <- maxLik::maxLik(apollo_logLike, grad=grad, start=beta_var_val,
                          method=estimationRoutine, finalHessian=FALSE,
                          control=maxLik_settings,
                          constraints=constraints,
                          countIter=TRUE, writeIter=writeIter, sumLL=FALSE, getNIter=FALSE)
    ### Store number of iterations from first (usually unscaled) optimisation
    model$nIter <- environment(apollo_logLike)$nIter
  }
  
  # Checks if main estimation was successful
  successfulEstimation <- FALSE
  if(exists("model")){
    if(estimationRoutine=="bfgs" & model$code==0) successfulEstimation <- TRUE
    if(estimationRoutine=="bgw"){
      apollo_print(model$message)
      if(model$code %in% c(3,4,5,6)) successfulEstimation <- TRUE
    }
    if(estimationRoutine=="bhhh" & (model$code %in% c(2,8)) ) successfulEstimation <- TRUE
    if(estimationRoutine=="nr" && model$code<=2) successfulEstimation <- TRUE
    if(is.null(model$maximum) || !is.finite(model$maximum)){
      successfulEstimation <- FALSE
      model$message       <- "Not converged"
      apollo_print(paste0("The estimation has led to a non-finite log-likelihood. ",
                          "Although the R optimiser indicates convergence, Apollo ",
                          "does not consider estimation to have been successful."), type="w")
    } 
    model$successfulEstimation <- successfulEstimation
    model$message=trimws(gsub("*","",model$message,fixed=TRUE))
  }
  
  # ################################## #
  #### Classical scaled estimation  ####
  # ################################## #
  
  ### Second checkpoint, for scaled estimation
  time2b <- Sys.time()

  ### Estimation using scaling at final estimates (optional convergence test)
  if(successfulEstimation && scaleAfterConvergence && all(model$estimate!=0)){
    if(!silent) apollo_print(paste0('Additional convergence test using scaled estimation. Parameters will be scaled by ', 
                                    'their current estimates and additional iterations will be performed.'))
    # De-scale non-fixed parameters
    b <- model$estimate
    b[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling*b[names(apollo_inputs$apollo_scaling)]
    # De-scale constraints if defined
    hasEqConst   <- length(constraints)==2 && all(names(constraints) %in% c('eqA', 'eqB'))
    hasIneqConst <- length(constraints)==2 && all(names(constraints) %in% c('ineqA', 'ineqB'))
    if(hasEqConst) for(a in colnames(constraints$eqA)) constraints$eqA[,a] <- constraints$eqA[,a]/apollo_inputs$apollo_scaling[a]
    if(hasIneqConst) for(a in colnames(constraints$ineqA)) constraints$ineqA[,a] <- constraints$ineqA[,a]/apollo_inputs$apollo_scaling[a]
    # Update scaling inside apollo_logLike and in this environment
    apollo_scaling_backup=apollo_inputs$apollo_scaling
    apollo_inputs$apollo_scaling <- abs(b)
    if(apollo_control$nCores==1) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- abs(b) else {
      ns <- abs(b)
      parallel::clusterExport(environment(apollo_logLike)$cl, 'ns', envir=environment())
      parallel::clusterEvalQ(environment(apollo_logLike)$cl, {apollo_inputs$apollo_scaling <- ns; rm(ns)})
    }
    # Update scaling of constraints, if needed
    if(hasEqConst) for(a in colnames(constraints$eqA)) constraints$eqA[,a] <- constraints$eqA[,a]*apollo_inputs$apollo_scaling[a]
    if(hasIneqConst) for(a in colnames(constraints$ineqA)) constraints$ineqA[,a] <- constraints$ineqA[,a]*apollo_inputs$apollo_scaling[a]
    # Check that new starting LL is the same than convergence one
    test1 <- model$maximum
    test2 <- tryCatch(
      apollo_logLike(b/abs(b), countIter=FALSE, sumLL=TRUE, writeIter=FALSE, getNIter=FALSE),
      error=function(e) NA)
    # Estimate again only if new starting LL matches maximum LL
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(test){
      if(debug) apollo_print("LL before and after re-scaling are the same, proceeding with re-scaled estimation.")
      model_prescaling=model
      if(!is.null(maxLik_settings$printLevel) && maxLik_settings$printLevel==3) maxLik_settings$printLevel <- 2
      if(estimationRoutine=="bgw"){
        f <- function(b) apollo_logLike(b, logP=FALSE, countIter=FALSE, writeIter=FALSE)
        model <- tryCatch(bgw::bgw_mle(calcR=f, betaStart=b/abs(b), calcJ=grad, bgw_settings=bgw_settings), 
                          error=function(e) NA)
        if(is.list(model)) model$nIter <- model$iterations else rm(model)
        rm(f)
      } else {
        model <- tryCatch(maxLik::maxLik(apollo_logLike, grad=grad, start=b/abs(b),
                                method=estimationRoutine, finalHessian=FALSE,
                                control=maxLik_settings,
                                constraints=constraints,
                                countIter=TRUE, writeIter=writeIter, sumLL=FALSE, 
                                getNIter=FALSE), error=function(e) NA)
        ### Store number of iterations from post-scaling estimation
        if(is.list(model)) model$nIter <- environment(apollo_logLike)$nIter - model_prescaling$nIter else rm(model)
      }
      # Checks if estimation with re-scaled parameters was successful
      successfulEstimation_scaled <- FALSE
      if(exists("model")){
        if(estimationRoutine=="bfgs" & model$code==0) successfulEstimation_scaled <- TRUE
        if(estimationRoutine=="bgw"){
          apollo_print(model$message)
          if(model$code %in% c(3,4,5,6)){
            successfulEstimation_scaled <- TRUE
          }
        }
        if(estimationRoutine=="bhhh" & (model$code %in% c(2,8)) ) successfulEstimation_scaled <- TRUE
        if(estimationRoutine=="nr" && model$code<=2) successfulEstimation_scaled <- TRUE
        model$message=trimws(gsub("*","",model$message,fixed=TRUE))
      }
      if(!successfulEstimation_scaled){
        apollo_print("The estimation of the scaled model failed, and the unscaled version will be returned instead.")
        model <- model_prescaling
        apollo_inputs$apollo_scaling=apollo_scaling_backup
        if(apollo_control$nCores==1) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- apollo_scaling_backup else {
          ns <- apollo_scaling_backup
          parallel::clusterExport(environment(apollo_logLike)$cl, 'ns', envir=environment())
          parallel::clusterEvalQ(environment(apollo_logLike)$cl, {apollo_inputs$apollo_scaling <- ns; rm(ns)})
          rm(ns)
        }
      }
    } else {
      apollo_print("The scaling of the model led to differences in the results, and was thus not used. This is unexpected. You may want to contact the developers.", pause=5, type="w")
      apollo_inputs$apollo_scaling=apollo_scaling_backup
      if(apollo_control$nCores==1) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- apollo_scaling_backup else {
        ns <- apollo_scaling_backup
        parallel::clusterExport(environment(apollo_logLike)$cl, 'ns', envir=environment())
        parallel::clusterEvalQ(environment(apollo_logLike)$cl, {apollo_inputs$apollo_scaling <- ns; rm(ns)})
        rm(ns)
      }
    }
  }
  
  ### Additional checkpoint, for after scaled estimation
  time2c <- Sys.time()
  
  ### reinstate panelProd
  if(environment(apollo_logLike)$singleCore){
    environment(apollo_logLike)$apollo_inputs$apollo_control$overridePanel=FALSE
  } else {
    parallel::clusterEvalQ(environment(apollo_logLike)$cl,
                           apollo_inputs$apollo_control$overridePanel <- FALSE)
  }
  
  ### temp BHHH
  ### anything marked ##8May## was commented out with now directly descaling the BHHH matrix
  if(successfulEstimation){
    #if(!silent) apollo_print('Computing score matrix...')
    # descale
    ##8May## btemp <- model$estimate[names(apollo_inputs$apollo_scaling)]*apollo_inputs$apollo_scaling
    ##8May## apollo_scaling_backup=apollo_inputs$apollo_scaling
    ##8May## apollo_inputs$apollo_scaling <- setNames(rep(1, length(btemp)), names(btemp))
    ##8May## if(apollo_control$nCores==1) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- setNames(rep(1, length(btemp)), names(btemp)) else {
    ##8May##   ns <- setNames(rep(1, length(btemp)), names(btemp))
    ##8May##   parallel::clusterExport(environment(apollo_logLike)$cl, 'ns', envir=environment())
    ##8May##   parallel::clusterEvalQ(environment(apollo_logLike)$cl, {apollo_inputs$apollo_scaling <- ns; rm(ns)})
    ##8May## }
    ##8May## #if(apollo_inputs$apollo_control$analyticGrad=='TRUE'){
    ##8May## if(!is.null(grad)){
    ##8May##   score <- grad(btemp)
    ##8May## } else { 
    ##8May##   score <- numDeriv::jacobian(apollo_logLike, btemp)
    ##8May## }
    ##8May## colnames(score)=names(btemp)
    ##8May## model$score=score
    ##8May## if(is.matrix(score)) model$BHHH_matrix   <- var(score)*nrow(score)
    ##8May## # restore scaling
    ##8May## apollo_inputs$apollo_scaling=apollo_scaling_backup
    ##8May## if(apollo_control$nCores==1) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- apollo_scaling_backup else {
    ##8May##   ns <- apollo_scaling_backup
    ##8May##   parallel::clusterExport(environment(apollo_logLike)$cl, 'ns', envir=environment())
    ##8May##   parallel::clusterEvalQ(environment(apollo_logLike)$cl, {apollo_inputs$apollo_scaling <- ns; rm(ns)})
    ##8May## }
    
    ##8May## NEW VERSION
    if(!is.null(model$varcovBGW)){
      model$BHHH_matrix=solve(model$varcovBGW)
    } else{
      btemp <- model$estimate
      if(!is.null(grad)){
        score <- grad(btemp)
      } else { 
        score <- numDeriv::jacobian(apollo_logLike, btemp)
      }
      colnames(score)=names(btemp)
      if(is.matrix(score)) model$BHHH_matrix   <- var(score)*nrow(score)
      rm(score)
    }
    # descale BHHH matrix if needed
    if(any(apollo_inputs$apollo_scaling!=1)){
      for(j in names(apollo_inputs$apollo_scaling)){
        scale=apollo_inputs$apollo_scaling[j]
        model$BHHH_matrix[j,]=model$BHHH_matrix[j,]/scale
        model$BHHH_matrix[,j]=model$BHHH_matrix[,j]/scale
      }
      if(exists("scale", inherits=FALSE)) rm(scale)
    }
    ##8May##
    model$BHHHvarcov  <- tryCatch(solve(model$BHHH_matrix), error=function(e) return(matrix(NA, nrow=length(btemp), ncol=length(btemp), 
                                                                                            dimnames=list(names(btemp), names(btemp)))))
    model$BHHHse      <- sqrt(diag(model$BHHHvarcov))
    model$BHHHcorrmat <- tryCatch(model$BHHHvarcov/(model$BHHHse%*%t(model$BHHHse)), error=function(e) return(matrix(NA, nrow=length(btemp), ncol=length(btemp), 
                                                                                                                     dimnames=list(names(btemp), names(btemp)))))
  
    ### add in fixed params    
    model$BHHHse    <- c(model$BHHHse, apollo_beta[apollo_fixed]*NA)[names(apollo_beta)]
  }
  ###
  
  ### Print estimated parameters
  if(exists("model") & !silent){
    cat("\n")
    if(!is.null(model$estimate)){
      tmp <- c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
      tmp[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling*tmp[names(apollo_inputs$apollo_scaling)]
      if(is.null(model$BHHHse)||all(is.na(model$BHHHse))){
        apollo_print("Estimated parameters:")
        tmp <- matrix(tmp, nrow=length(tmp), ncol=1, dimnames=list(names(tmp), 'Estimate'))
      }else{
        apollo_print("Estimated parameters with approximate standard errors from BHHH matrix:")
        tmp <- matrix(cbind(tmp,model$BHHHse,tmp/model$BHHHse), nrow=length(tmp), ncol=3, dimnames=list(names(tmp), c('Estimate','BHHH se','BHH t-ratio')))
      }
      apollo_print(tmp); rm(tmp)
    } else apollo_print("Estimated parameters: Not available.")
    apollo_print("\n")
    if(!is.null(model$maximum)){
      apollo_print(paste0("Final LL: " , round(model$maximum, 4)))
    } else apollo_print("Final LL: Not available.")
    apollo_print("\n")
  }
  
  ### If estimation failed, continue only if model exists
  if(!successfulEstimation){
    if(exists("model")){
      if(!is.null(model$message) && grepl("Function evaluation limit", model$message)){
        if(!silent) apollo_print("Function evaluation limit exceeded. No covariance matrix will be computed. You may wish to use the current estimates as starting values for a new estimation with a higher limit for functional evaluations.",  pause=3, type="w")
      }else if(!is.null(model$message) && grepl("limit", model$message)){
        if(!silent) apollo_print("Iteration limit exceeded. No covariance matrix will be computed. You may wish to use the current estimates as starting values for a new estimation with a higher limit for iterations.",  pause=3, type="w")
      }else{
        if(!silent) apollo_print("Estimation failed. No covariance matrix to compute.",  pause=3, type="w")
      }
      if(estimationRoutine %in% c("bfgs", "bhhh") && model$code==100){
        stop("CALCULATION ISSUE - Estimation failed, no estimated model to return.")
      }
      estimate_settings$hessianRoutine <- 'none'
    } else stop("CALCULATION ISSUE - Estimation failed, no estimated model to return.")
  }
  

  
  ### Third checkpoint
  time3 <- Sys.time()
  
  # ######################################## #
  #### Fit indices and packing up "model" ####
  # ######################################## #
  
  ### Calculate probabilities to identify outliers
  indLL <- apollo_logLike(model$estimate)
  P     <- exp(indLL)
  if(apollo_control$panelData){
    model$avgLL <- setNames(indLL/nObsPerIndiv, names(nObsPerIndiv))
    model$avgCP <- setNames(P^(1/nObsPerIndiv), names(nObsPerIndiv))
  } else {
    model$avgLL <- setNames(indLL, indiv)
    model$avgCP <- setNames(P, indiv)
  }
  
  if(exists("model_prescaling")){
    model$nIterPrescaling <- model_prescaling$nIter
    model$nIterPostscaling <- model$nIter
  } else {
    model$nIterPrescaling <- model$nIter
    model$nIterPostscaling <- 0
  }
  model$nIter <- model$nIterPrescaling + model$nIterPostscaling

  if(exists("model_prescaling")){
    model$model_prescaling <- model_prescaling
  } else{
    model$model_prescaling <- NULL
  }
  
  ### Restore fixed parameters to model$estimate
  temp           = c(model$estimate, apollo_beta[apollo_fixed])
  model$estimate = temp[names(apollo_beta)]
  
  ### Restore apollo_inputs (with draws and database, ~DOUBLING OF MEMORY)
  if(apollo_control$nCores>1 & is.null(apollo_inputs$database)){
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs")
    ns <- apollo_inputs$apollo_scaling
    apollo_inputs <- tryCatch(readRDS(file=fileName), error=function(e) NULL)
    apollo_inputs$apollo_scaling <- ns
    rm(ns, fileName)
  }
  
  ### Save and print overview
  model$componentReport <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="report"), 
                                    error=function(e) NULL)
  test <- is.list(model$componentReport) && !is.null(names(model$componentReport)) && any(c('data', 'param') %in% names(model$componentReport))
  if(test) model$componentReport <- list(model=model$componentReport)
  if(!silent){
    test <- FALSE
    for(r in model$componentReport) if(is.list(r) && !is.null(r$param) && length(r$param)>0){
      test <- TRUE
      for(j in r$param) cat(j, '\n', sep='')
    }; rm(r)
    if(test) cat('\n')
  }
  
  ### Calculate Zero LL
  if(!silent) apollo_print("Calculating log-likelihood at equal shares (LL(0)) for applicable models... ")
  model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
                        error=function(e) return(NA))
  if(is.list(model$LL0)){
    model$LL0=model$LL0[c("model",names(model$LL0)[names(model$LL0)!="model"])]
    for(s in 1:length(model$LL0)){
      if(is.list(model$LL0[[s]])) model$LL0[[s]]=model$LL0[[s]][["model"]]
    }
    if(!workInLogs) model$LL0=sapply(model$LL0,function(x) sum(log(x)))
    if( workInLogs) model$LL0=sapply(model$LL0,sum)
  } else {
    model$LL0 <- ifelse( workInLogs, sum(model$LL0), sum(log(model$LL0)) )
  }
  
  ### Calculate shares LL (constants only)
  if(apollo_control$calculateLLC){
    if(!silent) apollo_print("Calculating log-likelihood at observed shares from estimation data (LL(c)) for applicable models... ")
    model$LLC <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="shares_LL"),
                          error=function(e) return(NA))
    if(is.list(model$LLC)){
      model$LLC = model$LLC[c("model",names(model$LLC)[names(model$LLC)!="model"])]
      for(s in 1:length(model$LLC)) if(is.list(model$LLC[[s]])) model$LLC[[s]] = model$LLC[[s]][["model"]]
      if(workInLogs) model$LLC <- sapply(model$LLC,sum) else model$LLC <- sapply(model$LLC,function(x) sum(log(x)))
    } else model$LLC <- ifelse(workInLogs, sum(model$LLC), sum(log(model$LLC)) )
  } else model$LLC=NA
  
  ### Get LL at optimum for each component
  if(!silent) apollo_print("Calculating LL of each model component... ")
  LLout <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="output"),
                    error=function(e){ apollo_print("Could not complete validation using estimated parameters."); return(NA) }
  )
  if(!anyNA(LLout) && is.list(LLout)){
    LLout <- LLout[c("model",names(LLout)[names(LLout)!="model"])]
    for(s in 1:length(LLout)) if(is.list(LLout[[s]])) LLout[[s]]=LLout[[s]][["model"]]
    if(!workInLogs) LLout <- lapply(LLout, log)
    LLout <-lapply(LLout,sum)
    model$LLout <- unlist(LLout)
  } else{
    model$LLout <- list(NA)
    if(!silent) apollo_print("LL could not be calculated for all components.")
  }
  
  
  ### Calculate Rho2, AIC, and BIC
  nFreeParams <- length(apollo_beta) - length(apollo_fixed)
  test <- !is.null(modelTypeList) && !anyNA(model$LL0[1])
  test <- test && all(tolower(modelTypeList) %in% c("mnl", "nl", "cnl", "el", "dft", "lc", "rrm", "ol", "op"))
  if(test & !silent) apollo_print("Calculating other model fit measures")
  if(test) model$rho2_0 <- 1-(model$maximum/model$LL0[1]) else model$rho2_0 <- NA
  if(test) model$adjRho2_0 <- 1-((model$maximum-nFreeParams)/model$LL0[1]) else model$adjRho2_0 <- NA
  test <- test && is.numeric(model$LLC) && !anyNA(model$LLC[1])
  if(test) model$rho2_C <- 1-(model$maximum/model$LLC[1]) else model$rho2_C <- NA
  #if(test) model$adjRho2_C <- 1-((model$maximum-nFreeParams)/model$LLC[1]) else model$adjRho2_C <- NA
  if(test) model$adjRho2_C <- 1-((model$maximum-nFreeParams+sum(countAlt-1))/model$LLC[1]) else model$adjRho2_C <- NA
  
  
  model$AIC <- -2*model$maximum + 2*nFreeParams
  model$BIC <- -2*model$maximum + nFreeParams*log(ifelse(!is.null(nObsTot), nObsTot, nObs))
  model$nFreeParams <- nFreeParams
  
  ### Store functions and other relevant things inside "model"
  model$apollo_randCoeff     <- apollo_randCoeff_ORIG
  model$apollo_lcPars        <- apollo_lcPars_ORIG
  model$apollo_probabilities <- apollo_probabilities_ORIG
  model$apollo_fixed         <- apollo_fixed
  model$modelTypeList        <- environment(apollo_logLike)$mType
  model$nObsTot              <- environment(apollo_logLike)$nObsTot
  model$apollo_probabilities_mod <- apollo_probabilities
  environment(model$apollo_probabilities_mod) <- baseenv()
  model$apollo_lcPars_mod <- apollo_inputs$apollo_lcPars
  if(is.function(model$apollo_lcPars_mod)) environment(model$apollo_lcPars_mod) <- baseenv()
  model$apollo_randCoeff_mod <- apollo_inputs$apollo_randCoeff
  if(is.function(model$apollo_randCoeff_mod)) environment(model$apollo_randCoeff_mod) <- baseenv()
  model$se                <- setNames(rep(NA, length(apollo_beta)), names(apollo_beta))
  model$robse             <- model$se
  tmp                     <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)]
  model$varcov            <- matrix(NA, nrow=length(tmp), ncol=length(tmp), dimnames=list(tmp, tmp))
  model$corrmat           <- model$varcov
  model$robvarcov         <- model$varcov
  model$bootstrapSE       <- bootstrapSE
  model$apollo_beta       <- temp_start
  model$LLStart           <- sum(testLL)
  model$startTime         <- time1
  model$apollo_control    <- apollo_control
  model$nObs              <- nObs
  model$nIndivs           <- length(indiv)
  model$apollo_draws      <- apollo_inputs$apollo_draws
  model$estimationRoutine <- estimationRoutine
  model$scaling           <- apollo_inputs$apollo_scaling
  model$manualScaling     <- apollo_inputs$manualScaling
  model$timePre   <- as.numeric(difftime(time2, time1, units='secs'))
  model$timeEst   <- as.numeric(difftime(time3, time2, units='secs'))
  model$timeEstPrescaling <- as.numeric(difftime(time2b,time2,units='secs'))+as.numeric(difftime(time3,time2c,units='secs'))
  model$timeEstPostscaling <- as.numeric(difftime(time2c,time2b,units='secs'))
  model$timeTaken <- as.numeric(difftime(Sys.time(),time1,units='secs'))
  model$timePost  <- NA
  rm(tmp)
  
  ### Save constraints (NULL if none or given in maxLik format)
  test <- exists("apollo_constraints", envir=environment(), inherits=FALSE)
  if(test) model$apollo_constraints <- apollo_constraints else model$apollo_constraints <- NULL
  rm(test)
  
  ### If getting to this point has taken more than 10 minutes, write model object to disk
  test <- model$timeTaken >= 600
  if(test){
    modelName <- paste0(apollo_inputs$apollo_control$outputDirectory, 
                        apollo_control$modelName)
    fileName  <- paste0(modelName, "_model.rds")
    # Rename existing _model file to OLDn_model
    if(file.exists(fileName)){
      n <- 1
      while( file.exists( paste0(modelName, "_OLD", n, "_model.rds") ) ) n <- n + 1
      fileNameOld <- paste0(modelName, "_OLD", n, "_model.rds")
      if(file.exists(fileName)) file.rename(from=fileName, to=fileNameOld)
      rm(n, fileNameOld)
    }
    ### De-scale parameters
    model2 <- model
    model2$estimate[names(model2$scaling)] <- model2$estimate[names(model2$scaling)]*model2$scaling
    # Try writing it
    wrote <- tryCatch({saveRDS(model2, file=fileName); TRUE}, error=function(e) FALSE)
    if(wrote){
      txt <- paste0("Your model took more than 10 minutes to estimate, so it ", 
                    "was saved to file ", fileName, " before calculating its ", 
                    "covariance matrix. If calculation of the covariance ", 
                    "matrix fails or is stopped before finishing, you can ", 
                    "load the model up to this point using apollo_loadModel.")
      if(!is.null(model$BHHHse)) txt <- paste0(txt," You may also want to inspect the approximate BHHH standard errors shown above to determine whether you wish to continue this process.")
    } else txt <- paste0("Intermediate results of your model (up to this ", 
                         "point) could not be saved to file ", fileName, ".")
    if(!silent) apollo_print(txt, type="i")
    rm(modelName, fileName, wrote, txt, model2)
  }
  
  
  
  # ############################ #
  #### Hessian & covar matrix ####
  # ############################ #
  
  ### Create a de-scaled copy of model$estimate
  b <- model$estimate # c(model$estimate, apollo_beta[apollo_fixed])[names(apollo_beta)]
  b[names(apollo_inputs$apollo_scaling)] <- b[names(apollo_inputs$apollo_scaling)]*apollo_inputs$apollo_scaling
  
  ### Calculate varcov
  if((estimationRoutine!="bgw")&&!is.null(model$BHHH_matrix)){
    varcov <- apollo_varcov(apollo_beta=b, apollo_fixed, 
                            varcov_settings=list(BHHH_matrix       = model$BHHH_matrix,
                                                 hessianRoutine    = estimate_settings$hessianRoutine,
                                                 scaleBeta         = estimate_settings$scaleHessian, 
                                                 numDeriv_settings = estimate_settings$numDeriv_settings,
                                                 apollo_logLike    = apollo_logLike,
                                                 apollo_grad       = grad))
  }else{
    varcov <- apollo_varcov(apollo_beta=b, apollo_fixed, 
                            varcov_settings=list(hessianRoutine    = estimate_settings$hessianRoutine,
                                                 scaleBeta         = estimate_settings$scaleHessian, 
                                                 numDeriv_settings = estimate_settings$numDeriv_settings,
                                                 apollo_logLike    = apollo_logLike,
                                                 apollo_grad       = grad))
  }
  if(is.list(varcov)) for(i in names(varcov)) if(i!="apollo_beta") model[[i]] <- varcov[[i]]
  rm(varcov)
  
  ### Close cluster and delete apollo_logLike
  if(!environment(apollo_logLike)$singleCore){
    parallel::stopCluster(environment(apollo_logLike)$cl)
  }; rm(apollo_logLike)
  
  ### Restore apollo_inputs including database and draws into global environment
  test <- apollo_control$nCores>1
  test <- test && exists("apollo_inputs", envir=.GlobalEnv)
  test <- test && is.list(.GlobalEnv$apollo_inputs)
  test <- test && is.null(.GlobalEnv$apollo_inputs$database)
  if(test){
    ns <- apollo_inputs$apollo_scaling
    if(debug) cat('Restoring data to main thread...')
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName,"_inputs")
    apollo_inputs <- tryCatch(readRDS(file=fileName), error=function(e) NULL)
    if(!is.null(apollo_inputs)){
      tmp <- .GlobalEnv
      assign("apollo_inputs", apollo_inputs, envir=tmp)
      unlink(fileName)
      # copy original randCoeff and lcPars to apollo_inputs in global environment
      assign("fTemp", list(rc=apollo_randCoeff_ORIG, lp=apollo_lcPars_ORIG), envir=tmp)
      eval(quote(apollo_inputs$apollo_randCoeff <- fTemp$rc), envir=tmp)
      eval(quote(apollo_inputs$apollo_lcPars    <- fTemp$lp), envir=tmp)
      eval(quote(rm(fTemp))                                 , envir=tmp)
      if(debug) cat(' Done. ',sum(gc()[,2]),'MB of RAM in use\n',sep='') # do not change to apollo_print
    } else if(debug) cat(' Failed.\n')
    rm(fileName)
    apollo_inputs$apollo_scaling <- ns; rm(ns)
  }
  
  ### Calculate bootstrap s.e. if required
  if(bootstrapSE>0){
    if(!silent) apollo_print("\nStarting bootstrap calculation of standard errors.")
    # If using multiple cores, save a copy of apollo_inputs
    fileName <- paste0(tempdir(), "\\", apollo_inputs$apollo_control$modelName, "_inputs_extra")
    if(apollo_inputs$apollo_control$nCores>1) saveRDS(apollo_inputs, file=fileName)
    # Prepare inputs and call apollo_bootstrap
    tmp <- list(estimationRoutine=estimationRoutine, maxIterations=maxIterations,
                writeIter=FALSE, hessianRoutine="none", printLevel=printLevel,
                maxLik_settings=maxLik_settings, silent=silent)
    model$bootvarcov <- apollo_bootstrap(apollo_beta=b, apollo_fixed,
                                         apollo_probabilities_ORIG, apollo_inputs,
                                         estimate_settings=tmp,
                                         bootstrap_settings=list(nRep=bootstrapSE, seed=bootstrapSeed,
                                                                 calledByEstimate=TRUE))$varcov
    model$bootse <- sqrt(diag(model$bootvarcov))
    bVar         <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    dummyVCM     <- matrix(NA, nrow=length(bVar), ncol=length(bVar), dimnames=list(names(bVar), names(bVar)))
    model$bootcorrmat <- tryCatch(model$bootvarcov/(model$bootse%*%t(model$bootse)), error=function(e) return(dummyVCM))
    # update number of bootstrap repetitions
    bootstrapSE <- tryCatch(nrow(read.csv(paste0(apollo_inputs$apollo_control$outputDirectory,apollo_inputs$apollo_control$modelName, '_bootstrap_params.csv'))),
                            error=function(e) bootstrapSE)
    
    if(length(apollo_fixed)>0){
      model$bootse <- c(model$bootse, rep(NA, length(apollo_fixed)))
      names(model$bootse) <- c(colnames(model$bootvarcov), apollo_fixed)
      model$bootse <- model$bootse[names(apollo_beta)]
    }
    # If using multiple cores, restore apollo_inputs
    if(apollo_inputs$apollo_control$nCores>1) apollo_inputs <- tryCatch(readRDS(file=fileName), error=function(e) NULL)
    if(is.null(apollo_inputs)) stop('INTERNAL ISSUE - apollo_inputs could not be restored from disk after bootstrap')
    rm(bVar, dummyVCM)
  }
  
  # #################### #
  #### Prepare output ####
  # #################### #
  
  # ### Restore fixed parameters to model$estimate
  # temp           = c(model$estimate, apollo_beta[apollo_fixed])
  # model$estimate = temp[names(apollo_beta)]
  
  # ### Save and print overview
  # model$componentReport <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="report"), 
  #                                   error=function(e) NULL)
  # test <- is.list(model$componentReport) && !is.null(names(model$componentReport)) && any(c('data', 'param') %in% names(model$componentReport))
  # if(test) model$componentReport <- list(model=model$componentReport)
  # if(!silent){
  #   test <- FALSE
  #   for(r in model$componentReport) if(is.list(r) && !is.null(r$param) && length(r$param)>0){
  #     test <- TRUE
  #     for(j in r$param) cat(j, '\n', sep='')
  #   }; rm(r)
  #   if(test) cat('\n')
  # }
  
  # ### Calculate Zero LL
  # if(!silent) apollo_print("Calculating log-likelihood at equal shares (LL(0)) for applicable models... ")
  # model$LL0 <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"),
  #                       error=function(e) return(NA))
  # if(is.list(model$LL0)){
  #   model$LL0=model$LL0[c("model",names(model$LL0)[names(model$LL0)!="model"])]
  #   for(s in 1:length(model$LL0)){
  #     if(is.list(model$LL0[[s]])) model$LL0[[s]]=model$LL0[[s]][["model"]]
  #   }
  #   if(!workInLogs) model$LL0=sapply(model$LL0,function(x) sum(log(x)))
  #   if( workInLogs) model$LL0=sapply(model$LL0,sum)
  # } else {
  #   model$LL0 <- ifelse( workInLogs, sum(model$LL0), sum(log(model$LL0)) )
  # }
  
  # ### Calculate shares LL (constants only)
  # if(apollo_control$calculateLLC){
  #   if(!silent) apollo_print("Calculating log-likelihood at observed shares from estimation data (LL(c)) for applicable models... ")
  #   model$LLC <- tryCatch(apollo_probabilities(apollo_beta, apollo_inputs, functionality="shares_LL"),
  #                         error=function(e) return(NA))
  #   if(is.list(model$LLC)){
  #     model$LLC = model$LLC[c("model",names(model$LLC)[names(model$LLC)!="model"])]
  #     for(s in 1:length(model$LLC)) if(is.list(model$LLC[[s]])) model$LLC[[s]] = model$LLC[[s]][["model"]]
  #     if(workInLogs) model$LLC <- sapply(model$LLC,sum) else model$LLC <- sapply(model$LLC,function(x) sum(log(x)))
  #   } else model$LLC <- ifelse(workInLogs, sum(model$LLC), sum(log(model$LLC)) )
  # } else model$LLC=NA
  
  # ### Get LL at optimum for each component
  # if(!silent) apollo_print("Calculating LL of each model component... ")
  # LLout <- tryCatch(apollo_probabilities(model$estimate, apollo_inputs, functionality="output"),
  #                   error=function(e){ apollo_print("Could not complete validation using estimated parameters."); return(NA) }
  # )
  # if(!anyNA(LLout) && is.list(LLout)){
  #   LLout <- LLout[c("model",names(LLout)[names(LLout)!="model"])]
  #   for(s in 1:length(LLout)) if(is.list(LLout[[s]])) LLout[[s]]=LLout[[s]][["model"]]
  #   if(!workInLogs) LLout <- lapply(LLout, log)
  #   LLout <-lapply(LLout,sum)
  #   model$LLout <- unlist(LLout)
  # } else{
  #   model$LLout <- list(NA)
  #   if(!silent) apollo_print("LL could not be calculated for all components.")
  # }
  
  # ### Calculate Rho2, AIC, and BIC
  # nFreeParams <- length(apollo_beta) - length(apollo_fixed)
  # test <- !is.null(model$modelTypeList) && !anyNA(model$LL0[1])
  # test <- test && all(tolower(model$modelTypeList) %in% c("mnl", "nl", "cnl", "el", "dft", "lc", "rrm", "ol", "op"))
  # if(test & !silent) apollo_print("Calculating other model fit measures")
  # if(test) model$rho2_0 <- 1-(model$maximum/model$LL0[1]) else model$rho2_0 <- NA
  # if(test) model$adjRho2_0 <- 1-((model$maximum-nFreeParams)/model$LL0[1]) else model$adjRho2_0 <- NA
  # test <- test && is.numeric(model$LLC) && !anyNA(model$LLC[1])
  # if(test) model$rho2_C <- 1-(model$maximum/model$LLC[1]) else model$rho2_C <- NA
  # if(test) model$adjRho2_C <- 1-((model$maximum-nFreeParams)/model$LLC[1]) else model$adjRho2_C <- NA
  # model$AIC <- -2*model$maximum + 2*nFreeParams
  # model$BIC <- -2*model$maximum + nFreeParams*log(ifelse(!is.null(model$nObsTot), model$nObsTot, model$nObs))
  # model$nFreeParams <- nFreeParams
  
  # ### use pre-scaling values as starting values in output
  # model$bootstrapSE          <- bootstrapSE
  # model$apollo_beta          <- temp_start
  # model$LLStart              <- sum(testLL)
  # model$startTime            <- time1
  # model$apollo_control       <- apollo_control
  # model$nObs                 <- nObs
  # model$nIndivs              <- length(indiv)
  # model$apollo_draws         <- apollo_inputs$apollo_draws
  # model$estimationRoutine    <- estimationRoutine
  # model$scaling              <- apollo_inputs$apollo_scaling
  # model$manualScaling        <- apollo_inputs$manualScaling
  
  # ### Save constraints (NULL if none or given in maxLik format)
  # test <- exists("apollo_constraints", envir=environment(), inherits=FALSE)
  # if(test) model$apollo_constraints <- apollo_constraints else model$apollo_constraints <- NULL
  # rm(test)
  
  ### De-scale parameters
  model$estimate[names(model$scaling)] <- model$estimate[names(model$scaling)]*model$scaling
  
  ### Fourth checkpoint
  time4 <- Sys.time()
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  # model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))
  # model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))
  model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))

  ### assign apollo class to model
  class(model)<-c("apollo",class(model))  
  return(model)
}
