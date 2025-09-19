#' Conducts all the pre-estimation steps for a model
#'
#' Conducts all the pre-estimation steps for a model using the likelihood function defined by \code{apollo_probabilities}.
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
#'                                  \item \strong{\code{estimationRoutine}}: Character. Estimation method. Can take values "bfgs", "bgw", "bhhh", or "nr". Used only if \code{apollo_control$HB} is FALSE. Default is "bgw".
#'                                  \item \strong{\code{hessianRoutine}}: Character. Name of routine used to calculate the Hessian of the log-likelihood function after estimation. Valid values are \code{"analytic"} (default), \code{"numDeriv"} (to use the numeric routine in package numDeric), \code{"maxLik"} (to use the numeric routine in packahe maxLik), and \code{"none"} to avoid calculating the Hessian and the covariance matrix. Only used if \code{apollo_control$HB=FALSE}.
#'                                  \item \strong{\code{maxIterations}}: Numeric. Maximum number of iterations of the estimation routine before stopping. Used only if \code{apollo_control$HB} is FALSE. Default is 200.
#'                                  \item \strong{\code{maxLik_settings}}: List. Additional settings for maxLik. See argument \code{control} in \link[maxLik]{maxBFGS}, \link[maxLik]{maxBHHH} and \link[maxLik]{maxNM} for more details. Only used for maximum likelihood estimation.
#'                                  \item \strong{\code{numDeriv_method}}: Character. Method used for numerical differentiation when calculating the covariance matrix. Can be \code{"Richardson"} or \code{"simple"}, Only used if analytic gradients are available. See argument \code{method} in \link[numDeriv]{grad} for more details.
#'                                  \item \strong{\code{numDeriv_settings}}: List. Additional arguments to the method used by numDeriv to calculate the Hessian. See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                                  \item \strong{\code{printLevel}}: Higher values render more verbous outputs. Can take values 0, 1, 2 or 3. Ignored if apollo_control$HB is TRUE. Default is 3.
#'                                  \item \strong{\code{scaleAfterConvergence}}: Logical. Used to increase numerical precision of convergence. If TRUE, parameters are scaled to 1 after convergence, and the estimation is repeated from this new starting values. Results are reported scaled back, so it is a transparent process for the user. Default is FALSE.
#'                                  \item \strong{\code{scaleHessian}}: Logical. If TRUE, parameters are scaled to 1 for Hessian estimation. Default is TRUE.
#'                                  \item \strong{\code{scaling}}: Named vector. Names of elements should match those in \code{apollo_beta}. Optional scaling for parameters. If provided, for each parameter \code{i}, \code{(apollo_beta[i]/scaling[i])} is optimised, but \code{scaling[i]*(apollo_beta[i]/scaling[i])} is used during estimation. For example, if parameter b3=10, while b1 and b2 are close to 1, then setting \code{scaling = c(b3=10)} can help estimation, specially the calculation of the Hessian. Reports will still be based on the non-scaled parameters.
#'                                  \item \strong{\code{silent}}: Logical. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                  \item \strong{\code{validateGrad}}: Logical. If TRUE, the analytical gradient (if used) is compared to the numerical one. Default is FALSE.
#'                                  \item \strong{\code{writeIter}}: Logical. Writes value of the parameters in each iteration to a csv file. Works only if \code{estimation_routine=="bfgs"|"bgw"}. Default is TRUE.
#'                                 }
#' @return model object
#' @export
#' @importFrom numDeriv hessian grad
#' @importFrom maxLik maxLik
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov runif
#' @importFrom bgw bgw_mle
apollo_preEstimate  <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
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
  default <- list(estimationRoutine="bgw", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="analytic", printLevel=3L, constraints=NULL, 
                  maxLik_settings=NULL, numDeriv_method="Richardson", 
                  numDeriv_settings=list(), scaling=NA, 
                  bootstrapSE=0, bootstrapSeed=24, silent=FALSE, 
                  scaleHessian=TRUE, scaleAfterConvergence=FALSE,
                  validateGrad=FALSE,
                  bgw_settings=list())
  prtLvlMan <- is.list(estimate_settings) && !is.null(estimate_settings$printLevel)
  if(length(estimate_settings)==1 && is.na(estimate_settings)) estimate_settings <- default
  
  ### if using HB in apollo_control, set estimationRoutine to HB to avoid using BGW defaults with e.g. scaling
  if(apollo_inputs$apollo_control$HB) estimate_settings$estimationRoutine="HB"
  
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
  
  ### do not use BGW if maxIterations==0
  if(estimate_settings$maxIterations==0) estimate_settings$estimationRoutine="bfgs"
  
  ### Set missing BGW settings to default values
  default <- list(maxIterations    = estimate_settings$maxIterations,
                  maxFunctionEvals = max(1,3*estimate_settings$maxIterations),
                  silent           = estimate_settings$silent,
                  outputDirectory  = apollo_inputs$apollo_control$outputDirectory,
                  modelName        = apollo_inputs$apollo_control$modelName,
                  printLevel       = estimate_settings$printLevel, 
                  printNonDefaultSettings = FALSE, 
                  printStartingValues     = FALSE,
                  printFinalResults       = FALSE,
                  writeIter               = estimate_settings$writeIter,
                  writeIterMode           = "overWrite",
                  userScaleVector         = NULL)
  
  if(length(estimate_settings$bgw_settings)==1 && is.na(estimate_settings$bgw_settings)) estimate_settings$bgw_settings <- default
  tmp <- names(default)[!(names(default) %in% names(estimate_settings$bgw_settings))] # options missing in estimate_settings
  for(i in tmp) estimate_settings$bgw_settings[[i]] <- default[[i]]
  rm(default, tmp)
  #estimate_settings$bgw_settings$writeIter = FALSE ### FORCING BGW TO NOT WRITE ITERATIONS, EVEN IF REQUESTED
  
  
  if(exists("i")) rm(i)
  apollo_inputs$apollo_scaling <- estimate_settings$scaling
  apollo_inputs$scaling <- NULL
  
  ## change estimationRoutine to its lower case equivalent also inside estimate_settings which is used in some places
  estimate_settings$estimationRoutine = tolower( estimate_settings[["estimationRoutine"]] )
  
  ### if using BGW, put scaling into BGW, and avoid Apollo scaling, remembering that BGW uses inverse Apollo scaling!
  if(estimate_settings$estimationRoutine=="bgw" && !all(is.na(apollo_inputs$apollo_scaling))){
    if(!is.null(estimate_settings$bgw_settings$scalingMethod) && estimate_settings$bgw_settings$scalingMethod!="userScaling") stop("SYNTAX ISSUE - bgw_settings$scalingMethod needs to be set to 'userScaling' when providing a scaling vector!")
    estimate_settings$bgw_settings$scalingMethod = "userScaling"
    estimate_settings$bgw_settings$userScaleVector = 1/apollo_inputs$apollo_scaling
    apollo_inputs$apollo_scaling=setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                                          names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  }
  
  ### Stop if apollo_fixed passed to apollo_estimate is different from the one passed via apollo_inputs
  if(!all(apollo_fixed==apollo_inputs$apollo_fixed)) stop("INPUT ISSUE - The vector 'apollo_fixed' passed to 'apollo_estimate' differs from the one in 'apollo_inputs'. Please run 'apollo_inputs = apollo_validateInputs()' before proceeding with estimation!") 
  
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
  if( !(estimationRoutine %in% c("bfgs", "bgw", "bhhh", "nr", "hb")) ) stop("SYNTAX ISSUE - Invalid estimationRoutine. Use 'bfgs', 'bgw', 'bhhh', 'hb', or 'nr'.")
  if( !(estimate_settings$hessianRoutine %in% c('analytic', 'numDeriv', 'maxLik', 'none')) ) stop("SYNTAX ISSUE - Invalid hessianRoutine. Use 'analytic', 'numDeriv', 'maxLik' or 'none'.")
  
  # Check apollo_beta and apollo_fixed
  if(!is.numeric(apollo_beta) | !is.vector(apollo_beta) | is.null(names(apollo_beta))) stop("INPUT ISSUE - The \"apollo_beta\" argument needs to be a named vector")
  if(length(apollo_fixed)>0 && !is.character(apollo_fixed)) stop("INPUT ISSUE - 'apollo_fixed' is not an empty vector nor a vector of names.")
  if(length(unique(names(apollo_beta)))<length(apollo_beta)) stop("INPUT ISSUE - The \"apollo_beta\" argument contains duplicate elements")
  if(length(unique(apollo_fixed))<length(apollo_fixed)) stop("INPUT ISSUE - The \"apollo_fixed\" argument contains duplicate elements")
  if(!all(apollo_fixed %in% names(apollo_beta))) stop("INPUT ISSUE - Some parameters included in 'apollo_fixed' are not included in 'apollo_beta'.")
  
  # Check numDeriv_settings
  test <- estimate_settings[["numDeriv_method"]] %in% c("Richardson", "simple")
  if(!test) stop("SYNTAX ISSUE - Setting 'numDeriv_method' must be one of 'Richardson' or 'simple'.")
  
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
  
  if(apollo_inputs$apollo_control$nCores>1) stop('SYNTAX ISSUE - This function can only be run for single core')
  if(apollo_control$HB) stop('SYNTAX ISSUE - This function cannot be run for HB')
  
  ### Check if any of the fixed parameters are not zero
  if(!apollo_control$noValidation){
    if(length(apollo_fixed)>0){
      fixed_params=apollo_beta[apollo_fixed]
      if(any(!(fixed_params%in%c(0,1)))){
        tmp=names(fixed_params[!(fixed_params%in%c(0,1))])
        one=(length(tmp)==1)
        txt <- paste0('Element', ifelse(one,' ', 's '), paste0(tmp, collapse=', '), 
                      ' in \'apollo_fixed\' ', ifelse(one,'is ', 'are '),' constrained to ', ifelse(one,'a value ', 'values '),' other than zero or one. This may be intentional. ',
                      'If not, stop this function by pressing the "Escape" key and adjust the starting values accordingly.')
        apollo_print(txt, pause=0, type="w")
      }
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
  
  ### check if scaling was inserted
  tmp <- as.character(body(L$apollo_probabilities))
  tmp <- grep("apollo_inputs$apollo_scaling", tmp, fixed=TRUE)
  scaling_failure=FALSE
  if(length(tmp)==0) scaling_failure=TRUE
  
  ### Update functions if modification was successful, or fall back if not
  #if(L$success){
  apollo_probabilities           <- L$apollo_probabilities
  apollo_inputs$apollo_randCoeff <- L$apollo_randCoeff
  apollo_inputs$apollo_lcPars    <- L$apollo_lcPars
  apollo_inputs$apollo_scaling   <- L$apollo_scaling
  apollo_inputs$manualScaling    <- L$manualScaling
  if(scaling_failure){
    apollo_inputs$apollo_scaling=setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                                          names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
    estimate_settings$scaleHessian=FALSE
    scaleAfterConvergence=FALSE
  } 
  #} else {
  # apollo_inputs$apollo_scaling   <- L$apollo_scaling
  # apollo_inputs$manualScaling    <- all(L$apollo_scaling==1)
  # estimate_settings$scaleHessian <- FALSE
  # estimate_settings$scaleAfterConvergence <- FALSE
  # scaleAfterConvergence          <- estimate_settings[['scaleAfterConvergence']]
  #}; 
  rm(L)
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  
  # Scale starting parameters (this should always happen except with BGW)
  if(estimate_settings$estimationRoutine!="bgw" && length(apollo_inputs$apollo_scaling)>0 && !anyNA(apollo_inputs$apollo_scaling)){
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
      apollo_inputs$apollo_control$analyticHessian <- FALSE 
      environment(apollo_logLike)$analyticGrad  <- FALSE
      environment(apollo_logLike)$nIter         <- 1     # fix for iteration counting
    }
  } else grad <- NULL
  
  ### Create hessian function
  if(!is.null(apollo_inputs$apollo_control$analyticHessian) && apollo_inputs$apollo_control$analyticHessian){
    # apollo_makeHessian will create gradient function ONLY IF all components have analytical hessians.
    hessian <- apollo_makeHessian(apollo_beta, apollo_fixed, apollo_logLike)
    if(is.null(hessian)){
      if(!silent) apollo_print(paste0("Analytical hessians could not be calculated for all components, ", 
                                      "numerical hessian will be used."))
      apollo_inputs$apollo_control$analyticHessian <- FALSE 
      environment(apollo_logLike)$analyticHessian  <- FALSE
    }
  } else hessian <- NULL
  
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
  
  ### Scale constraints if necessary
  if(hasEqConst) for(a in colnames(constraints$eqA)) constraints$eqA[,a] <- constraints$eqA[,a]*apollo_inputs$apollo_scaling[a]
  if(hasIneqConst) for(a in colnames(constraints$ineqA)) constraints$ineqA[,a] <- constraints$ineqA[,a]*apollo_inputs$apollo_scaling[a]
  
  rm(apollo_beta) ## apollo_beta is not included in env that is returned
  return(environment())
}
