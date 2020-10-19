#' Uses EM for models with continuous random coefficients
#'
#' Uses the EM algorithm for estimating a model with continuous random coefficients.
#'
#' This function uses the EM algorithm for estimating a model with continuous random coefficients. It is only suitable for models where all parameters are random, with a full covariance matrix.
#' All random parameters need to be based on underlying Normals with a full covariance matrix, but any transform thereof can be used.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters. These need to be provided in the following order. With K random parameters, K means for the underlying Normals, followed by the elements of the lower triangle of the Cholesky matrix, by row.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param mixEM_settings List. Options controlling the EM process.
#'                                 \itemize{
#'                                   \item \strong{stoppingCriterion}: Numeric. Convergence criterion. The EM process will stop when improvements in the log-likelihood fall below this value. Default is 10^-5.
#'                                   \item \strong{EMmaxIterations}: Numeric. Maximum number of iterations of the EM algorithm before stopping. Default is 100.
#'                                                                 Used only if \code{apollo_control$HB} is FALSE. Default is 200.
#'                                   \item \strong{postEM}: Numeric scalar. Determines the number of tasks performed by this function 
#'                                               after the EM algorithm has converged. Can take values \code{0}, \code{1} 
#'                                               or \code{2} only. If value is \code{0}, only the EM algorithm will be 
#'                                               performed, and the results will be a model object without a covariance 
#'                                               matrix (i.e. estimates only.). If value is \code{1}, after the EM 
#'                                               algorithm the covariance matrix of the model will be calculated as well, 
#'                                               and the result will be a model object with a covariance matrix. If value 
#'                                               is \code{2}, after the EM algorithm, the estimated parameter values will 
#'                                               be used as starting value for a maximum likelihood estimation process, 
#'                                               which will render a model object with a covariance matrix. Performing 
#'                                               maximum likelihood estimation after the EM algorithm is useful, as there 
#'                                               may be room for further improvement. Default is \code{2}.
#'                                   \item \strong{silent}: Boolean. If TRUE, no information is printed to the console during estimation. Default is FALSE.
#'                                   \item \strong{transforms}: List. Optional argument, with one entry per parameter, showing the inverse transform to return from beta to the underlying Normal. E.g. if the first parameter is specified as negative logormal inside apollo_randCoeff, then the entry in transforms should be transforms[[1]]=function(x) log(-x)
#'                                 }
#' @param estimate_settings List. Options controlling the estimation process within each EM iteration. See \link{apollo_estimate} for details.
#' @return model object
#' @export
apollo_mixEM=function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, mixEM_settings=NA, estimate_settings=NA){
  
  # # # # # # # # # # #
  #### Initialise  ####
  # # # # # # # # # # #
  
  ### First checkpoint
  time1 <- Sys.time()
  
  if(is.na(mixEM_settings) || is.null(mixEM_settings$transforms)) cat("The list \'transforms\' has not been provided in \'mixEM_settings\'. Apollo will assume that all parameters use untransformed Normal distributions.\n\n")
  ### Set missing settings to default values
  default <- list(EMstoppingCriterion=10^-5, EMmaxIterations=100, postEM=2, silent=FALSE, transforms=NA)
  if(length(mixEM_settings)==1 && is.na(mixEM_settings)) mixEM_settings <- default
  tmp <- names(default)[!(names(default) %in% names(mixEM_settings))] # options missing in estimate_settings
  for(i in tmp) mixEM_settings[[i]] <- default[[i]]
  test <- is.list(mixEM_settings) && !is.null(mixEM_settings$postEM) && (mixEM_settings$postEM %in% 0:2)
  if(!test) stop('Setting "postEM" inside argument "mixEM_settings" can only take values 0, 1 or 2.')
  
  
  ### Extract variables from apollo_input
  stopping_criterion = mixEM_settings[["EMstoppingCriterion"]]
  EMmaxIterations    = mixEM_settings[["EMmaxIterations"]]
  calculate_SE       = mixEM_settings[["calculateSE"]]
  silent             = mixEM_settings[["silent"]]
  transforms         = mixEM_settings[["transforms"]]
  
  ### Load estimate_settings defaults
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=FALSE, 
                  hessianRoutine="analytic", printLevel=3L, constraints=NULL, maxLik_settings=NULL, 
                  numDeriv_settings=list(), apollo_scaling=NA, scaling=NA, bootstrapSE=0, bootstrapSeed=24, silent=FALSE)
  test <- length(estimate_settings)==1 && is.na(estimate_settings)
  if(test) estimate_settings <- default else{
    test <- !(names(default) %in% names(estimate_settings))
    if(any(test)) for(i in names(default)[test]) estimate_settings[[i]] <- default[[i]]
    rm(i)
  }; rm(default, test)
  if(is.null(estimate_settings$maxLik_settings)){
    estimate_settings$maxLik_settings <- list(printLevel=3, iterlim=200)
    if(!is.null(estimate_settings$printLevel)) estimate_settings$maxLik_settings$printLevel <- estimate_settings$printLevel
    if(!is.null(estimate_settings$maxIterations)) estimate_settings$maxLik_settings$iterlim <- estimate_settings$maxIterations
  }
  estimate_settings$apollo_scaling <- estimate_settings$scaling
  
  ### Create an unmodified copy of apollo_probabilities
  apollo_probabilities_ORIG <- apollo_probabilities
  if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
    apollo_randCoeff_ORIG <- apollo_inputs$apollo_randCoeff
  }
  if(is.function(apollo_inputs$apollo_lcPars)) apollo_lcPars_ORIG <- apollo_inputs$apollo_lcPars
  
  # Constraints are not supported. It would require altering the constraints based on which parameters are fixed
  if(!is.null(estimate_settings$constraints)) stop('Constraints are not supported')
  
  ### Validate inputs
  if(apollo_inputs$apollo_control$mixing==FALSE) stop("The apollo_mixEM function cannot be used for models that do not include continuous random parameters!")
  if(apollo_inputs$apollo_control$HB==TRUE) stop("The apollo_mixEM function cannot be used with Bayesian estimation!")
  if(!is.na(apollo_inputs$apollo_lcPars)) stop("The apollo_mixEM function cannot be used for models that use latent classes!")
  if(apollo_inputs$apollo_draws$intraNDraws!=0) stop("The apollo_mixEM function cannot be used for models using intra-individual draws!")
  if(!is.null(apollo_inputs$apollo_draws$interUnifDraws)) stop("The apollo_mixEM function cannot be used for models using uniform draws!")
  if(!apollo_inputs$apollo_draws$interDrawsType%in%c("halton","mlhs","pmc","sobol","sobolOwen","sobolFaureTezuka","sobolOwenFaureTezuka")) stop("The apollo_mixEM function cannot be used for models importing user defined draws!")
  
  apollo_print("The use of apollo_mixEM has a number of requirements. No checks are run for these, so the user needs to ensure these conditions are met by their model:")
  apollo_print("1:This function is only suitable for single component models, i.e. no use of apollo_combineModels or manual multiplication of model components.")
  apollo_print("2:All parameters need to be random, and a full covariance matrix needs to be estimated/specified in apollo_randCoeff.")
  apollo_print("3:All random parameters need to be based on Normal distributions or transformations thereof.")
  apollo_print("4:With K random parameters, the order of the elements in \'apollo_beta\' needs to be as follows: K means for the underlying Normals, followed by the elements of the lower triangle of the Cholesky matrix, by row.")
  apollo_print('\n')
  
  apollo_print("Initialising EM algorithm")
  if(!silent){
    apollo_print("Validating inputs of likelihood function (apollo_probabilities)")
    apollo_print('\n')
  }
  
  ### Validate model
  apollo_inputs$apollo_control$analyticGrad=FALSE
  database <- apollo_inputs$database # neccesary to circumvent deletion of database across the call stack in apollo_estimate
  draws    <- apollo_inputs$draws    # neccesary to circumvent deletion of draws across the call stack in apollo_estimate
  apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                  estimate_settings=list(maxIterations=0, writeIter=FALSE,hessianRoutine="none",silent=TRUE))
  apollo_inputs$database <- database; rm(database)
  apollo_inputs$draws    <- draws;    rm(draws)
  apollo_inputs$apollo_scaling <- estimate_settings$scaling
  apollo_inputs$scaling <- NULL  
  
  # # # # # # # # # # # # # #
  #### Prepare functions ####
  # # # # # # # # # # # # # #
  
  ### Validate scaling
  test <- is.null(apollo_inputs$apollo_scaling)
  test <- test || (length(apollo_inputs$apollo_scaling)==1 && is.na(apollo_inputs$apollo_scaling))
  if(test) apollo_inputs$apollo_scaling <- setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                                                    names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  if(any(!(names(apollo_inputs$apollo_scaling) %in% names(apollo_beta)))) stop("Some parameters included in 'scaling' are not included in 'apollo_beta'")
  if(any((names(apollo_inputs$apollo_scaling) %in% apollo_fixed))) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
  if(any(apollo_inputs$apollo_scaling<0)){
    apollo_inputs$apollo_scaling <- abs(apollo_inputs$apollo_scaling)
    txt <- 'WARNING: Some values in "scaling" were negative, they were replaced by their absolute value.'
    if(!silent) apollo_print(txt) else warning('Some negative values in "scaling" were replaced by their absolute value')
  }
  if(any(apollo_inputs$apollo_scaling<=0)) stop('All terms in "scaling" should be strictly positive!')
  txt <- "During estimation, parameters will be scaled using the values in estimate_settings$scaling"
  if(!all(apollo_inputs$apollo_scaling==1)){ if(!silent) apollo_print(txt) else warning(txt)}
  rm(txt)
  
  ### Scale apollo_probabilities, apollo_randCoeff and apollo_lcPars, and names components in apollo_probabilities
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
  test <- is.function(apollo_inputs$apollo_lcPars)
  if(test) apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
  
  # Scale starting parameters
  if(length(apollo_inputs$apollo_scaling)>0 && !is.na(apollo_inputs$apollo_scaling)){
    r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- apollo_beta[r]/apollo_inputs$apollo_scaling[r]
  }
  
  ### Create temporary model object
  model=list()
  
  ### Second checkpoint
  time2 <- Sys.time()

  ### Calculate initial likelihood at the level of each draw
  Ln=apollo_probabilities(apollo_beta, apollo_inputs, functionality="conditionals")
  cat("\n\nInitial LL:",sum(log(rowMeans(Ln))),"\n\n")
  LLStart <- tryCatch(sum(log(rowMeans(Ln))), error=function(e) NA)
  
  # # # # # # # # # # # # # # # #
  #### Main estimation loop  ####
  # # # # # # # # # # # # # # # #
  
  ### Loop over repeated EM iterations until convergence has been reached
  iteration=1
  stop=0
  apollo_inputs$silent=TRUE
  while((stop==0) & (iteration<= EMmaxIterations)){
    cat("Starting iteration: ",iteration,"\n",sep="")
    cat("Current LL: ",sum(log(rowMeans(Ln))),"\n",sep="")
    
    ### Calculate weight for each individual and for each draw 
    wn=Ln/(rowMeans(Ln))
    
    ### Copy current parameter values into temporary model object
    model$estimate=apollo_beta
    
    ### Produce draws for random coefficients with current vector of parameters
    d=apollo_unconditionals(model,apollo_probabilities,apollo_inputs)
    
    ### Translate draws back to Normal from negative Lognormal
    if(!(length(transforms)==1 && is.na(transforms))){
      for(s in 1:length(d)){
        if(is.function(transforms[[s]])) d[[s]]=transforms[[s]](d[[s]])
      }
    }
    
    ### Apply weights to individual draws and turn into a matrix
    dwn=lapply(d,"*",wn)
    dwn=lapply(dwn,as.vector)
    dwn=do.call(cbind,dwn)
    
    ### Calculate means for weighted draws
    mu=colMeans(dwn)
    
    ### Calculate weighted covariance matrix
    K = length(mu)   # n of coefficients
    R = ncol(d[[1]]) # n of draws
    N = nrow(d[[1]]) # n of individuals
    tmp=matrix(0, nrow=N, ncol=K)
    Omega=matrix(0, nrow=K, ncol=K)
    for(r in 1:R){
      for(k in 1:K) tmp[,k] = d[[k]][,r] - mu[k]
      for(n in 1:N) Omega = Omega + wn[n,r]*(tmp[n,] %*% t(tmp[n,]))
    } 
    
    ### Compute Cholesky of average weighted covariance matrix
    cholesky = chol(Omega/(N*R))
    cholesky = cholesky[upper.tri((cholesky),diag=TRUE)]
    
    ### Update vector of model parameters on the basis of calculated mu and Omega
    apollo_beta[1:4]=mu
    apollo_beta[5:14]=cholesky
    
    ### Calculate likelihood with new parameters
    Lnew=((apollo_probabilities(apollo_beta, apollo_inputs, functionality="conditionals")))
    
    ### Compute improvement
    change=sum(log(rowMeans(Lnew)))-sum(log(rowMeans(Ln)))
    cat("New LL: ",sum(log(rowMeans(Lnew))),"\n",sep="")
    cat("Improvement: ",change,"\n\n",sep="")
    Ln=Lnew
    
    ### Determine whether convergence has been reached
    if(change<stopping_criterion) stop=1
    iteration=iteration+1
  }
  
  if(iteration>=EMmaxIterations){
    apollo_print("EM algorithm stopped: maximum number of iterations reached!") 
    apollo_print(paste("No covariance matrix will be computed as convergence has not been reached.",
                       "You may call apollo_addCovariance(model,apollo_inputs) to compute a covariance",
                       "matrix at the current estimates.")) 
  } else apollo_print("EM algorithm stopped: improvements in LL smaller than convergence criterion.")
  
  apollo_inputs = apollo_validateInputs(silent=TRUE)
  
  ### Third checkpoint
  time3 <- Sys.time()

  # De-scale parameters
  if(length(apollo_inputs$apollo_scaling)>0 && !is.na(apollo_inputs$apollo_scaling)){
    r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- apollo_inputs$apollo_scaling[r]*apollo_beta[r]
  }
  
  #### POST-EM TASKS
  apollo_inputs$EM = FALSE
  if(mixEM_settings$postEM>0 && iteration<EMmaxIterations){
    if(mixEM_settings$postEM==1) apollo_print("Computing covariance matrix...")
    if(mixEM_settings$postEM>1) apollo_print("Continuing with classical estimation...")
    ### Reinstate original functions inside apollo_inputs
    test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
    if(test) apollo_inputs$apollo_randCoeff <- apollo_randCoeff_ORIG
    if(is.function(apollo_inputs$apollo_lcPars)) apollo_inputs$apollo_lcPars <- apollo_lcPars_ORIG
    ### Avoid diagnostics and validation
    apollo_inputs$apollo_control$noValidation  <- TRUE
    apollo_inputs$apollo_control$noDiagnostics <- TRUE
    ### Set maxIter and writeIter
    if(!is.list(estimate_settings)) estimate_settings = list()
    if(mixEM_settings$postEM<2) estimate_settings$maxIterations <- 0
    estimate_settings$writeIter <- FALSE
    estimate_settings$silent    <- ifelse(mixEM_settings$postEM<2, TRUE, FALSE)
    ### Calculate S.E.
    model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities_ORIG, apollo_inputs, estimate_settings)
  } else model=apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, 
                               apollo_inputs, estimate_settings=list(maxIterations=0,hessianRoutine="none",silent=TRUE))  
  if(mixEM_settings$postEM>1) classicalIter=model$nIter
  time4 <- Sys.time()
  
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  ### combine EM and classical times
  model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))+ifelse(!is.null(model$timePre),model$timePre,0)
  model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))+ifelse(!is.null(model$timeEst),model$timeEst,0)
  ##model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))
  model$timePost  <- model$timeTaken-model$timePre-model$timeEst
  model$estimationRoutine = paste0('EM algorithm')
  if(mixEM_settings$postEM>1) model$estimationRoutine <- paste0('EM algorithm (', estimate_settings$estimationRoutine, 
                                                               ') -> Maximum likelihood (', 
                                                               estimate_settings$estimationRoutine, ')')
  model$nIter = iteration
  if(mixEM_settings$postEM>1) model$nIter = paste0(model$nIter," (EM) & ",classicalIter, " (",estimate_settings$estimationRoutine,")")
  if(iteration>=EMmaxIterations) model$nIter = paste0(iteration," (convergence not reached)")
  model$LLStart <- LLStart
  
  return(model)
}