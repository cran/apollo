#' Uses EM for latent class model
#'
#' Uses the EM algorithm for estimating a latent class model.
#'
#' This function uses the EM algorithm for estimating a Latent Class model. It is only suitable for models without 
#' continuous mixing. All parameters that vary across classes need to be included in the \code{apollo_lcPars} 
#' function which is used by \code{apollo_lcEM}.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not 
#'                     change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                             \itemize{
#'                               \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                               \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                               \item functionality: Character. Can be either "estimate" (default), "prediction", 
#'                                     "validate", "conditionals", "zero_LL", or "raw".
#'                             }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param lcEM_settings List. Options controlling the EM process.
#'                      \itemize{
#'                        \item \strong{stoppingCriterion}: Numeric. Convergence criterion. The EM process will stop when 
#'                                                          improvements in the log-likelihood fall below this value. 
#'                                                          Default is 10^-5.
#'                        \item \strong{EMmaxIterations}: Numeric. Maximum number of iterations of the EM algorithm before 
#'                                                        stopping. Default is 100. Used only if \code{apollo_control$HB} 
#'                                                        is FALSE. Default is 200.
#'                        \item \strong{postEM}: Numeric scalar. Determines the number of tasks performed by this function 
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
#'                        \item \strong{silent}: Boolean. If TRUE, no information is printed to the console during 
#'                                               estimation. Default is FALSE.
#'                      }
#' @param estimate_settings List. Options controlling the estimation process within each EM iteration. See 
#'                          \link{apollo_estimate} for details.
#' @return model object
#' @importFrom stats rnorm
#' @importFrom utils stack
#' @importFrom parallel clusterCall stopCluster
#' @export
apollo_lcEM=function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, lcEM_settings=NA, estimate_settings=NA){
  
  apollo_print("The use of apollo_lcEM has a number of requirements. No checks are run for these, so the user needs to ensure these conditions are met by their model:")
  apollo_print("1:This function is only suitable for single component models, i.e. no use of apollo_combineModels or manual multiplication of model components.")
  apollo_print("2:Any parameters that vary across classes need to be included in the definition of random parameters in apollo_lcPars.")
  apollo_print("3:The entries in the lists in apollo_lcPars need to be individual parameters, not functions thereof.")
  apollo_print('\n')

  ### First checkpoint
  time1 <- Sys.time()
  
  ### Set missing lcEM_settings to default values
  default <- list(EMstoppingCriterion=10^-5, EMmaxIterations=100, postEM=2, silent=FALSE) # calculateSE=TRUE
  if(length(lcEM_settings)==1 && is.na(lcEM_settings)) lcEM_settings <- default
  tmp <- names(default)[!(names(default) %in% names(lcEM_settings))] # options missing in lcEM_settings
  for(i in tmp) lcEM_settings[[i]] <- default[[i]]
  test <- is.list(lcEM_settings) && !is.null(lcEM_settings$postEM) && (lcEM_settings$postEM %in% 0:2)
  if(!test) stop('Setting "postEM" inside argument "lcEM_settings" can only take values 0, 1 or 2.')
  rm(i,tmp)
  
  ### Extract variables from apollo_input
  database           = apollo_inputs[["database"]]
  apollo_lcPars      = apollo_inputs[["apollo_lcPars"]]
  stopping_criterion = lcEM_settings[["EMstoppingCriterion"]]
  EMmaxIterations    = lcEM_settings[["EMmaxIterations"]]
  silent             = lcEM_settings[["silent"]]
  
  if(!silent){
    apollo_print("Validating inputs of likelihood function (apollo_probabilities)")
    apollo_print('\n')
  } 

  ### Load estimate_settings defaults
  default <- list(estimationRoutine="bfgs", maxIterations=200, writeIter=FALSE, 
                  hessianRoutine="numDeriv", printLevel=3L, constraints=NULL, maxLik_settings=NULL, 
                  numDeriv_settings=list(), scaling=NA, bootstrapSE=0, bootstrapSeed=24, silent=FALSE)
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
  apollo_inputs$apollo_scaling <- estimate_settings$scaling
  apollo_inputs$scaling <- NULL
  if(anyNA(apollo_inputs$apollo_scaling)){
    tmp <- abs(apollo_beta[!(names(apollo_beta) %in% apollo_fixed)])
    apollo_inputs$apollo_scaling <- setNames(rep(1, length(tmp)), names(tmp))
    rm(tmp)
  }
  
  ### Test apollo_probabilities
  #tmp = apollo_inputs$apollo_control$nCores
  #apollo_inputs$apollo_control$nCores = 1
  apollo_inputs$class_specific        = 0
  database <- apollo_inputs$database # necessary to circumvent deletion of database across the call stack in apollo_estimate
  draws    <- apollo_inputs$draws    # necessary to circumvent deletion of draws across the call stack in apollo_estimate
  apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                  estimate_settings=list(maxIterations=0, writeIter=FALSE, hessianRoutine="none",silent=TRUE))
  apollo_inputs$database <- database
  apollo_inputs$draws    <- draws;    rm(draws)
  #apollo_inputs$apollo_control$nCores = tmp
  #rm(tmp)
  
  ### Create an unmodified copy of apollo_probabilities, apollo_randCoeff and apollo_lcPars
  apollo_probabilities_ORIG <- apollo_probabilities
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_randCoeff_ORIG <- apollo_inputs$apollo_randCoeff
  if(is.function(apollo_inputs$apollo_lcPars)) apollo_lcPars_ORIG <- apollo_inputs$apollo_lcPars
  
  ### Checks
  # Using constraints would require altering the constraints based on which parameters are fixed
  if(!is.null(estimate_settings$constraints)) stop('Constraints are not supported')
  if(apollo_inputs$apollo_control$mixing==TRUE) stop("The apollo_lcEM function cannot be used for models that include continuous random parameters!")
  if(apollo_inputs$apollo_control$HB==TRUE) stop("The apollo_lcEM function cannot be used with Bayesian estimation!")
  
  ### Commence EM
  apollo_print("Initialising EM algorithm")

  ### Add weighting to apollo_control
  apollo_inputs$apollo_control$weights="weights"
  
  ### Keep backup of vector of fixed parameters as this changes throughout
  apollo_fixed_base = apollo_fixed
  
  ### get number of classes
  lcPars <- apollo_lcPars
  environment(lcPars) <- list2env( c(apollo_inputs$database, as.list(apollo_beta)), hash=TRUE)
  classes = length(lcPars(apollo_beta,apollo_inputs)[["pi_values"]])
  
  apollo_EMClassifyParams <- function(apollo_beta, apollo_inputs){
    # Validate apollo_beta & apollo_lcPars
    if(is.null(names(apollo_beta))) stop("apollo_beta must be a named vector.")
    if(is.null(apollo_inputs$apollo_lcPars)) stop("apollo_lcPars is missing from apollo_inputs.")
    if(!is.function(apollo_inputs$apollo_lcPars)) stop("apollo_lcPars is not a functions.")
    
    # Run apollo_lcPars with dummy values to identify them later
    b <- c(0,0)
    while(length(b)!=length(unique(b))){
      b <- setNames(round(rnorm(length(apollo_beta), mean=1000, sd=200),0), names(apollo_beta))
    }
    lcPars <- apollo_inputs$apollo_lcPars
    environment(lcPars) <- list2env(c(as.list(b), 
                                      apollo_inputs$database#, 
                                      #apollo_inputs$draws
    ), hash=TRUE)
    lcPars <- lcPars(apollo_beta, apollo_inputs)
    
    # Validate return of apollo_lcPars, get number of classes and remove the pi_values (NAMED ENFORCED)
    if(!is.list(lcPars)) stop("The apollo_lcPars function should return a list.")
    if(!all(sapply(lcPars, is.list))) stop("All elements inside the list returned by apollo_lcPars should be lists.")
    S <- max(sapply(lcPars, length))
    if(!all(sapply(lcPars, length)==S)) stop("All elements inside the list returned by apollo_lcPars should have the same length.")
    if(S==1) stop("At least two classes must be defined in apollo_lcPars.")
    if(is.null(lcPars$pi_values)) stop("The list returned by apollo_lcPars should contain an element called 'pi_values'.")
    lcPars <- lcPars[names(lcPars)[names(lcPars)!="pi_values"]]
    
    # Scale b if necessary
    test <- !is.null(apollo_inputs$apollo_scaling) && !is.null(names(apollo_inputs$apollo_scaling))
    test <- test && length(grep('apollo_scaling', capture.output(apollo_beta)))>0
    if(test) for(i in names(apollo_inputs$apollo_scaling)) b[i] <- b[i]*apollo_inputs$apollo_scaling[i]
    
    # Recover names of parameters used in utilities
    ans <- list()
    for(s in 1:S){
      tmp <- c()
      for(i in lcPars){
        test <- (i[[s]] %in% b) || (length(i[[s]])==1 && is.numeric(i[[s]])) # is param or value
        if(!test) stop(paste("A value inside the return of apollo_lcPars (except for pi_values)",
                             "is not inside apollo_beta nor a fixed value. E.g. The second ",
                             "element of lcpars[[1]] = list(b1, b1+1) is not acceptable."))
        tmp <- c(tmp, names(b)[which(b==i[[s]])])
      } 
      ans[[paste0("class_",s)]] <- tmp
    }
    rm(b, tmp)
    
    # Put repeated values in a separate list, and remove them from the existing ones
    gen <- c()
    for(s in 1:(S-1)) for(i in ans[[paste0("class_",s)]]) {
      for(s2 in (s+1):S) if(i %in% ans[[paste0("class_",s2)]]) gen <- c(gen, i)
    }
    if(length(gen)>0) for(g in gen) for(s in paste0("class_",1:S)){
      ans[[s]] <- ans[[s]][ans[[s]]!=g] # This could lead to some elements being c()
    }
    ans$generic <- gen
    #rm(gen, i, g, s, s2)
    
    # Recover names used in class allocation (i.e. those used inside apollo_lcPars but not already in ans)
    lcPars <- body(apollo_inputs$apollo_lcPars)
    tmp    <- all.vars(lcPars)
    tmp    <- tmp[tmp %in% names(apollo_beta)]
    tmp    <- tmp[!(tmp %in% unlist(ans))]
    ans$class_alloc <- tmp
    rm(tmp)
    
    # Put any remaining element in apollo_beta into ans$generic
    tmp <- names(apollo_beta)
    tmp <- tmp[!(tmp %in% unlist(ans))]
    if(length(tmp)>0) ans$generic <- c(ans$generic, tmp)
    rm(tmp)
    
    return(ans)
  }
  
  apollo_beta_list = apollo_EMClassifyParams(apollo_beta, apollo_inputs)
  
  ### new line 18 Nov
  #if(all(apollo_beta_list$generic%in%apollo_fixed)) apollo_beta_list$generic=c()
  
  ### Stop if there are generic parameters
  test <- length(apollo_beta_list$generic)>0 && !all(apollo_beta_list$generic %in% apollo_fixed)
  if(test) stop('The EM algorithm for latent classes cannot handle generic parameters across classes (', 
                paste0(apollo_beta_list$generic, collapse=', '),')')
  
  ### DEFINE MODEL AND LIKELIHOOD FUNCTION FOR CLASS ALLOCATION
  apollo_probabilities_class = function(apollo_beta, apollo_inputs, functionality="estimate"){
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P = list()
    h = apollo_inputs$h
    ###new line 16 Oct
    #lcPars = apollo_inputs$apollo_lcPars(apollo_beta, apollo_inputs)
    #log_pi_values = lapply(lcPars$pi_values,log)
    log_pi_values = lapply(get('pi_values'),log)
    P[["model"]]  = exp(Reduce('+', mapply('*',h,log_pi_values,SIMPLIFY = FALSE)))
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  
  ### Set flag for EM
  apollo_inputs$EM = TRUE
  
  ### Create logLike for apollo_probabilities_class
  classLL <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities_class, apollo_inputs, 
                                apollo_estSet=estimate_settings)
  
  ### Create logLike for apollo_probabilities
  withinLL <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                                     apollo_estSet=estimate_settings)
  
  ### Create function to update h and class_specific
  updateHCS <- function(h, cs, w=FALSE, LL){
    if(!is.logical(w)) stop('Argument "w" must be a scalar Logical')
    ### Update h & class_specific inside apollo_inputs
    if(anyNA(environment(LL)$cl)){
      # Single core
      if(!is.null(h) ) environment(LL)$apollo_inputs$h                <- h
      if(!is.null(cs)) environment(LL)$apollo_inputs$class_specific   <- cs
      if(w){
        nObsPerIndiv <- environment(LL)$apollo_inputs$database[,environment(LL)$apollo_inputs$apollo_control$indivID]
        nObsPerIndiv <- as.vector(table(nObsPerIndiv))
        environment(LL)$apollo_inputs$database$weights <- rep(h[[cs]], times=nObsPerIndiv)
      } 
    } else {
      # Multi-core
      if(is.list(h)){
        # split h into lists with only the relevant data for each core (assumes length(h[[i]])==nInd)
        nInd<- length(unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID]))
        ind <- parallel::clusterEvalQ(environment(LL)$cl, length(unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])))
        ind <- c(0, unlist(ind))
        for(i in 3:length(ind)) ind[i] <- ind[i-1] + ind[i]
        h2  <- vector(mode='list', length=length(ind)-1)
        for(nc in 1:length(h2)) h2[[nc]] <- lapply(h, apollo_keepRows, 1:nInd %in% (ind[nc]+1):ind[nc+1])
        h   <- h2; rm(h2)
      } else h <- vector(mode='list', length=length(environment(LL)$cl))
      # copy elements into each core
      parallel::parLapply(environment(LL)$cl, h, function(hh, w, cs){
        apollo_inputs   <- get('apollo_inputs', envir=globalenv())
        if(!is.null(cs)) apollo_inputs$class_specific <- cs
        if( is.list(hh)) apollo_inputs$h              <- hh
        nObsPerIndiv    <- as.vector(table(apollo_inputs$database[,apollo_inputs$apollo_control$indivID]))
        if(w & !is.null(cs)) apollo_inputs$database$weights <- rep(hh[[cs]], times=nObsPerIndiv)
        tmp <- globalenv()
        assign("apollo_inputs", apollo_inputs, tmp)
        return(TRUE)
      }, w=w, cs=cs)
    }
  }
  
  ### Scale apollo_probabilities
  test <- !is.null(apollo_inputs$apollo_scaling)
  if(test){
    ### new 13 Oct  
    r <- names(apollo_beta) %in% names(apollo_inputs$apollo_scaling)
    r <- names(apollo_beta)[r]
    apollo_beta[r] <- apollo_beta[r]/apollo_inputs$apollo_scaling[r]
    ### end
    apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
    test2 <- is.function(apollo_inputs$apollo_lcPars)
    if(test2) apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
    test2 <- is.function(apollo_inputs$apollo_randCoeff)
    if(test2) apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, apollo_inputs$apollo_scaling)
  } 
  
  ### Second checkpoint
  time2 <- Sys.time()
  
  ### Step 1
  
  ### Create temporary model object
  model=list()
  
  ### Loop over repeated EM iterations until convergence has been reached
  iteration        = 1
  stop             = 0
  LLStart          = -Inf
  while((stop==0) & (iteration<= EMmaxIterations)){
    cat("Starting iteration: ",iteration,"\n",sep="")
    
    # ########## #
    ### Step 2 ###
    # ########## #
    
    ### Calculate model likelihood and class specific likelihoods
    L = apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")
    if(iteration==1) LLStart <- tryCatch(sum(log(L$model)), error=function(e) NA)
    
    # ########## #
    ### Step 3 ###
    # ########## #
    
    ### Calculate class specific conditional likelihoods
    model$estimate = apollo_beta
    apollo_inputs$silent=TRUE
    pi             = apollo_lcUnconditionals(model, apollo_probabilities, apollo_inputs)[["pi_values"]]
    h              = list()
    for(s in 1:classes){
      h[[s]] = as.vector( pi[[s]]*L[[s]]/L[[(classes+1)]] )
    }
    
    ### Calculate current log-likelihood for LC model
    Lcurrent = sum( log(L[[(classes+1)]]) )
    cat("Current LL : ",Lcurrent,"\n",sep="")
    
    # ########## #
    ### Step 4 ###
    # ########## #  
    
    ### Set fixed parameters (only estimating parameters for class allocation model)
    fix_temp                  = apollo_beta_list
    fix_temp[["class_alloc"]] = NULL
    apollo_fixed              = unique( c(apollo_fixed_base, stack(fix_temp)[,1]) )
    
    ### Update h & class_specific inside apollo_inputs
    updateHCS(h=h, cs=0, w=FALSE, LL=classLL)
    
    ### Estimate class allocation model
    theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    theta_fix_val <- apollo_beta[apollo_fixed]
    environment(classLL)$bFixedVal <- theta_fix_val
    # Estimate
    model <- maxLik::maxLik(classLL, start=theta_var_val,
                            method=estimate_settings$estimationRoutine, 
                            control=estimate_settings$maxLik_settings, 
                            constraints=estimate_settings$constraints,
                            print.level=0, finalHessian=FALSE,
                            countIter=FALSE, writeIter=FALSE, sumLL=FALSE)
    # Restore fixed parameters
    temp           = c(model$estimate, apollo_beta[apollo_fixed])
    model$estimate = temp[names(apollo_beta)]
    
    ### Update overall parameters
    apollo_beta=model$estimate
    
    # ########## #
    ### Step 5 ###
    # ########## #
    
    
    ### Update coefficients in class specific models by estimating class specific models 
    ### using posterior class allocation probabilities as weights
    
 
    for(s in 1:classes){
      
      
      ### Set fixed parameters (only estimating parameters for class 1)
      fix_temp = apollo_beta_list
      fix_temp[[paste0("class_",s)]] = NULL
      apollo_fixed = unique(c(apollo_fixed_base,(stack(fix_temp)[,1])))
      
      ### Set class index to use inside apollo_probabilities_within_class
      updateHCS(h=h, cs=s, w=TRUE, LL=withinLL)
      
      # Create logLike function
      theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
      theta_fix_val <- apollo_beta[apollo_fixed]
      environment(withinLL)$bFixedVal <- theta_fix_val
      # Estimate
      model <- maxLik::maxLik(withinLL, start=theta_var_val,
                              method=estimate_settings$estimationRoutine, 
                              control=estimate_settings$maxLik_settings, 
                              constraints=estimate_settings$constraints,
                              print.level=0, finalHessian=FALSE,
                              countIter=FALSE, writeIter=FALSE, sumLL=FALSE)
      # Restore fixed parameters
      temp           = c(model$estimate, apollo_beta[apollo_fixed])
      model$estimate = temp[names(apollo_beta)]
      
      ### Update overall parameters
      apollo_beta=model$estimate
    }
    
    ### Update generic coefficients if present
    test <- !is.null(apollo_beta_list$generic) && (length(apollo_beta_list$generic)!=0)
    test <- test && !all(apollo_beta_list$generic%in%apollo_fixed_base)
    if(test){
      
      ### Set class index to use inside apollo_probabilities_within_class
      updateHCS(h=NULL, cs=0, w=FALSE, LL=withinLL)
      
      ### Set fixed parameters (only estimating generic parameters)
      fix_temp=apollo_beta_list
      fix_temp[["generic"]]=NULL
      apollo_fixed=unique(c(apollo_fixed_base,(stack(fix_temp)[,1])))
      # Create logLike function
      theta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
      theta_fix_val <- apollo_beta[apollo_fixed]
      environment(withinLL)$bFixedVal <- theta_fix_val
      # Estimate
      model <- maxLik::maxLik(withinLL, start=theta_var_val,
                              method=estimate_settings$estimationRoutine, 
                              control=estimate_settings$maxLik_settings, 
                              constraints=estimate_settings$constraints,
                              print.level=0, finalHessian=FALSE,
                              countIter=FALSE, writeIter=FALSE, sumLL=FALSE)
      # Restore fixed parameters
      temp           = c(model$estimate, apollo_beta[apollo_fixed])
      model$estimate = temp[names(apollo_beta)]
      
      ### Update overall parameters
      apollo_beta=model$estimate
      
    }
    
    # ########## #
    ### Step 6 ###
    # ########## #
    
    ### Calculate new log-likelihood and compute improvement
    apollo_inputs$class_specific = 0
    Lnew   = sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="output")[[(classes+1)]]))
    change = Lnew - Lcurrent
    cat("New LL     : ",Lnew,"\n",sep="")
    cat("Improvement: ",change,"\n\n",sep="")
    
    ### Determine whether convergence has been reached
    if(change<stopping_criterion) stop = 1
    iteration = iteration + 1
  }
  
  
  if(iteration>=EMmaxIterations){
    apollo_print("EM algorithm stopped: maximum number of iterations reached!") 
    apollo_print(paste("No covariance matrix will be computed as convergence has not been reached.",
                      "You may call apollo_addCovariance(model,apollo_inputs) to compute a covariance",
                      "matrix at the current estimates.")) 
  } else apollo_print("EM algorithm stopped: improvements in LL smaller than convergence criterion.")
  
  ### Close cluster
  if(!anyNA(environment(classLL)$cl)) parallel::stopCluster(environment(classLL)$cl)
  if(!anyNA(environment(withinLL)$cl)) parallel::stopCluster(environment(withinLL )$cl)
  
  ### Update apollo_inputs
  tmp <- apollo_inputs$apollo_scaling
  apollo_inputs = apollo_validateInputs(silent=TRUE)
  apollo_inputs$apollo_scaling <- tmp
  rm(tmp)
  
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
  ### Reinstate original functions inside apollo_inputs
  test <- apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)
  if(test) apollo_inputs$apollo_randCoeff <- apollo_randCoeff_ORIG
  if(is.function(apollo_inputs$apollo_lcPars)) apollo_inputs$apollo_lcPars <- apollo_lcPars_ORIG
  ### Avoid diagnostics and validation
  apollo_inputs$apollo_control$noValidation  <- TRUE
  apollo_inputs$apollo_control$noDiagnostics <- TRUE
  ### Get a model object. Call to apollo_estimate is different depending on needing a covariance matrix or not
  if(lcEM_settings$postEM>0 && iteration<EMmaxIterations){
    if(lcEM_settings$postEM==1) apollo_print("Computing covariance matrix...")
    if(lcEM_settings$postEM>1) apollo_print("Continuing with classical estimation...")
    ### Set maxIter and writeIter
    if(!is.list(estimate_settings)) estimate_settings = list()
    if(lcEM_settings$postEM<2) estimate_settings$maxIterations <- 0
    estimate_settings$writeIter <- FALSE
    estimate_settings$silent    <- ifelse(lcEM_settings$postEM<2, TRUE, FALSE)
    ### Calculate S.E.
    model = apollo_estimate(apollo_beta, apollo_fixed_base, apollo_probabilities_ORIG, apollo_inputs, estimate_settings)
  #fixed 22 Nov
  #  } else model=apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities_ORIG, 
  #                             apollo_inputs, estimate_settings=list(maxIterations=0,hessianRoutine="none",silent=TRUE))  
    } else model=apollo_estimate(apollo_beta, apollo_fixed_base, apollo_probabilities_ORIG, 
                                 apollo_inputs, estimate_settings=list(maxIterations=0,hessianRoutine="none",silent=TRUE))  
    if(lcEM_settings$postEM>1) classicalIter=model$nIter
  time4 <- Sys.time()
  
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  ### combine EM and classical times
  model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))+ifelse(!is.null(model$timePre),model$timePre,0)
  model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))+ifelse(!is.null(model$timeEst),model$timeEst,0)
  ##model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))
  model$timePost  <- model$timeTaken-model$timePre-model$timeEst
  model$estimationRoutine = paste0('EM algorithm')
  if(lcEM_settings$postEM>1) model$estimationRoutine <- paste0('EM algorithm (', estimate_settings$estimationRoutine, 
                                                               ') -> Maximum likelihood (', 
                                                               estimate_settings$estimationRoutine, ')')
  model$nIter = iteration
  if(lcEM_settings$postEM>1) model$nIter = paste0(model$nIter," (EM) & ",classicalIter, " (",estimate_settings$estimationRoutine,")")
  if(iteration>=EMmaxIterations) model$nIter = paste0(iteration," (convergence not reached)")
  model$LLStart <- LLStart
  
  return(model)
}
