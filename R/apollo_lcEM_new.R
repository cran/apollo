#' Uses EM for latent class model
#'
#' Uses the EM algorithm for estimating a latent class model.
#'
#' This function uses the EM algorithm for estimating a Latent Class model. It is only suitable for models without 
#' continuous mixing. All parameters need to vary across classes and need to be included in the \code{apollo_lcPars} 
#' function which is used by \code{apollo_lcEM}.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not 
#'                     change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either 
#'                            \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param lcEM_settings List. Options controlling the EM process.
#'                      \itemize{
#'                        \item \strong{EMmaxIterations}: Numeric. Maximum number of iterations of the EM algorithm before 
#'                                                        stopping. Default is 100. 
#'                        \item \strong{postEM}: Numeric scalar. Determines the tasks performed by this function 
#'                                               after the EM algorithm has converged. Can take values \code{0}, \code{1} 
#'                                               or \code{2} only. If value is \code{0}, only the EM algorithm will be 
#'                                               performed, and the results will be a model object without a covariance 
#'                                               matrix (i.e. estimates only). If value is \code{1}, after the EM 
#'                                               algorithm, the covariance matrix of the model will be calculated as well, 
#'                                               and the result will be a model object with a covariance matrix. If value 
#'                                               is \code{2}, after the EM algorithm, the estimated parameter values will 
#'                                               be used as starting value for a maximum likelihood estimation process, 
#'                                               which will render a model object with a covariance matrix. Performing 
#'                                               maximum likelihood estimation after the EM algorithm is useful, as there 
#'                                               may be room for further improvement. Default is \code{2}.
#'                        \item \strong{silent}: Boolean. If TRUE, no information is printed to the console during 
#'                                               estimation. Default is FALSE.
#'                        \item \strong{stoppingCriterion}: Numeric. Convergence criterion. The EM process will stop when 
#'                                                          improvements in the log-likelihood fall below this value. 
#'                                                          Default is 10^-5.
#'                      }
#' @param estimate_settings List. Options controlling the estimation process within each EM iteration. See 
#'                          \link{apollo_estimate} for details.
#' @return model object
#' @importFrom stats rnorm
#' @importFrom utils stack
#' @importFrom parallel clusterCall stopCluster
apollo_lcEM_new=function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, lcEM_settings=NA, estimate_settings=NA){
  
  apollo_print("The use of apollo_lcEM has a number of requirements. No checks are run for these, so the user needs to ensure these conditions are met by their model:", type="i")
  apollo_print("1: This function is only suitable for single component models, i.e. no use of apollo_combineModels or manual multiplication of model components.")
  apollo_print("2: Any parameters that vary across classes need to be included in the definition of random parameters in apollo_lcPars.")
  apollo_print("3: The entries in the lists in apollo_lcPars need to be individual parameters, not functions thereof.")
  apollo_print('\n')

  ### First checkpoint
  time1 <- Sys.time()
  
  ### Set missing lcEM_settings to default values
  default <- list(EMstoppingCriterion=10^-5, EMmaxIterations=100, postEM=2, silent=FALSE) # calculateSE=TRUE
  if(length(lcEM_settings)==1 && is.na(lcEM_settings)) lcEM_settings <- default
  tmp <- names(default)[!(names(default) %in% names(lcEM_settings))] # options missing in lcEM_settings
  for(i in tmp) lcEM_settings[[i]] <- default[[i]]
  test <- is.list(lcEM_settings) && !is.null(lcEM_settings$postEM) && (lcEM_settings$postEM %in% 0:2)
  if(!test) stop('SYNTAX ISSUE - Setting "postEM" inside argument "lcEM_settings" can only take values 0, 1 or 2.')
  rm(i,tmp)
  
  ### Extract variables from apollo_input
  database           = apollo_inputs[["database"]]
  apollo_lcPars      = apollo_inputs[["apollo_lcPars"]]
  stopping_criterion = lcEM_settings[["EMstoppingCriterion"]]
  EMmaxIterations    = lcEM_settings[["EMmaxIterations"]]
  silent             = lcEM_settings[["silent"]]
  # set analytic gradient to FALSE for EM part if apollo_lcPars does not use mnl or classAlloc
  test <- any(grepl("apollo_mnl",deparse(apollo_lcPars))|grepl("apollo_classAlloc",deparse(apollo_lcPars)))
  if(!test) apollo_inputs$apollo_control$analyticGrad = FALSE
  
  if(!silent){
    apollo_print("Validating inputs of likelihood function (apollo_probabilities)")
    apollo_print('\n')
  } 

  ### Load estimate_settings defaults
  default <- list(estimationRoutine="bgw", maxIterations=200, writeIter=TRUE, 
                  hessianRoutine="analytic", printLevel=3L, constraints=NULL, 
                  maxLik_settings=NULL, numDeriv_method="Richardson", 
                  numDeriv_settings=list(), scaling=NA, 
                  bootstrapSE=0, bootstrapSeed=24, silent=FALSE, 
                  scaleHessian=TRUE, scaleAfterConvergence=FALSE,
                  validateGrad=FALSE,
                  bgw_settings=list())
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
                  estimate_settings=list(maxIterations=1, writeIter=FALSE, hessianRoutine="none",silent=TRUE))
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
  if(!is.null(estimate_settings$constraints)) stop('INCORRECT FUNCTION/SETTING USE - Constraints are not supported with the EM algorithm!')
  if(apollo_inputs$apollo_control$mixing==TRUE) stop("INCORRECT FUNCTION/SETTING USE - he apollo_lcEM function cannot be used for models that include continuous random parameters!")
  if(apollo_inputs$apollo_control$HB==TRUE) stop("INCORRECT FUNCTION/SETTING USE - The apollo_lcEM function cannot be used with Bayesian estimation!")

  ### Commence EM
  apollo_print("Initialising EM algorithm")

  ### Second checkpoint
  time2 <- Sys.time()

  ### Create functions required for EM
  make_modified_copy <- function(fun,search_lines, repl_lines) {
    strip_ws <- function(x) gsub("\\s+", "", x)
    code <- deparse(body(fun))              # turn body into character vector
    # locate the two needles
    code_key   <- strip_ws(code)
    needle_key <- strip_ws(search_lines)
    
    idx1 <- match(needle_key[1], code_key)
    idx2 <- match(needle_key[2], code_key)
    if(is.na(idx1) || is.na(idx2) || idx2 != idx1 + 1){
      txt=paste0("Original apollo_probabilities function did not contain the following two lines next to each other:\n",
                 search_lines[1], "\n", search_lines[2], "\n",
                 "Please check the function and try again.")
      stop(apollo_print(txt))
    }
    # replace the two-line chunk with your five-line block
    new_code <- c(
      code[1:(idx1 - 1)],   # everything before the first target line
      repl_lines,           # your five replacement lines
      code[(idx2 + 1):length(code)]   # everything after the second target line
    )
    
    # stitch back together → language object → function
    new_body <- parse(text = paste(new_code, collapse = "\n"))[[1]]
    new_fun  <- fun                         # copy the full function object
    body(new_fun) <- new_body                   # overwrite its body
    environment(new_fun) <- environment(fun)    # keep same environment
    new_fun
  }
  
  apollo_probabilities_EM <- make_modified_copy(
    fun          = apollo_probabilities,
    search_lines = c("lc_settings = list(inClassProb = P, classProb=pi_values)",
                     "P[[\"model\"]] = apollo_lc(lc_settings, apollo_inputs, functionality)"),
    repl_lines   = c(
        'P[["part_1"]] = 0',
        'P[["part_2"]] = 0',
        'classProb <- apollo_firstRow(pi_values, apollo_inputs)',
        'for(ss in 1:length(apollo_inputs$pi)){',
        ' P[["part_1"]] = P[["part_1"]] + apollo_inputs$pi[[ss]] * log(P[[ss]])',
        ' P[["part_2"]] = P[["part_2"]] + apollo_inputs$pi[[ss]] * log(classProb[[ss]])',
        '}',
        ' P[["part_1"]] = exp(P[["part_1"]])',
        ' P[["part_2"]] = exp(P[["part_2"]])',
        'P <- apollo_combineModels(P, apollo_inputs, functionality,components=c("part_1","part_2"))'
    )
  )
  
  posteriors=function(apollo_beta){
    tmp=apollo_conditionals(model=list(estimate=apollo_beta),apollo_probabilities,apollo_inputs)
    pi=list()
    for(s in 2:length(tmp)) pi[[s-1]]=tmp[,s]
    return(pi)
  }
  
  ###
  LLStart          = -Inf
  for(iteration in 1:EMmaxIterations){
    cat("\n\n Starting iteration: ",iteration,"out of ",EMmaxIterations,"\n")
    if(iteration==1) LLStart <- tryCatch(sum(log(apollo_probabilities(apollo_beta, apollo_inputs,functionality="estimate"))), error=function(e) NA)
    pi=posteriors(apollo_beta)
    apollo_inputs$pi=pi
    model = apollo_estimate(apollo_beta, apollo_fixed, 
                            apollo_probabilities_EM, apollo_inputs, estimate_settings=list(silent=TRUE))
    cat("\n    LL start (EM):   ",model$LLStart)
    cat("\n    LL end   (EM):   ",model$maximum)
    Lcurrent=sum(log(apollo_probabilities(apollo_beta,apollo_inputs, functionality="estimate")))
    cat("\n    LL start (true): ",Lcurrent)
    Lnew=sum(log(apollo_probabilities(model$estimate,apollo_inputs, functionality="estimate")))
    cat("\n    LL end   (true): ",Lnew)
    apollo_beta=model$estimate
    change = Lnew - Lcurrent
    cat("\n    Improvement:     ",change)

    ### Determine whether convergence has been reached
    if(change<stopping_criterion) break
  }
  
  
  if(iteration>=EMmaxIterations){
    apollo_print("\n\nEM algorithm stopped: maximum number of iterations reached!") 
    apollo_print(paste("No covariance matrix will be computed as convergence has not been reached.",
                      "You may call apollo_addCovariance(model,apollo_inputs) to compute a covariance",
                      "matrix at the current estimates.")) 
  } else apollo_print("\n\nEM algorithm stopped: improvements in LL smaller than convergence criterion.")
  
  ### Third checkpoint
  time3 <- Sys.time()
  

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
  if(lcEM_settings$postEM>0){#} && iteration<EMmaxIterations){
    if(lcEM_settings$postEM==1) apollo_print("Computing covariance matrix...")
    if(lcEM_settings$postEM>1) apollo_print("Continuing with classical estimation...")
    ### Set maxIter and writeIter
    if(!is.list(estimate_settings)) estimate_settings = list()
    if(lcEM_settings$postEM<2) estimate_settings$maxIterations <- 0
    estimate_settings$writeIter <- FALSE
    estimate_settings$silent    <- ifelse(lcEM_settings$postEM<2, TRUE, FALSE)
    ### Calculate S.E.
    model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities_ORIG, apollo_inputs, estimate_settings)
  #fixed 22 Nov
  #  } else model=apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities_ORIG, 
  #                             apollo_inputs, estimate_settings=list(maxIterations=0,hessianRoutine="none",silent=TRUE))  
    } else model=apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities_ORIG, 
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
