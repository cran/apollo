#' Calculates variance-covariance matrix of an Apollo model
#' 
#' Calculates the Hessian, variance-covariance matrix and standard errors of an Apollo model as defined by its likelihood function 
#' and \code{apollo_inputs} list of settings. Performs automatic scaling for increased numeric stability.
#' 
#' It calculates the Hessian, variance-covariance, and standard errors at \code{apollo_beta} values of an 
#' estimated model. At least one of the following settings must be provided (ordered by speed of computation): \code{apollo_grad}, 
#' \code{apollo_logLike}, or (\code{apollo_probabilities} and \code{apollo_inputs}). If more than one is provided, 
#' then the priority is: \code{apollo_grad}, \code{apollo_logLike}, (\code{apollo_probabilities} and \code{apollo_inputs}).
#' 
#' @param apollo_beta Named numeric vector. Names and values of parameters at which to calculate the covariance matrix.
#'                    Values \strong{must not be scaled}, and they must include any fixed parameter.
#' @param apollo_fixed Character vector. Names of fixed parameters.
#' @param varcov_settings List of settings defining the behaviour of this function. It must contain at least one of
#'                        the following: \code{apollo_logLike}, \code{apollo_grad} or \code{apollo_inputs} together with \code{apollo_probabilities}.
#'                        \itemize{
#'                          \item \strong{\code{apollo_grad}}: Function to calculate the gradient of the model, as returned by \link{apollo_makeGrad}.
#'                          \item \strong{\code{apollo_inputs}}: List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#'                          \item \strong{\code{apollo_logLike}}: Function to calculate the log-likelihood of the model, as returned by \link{apollo_makeLogLike}.
#'                          \item \strong{\code{apollo_probabilities}}: apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either 
#'                            \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#'                          \item \strong{\code{BHHH_matrix}}: Matrix. Optional input, providing the BHHH matrix so it does not get recalculated.
#'                          \item \strong{\code{hessianRoutine}}: Character. Name of routine used to calculate the Hessian. Valid values are \code{"analytic"}, \code{"numDeriv"}, \code{"maxLik"} or \code{"none"} to avoid estimating the Hessian and covariance matrix.
#'                          \item \strong{\code{numDeriv_settings}}: List. Additional arguments to the Richardson method used by numDeriv to calculate the Hessian. See argument \code{method.args} in \link[numDeriv]{grad} for more details.
#'                          \item \strong{\code{scaleBeta}}: Logical. If TRUE (default), parameters are scaled by their own value before calculating the Hessian to increase numerical stability. However, the output is de-scaled, so they are in the same scale as the \code{apollo_beta} argument.
#'                        }
#' @return List with the following elements
#'         \itemize{
#'           \item \strong{\code{apollo_beta}}: Named numerical vector. Parameter estimates (\code{model$estimate}, not scaled).
#'           \item \strong{\code{corrmat}}: Numerical matrix. Correlation between parameter estimates.
#'           \item \strong{\code{hessian}}: Numerical matrix. Hessian of the model at parameter estimates (\code{model$estimate}).
#'           \item \strong{\code{hessianScaling}}: Named numeric vector. Scales used on the paramaters to calculate the Hessian (non-fixed only).
#'           \item \strong{\code{methodsAttempted}}: Character vector. Name of methods attempted to calculate the Hessian.
#'           \item \strong{\code{methodUsed}}: Character. Name of method used to calculate the Hessian.
#'           \item \strong{\code{robcorrmat}}: Numerical matrix. Robust correlation between parameter estimates.
#'           \item \strong{\code{robse}}: Named numerical vector. Robust standard errors of parameter estimates.
#'           \item \strong{\code{robvarcov}}: Numerical matrix. Robust variance-covariance matrix.
#'           \item \strong{\code{se}}: Named numerical vector. Standard errors of parameter estimates.
#'           \item \strong{\code{varcov}}: Numerical matrix. Variance-covariance matrix.
#'         }
#' @importFrom stats var
#' @export
apollo_varcov <- function(apollo_beta, apollo_fixed, varcov_settings){
  # # # # # # #  # # # # #
  #### Initialisation ####
  # # # # # # #  # # # # #
  
  ### Load default settings and perform checks on them
  if(!is.list(varcov_settings)) stop('SYNTAX ISSUE - Argument "varcov_settings" must be a list.')
  default <- list(hessianRoutine='analytic', scaleBeta=TRUE, numDeriv_settings=list())
  for(i in names(default)) if(is.null(varcov_settings[[i]])) varcov_settings[[i]] <- default[[i]]
  # At least one of apollo_grad, apollo_logLike or apollo_inputs must be provided
  test <- is.function(varcov_settings$apollo_grad) || is.function(varcov_settings$apollo_logLike)
  test <- test || (is.function(varcov_settings$apollo_probabilities) && is.list(varcov_settings$apollo_inputs))
  if(!test) stop('SYNTAX ISSUE - Argument varcov_settings must contain at least one of the following: "apollo_logLike", "apollo_grad" or ("apollo_probabilities" and "apollo_inputs").')
  # Routine can only be 'analytic','numDeriv', 'maxLik', or 'none'
  test <- varcov_settings$hessianRoutine %in% c('analytic','numDeriv', 'maxLik', 'none')
  if(!test) stop('SYNTAX ISSUE - Setting "hessianRoutine" must take one of the following values: "analytic", "numDeriv", "maxLik" or "none".')
  # If numeric routine is requested, apollo_logLike or apollo_inputs must be available
  test <- varcov_settings$hessianRoutine %in% c('numDeriv', 'maxLik')
  test <- test && !is.list(varcov_settings$apollo_inputs) && !is.function(varcov_settings$apollo_logLike)
  if(test) stop('SYNTAX ISSUE - Cannot use numeric routine ("numDeriv" or "maxLik") if only the analytical gradient is given ("apollo_grad").')
  # If apollo_control$analyticGrad is FALSE, then no analytic routine is possible, and it defaults to numeric.
  aGrad <- FALSE # value of apollo_control$analyticGrad (wherever that is stored)
  if(is.function(varcov_settings$apollo_grad)) aGrad <- TRUE else {
    if(is.function(varcov_settings$apollo_logLike)){
      aGrad <- environment(varcov_settings$apollo_logLike)$analyticGrad
    } else if(is.list(varcov_settings$apollo_inputs)) aGrad <- varcov_settings$apollo_inputs$apollo_control$analyticGrad
  }
  if(!aGrad && varcov_settings$hessianRoutine=='analytic') varcov_settings$hessianRoutine <- 'numDeriv'
  
  ### Extract relevant values
  apollo_probabilities <- varcov_settings$apollo_probabilities
  apollo_inputs  <- varcov_settings$apollo_inputs
  apollo_logLike <- varcov_settings$apollo_logLike
  apollo_grad    <- varcov_settings$apollo_grad
  hessianRoutine <- varcov_settings$hessianRoutine
  dummyVCM       <- matrix(NA, nrow=length(apollo_beta), ncol=length(apollo_beta), 
                           dimnames=list(names(apollo_beta), names(apollo_beta)))
  test <- is.function(apollo_grad) && !is.null(environment(apollo_grad)$silent)
  if(test) silent <- environment(apollo_grad)$silent else {
    test <- is.function(apollo_logLike) && !is.null(environment(apollo_logLike)$silent)
    if(test) silent <- environment(apollo_logLike)$silent else {
      test <- is.list(apollo_inputs) && !is.null(apollo_inputs$silent)
      if(test) silent <- apollo_inputs$silent else silent <- FALSE
    }
  }
  setSMulti <- function(s){
    genv <- globalenv()
    x <- get('apollo_inputs', envir=genv, inherits=FALSE)
    x$apollo_scaling <- s
    assign('apollo_inputs', x, envir=genv) }
  
  ### If hessianRoutine=='none', return dummyVCM
  if(hessianRoutine=="none"){
    varInd <- which(!(names(apollo_beta) %in% apollo_fixed))
    if(length(varInd)==0) stop('SYNTAX ISSUE - All parameters are fixed. Covariance matrix does not exist.')
    L <- list(hessian     = dummyVCM[varInd, varInd, drop=FALSE], 
              varcov      = dummyVCM[varInd, varInd, drop=FALSE], 
              se          = diag(dummyVCM), 
              corrmat     = dummyVCM[varInd, varInd, drop=FALSE], 
              robvarcov   = dummyVCM[varInd, varInd, drop=FALSE], 
              robse       = diag(dummyVCM), 
              robcorrmat  = dummyVCM[varInd, varInd, drop=FALSE],  
              hessianMethodUsed       = hessianRoutine,
              hessianMethodsAttempted = hessianRoutine,
              apollo_beta             = apollo_beta, 
              hessianScaling          = setNames(rep(1, length(apollo_beta)), names(apollo_beta)),
              eigValue                = NA)
    return(L)
  }
  
  ### If apollo_inputs is provided, then try to build apollo_logLike
  if(!is.function(apollo_grad) && !is.function(apollo_logLike) && is.function(apollo_probabilities) && is.list(apollo_inputs)){
    apollo_logLike <- apollo_makeLogLike(apollo_beta          = apollo_beta, 
                                         apollo_fixed         = apollo_fixed, 
                                         apollo_probabilities = apollo_probabilities, 
                                         apollo_inputs        = apollo_inputs,
                                         apollo_estSet        = list(estimationRoutine='bhhh'))
    apollo_grad <- NULL
  }
  if(!is.function(apollo_logLike)) stop('CALCULATION ISSUE - Creation of log-likelihood function (apollo_logLike) failed. Varcov matrix cannot be calculated.')
  
  ### If apollo_grad is not provided, but analytic gradient is requested, then construct analytical gradient
  test <- !is.function(apollo_grad) && hessianRoutine=='analytic'
  if(test) apollo_grad <- apollo_makeGrad(apollo_beta    = apollo_beta, 
                                          apollo_fixed   = apollo_fixed, 
                                          apollo_logLike = apollo_logLike, 
                                          validateGrad   = TRUE)
  
  ### Scale if requested
  if(varcov_settings$scaleBeta){
    # Extract old scaling
    oneCore <- environment(apollo_logLike)$singleCore
    if( oneCore) s <- environment(apollo_logLike)$apollo_inputs$apollo_scaling
    if(!oneCore) s <- parallel::clusterEvalQ(cl=environment(apollo_logLike)$cl,
                                             apollo_inputs$apollo_scaling)[[1]]
    oldScaling <- s
    # Calculate LL before scaling
    b <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    b[names(s)] <- b[names(s)]/s
    test1 <- apollo_logLike(b, countIter=FALSE, sumLL=TRUE, writeIter=FALSE, getNIter=FALSE)
    # Calculate new (hessian) scaling and turn scales equal to 0 into 1 to avoid divisions by zero.
    scaling <- abs(apollo_beta[!(names(apollo_beta) %in% apollo_fixed)])
    if(any(scaling==0)) scaling[scaling==0] <- 1
    # Update apollo_inputs$apollo_scaling
    if( oneCore) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- scaling
    if(!oneCore) parallel::clusterCall(cl=environment(apollo_logLike)$cl, fun=setSMulti, s=scaling)
    # Scale apollo_beta to hessian scaling
    apollo_beta[names(scaling)] <- apollo_beta[names(scaling)]/scaling
    # Calculate LL after scaling and make sure it is the same as before
    b     <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    test2 <- apollo_logLike(b, countIter=FALSE, sumLL=TRUE, writeIter=FALSE, getNIter=FALSE)
    test <- is.numeric(test1) && is.numeric(test2) && !any(is.nan(test1)) && !any(is.nan(test2))
    test <- test && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(is.na(test)||!test){
      # If scaling didn't work, undo changes
      if(!silent)apollo_print("Parameters could not be scaled for the Hessian calculation.")
      apollo_beta[names(scaling)] <- apollo_beta[names(scaling)]*oldScaling
      if(environment(apollo_logLike)$singleCore) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- scaling
      if(!environment(apollo_logLike)$singleCore) parallel::clusterCall(cl=environment(apollo_logLike)$cl, 
                                                                        fun=setSMulti, s=s)
    }
  ### NEW 9 May
  }else{
    oneCore <- environment(apollo_logLike)$singleCore
    if( oneCore) s <- environment(apollo_logLike)$apollo_inputs$apollo_scaling
    if(!oneCore) s <- parallel::clusterEvalQ(cl=environment(apollo_logLike)$cl,
                                             apollo_inputs$apollo_scaling)[[1]]
    apollo_beta[names(s)] <- apollo_beta[names(s)]/s
  }
  
  
  # # # # # # # # # # # # # #
  #### Calculate Hessian ####
  # # # # # # # # # # # # # # 
  
  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  H <- NULL
  methodUsed <- NULL
  
  ### Attempt to calculate Hessian as as the Jacobian of the analytical gradient
  if(is.null(H) && is.function(apollo_grad) && hessianRoutine=='analytic'){
    methodUsed <- c(methodUsed, 'numerical jacobian of LL analytical gradient')
    if(!silent) apollo_print("Computing covariance matrix using analytical gradient.")
    i <- 0; k <- length(beta_var_val); I <- 1+8*k; di <- ceiling(I/20);
    sumGradLL <- function(theta){
      ans <- colSums( apollo_grad(theta) )
      i <<- i+1
      if(!silent && i%%di==0) if(i%%(5*di)==0) cat(i/(5*di)*25,'%',sep='') else cat('.')
      return(ans)
    }
    if(!silent) cat(" 0%")
    H <- tryCatch(numDeriv::jacobian(func=sumGradLL, x=beta_var_val, method.args=varcov_settings$numDeriv_settings),
                  error = function(e) return(NA))
    if(!silent) cat('100%')
    if(anyNA(H)){
      if(!silent && length(H)==1) cat(' Failed')
      if(!silent && is.matrix(H)) cat(' (', sum(is.na(H))," NA values)", sep='')
      H <- NULL
      hessianRoutine <- "numDeriv"
    }; if(!silent) cat('\n')
    rm(i, k, I, di)
  }
  
  ### Attempt to calculate Hessian numerically using numDeriv
  if(is.null(H) && is.function(apollo_logLike) && hessianRoutine=="numDeriv"){
    methodUsed <- c(methodUsed, 'numerical second derivative of LL (using numDeriv)')
    if(!silent) apollo_print("Computing covariance matrix using numerical methods (numDeriv).")
    i <- 0; k <- length(beta_var_val); I <- 2+8*( k*(k+1)/2 ); di <- ceiling(I/20)
    sumLogLike <- function(theta){
      ans <- apollo_logLike(theta, countIter=FALSE, writeIter=FALSE, sumLL=TRUE)
      i <<- i+1
      if(!silent && i%%di==0) if(i%%(5*di)==0) cat(i/(5*di)*25,'%',sep='') else cat('.')
      return(ans)
    }
    if(!silent) cat(' 0%')
    H <- tryCatch(numDeriv::hessian(func=sumLogLike, x=beta_var_val, method.args=varcov_settings$numDeriv_settings),
                  error = function(e) return(NA))
    if(!silent) cat('100%')
    if(anyNA(H)){
      if(!silent && length(H)==1) cat(' Failed')
      if(!silent && is.matrix(H)) cat(' (', sum(is.na(H))," NA values)", sep='')
      H <- NULL
      hessianRoutine <- "maxLik"
    }; if(!silent) cat('\n')
    rm(i, k, I, di)
  }
  
  ### Attempt to calculate Hessian numerically using maxLik
  if(is.null(H) && is.function(apollo_logLike) && hessianRoutine=="maxLik"){
    methodUsed <- c(methodUsed, 'numerical second derivative of LL (using maxLik)')
    if(!silent) apollo_print("Computing covariance matrix using numerical methods (maxLik). This may take a while, no progress bar displayed.")
    H <- tryCatch(maxLik::maxLik(apollo_logLike, start=beta_var_val, method="BFGS", print.level=0,
                                 finalHessian=TRUE, iterlim=2, countIter=FALSE, writeIter=FALSE, sumLL=FALSE)$hessian,
                  error=function(e) return(NA))
    if(anyNA(H)){
      if(!silent && length(H)==1) apollo_print("Hessian calculation failed using numerical methods (maxLik).")
      if(!silent && is.matrix(H)) apollo_print(paste0('(', sum(is.na(H))," NA values)"))
    }
  }
  
  ### Check result
  if(!silent && (is.null(H) || (length(H)==1 && is.na(H))) ){
    apollo_print("Hessian could not be calculated.", type="w")
    methodUsed <- c(methodUsed, 'Failed')
  } 
  if(is.matrix(H)){
    # Add names to Hessian
    colnames(H) <- names(beta_var_val)
    rownames(H) <- names(beta_var_val)
    # De-scale Hessian
    if(varcov_settings$scaleBeta){
      for(i in names(scaling)){
        H[i,] <- H[i,]/scaling[i]
        H[,i] <- H[,i]/scaling[i]
      }
    ### new bit 9 May 2023
      }else{
        oneCore <- environment(apollo_logLike)$singleCore
        if( oneCore) s <- environment(apollo_logLike)$apollo_inputs$apollo_scaling
        if(!oneCore) s <- parallel::clusterEvalQ(cl=environment(apollo_logLike)$cl,
                                                 apollo_inputs$apollo_scaling)[[1]]
        
      if(any(s!=1)){
        for(i in names(s)){
          scale=s[i]
          H[i,] <- H[i,]/scale[i]
          H[,i] <- H[,i]/scale[i]
        }
      }
      }
    ####
    # Check invertibility
    invertible <- tryCatch(is.matrix(solve(H)), error=function(e) FALSE)
    if(!invertible){
      if(!silent) apollo_print('Singular Hessian, cannot calculate s.e.', type="w")
      tryCatch({
        outD <- tryCatch(environment(apollo_logLike)$apollo_inputs$apollo_control$outputDirectory,
                         error=function(e) ".")
        if(!(substr(outD, nchar(outD), nchar(outD)) %in% c('/','\\'))) outD <- paste0(outD,'/')
        fileName <- paste0(outD, environment(apollo_logLike)$modelName, '_hessian.csv')
        utils::write.csv(H, fileName)
        if(!silent) apollo_print(paste0("Hessian written to ", fileName))
        rm(test, fileName)
      }, error=function(e) if(!silent) apollo_print("Could not write hessian to a file."))
    }
    # Calculate eigenvalues
    eigValue <- tryCatch(eigen(H, only.values=TRUE)$values, error=function(e) return(NA))
    if(!silent && is.complex(eigValue)) apollo_print('Some eigenvalues of the Hessian are complex, indicating that the Hessian is not symmetrical.', type="w")
    if(!silent && !anyNA(eigValue) && !is.complex(eigValue)){
      if(any(eigValue>0)) apollo_print("Some eigenvalues of the Hessian are positive, indicating convergence to a saddle point!", type="w") 
      if(all(eigValue<0)) apollo_print(paste0("Negative definite Hessian with maximum eigenvalue: ",round(max(eigValue),6)))
    }
  }
  
  # # # # # # # # # # # # # # # # #
  #### Calculate varcov & s.e. ####
  # # # # # # # # # # # # # # # # #
  
  ### De-scale model$estimate from the hessian scaling
  if(varcov_settings$scaleBeta){
    apollo_beta[names(scaling)] <- apollo_beta[names(scaling)]*scaling
    ### new bit 9 May 2023
  }else{
    oneCore <- environment(apollo_logLike)$singleCore
    if( oneCore) s <- environment(apollo_logLike)$apollo_inputs$apollo_scaling
    if(!oneCore) s <- parallel::clusterEvalQ(cl=environment(apollo_logLike)$cl,
                                             apollo_inputs$apollo_scaling)[[1]]
    if(any(s!=1)) apollo_beta[names(s)] <- apollo_beta[names(s)]*s
  }
  ####
  
  ### Create dummy matrix without fixed params
  varInd <- which(!(names(apollo_beta) %in% apollo_fixed))
  dummyVCM_small <- dummyVCM[varInd, varInd]
  
  ### Calculate variance-covariance matrix, s.e. and correlation matrix
  if(is.matrix(H) && invertible){
    varcov  <- solve(-H)
    varcov  <- (varcov + t(varcov))/2
    se      <- sqrt(diag(varcov))
    corrmat <- tryCatch(varcov/(se%*%t(se)), error=function(e) return(dummyVCM_small))
  } else {
    varcov  <- dummyVCM_small
    se      <- diag(dummyVCM)
    corrmat <- dummyVCM_small
  }
  
  ### Calculate scores
  bVar       <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  newScaling <- setNames(rep(1, length(bVar)), names(bVar))
  if(is.null(varcov_settings$BHHH_matrix)){
    if(!silent) apollo_print('Computing score matrix...')
    # Remove hessian scaling
    # If using analytic
    if(hessianRoutine=='analytic'){
      # Remove hessian scaling
      ### 9 May - should always do this, not just if scaleBeta, even if scaling is just 1
      ###if(varcov_settings$scaleBeta){
      singleCore <- environment(apollo_grad)$singleCore
      if( singleCore) environment(apollo_grad)$apollo_inputs$apollo_scaling <- newScaling
      if(!singleCore) parallel::clusterCall(cl=environment(apollo_grad)$cl, fun=setSMulti, s=newScaling)
      ###}
      # Calculate scores
      score <- apollo_grad(bVar)
    } else { # If using numeric
      # Remove hessian scaling
      ### 9 May - should always do this, not just if scaleBeta, even if scaling is just 1
      ###if(varcov_settings$scaleBeta){
      singleCore <- environment(apollo_logLike)$singleCore
      if( singleCore) environment(apollo_logLike)$apollo_inputs$apollo_scaling <- newScaling
      if(!singleCore) parallel::clusterCall(cl=environment(apollo_logLike)$cl, fun=setSMulti, s=newScaling)
      ###}
      # Calculate scores
      score <- numDeriv::jacobian(apollo_logLike, bVar)
    }
    BHHH_matrix=var(score)*nrow(score)
  }else{
    BHHH_matrix=varcov_settings$BHHH_matrix
  }

  
  ### Calculate robust variance-covariance matrix, s.e. and correlation matrix
  #if(is.matrix(H) && is.matrix(score)){
  if(is.matrix(H) && is.matrix(BHHH_matrix)){
    #bread <- solve(-H)
    #meat  <- split(score, rep(1:nrow(score), ncol(score)))
    #meat  <- Reduce('+', lapply(meat, function(r) r%*%t(r))) # it should include /nrow(score) but with it, it doesn't work.
    bread  <- varcov
    #meat   <- var(score)*nrow(score) # not sure why the *nrow(score) but without it, it doesn't work.
    meat   <- BHHH_matrix
    robvarcov     <- bread %*% meat %*% bread
    robse         <- sqrt(diag(robvarcov))
    robcorrmat    <- tryCatch(robvarcov/(robse%*%t(robse)) , error=function(e) return(dummyVCM_small))
  } else {
    robvarcov  <- dummyVCM_small
    robse      <- diag(dummyVCM)
    robcorrmat <- dummyVCM_small
  }
  
  ### Restore fixed parameters to s.e. and rob s.e. only
  se    <- c(se, apollo_beta[apollo_fixed]*NA)[names(apollo_beta)]
  robse <- c(robse, apollo_beta[apollo_fixed]*NA)[names(apollo_beta)]
  
  
  # # # # # # # # # # # # #
  #### Pack and return ####
  # # # # # # # # # # # # #
  
  L <- list(hessian     = H,
            varcov      = varcov,
            se          = se,
            corrmat     = corrmat,
            robvarcov   = robvarcov, 
            robse       = robse, 
            robcorrmat  = robcorrmat,
            apollo_beta = apollo_beta, 
            hessianMethodUsed  = ifelse(!is.null(methodUsed), methodUsed[length(methodUsed)], 'Failed'),
            hessianMethodsAttempted = methodUsed)
  if(varcov_settings$scaleBeta) L$hessianScaling <- scaling else {
    L$hessianScaling <- environment(apollo_logLike)$apollo_inputs$apollo_scaling
  }
  if(is.matrix(H)){
    L$hessianEigenValue = eigValue
    if(!is.complex(eigValue)) L$eigValue = round(max(eigValue),6)
  }
  return(L)
}