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
#' @importFrom RSGHB doHB
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats sd cor cov runif
#' @export
apollo_estimateHB <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings=NA){
  
  ### First checkpoint
  time1 <- Sys.time()

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
  
  # Checks for scaling
  if(length(apollo_inputs$apollo_scaling)>0 && !is.na(apollo_inputs$apollo_scaling)){
    if(!is.null(apollo_HB$gVarNamesFixed)){
      r <- ( names(apollo_beta) %in% names(apollo_inputs$apollo_scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
      r <- names(apollo_beta)[r]
      apollo_HB$FC[r] <- 1/apollo_inputs$apollo_scaling[r]*apollo_HB$FC[r]
      rm(r)
    }
    if(!is.null(apollo_HB$gVarNamesNormal)){
      r <- ( names(apollo_beta) %in% names(apollo_inputs$apollo_scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
      r <- names(apollo_beta)[r]
      hb<- apollo_HB
      dists_normal= names(hb$gDIST[r][hb$gDIST[r]==1])
      dists_lnp   = names(hb$gDIST[r][hb$gDIST[r]==2])
      dists_lnn   = names(hb$gDIST[r][hb$gDIST[r]==3])
      dists_cnp   = names(hb$gDIST[r][hb$gDIST[r]==4])
      dists_cnn   = names(hb$gDIST[r][hb$gDIST[r]==5])
      dists_sb    = names(hb$gDIST[r][hb$gDIST[r]==6])
      s <- apollo_inputs$apollo_scaling
      if(length(dists_normal)>0) apollo_HB$svN[dists_normal] <- 1/s[dists_normal]*hb$svN[dists_normal]
      if(length(dists_lnp)>0) apollo_HB$svN[dists_lnp] <- -log(s[dists_lnp]) + hb$svN[dists_lnp]
      if(length(dists_lnn)>0) apollo_HB$svN[dists_lnn] <- -log(s[dists_lnn]) + hb$svN[dists_lnn]
      if(length(dists_cnp)>0) apollo_HB$svN[dists_cnp] <- 1/s[dists_cnp]*hb$svN[dists_cnp]
      if(length(dists_cnn)>0) apollo_HB$svN[dists_cnn] <- 1/s[dists_cnn]*hb$svN[dists_cnn]
      if(length(dists_sb)>0){
        names(apollo_HB$gMINCOEF)=names(hb$svN)
        names(apollo_HB$gMAXCOEF)=names(hb$svN)
        apollo_HB$gMINCOEF[dists_sb] <- 1/s[dists_sb]*hb$gMINCOEF[dists_sb]
        apollo_HB$gMAXCOEF[dists_sb] <- 1/s[dists_sb]*hb$gMAXCOEF[dists_sb]
      }
      rm(r, hb, dists_normal, dists_lnp, dists_lnn, dists_cnp, dists_cnn, dists_sb)
    }
  }
  
  ### Insert component names
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  
  ### Start clock
  starttime <- Sys.time()
  
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  if(!apollo_control$noValidation){
    ### Validation using HB estimation
    
    if(!silent) cat("Testing probability function (apollo_probabilities)\n")
    
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
    testLL = apollo_probabilities(apollo_test_beta, apollo_inputs, functionality="estimate")
    if(!workInLogs) testLL=log(testLL)
    # Maybe here we could return the value of the likelihood and print and error with cat, instead of simply stopping
    if(anyNA(testLL)) stop('Log-likelihood calculation fails at starting values!')
    ### Test for unused parameters
    #apollo_beta_base=apollo_beta+0.001
    apollo_beta_base=apollo_test_beta+(!(names(apollo_beta)%in%(apollo_fixed)))*0.001*runif(length(apollo_beta))
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
    
  }
  
  # ################################## #
  #### HB estimation                ####
  # ################################## #
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
  
  ### Second checkpoint
  time2 <- Sys.time()

  model <- RSGHB::doHB(apollo_HB_likelihood, database, apollo_HB)
  
  ### Third checkpoint
  time3 <- Sys.time()

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
    colnames(model$posterior_mean)=c("Mean","SD")
    
    ### create matrix of draws from distributions
    
    draws=10000
    covMat=random_cov_mean
    meanA=random_mean[,1]
    pars = length(meanA)
    covMat=as.matrix(covMat)
    Ndraws=mvtnorm::rmvnorm(draws,meanA,covMat,method="chol")
    for(i in 1:pars){
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
    }
    
  }
  
  if(length(scaling)>0 && !is.na(scaling)){
    for(s in 1:length(scaling)){
      ss=names(scaling)[s]
      if(ss%in%colnames(model$C)) model$C[,ss]=scaling[s]*model$C[,ss]
      if(ss%in%colnames(model$Csd)) model$Csd[,ss]=scaling[s]*model$Csd[,ss]
      if(ss%in%colnames(model$F)) model$F[,ss]=scaling[s]*model$F[,ss]
      if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){if(ss%in%colnames(Ndraws)) Ndraws[,ss]=scaling[s]*Ndraws[,ss]}
      if(ss%in%rownames(model$chain_non_random)) model$chain_non_random[ss,]=scaling[s]*model$chain_non_random[ss,]
      if(ss%in%rownames(model$posterior_mean)) model$posterior_mean[ss,]=scaling[s]*model$posterior_mean[ss,]
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
  if(exists("apollo_HBcensor", envir=globalenv()) && !apollo_inputs$silent){
    apollo_print(paste0('WARNING: RSGHB has censored the probabilities. ',
                        'Please note that in at least some iterations RSGHB has ',
                        'avoided numerical issues by left censoring the ',
                        'probabilities. This has the side effect of zero or ',
                        'negative probabilities not leading to failures!'))
    rm("apollo_HBcensor", envir=globalenv())
  }
  
  ### Fourth checkpoint
  time4 <- Sys.time()
 
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))
  model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))
  model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))
  
  return(model)
}