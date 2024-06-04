#' Estimates model using Bayesian estimation
#'
#' Estimates a model using Bayesian estimation on the likelihood function defined by \code{apollo_probabilities}.
#'
#' This is a sub function of \link{apollo_estimate} which is called when using Bayesian estimation. 
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
#' @param estimate_settings List. Options controlling the estimation process, as used for in \link{apollo_estimate}.
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
  scaling           = apollo_inputs$apollo_scaling    #estimate_settings[["scaling"]]
  bootstrapSE       = estimate_settings[["bootstrapSE"]]
  bootstrapSeed     = estimate_settings[["bootstrapSeed"]]
  
  ### added 1 Feb
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                       apollo_probabilities, apollo_inputs, 
                                       apollo_estSet=estimate_settings, cleanMemory=TRUE)
  
  # Checks for scaling
  if(!is.null(apollo_inputs$apollo_scaling) && !anyNA(apollo_inputs$apollo_scaling)){
    if(!is.null(apollo_HB$gVarNamesFixed)){
      r <- ( names(apollo_beta) %in% names(apollo_inputs$apollo_scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
      r <- names(apollo_beta)[r]
      apollo_HB$FC[r] <- 1/apollo_inputs$apollo_scaling[r]*apollo_HB$FC[r]
      rm(r)
    }
    if(!is.null(apollo_HB$gVarNamesNormal)){
      r <- ( names(apollo_beta) %in% names(apollo_inputs$apollo_scaling) ) & ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
      r <- names(apollo_beta)[r]
      dists_normal= names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==1])
      dists_lnp   = names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==2])
      dists_lnn   = names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==3])
      dists_cnp   = names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==4])
      dists_cnn   = names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==5])
      dists_sb    = names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==6])
      s <- apollo_inputs$apollo_scaling
      if(length(dists_normal)>0) apollo_HB$svN[dists_normal] <- 1/s[dists_normal]*apollo_HB$svN[dists_normal]
      if(length(dists_lnp)>0) apollo_HB$svN[dists_lnp] <- -log(s[dists_lnp]) + apollo_HB$svN[dists_lnp]
      if(length(dists_lnn)>0) apollo_HB$svN[dists_lnn] <- -log(s[dists_lnn]) + apollo_HB$svN[dists_lnn]
      if(length(dists_cnp)>0) apollo_HB$svN[dists_cnp] <- 1/s[dists_cnp]*apollo_HB$svN[dists_cnp]
      if(length(dists_cnn)>0) apollo_HB$svN[dists_cnn] <- 1/s[dists_cnn]*apollo_HB$svN[dists_cnn]
      if(length(dists_sb)>0){
        names(apollo_HB$gMINCOEF)=names(apollo_HB$svN)  
        names(apollo_HB$gMAXCOEF)=names(apollo_HB$svN)
        apollo_HB$gMINCOEF[dists_sb] <- 1/s[dists_sb]*apollo_HB$gMINCOEF[dists_sb]
        apollo_HB$gMAXCOEF[dists_sb] <- 1/s[dists_sb]*apollo_HB$gMAXCOEF[dists_sb]
      }
      rm(r, dists_normal, dists_lnp, dists_lnn, dists_cnp, dists_cnn, dists_sb)
    }
  }
  
  ### Insert component names
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  
  ### Start clock
  starttime <- Sys.time()
  
  
  # ####################################### #
  #### Validation of likelihood function ####
  # ####################################### #
  apollo_test_beta=apollo_beta
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
      rm(dists_normal, dists_lnp, dists_lnn, dists_cnp, dists_cnn, dists_sb)
    }
    apollo_probabilities(apollo_test_beta, apollo_inputs, functionality="validate")
    testLL = apollo_probabilities(apollo_test_beta, apollo_inputs, functionality="estimate")
    if(!workInLogs) testLL=log(testLL)
    # Maybe here we could return the value of the likelihood and print and error with cat, instead of simply stopping
    if(anyNA(testLL)) stop('CALCULATION ISSUE - Log-likelihood calculation fails at starting values!')
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
        test1_LL=sum(log( ifelse(!is.finite(test1_LL) | test1_LL<=0, .1, test1_LL) ))
        test2_LL=sum(log( ifelse(!is.finite(test2_LL) | test2_LL<=0, .1, test2_LL) ))
      }
      if(is.na(test1_LL)) test1_LL <- base_LL + 1 # Avoids errors if test1_LL is NA
      if(is.na(test2_LL)) test2_LL <- base_LL + 2 # Avoids errors if test2_LL is NA
      if(base_LL==test1_LL & base_LL==test2_LL) stop("SPECIFICATION ISSUE - Parameter ",p," does not influence the log-likelihood of your model!")
    }
    
  }
  
  # ################################## #
  #### HB estimation                ####
  # ################################## #
  ### Check that apollo_fixed hasn't changed since calling apollo_validateInputs
  tmp <- tryCatch( get("apollo_fixed", envir=globalenv()), error=function(e) 1 )
  if( length(tmp)>0 && any(tmp %in% c(apollo_HB$gVarNamesFixed, apollo_HB$gVarNamesNormal)) ) stop("INTERNAL ISSUE - apollo_fixed seems to have changed since calling apollo_inputs.")
  
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
  
  # ### Change ids to numeric ones if necessary
  # if(!is.numeric(database[,apollo_control$indivID])){
  #   originalID <- database[,apollo_control$indivID]
  #   idDict <- unique(originalID)
  #   idDict <- stats::setNames(1:length(idDict), idDict)
  #   database[,apollo_control$indivID] <- idDict[originalID]
  #   apollo_inputs$database[,apollo_control$indivID] <- idDict[originalID]
  # }
  
  ### Actual estimation with RSGHB
  ## removed 2 May
  currentWD <- getwd()  
  if(dir.exists(apollo_inputs$apollo_control$outputDirectory)) setwd(apollo_inputs$apollo_control$outputDirectory)
  model <- RSGHB::doHB(apollo_HB_likelihood, database, apollo_HB)
  setwd(currentWD)
  
  ### Rename RSGHB components in model object
  if(!is.null(model$A)) names(model)[which(names(model)=="A")]="HB_iterations_means"
  if(!is.null(model$B)) names(model)[which(names(model)=="B")]="HB_indiv_draws_means"
  if(!is.null(model$Bsd)) names(model)[which(names(model)=="Bsd")]="HB_indiv_draws_sd"
  if(!is.null(model$C)) names(model)[which(names(model)=="C")]="HB_posterior_means"
  if(!is.null(model$Csd)) names(model)[which(names(model)=="Csd")]="HB_posterior_sd"
  if(!is.null(model$D)) names(model)[which(names(model)=="D")]="HB_iterations_covar"
  if(!is.null(model$F)) names(model)[which(names(model)=="F")]="HB_iterations_non_random"
  if(!is.null(model$params.fixed)) names(model)[which(names(model)=="params.fixed")]="HB_names_nonrandom_params"
  if(!is.null(model$params.vary)) names(model)[which(names(model)=="params.vary")]="HB_names_random_params"
  if(!is.null(model$iter.detail)) names(model)[which(names(model)=="iter.detail")]="HB_iterations_detail"
  if(!is.null(model$cmcLLout)) model$cmcLLout=NULL
  if(!is.null(model$cmcRLHout)) model$cmcRLHout=NULL
  if(!is.null(model$distributions)) model$distributions=NULL
  #if(!is.null(model$gNCREP)) model$gNCREP=NULL
  #if(!is.null(model$gNEREP)) model$gNEREP=NULL
  if(!is.null(model$gNOBS)) model$gNOBS=NULL
  if(!is.null(model$gNP)) model$gNP=NULL
  if(!is.null(model$gSeed)) model$gSeed=NULL
  if(!is.null(model$gSIGDIG)) model$gSIGDIG=NULL
  if(!is.null(model$llf)) model$llf=NULL
  if(!is.null(model$ll0)) model$llo=NULL
  
  
  
  
  # ### Change ids back to the original ones
  # if(!is.numeric(originalID)){
  #   database[,apollo_control$indivID]               <- originalID
  #   apollo_inputs$database[,apollo_control$indivID] <- originalID
  # }
  
  ### Third checkpoint
  time3 <- Sys.time()

  ### added 1 Feb
  model$modelTypeList        <- environment(apollo_logLike)$mType
  model$nObsTot              <- environment(apollo_logLike)$nObsTot
  
  model$apollo_HB   <- apollo_HB
  ### use pre-scaling values as starting values in output 
  model$apollo_beta <- apollo_test_beta
  model$LLStart <- sum(testLL)
  ### model$LL0         <- sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL")))
  #if(workInLogs) model$LL0 <- sum((apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL"))) else model$LL0 <- sum(log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="zero_LL")[["model"]]))
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
    if(workInLogs) model$LL0=sapply(model$LL0,sum)
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
    } else {
      model$LLC <- ifelse(workInLogs, sum(model$LLC), sum(log(model$LLC)) )
    }
  } else {
    model$LLC=NA
  }
  
  model$startTime   <- starttime
  model$apollo_control <- apollo_control
  model$nObs        <- nrow(database)
  model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
  endtime           <- Sys.time()
  model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
  model$apollo_fixed <- apollo_fixed
  model$estimationRoutine <- "Hierarchical Bayes"
  
  if(!is.null(model$HB_iterations_non_random)){
    tmp <- coda::geweke.diag(model$HB_iterations_non_random[,2:(ncol(model$HB_iterations_non_random))], frac1=0.1, frac2=0.5)[[1]]
    names(tmp) <- model$HB_names_nonrandom_params
    model$HB_Geweke_test_non_random=tmp
    rm(tmp)
  }
  if(!is.null(model$HB_iterations_means)){
    tmp <- coda::geweke.diag(model$HB_iterations_means[,2:(ncol(model$HB_iterations_means))], frac1=0.1, frac2=0.5)[[1]]
    model$HB_Geweke_test_means=tmp
    rm(tmp)
  }
  if(!is.null(model$HB_iterations_covar)){
    # This assumes the matrix is square
    tmp <- c()
    for(i in 1:dim(model$HB_iterations_covar)[1]) for(j in 1:i){
      if(i==1 & j==1) Dmatrix <- as.matrix(model$HB_iterations_covar[i,j,]) else Dmatrix <- cbind(Dmatrix, as.vector(model$HB_iterations_covar[i,j,]))
      tmp <- c(tmp, paste(colnames(model$HB_iterations_means)[i+1],colnames(model$HB_iterations_means)[j+1], sep="_"))
    }
    colnames(Dmatrix) <- tmp
    tmp <- coda::geweke.diag(Dmatrix, frac1=0.1, frac2=0.5)[[1]]
    model$HB_Geweke_test_covar=tmp
    rm(tmp)
  }
  
  if(length(apollo_HB$gVarNamesFixed)>0 | length(model$apollo_fixed)>0){
    if(length(apollo_HB$gVarNamesFixed)>0){
      non_random=matrix(0,nrow=length(apollo_HB$gVarNamesFixed),2)
      non_random[,1]=colMeans(model$HB_iterations_non_random)[2:ncol(model$HB_iterations_non_random)]
      non_random[,2]=apply(model$HB_iterations_non_random,FUN=stats::sd,2)[2:ncol(model$HB_iterations_non_random)]
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
    model$HB_chains_non_random=non_random[originalOrder,,drop=FALSE]
    rm(non_random)
  }
  
  apollo_HB$gVarNamesFixed <- model$HB_names_nonrandom_params
  apollo_HB$gVarNamesNormal <- model$HB_names_random_params
  if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
    random_mean     = matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
    random_mean[,1] = colMeans(model$HB_iterations_means)[2:ncol(model$HB_iterations_means)]
    random_mean[,2] = apply(model$HB_iterations_means,FUN=stats::sd,2)[2:ncol(model$HB_iterations_means)]
    rownames(random_mean)=apollo_HB$gVarNamesNormal
    colnames(random_mean)=c("Mean","SD")
    model$HB_chains_normals_means=random_mean
    
    random_cov_mean           = apply(model$HB_iterations_covar,FUN=mean,c(1,2))
    random_cov_sd             = apply(model$HB_iterations_covar,FUN=stats::sd,c(1,2))
    rownames(random_cov_mean) = apollo_HB$gVarNamesNormal
    colnames(random_cov_mean) = apollo_HB$gVarNamesNormal
    model$HB_chains_normals_mean_of_covar=random_cov_mean
    
    rownames(random_cov_sd) = apollo_HB$gVarNamesNormal
    colnames(random_cov_sd) = apollo_HB$gVarNamesNormal
    model$HB_chains_normals_sd_of_covar=random_cov_sd
    
    posterior=matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
    posterior[,1]=colMeans(model$HB_posterior_means)[3:ncol(model$HB_posterior_means)]
    posterior[,2]=apply(model$HB_posterior_means,FUN=stats::sd,2)[3:ncol(model$HB_posterior_means)]
    rownames(posterior)=apollo_HB$gVarNamesNormal
    model$HB_posterior_means_summary=posterior
    colnames(model$HB_posterior_means_summary)=c("Mean","SD")
    
    ### create matrix of draws from distributions
    
    draws=10000
    covMat=random_cov_mean
    meanA=random_mean[,1]
    pars = length(meanA)
    covMat=as.matrix(covMat)
    Ndraws=mvtnorm::rmvnorm(draws,meanA,covMat,method="chol")
    ### add names if only single random coeff (bug fix 3 Sept)
    if(length(meanA)==1) colnames(Ndraws)=rownames(random_mean)
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
  
  if(length(scaling)>0 && !anyNA(scaling)){
    for(s in 1:length(scaling)){
      ss=names(scaling)[s]
      if(ss%in%colnames(model$HB_posterior_means)) model$HB_posterior_means[,ss]=scaling[s]*model$HB_posterior_means[,ss]
      if(ss%in%colnames(model$HB_posterior_sd)) model$HB_posterior_sd[,ss]=scaling[s]*model$HB_posterior_sd[,ss]
      if(ss%in%colnames(model$HB_iterations_non_random)) model$HB_iterations_non_random[,ss]=scaling[s]*model$HB_iterations_non_random[,ss]
      if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){if(ss%in%colnames(Ndraws)) Ndraws[,ss]=scaling[s]*Ndraws[,ss]}
      if(ss%in%rownames(model$HB_chains_non_random)) model$HB_chains_non_random[ss,]=scaling[s]*model$HB_chains_non_random[ss,]
      if(ss%in%rownames(model$HB_posterior_means_summary)) model$HB_posterior_means_summary[ss,]=scaling[s]*model$HB_posterior_means_summary[ss,]
    }
    model$scaling <- scaling
  }
  
  if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
    model$HB_random_params_mean_sd=cbind(colMeans(Ndraws),apply(Ndraws,2,sd))
    colnames(model$HB_random_params_mean_sd)=c("Mean","SD")
    if(length(apollo_HB$gVarNamesNormal)>1){
      model$HB_random_params_covar=cov(Ndraws)
      model$HB_random_params_corr=cor(Ndraws)
    }
  }
  
  ### produce model$estimate
  
  panelData <- apollo_control$panelData
  indivID   <- database[,apollo_control$indivID]
  nObs  <- nrow(apollo_inputs$database)
  nIndiv <- length(unique(indivID))
  indiv <- unique(apollo_inputs$database[,apollo_inputs$apollo_control$indivID])
  nObsPerIndiv <- setNames(sapply(as.list(indiv),function(x) sum(indivID==x)),indiv)
  
  if(is.null(model$HB_chains_non_random)){
    fc1 <- NULL
  }else{
    ###fc1 <- stats::setNames(as.list(model$HB_chains_non_random[,1]),names(model$HB_chains_non_random[,1]))
    ### cgabge 7 Feb 2023
    fc1 <- stats::setNames(as.list(model$HB_chains_non_random[,1]),rownames(model$HB_chains_non_random))
  }
  
  if(is.null(model$HB_posterior_means)){
    b1 <- NULL
  }else{
    M=model$HB_posterior_means[,-c(1,2),drop=FALSE]
    M1 <- matrix(0, nrow=nObs, ncol=ncol(M))
    r1 <- 1
    for(i in 1:nIndiv){
      r2 <- r1 + nObsPerIndiv[i] - 1
      M1[r1:r2,] <- matrix(as.vector(M[i,]), nrow=r2-r1+1, ncol=ncol(M), byrow=TRUE)
      r1 <- r2 + 1
    }
    b1  <- stats::setNames(as.data.frame(M1), colnames(M))
  } 
  model$estimate=c(fc1,b1)
  
  # Report if the probs have been censored
  if(exists("apollo_HBcensor", envir=globalenv()) && !apollo_inputs$silent){
    apollo_print(paste0('RSGHB has censored the probabilities. ',
                        'Please note that in at least some iterations RSGHB has ',
                        'avoided numerical issues by left censoring the ',
                        'probabilities. This has the side effect of zero or ',
                        'negative probabilities not leading to failures!'), type="w")
    rm("apollo_HBcensor", envir=globalenv())
  }
  
  ### New calculations for log-likelihood (29 Jan 2022)
  
  if(!silent) apollo_print("Calculating LL of each model component... ")
  if(!silent && !is.null(b1)) apollo_print("...as you have used Bayesian estimation, sampling will be used, which may take a while... ")
  
  if(is.null(b1)){
     Lout <- tryCatch(apollo_probabilities(fc1, apollo_inputs, functionality="output"),
                      error=function(e){ apollo_print("Could not complete validation using estimated parameters."); return(NA) })
     for(s in 1:length(Lout)){
       Lout[[s]]=exp(rowsum(log(Lout[[s]]), group=indivID, reorder=FALSE))
     }
    }else{
      LLdraws=500
      Lout=NULL
      for(s in 1:LLdraws){
        pars=Ndraws[s,]
        tmp <- tryCatch(apollo_probabilities(c(fc1,pars), apollo_inputs, functionality="output"),
                      error=function(e){ apollo_print("Could not complete validation using estimated parameters."); return(NA) }) 
        for(s in 1:length(tmp)){
          tmp[[s]]=exp(rowsum(log(tmp[[s]]), group=indivID, reorder=FALSE))
        }
        if(is.null(Lout)){
          Lout=tmp
          for(s in 1:length(Lout)) Lout[[s]]=Lout[[s]]/LLdraws ## for averaging across draws
        }else{
          for(s in 1:length(Lout)){
            Lout[[s]]=Lout[[s]]+tmp[[s]]/LLdraws ## for averaging across draws
          }
        }
      }
    }
  
  if(!anyNA(Lout) && is.list(Lout)){
    Lout <- Lout[c("model",names(Lout)[names(Lout)!="model"])]
    LLout <- Lout
    for(s in 1:length(LLout)){
      if(is.list(LLout[[s]])) LLout[[s]]=LLout[[s]][["model"]]
    }
    if(!workInLogs) LLout <- lapply(LLout, log)
    ind_L=Lout[["model"]]
    ind_LL=log(ind_L)
    model$maximum=sum(ind_LL)
    if(apollo_control$panelData){
      model$avgLL <- setNames(ind_LL/nObsPerIndiv, names(nObsPerIndiv))
      model$avgCP <- setNames(ind_L^(1/nObsPerIndiv), names(nObsPerIndiv))
    } else {
      model$avgLL <- setNames(ind_LL, indiv)
      model$avgCP <- setNames(ind_L , indiv)
    }
    LLout <-lapply(LLout,sum)
    #model$Pout  <- LLout
    model$LLout <- unlist(LLout)
  } else{
    #model$Pout  <- LLout
    model$LLout <- list(NA)
    if(!silent) apollo_print("LL could not be calculated for all components.")
  }
  
  ### Calculate Rho2, AIC, and BIC
  nRandom <- length(model$apollo_HB$gVarNamesNormal)
  if(nRandom>0){
    random_means = nRandom-sum(!is.na(model$apollo_HB$fixedA)) ### fixed means
    random_covar = (nRandom 
                    +model$apollo_HB$gFULLCV*(nRandom*(nRandom-1)/2) ### off-diagonal covar
                    -sum(!is.na(model$apollo_HB$fixedD))) ### fixed vars
  } else {
    random_means = 0
    random_covar = 0
  }
  nFreeParams <- length(model$apollo_HB$gVarNamesFixed) + random_means + random_covar
  test <- !is.null(model$modelTypeList) && !anyNA(model$LL0[1])
  test <- test && all(tolower(model$modelTypeList) %in% c("mnl", "nl", "cnl", "el", "dft", "lc", "rrm"))
  if(test & !silent) apollo_print("Calculating other model fit measures")
  if(test) model$rho2_0 <- 1-(model$maximum/model$LL0[1]) else model$rho2_0 <- NA
  if(test) model$adjRho2_0 <- 1-((model$maximum-nFreeParams)/model$LL0[1]) else model$adjRho2_0 <- NA
  test <- test && is.numeric(model$LLC) && !anyNA(model$LLC[1])
  if(test) model$rho2_C <- 1-(model$maximum/model$LLC[1]) else model$rho2_C <- NA
  if(test) model$adjRho2_C <- 1-((model$maximum-nFreeParams)/model$LLC[1]) else model$adjRho2_C <- NA
  model$AIC <- -2*model$maximum + 2*nFreeParams
  model$BIC <- -2*model$maximum + nFreeParams*log(ifelse(!is.null(model$nObsTot), model$nObsTot, model$nObs))
  model$nFreeParams <- nFreeParams
  model$HB_n_nonrandom <- nFreeParams-random_means-random_covar
  model$HB_n_random_means <- random_means
  model$HB_n_random_covar <- random_covar
  
  ### Fourth checkpoint
  time4 <- Sys.time()
 
  model$timeTaken <- as.numeric(difftime(time4,time1,units='secs'))
  model$timePre   <- as.numeric(difftime(time2,time1,units='secs'))
  model$timeEst   <- as.numeric(difftime(time3,time2,units='secs'))
  model$timePost  <- as.numeric(difftime(time4,time3,units='secs'))
  
  ### assign apollo class to model
  class(model)<-c("apollo",class(model))  
  return(model)
}