#' Estimates model
#'
#' Estimates a model using the likelihood function defined by \code{apollo_probabilities}.
#'
#' This is the main function of the cmcRcode package. The estimation process begins by checking the definition of
#' \code{apollo_probabilities} by estimating it at the starting values. Then it runs the function with argument \code{functionality="validate"}.
#' If the user requested more than one core for estimation (i.e. \code{apollo_control$nCores>1}), and no bayesian estimation is used
#' (i.e. \code{apollo_control$HB=FALSE}), then a cluster is created. Using a cluster at least doubles the requires RAM, as the database
#' must be copied into the cluster.
#' If all checks are passed, estimation begins. There is no limit to estimation time other than reaching the maximum number of
#' iterations. If bayesian estimation is used, estimation will finish once the predefined number of draws are completed.
#' This functions does not save results into a file nor prints them into the console, so if users want to see and store estimation the results,
#' they must make sure to call function \code{apollo_modeloutput} and \code{apollo_saveoutput} afterwards.
#'
#' @param apollo_beta Named numeric vector. Names and starting values for parameters
#' @param apollo_beta_fixed Character vector. Names of parameters to be fixed at starting values.
#' @param database data.frame. Data used by model.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive seven arguments:
#'                          \describe{
#'                            \item{apollo_beta}{Named numeric vector. Names and values of model parameters.}
#'                            \item{database}{Same as described above.}
#'                            \item{apollo_control}{Same as described below.}
#'                            \item{draws}{Same as described below. Optional. NA by default.}
#'                            \item{apollo_randcoeff}{Same as described below. NA by default.}
#'                            \item{functionality}{Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".}
#'                            \item{work_in_logs}{Same as described below. FALSE by default.}
#'                          }
#' @param apollo_control List. Options controlling the running of the code.
#'                    \describe{
#'                      \item{modelName}{Character. Name of the model. Used when saving the output to files.}
#'                      \item{modelDescr}{Character. Description of the model. Used in output files.}
#'                      \item{indivID}{Character. Name of column in the database with each decision maker's ID.}
#'                      \item{mixing}{Boolean. TRUE for models that include random parameters.}
#'                      \item{nCores}{Numeric>0. Number of cores to use in calculations of the model likelihood.}
#'                      \item{seed}{Numeric. Seed for random number generation.}
#'                      \item{HB}{Boolean. TRUE if using RSGHB for Bayesian estimation of model.}
#'                      \item{noValidation}{Boolean. TRUE if user does not wish model input to be validated before estimation - FALSE by default.}
#'                      \item{noDiagnostics}{Boolean. TRUE if user does not wish model diagnostics to be printed - FALSE by default.}
#'                      \item{panel}{Boolean. TRUE if using panel data (created automatically by \code{apollo_validatecontrol)}}
#'                    }
#' @param apollo_draws List of arguments describing the inter and intra individual draws.
#'                  \describe{
#'                    \item{inter_drawsType}{Character. Type of inter-individual draws ('MLHS', 'halton' or 'pmc').}
#'                    \item{inter_nDraws}{Numeric scalar (>=0). Number of inter-individual draws per individual. Set to 0 if not using them.}
#'                    \item{inter_unifDraws}{Character vector. Names of uniform-distributed inter-individual draws.}
#'                    \item{inter_normDraws}{Character vector. Names of normaly distributed inter-individual draws.}
#'                    \item{intra_drawsType}{Character. Type of intra-individual draws ('MLHS', 'halton' or 'pmc').}
#'                    \item{intra_nDraws}{Numeric scalar (>=0). Number of intra-individual draws per individual. Set to 0 if not using them.}
#'                    \item{intra_unifDraws}{Character vector. Names of uniform-distributed intra-individual draws.}
#'                    \item{intra_normDraws}{Character vector. Names of normaly distributed intra-individual draws.}
#'                  }
#' @param apollo_randcoeff Function. Used with mixing models. Constructs the random parameters of a mixing model. Receives two arguments:
#'                      \describe{
#'                        \item{apollo_beta}{Same as described above}
#'                        \item{draws}{Same as described above}
#'                      }
#'                      Takes value NA if omitted. Required for models with mixing.
#' @param apollo_lcpars Function. Used with latent class models. Constructs a list of parameters for each latent class. Receives two arguments:
#'                      \describe{
#'                        \item{apollo_beta}{Same as described above}
#'                        \item{randcoeff}{The output of function \code{apollo_randcoeff}, as described above.
#'                                         If the model does not contain mixing, then defaults to \code{NA}.}
#'                      }
#'                      Takes value NA if omitted. Required for models with mixing.
#' @param HB_control List. Contains options for bayesian estimation. Required if \code{apollo_control$HB} is TRUE.
#' @param estimation_routine Character. Estimation method. Can take values "bfgs" (recommended), "bhhh", or "nr".
#'                           Used only if \code{apollo_control$HB} is FALSE.
#' @param max_iterations Numeric. Maximum number of iterations of the search algorithm before stopping.
#'                       Used only if \code{apollo_control$HB} is FALSE.
#' @param work_in_logs Boolean. TRUE for higher numeric stability at the expense of computational time. Useful for panel models only. Default is FALSE.
#' @param numDerivHessian Boolean. Determines the R package used to calculate the Hessian after estimation. If TRUE, the code{numDeriv}
#'                        package is used. If FALSE, the code{maxLik} package is used. Only used if \code{apollo_control$HB=FALSE}. Default is TRUE.
#' @param printLevel Numeric scalar. Can take integer values 0, 1, 2 or 3. The higher the number the more information is printed during estimation.
#'                   Ignored if apollo_control$HB is TRUE. Default value is 3.
#' @param silent Boolean. If TRUE, no information is printed to console during estimation. Default is FALSE.
#' @return model object
#'
#' @examples
#' ### Set core controls
#' apollo_control = list(
#'   modelName ="MNL", # Make sure to use a new name for every model
#'   indivID   ="ID",  # Name of column in the database with each individual's ID
#'   mixing    = FALSE,# TRUE for models that include random parameters
#'   nCores    = 1     # Number of cores to use in estimation
#' )
#'
#' ### Load data
#' data(apollo_modeChoiceData)
#'
#' ### Model parameters
#' apollo_beta = c(asc_1=0, asc_2=0,
#'                 asc_3=0, asc_4=0,
#'                 tt   =0, tc   =0,
#'                 acc  =0)
#'
#' ### Name of parameters fixed to starting values.
#' apollo_beta_fixed = c("asc_2")
#'
#' ### Likelihood function (do not change the arguments)
#' ### b contains the parameters, x contains the explanatory variables
#' apollo_probabilities=function(b, x, functionality="estimate"){
#'   P <- list() ### Do not delete. Store probabilities here.
#'
#'   ### Enumerate alternatives and availability, and select choice variable.
#'   alternatives = c(car=1, bus=2, air=3, rail=4)
#'   avail        = list(car=x$av_car, bus=x$av_bus, air=x$av_air, rail=x$av_rail)
#'   choiceVar    = x$choice
#'
#'   ### List of utilities
#'   V = list()
#'   V[['car' ]] = b$asc_1 + b$tt*x$time_car  + b$tc*x$cost_car
#'   V[['bus' ]] = b$asc_2 + b$tt*x$time_bus  + b$tc*x$cost_bus  + b$acc*x$access_bus
#'   V[['air' ]] = b$asc_3 + b$tt*x$time_air  + b$tc*x$cost_air  + b$acc*x$access_air
#'   V[['rail']] = b$asc_4 + b$tt*x$time_rail + b$tc*x$cost_rail + b$acc*x$access_rail
#'
#'   ### Compute choice probabilities using MNL model
#'   P[['model']] = apollo_mnl(alternatives, avail, choiceVar, V, functionality)
#'
#'   return(P)
#' }
#'
#' ### Estimate model
#' model = apollo_estimate(apollo_beta, apollo_beta_fixed, database,
#'                         apollo_probabilities, apollo_control)
#'
#' ### Show output in screen
#' apollo_modeloutput(model)
#'
#' @export
apollo_estimate <- function(apollo_beta, apollo_beta_fixed, database,
                         apollo_probabilities, apollo_control, apollo_draws=NA,
                         apollo_randcoeff=NA, apollo_lcpars=NA,
                         HB_control=NA, estimation_routine="bfgs",
                         max_iterations=200,
                         work_in_logs=FALSE, numDerivHessian=TRUE,
                         printLevel=3L, silent=FALSE){

  # ################################## #
  #### initial processing & testing ####
  # ################################## #

  writeIter=FALSE
  apollo_control = apollo_validatecontrol(database,apollo_control)
  database = apollo_validatedata(database, apollo_control)
  if(apollo_control$mixing) draws=apollo_makeDraws(apollo_control, apollo_draws, database) else draws = NA
  assign("apollo_control", apollo_control, envir=environment(apollo_probabilities))
  if(apollo_control$HB) HB_control=apollo_validateHBcontrol(HB_control)

  estimation_routine <- tolower(estimation_routine)
  if( !(estimation_routine %in% c("bfgs","bhhh", "nr")) ) stop("Invalid estimation_routine. Use 'bfgs', 'bhhh' or 'nr'.")
  if((length(apollo_beta_fixed)>0) & any(!(apollo_beta_fixed %in% names(apollo_beta)))) stop("Some parameters included in 'apollo_beta_fixed' are not included in 'apollo_beta'")
  if(!is.integer(printLevel)) printLevel <- as.integer(round(printLevel,0))
  if(printLevel<0L) printLevel <- 0L
  if(3L<printLevel) printLevel <- 3L

  if(max_iterations<1) stop("Need at least one iteration!")
  max_iterations=round(max_iterations,0)
  if(work_in_logs!=TRUE) work_in_logs=FALSE
  if(estimation_routine!="bfgs" & writeIter==TRUE){
    writeIter = FALSE
    cat("witeIter set to FALSE. Writing parameters values\nat each iteration is only available\nfor BFGS estrimation method.\n")
  }

  if(!apollo_control$mixing | apollo_control$HB) {
    draws <- NA
    apollo_randcoeff <- NA
  }
  if(!apollo_control$HB) HB_control <- NA
  #if(!exists("lc_pars", inherits=FALSE)) lc_pars <- NA

  if(apollo_control$HB && anyNA(HB_control)) stop("Argument 'HB_control' must be provided when using Bayesian estimation (see ?apollo_validateHBcontrol).")

  if(apollo_control$mixing){
    if(anyNA(draws)) stop("Argument 'draws' must be provided when estimating mixing models. Use apollo_makeDraws.")
    if(!is.function(apollo_randcoeff)) stop("Argument 'apollo_randcoeff' must be provided when estimating mixing models.")
    if(!apollo_control$panel & dim(draws[[1]])[2]>1) warning('Inter-person draws are used without a panel structure. This is unusual.')
  }

  tempOutputFile <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
  tempOutputFile <- file.path(tempdir(),tempOutputFile)
  if(file.exists(tempOutputFile)){
    tryCatch( file.remove(tempOutputFile),
              error=function(e) cat("Could not delete old file ",tempOutputFile,".\n", sep=""))
  }

  starttime <- Sys.time()

  if(!silent & !apollo_control$noValidation) cat("\nTesting probability function (apollo_probabilities)\n")
  b <- as.list(apollo_beta)
  if(apollo_control$mixing) b <- c(b, apollo_randcoeff(apollo_beta, draws))
  apollo_probabilities(b, database, functionality="validate")
  testLL = apollo_probabilities(b, database, functionality="estimate")
  isVec <- is.vector(testLL[["model"]])
  isMat <- is.matrix(testLL[["model"]])
  isCub <- (is.array(testLL[["model"]]) && length(dim(testLL[["model"]]))==3)
  ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
  if(isCub) testLL = apollo_avgIntraDraws(testLL, apollo_control, functionality="estimate")
  ### Product across choices for same individual (comment out this line if not using panel data)
  if(apollo_control$panelData) testLL = apollo_panelProd(testLL,apollo_control,functionality="estimate", work_in_logs, database[,apollo_control$indivID])
  ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
  if(isMat || isCub) testLL = apollo_avgInterDraws(testLL, apollo_control, functionality="estimate", database[,apollo_control$indivID])
  ### Prepares output for function (do not change this line)
  testLL = apollo_prepareProb(testLL,apollo_control,functionality="estimate")
  if(anyNA(testLL)) stop('Log-likelihood calculation fails at starting values!')

  # ################################## #
  #### HB version                   ####
  # ################################## #

  if(apollo_control$HB){

    apollo_HB_likelihood=function(fc,b){

      if(is.null(HB_control$gVarNamesFixed)) fc1=NULL else {
        fc1=as.list(fc)
        names(fc1)=HB_control$gVarNamesFixed
      }
      if(is.null(HB_control$gVarNamesNormal)) b1=NULL else {
        b1=as.data.frame(b)
        names(b1)=HB_control$gVarNamesNormal
      }
      if(length(apollo_beta_fixed)==0) fp=NULL else {
        fp=as.list(apollo_beta[apollo_beta_fixed])
        names(fp)=apollo_beta_fixed
      }
      theta=c(fc1,b1,fp)

      P = apollo_probabilities(theta, database, functionality="estimate")
      isVec <- is.vector(P[["model"]])
      isMat <- is.matrix(P[["model"]])
      isCub <- (is.array(P[["model"]]) && length(dim(P[["model"]]))==3)
      ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
      if(isCub) P = apollo_avgIntraDraws(P, apollo_control, functionality="estimate")
      ### Product across choices for same individual (comment out this line if not using panel data)
      if(apollo_control$panelData) P = apollo_panelProd(P,apollo_control,functionality="estimate", work_in_logs, database[,apollo_control$indivID])
      ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
      if(isMat || isCub) P = apollo_avgInterDraws(P, apollo_control, functionality="estimate", database[,apollo_control$indivID])
      ### Prepares output for function (do not change this line)
      P = apollo_prepareProb(P,apollo_control,functionality="estimate")
      return(P)
    }

    model   = RSGHB::doHB(apollo_HB_likelihood, database, HB_control)

    model$apollo_beta <- apollo_beta
    if(work_in_logs) model$LLStart <- sum(testLL) else model$LLStart <- sum(log(testLL))
    b <- as.list(apollo_beta)
    if(apollo_control$mixing) b <- c(b, apollo_randcoeff(apollo_beta, draws))
    P = apollo_probabilities(b, database, functionality="zero_LL")
    isVec <- is.vector(P[["model"]])
    isMat <- is.matrix(P[["model"]])
    isCub <- (is.array(P[["model"]]) && length(dim(P[["model"]]))==3)
    ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
    if(isCub) P = apollo_avgIntraDraws(P, apollo_control, functionality="zero_LL")
    ### Product across choices for same individual (comment out this line if not using panel data)
    if(apollo_control$panelData) P = apollo_panelProd(P,apollo_control,functionality="zero_LL", work_in_logs, database[,apollo_control$indivID])
    ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
    if(isMat || isCub) P = apollo_avgInterDraws(P, apollo_control, functionality="zero_LL", database[,apollo_control$indivID])
    ### Prepares output for function (do not change this line)
    P = apollo_prepareProb(P,apollo_control,functionality="zero_LL")
    model$LL0         <- sum(log(P))
    model$startTime   <- starttime
    model$apollo_control <- apollo_control
    model$nObs        <- nrow(database)
    model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
    endtime           <- Sys.time()
    model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
    model$apollo_beta_fixed <- apollo_beta_fixed
    model$estimation_routine <- "Hierarchical Bayes"
    model$HB_control  <- HB_control

    return(model)
  }

  # ################################## #
  #### classical version            ####
  # ################################## #

  if(!apollo_control$HB){
    if(apollo_control$nCores==1) cl <- NA else {
      cl <- apollo_makeCluster(apollo_control, apollo_probabilities, database,
                            draws, apollo_randcoeff, apollo_lcpars, silent=silent)
      apollo_control$nCores <- length(cl)
    }

    beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_beta_fixed)]
    beta_fix_val <- apollo_beta[apollo_beta_fixed]

    apollo_logLike <- apollo_makeLogLike(beta_fix_val, database, apollo_probabilities,
                                   apollo_control, draws, apollo_randcoeff, apollo_lcpars,
                                   cl=cl, estimation_routine, work_in_logs)

    if(!silent) cat("\n\nStarting main estimation\n")
    initial <- round(apollo_beta,4)
    if(silent) printLevel=0
    model <- maxLik::maxLik(apollo_logLike, start=beta_var_val,
                            method=estimation_routine, finalHessian=FALSE,
                            control=list(printLevel=printLevel, iterlim=max_iterations),
                            countIter=TRUE, writeIter=writeIter,
                            sumLL=FALSE)

    succesfulEstimation <- FALSE
    if(exists("model")){
      if(estimation_routine=="bfgs" & model$code==0) succesfulEstimation <- TRUE
      if(estimation_routine=="bhhh" & (model$code %in% c(2,8)) ) succesfulEstimation <- TRUE
      if(estimation_routine=="nr" && model$code<=2) succesfulEstimation <- TRUE
    }

    if(!succesfulEstimation){
      cat("ERROR: Estimation failed. No covariance matrix to compute.\n")
      if(exists("model")){
        print(as.matrix(model$estimate, ncol=1))
        return(model)
      } else stop("No estimated model to return.\n")
    }

    ### Print coeffs
    print(as.matrix(round(model$estimate, 4)))

    {
      success <- FALSE
      nNA <- -1

      # If Hessian is to be calculated with numDeriv
      if(numDerivHessian){
        if(!silent) cat("Computing covariance matrix using numDeriv package.\n (this may take a while)\n")

        # Create closure to keep track of progress and estimate hessian using numDeriv
        sumLogLike <- function(k, silent){
          i <- 0
          I <- 2+8*( k*(k+1)/2 )
          step <- ceiling(I/20)

          function(theta){
            if(i==0 & !silent) cat('0%')
            tmp <- apollo_logLike(theta, countIter=FALSE, writeIter=FALSE, sumLL=TRUE)
            i <<- i+1
            if(i%%step==0 & !silent){
              if(i%%(5*step)==0) cat(i/(5*step)*25,'%',sep='') else cat('.')
            }
            if(i==I & !silent) cat('100%\n')
            return(tmp)
          }
        }

        H <- tryCatch(numDeriv::hessian(func=sumLogLike(length(model$estimate), silent),
                                        x=model$estimate), error = function(e) return(NA))

        if(length(H)==1 && anyNA(H)){
          if(!silent) cat("\nERROR: Hessian estimation using numDeriv failed.\n")
          numDerivHessian <- FALSE
        } else {
          success <- TRUE
          nNA <- sum(anyNA(H))
          if(nNA>0 & !silent) cat("\nSome NA values found in numDeriv Hessian.\n")
        }
      }

      # If Hessian is to be calculated with maxLik
      if(!numDerivHessian | nNA>0){

        if(!silent) cat("Computing covariance matrix using maxLik package.\n (this may take a while, no progress bar displayed)\n")
        model2 <- tryCatch(maxLik::maxLik(apollo_logLike, start=model$estimate, method=estimation_routine, print.level=0,
                                          finalHessian=TRUE, iterlim=10, countIter=FALSE, writeIter=FALSE, sumLL=FALSE),
                           error=function(e) return(NA))

        if(length(model2)==1 && anyNA(model2)){
          if(!silent) cat("\nERROR: Hessian estimation using maxLik failed.\n")
        } else {
          if(!numDerivHessian){
            H <- model2$hessian
            success <- TRUE
            if(!silent) cat("\nHessian estimated with maxLik will be used.\n")
          } else {
            if(sum(anyNA(model2$hessian))<nNA){
              H <- model2$hessian
              if(!silent) cat("\nHessian estimated with maxLik will be used.\n")
            }
          }
        }
      }

      if(!success) H <- NULL
    }

    model$hessian <- H
    if(is.null(model$hessian)){
      if(!silent) cat("ERROR: Hessian could not be estimated. Postprocessing aborted.\n")
      return(model)
    }
    if(!is.matrix(try(solve(model$hessian),silent=T))){
      if(!silent) cat('ERROR: Singular Hessian, cannot calculate s.e. Postprocessing aborted.\n')
      utils::write.csv(model$hessian, paste(apollo_control$modelName, "hessian.csv", sep="_"))
      if(!silent) cat("Hessian written to", paste(apollo_control$modelName, "hessian.csv", sep="_"), "file\n.")
      return(model)
    }
    dummyVCM           <- matrix(NA, nrow=length(model$estimate), ncol=length(model$estimate))
    rownames(dummyVCM) <- names(model$estimate)
    colnames(dummyVCM) <- names(model$estimate)
    model$varcov    <- tryCatch(stats::vcov(model), error=function(e) return(dummyVCM))
    model$robvarcov <- tryCatch(sandwich::sandwich(model), error=function(e) return(dummyVCM))
    model$se        <- sqrt(diag(model$varcov))
    model$robse     <- sqrt(diag(model$robvarcov))

    model$corrmat    <- tryCatch(model$varcov/(model$se%*%t(model$se)), error=function(e) return(dummyVCM))
    model$robcorrmat <- tryCatch(model$robvarcov/(model$robse%*%t(model$robse)) , error=function(e) return(dummyVCM))

    P <- exp(apollo_logLike(model$estimate))
    if(apollo_control$panelData){
      nObsPerIndiv <- as.vector(table(database[,apollo_control$indivID]))
    } else nObsPerIndiv <- 1
    avgCP <- P^(1/nObsPerIndiv)
    names(avgCP) <- unique(database[,apollo_control$indivID])

    b <- as.list(model$estimate)
    if(apollo_control$mixing) b <- c(b, apollo_randcoeff(model$estimate, draws))
    P <- apollo_probabilities(b, database, functionality="zero_LL")
    isVec <- is.vector(P[["model"]])
    isMat <- is.matrix(P[["model"]])
    isCub <- (is.array(P[["model"]]) && length(dim(P[["model"]]))==3)
    ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
    if(isCub) P = apollo_avgIntraDraws(P, apollo_control, functionality="zero_LL")
    ### Product across choices for same individual (comment out this line if not using panel data)
    if(apollo_control$panelData) P = apollo_panelProd(P,apollo_control,functionality="zero_LL", work_in_logs, database[,apollo_control$indivID])
    ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
    if(isCub || isMat) P = apollo_avgInterDraws(P, apollo_control, functionality="zero_LL", database[,apollo_control$indivID])
    ### Prepares output for function (do not change this line)
    P = apollo_prepareProb(P,apollo_control,functionality="zero_LL")
    LL0 = sum(log(P))

    temp=c(model$estimate,apollo_beta[apollo_beta_fixed])
    model$estimate=temp[names(apollo_beta)]

    temp=c(model$se,apollo_beta[apollo_beta_fixed])
    model$se=temp[names(apollo_beta)]

    temp=c(model$robse,apollo_beta[apollo_beta_fixed])
    model$robse=temp[names(apollo_beta)]

    model$se[apollo_beta_fixed]         <- NA
    model$robse[apollo_beta_fixed]      <- NA

    b <- as.list(model$estimate)
    if(apollo_control$mixing) b <- c(b, apollo_randcoeff(model$estimate, draws))
    Pout <- apollo_probabilities(b, database, functionality="output")
    isVec <- is.vector(Pout[["model"]])
    isMat <- is.matrix(Pout[["model"]])
    isCub <- (is.array(Pout[["model"]]) && length(dim(Pout[["model"]]))==3)
    ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
    if(isCub) Pout = apollo_avgIntraDraws(Pout, apollo_control, functionality="output")
    ### Product across choices for same individual (comment out this line if not using panel data)
    if(apollo_control$panelData) Pout = apollo_panelProd(Pout,apollo_control,functionality="output", work_in_logs, database[,apollo_control$indivID])
    ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
    if(isMat || isCub) Pout = apollo_avgInterDraws(Pout, apollo_control, functionality="output", database[,apollo_control$indivID])
    ### Prepares output for function (do not change this line)
    Pout = apollo_prepareProb(Pout,apollo_control,functionality="output")
    if(is.list(Pout)){
      LLout = c(rep(0,length(Pout)))
      if(work_in_logs){
        LLout[1] <- sum(Pout[["model"]])
      } else {
        LLout[1] <- sum(log(Pout[["model"]]))
      }
      j=1
      k=2
      while(j<=length(Pout)){
        if(names(Pout)[j]!="model"){
          if(work_in_logs){
            LLout[k] <- sum(Pout[[j]])
          } else {
            LLout[k] <- sum(log(Pout[[j]]))
          }
          k <- k+1
        }
        j=j+1
      }
      names(LLout) <- c("model",names(Pout)[which(names(Pout)!="model")])
      for(i in 1:length(LLout)) if(names(LLout)[i]=="") names(LLout)[i] <- paste("component_",i-1,sep="")
    }

    if(exists('cl') & apollo_control$nCores>1) parallel::stopCluster(cl)

    model$apollo_beta <- beta_var_val
    if(work_in_logs) model$LLStart <- sum(testLL) else model$LLStart <- sum(log(testLL))
    model$LL0         <- LL0
    model$Pout        <- Pout
    model$LLout       <- LLout
    model$avgCP       <- avgCP
    model$startTime   <- starttime
    model$nIter       <- ifelse(estimation_routine=="bfgs", apollo_logLike(NA, getNIter=TRUE), model$iterations)
    model$apollo_control <- apollo_control
    model$nObs        <- nrow(database)
    model$nIndivs     <- length(unique(database[,apollo_control$indivID]))
    if(apollo_control$mixing) model$apollo_draws <- draws$apollo_draws else model$apollo_draws <- NA
    model$apollo_randcoeff<-apollo_randcoeff
    model$apollo_lcpars  <- apollo_lcpars
    endtime           <- Sys.time()
    model$timeTaken   <- as.numeric(difftime(endtime,starttime,units='secs'))
    model$apollo_beta_fixed <- apollo_beta_fixed
    model$estimation_routine <- estimation_routine

    return(model)
  }
}
