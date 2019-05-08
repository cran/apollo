#' Searches for better starting values.
#'
#' Given a set of starting values and a range for them, searches for points with a better likelihood.
#'
#' This function implements a simplified version of the algorithm proposed by Bierlaire, Themans, & Zufferey (2010). The main difference
#' lies in it implementing only two out of three tests on the candidates described by the authors. The implemented algorithm has the
#' following steps.
#' \enumerate{
#'   \item Randomly draw nCandidates candidates from an interval given by the user.
#'   \item Label all candidates with a valid log-likelihood (LL) as active.
#'   \item Apply \code{bfgsIter} iterations of the BFGS algorithm to each active candidate.
#'   \item Apply the following tests to each active candidate:
#'   \enumerate{
#'     \item Has the BGFS search converged?
#'     \item Are the candidate parameters after BFGS closer than \code{dTest} from any other candidate with higher LL?
#'     \item Is the LL of the candidate after BFGS further than \code{distLL} from a candidate with better LL, and its gradient smaller than \code{gTest}?
#'   }
#'   \item Mark any candidates for which at least one test results in yes as inactive.
#'   \item Go back to step 3, unless only one candidate is active, or the maximum number of iterations (\code{maxStages}) has been reached.
#' }
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param searchStart_settings List containing options for the search of starting values. The following are valid options.
#'                                   \describe{
#'                                     \item{nCandidates}{Numeric scalar. Number of candidate sets of parameters to be used at the start. Should be an integer bigger than 1. Default is 100.}
#'                                     \item{smartStart}{Boolean. If TRUE, candidates are randomly generated with more chances in the directions the Hessian indicates improvement of the LL function. Default is FALSE.}
#'                                     \item{apolloBetaMin}{Vector. Minimum possible value of parameters when generating candidates. Ignored if smartStart is TRUE. Default is \code{apollo_beta - 0.1}.}
#'                                     \item{apolloBetaMax}{Vector. Maximum possible value of parameters when generating candidates. Ignored if smartStart is TRUE. Default is \code{apollo_beta + 0.1}.}
#'                                     \item{maxStages}{Numeric scalar. Maximum number of search stages. The algorithm will stop when there is only one candidate left, or if it reaches this number of stages. Default is 5.}
#'                                     \item{dTest}{Numeric scalar. Tolerance for test 1. A candidate is discarded if its distance in parameter space to a better one is smaller than \code{dTest}. Default is 1.}
#'                                     \item{gTest}{Numeric scalar. Tolerance for test 2. A candidate is discarded if the norm of its gradient is smaller than \code{gTest} AND its LL is further than \code{llTest} from a better candidate. Default is 10^(-3).}
#'                                     \item{llTest}{Numeric scalar. Tolerance for test 2. A candidate is discarded if the norm of its gradient is smaller than \code{gTest} AND its LL is further than \code{llTest} from a better candidate. Default is 3.}
#'                                     \item{bfgsIter}{Numeric scalar. Number od BFGS iterations to perform at each stage to each remaining candidate. Default is 20.}
#'                                   }
#' @return named vector of model parameters. These are the best values found.
#' @export
#' @importFrom numDeriv hessian
#' @importFrom maxLik maxLik
apollo_searchStart <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, searchStart_settings=NA){

 database        = apollo_inputs[["database"]]
 apollo_control  = apollo_inputs[["apollo_control"]]
 draws           = apollo_inputs[["draws"]]
 apollo_randCoeff= apollo_inputs[["apollo_randCoeff"]]
 apollo_lcPars   = apollo_inputs[["apollo_lcPars"]]
 workInLogs      = apollo_control$workInLogs

 default <- list(nCandidates=100, apolloBetaMin=apollo_beta - 0.1, apolloBetaMax=apollo_beta + 0.1,
                 smartStart=FALSE, maxStages=5, dTest=1, gTest=10^-3, llTest=3, bfgsIter=20)
 if(length(searchStart_settings)==1 && is.na(searchStart_settings)) searchStart_settings <- default
 tmp <- names(default)[!(names(default) %in% names(searchStart_settings))] # options missing in searchStart_settings
 for(i in tmp) searchStart_settings[[i]] <- default[[i]]

 nCandidates     = searchStart_settings[["nCandidates"]]
 apolloBetaMin = searchStart_settings[["apolloBetaMin"]]
 apolloBetaMax = searchStart_settings[["apolloBetaMax"]]
 smartStart      = searchStart_settings[["smartStart"]]
 maxStages       = searchStart_settings[["maxStages"]]
 dTest           = searchStart_settings[["dTest"]]
 gTest           = searchStart_settings[["gTest"]]
 llTest           = searchStart_settings[["llTest"]]
 bfgsIter        = searchStart_settings[["bfgsIter"]]

  ### checks
  if(nCandidates<2) stop("Argument 'nCandidates' should be at least 2.")
  if(maxStages<1) stop("Argument 'maxStages' should be at least 1.")
  if(anyNA(c(apolloBetaMin,apolloBetaMax)) & !smartStart) stop("Invalid 'apolloBetaMin' and/or 'apolloBetaMax' parameters.")
  cat("Testing probability function (apollo_probabilities)")
  apollo_control$noDiagnostics <- TRUE
  apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate")

  ### Separate beta into variable and fixed part
  beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  beta_fix_val <- apollo_beta[apollo_fixed]
  K <- length(beta_var_val)

  ### Create logLike function
  if(apollo_control$nCores==1) cl <- NA else {
    cat("Initializing cluster.\n")
    cl <- apollo_makeCluster(apollo_probabilities, apollo_inputs, silent=TRUE)
    apollo_control$nCores <- length(cl)
  }
  apollo_estSet  <- list(estimationRoutine="bfgs", maxIterations=2, writeIter=FALSE,
                         hessianRoutine="none", printLevel=0, silent=FALSE)
  apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, apollo_estSet, cl=cl)

  ### Create candidates.
  ### Only non-fixed parameters are considered.
  ### Candidates are stored in a list of matrices called 'candidates'
  ### Each element of the list is a matrix where each row is a candidate
  set.seed(2)
  candidates <- list()
  cat("Creating initial set of",nCandidates,"candidate values.\n")
  if(smartStart){
    # Create neighbors of starting value using Hessian (Bierlaire et al. 2007)
    cat(" (using Hessian, this might take a while).\n")
    H <- numDeriv::hessian(apollo_logLike, beta_var_val, sumLL=TRUE)
    eig <- eigen(H)
    w <- apply(eig$vectors, MARGIN=2, function(wi) wi/sqrt(sum(wi^2)))
    w <- cbind(w,-w)
    p <- rep(exp(0.5*eig$values),2)
    p <- p/sum(p)
    w <- w[, sample(ncol(w), nCandidates-1, replace=TRUE, prob=p)]
    w <- split(w, rep(1:ncol(w), each = nrow(w)))
    a <- as.list(stats::runif(length(w), min=0.75, max=1))
    z <- beta_var_val + mapply(function(wi,ai) wi*ai , w, a)
    z <- cbind(beta_var_val,z) # matrix where each column is a candidate
    candidates[[1]] <- t(z)
    rm(H, eig, w, p, a, z)
  } else{
    # Create random starting values using apolloBetaMin & apolloBetaMax.
    apolloBetaMin <- apolloBetaMin[names(beta_var_val)]
    apolloBetaMax <- apolloBetaMax[names(beta_var_val)]
    candidates[[1]] <- t(apollo_mlhs(nCandidates-1,K,1)) #matrix(stats::runif((nCandidates-1)*K), nrow=K, ncol=nCandidates-1)
    candidates[[1]] <- apolloBetaMin + (apolloBetaMax-apolloBetaMin)*candidates[[1]]
    candidates[[1]] <- cbind(beta_var_val, candidates[[1]])
    candidates[[1]] <- t(candidates[[1]])
  }

  ### Calculate LL of all candidates
  LL <- rep(0, nCandidates)
  cat("Calculating LL of candidates 0%")
  for(j in 1:nCandidates){
    beta_test        <- as.vector(candidates[[1]][j,])
    names(beta_test) <- names(beta_var_val)
    LL[j]            <- apollo_logLike(beta_test, sumLL=TRUE)
    if(j==ceiling(nCandidates/2)) cat("50%") else if(j%%ceiling(nCandidates/10)==0) cat(".")
  }
  cat("100%")

  ### Remove candidates with infinite or NA loglike
  if(any(!is.finite(LL))){
    removeRows <- which(!is.finite(LL))
    candidates[[1]] <- candidates[[1]][-removeRows,]
    LL <- LL[-removeRows]
    nCandidates <- nCandidates - length(removeRows)
    cat("\n", length(removeRows), " candidates removed due to non-finite starting LL. ", nCandidates, " remain.", sep="")
    if(nCandidates<1) stop("All initial candidates had non-finite starting LL.")
    if(nCandidates==1) return(candidates[[1]])
  }

  ### Write candidates
  fileName <- paste0(apollo_control$modelName,"_searchStart.csv")
  utils::write.csv(cbind(candidates[[1]], LL=LL, stage=1), fileName, row.names=FALSE)

  ### MAIN LOOP
  active    <- rep(TRUE, nCandidates)
  converged <- rep(FALSE, nCandidates)
  i <- 2
  while(sum(active)>1 & i<maxStages+1){
    activeRows <- which(active)

    cat("\n\nStage ", i-1, ", ", length(activeRows), " active candidates.", sep="")
    cat("\n Estimating", bfgsIter, "BFGS iteration(s) for each active candidate.")
    cat("\n Candidate......LLstart.....LLfinish.....GradientNorm...Converged")

    # Apply BFGS for each active candidate
    LL              <- cbind(LL, NA)
    candidates[[i]] <- matrix(NA, nrow=nCandidates, ncol=K)
    gradientNorm    <- rep(NA, nCandidates)
    for(j in activeRows){
      LLstart <- as.character(round(LL[j,i-1],0))
      cat("\n ",rep(" ",9-nchar(as.character(j))), j,
          " ", rep(" ",12-nchar(LLstart)),LLstart, sep="")

      candidateParam <- as.vector(candidates[[i-1]][j,])
      names(candidateParam) <- names(beta_var_val)
      model <- maxLik::maxLik(apollo_logLike, start=candidateParam,
                              method="bfgs", print.level=0,
                              finalHessian=FALSE, iterlim=bfgsIter)
      candidates[[i]][j,] <- model$estimate
      LL[j,i]             <- model$maximum
      gradientNorm[j]     <- sqrt(sum(model$gradient^2))
      converged[j]        <- ifelse(model$code==0, TRUE, FALSE)

      LLfinish <- as.character(round(LL[j,i],0))
      gradFin  <- as.character(round(gradientNorm[j],3))
      cat(" ", rep(" ",12-nchar(LLfinish)), LLfinish,
          " ", rep(" ",16-nchar(gradFin)), gradFin,
          " ", rep(" ",10), 1*converged[j], sep="")
    }
    rm(LLstart, candidateParam, model, LLfinish, gradFin)

    # Update active list
    for(j in activeRows){
      candParam    <- as.vector(candidates[[i]][j,])
      candLL       <- LL[j,i]
      betterLLRows <- activeRows[LL[activeRows,i]>=candLL & activeRows!=j]

      # Test 0: It hasn't converged
      if(converged[j]) active[j] <- FALSE

      # Apply tests only if it is not the best LL of the stage
      if(length(betterLLRows)>0 & !converged[j]){
        # Test 1: Distance in parameter space to a better one
        betterParams <- candidates[[i]][betterLLRows,]
        if(length(betterLLRows)==1) betterParams <- matrix(betterParams, nrow=1)
        distParam <- apply(betterParams, MARGIN=1, function(x) sqrt(sum((x-candParam)^2)) )
        failedT1Rows <- betterLLRows[ distParam<dTest ]

        # Test 2: small gradient norm and close to another
        distLL <- abs(as.vector(LL[betterLLRows,i] - candLL))
        failedT2Rows <- betterLLRows[ gradientNorm[betterLLRows]<gTest & distLL>=llTest ]

        if(length(failedT1Rows)>0 | length(failedT2Rows)>0){
          active[j] <- FALSE
          cat("\n Candidate", j, "dropped.")
          if(length(failedT1Rows)>0) cat("\n  Failed test 1 against", failedT1Rows)
          if(length(failedT2Rows)>0) cat("\n  Failed test 2 against", failedT2Rows)
        }
      }
    }
    rm(candParam, candLL, betterLLRows, betterParams, distParam, failedT1Rows, failedT2Rows)

    # Print best candidate to screen
    bestCandRow <- which.max(LL[,i])
    bestCandParam <- as.vector(candidates[[i]][j,])
    names(bestCandParam) <- names(beta_var_val)
    bestCandParam <- c(bestCandParam, apollo_beta[apollo_fixed])
    bestCandParam <- bestCandParam[names(apollo_beta)]
    cat("\n Best candidate so far (LL=", round(LL[bestCandRow,i],1), ")\n", sep="")
    print(as.matrix(round(bestCandParam,4)))

    # Update candidates file
    if(!file.exists(fileName)){
      tmp <- cbind(candidates[[1]], logLike=LL[,1], stage=1)
      tmp[activeRows, ] <- cbind(candidates[[i]][activeRows,], logLike=LL[activeRows,i], stage=i)
      tryCatch(utils::write.csv(tmp, fileName, row.names=FALSE), error=cat("\n Stage update could not be written to file", fileName))
    } else{
      carryOn <- FALSE
      tryCatch({tmp <- utils::read.csv(fileName); carryOn <- TRUE}, error=function(e) cat("\n File", fileName, "could not be accessed."))
      if(carryOn){
        tmp[activeRows,] <- cbind(candidates[[i]][activeRows,], logLike=LL[activeRows,i], stage=i)
        tryCatch(utils::write.csv(tmp, fileName, row.names=FALSE), error=function(e) cat("\n Stage update could not be written to file."))
      }
    }


    # Next iteration
    i <- i+1
  }

  invisible(return(bestCandParam))
}
