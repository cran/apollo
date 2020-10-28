#' Searches for better starting values.
#'
#' Given a set of starting values and a range for them, searches for points with a better likelihood.
#'
#' This function implements a simplified version of the algorithm proposed by Bierlaire, Themans, & Zufferey (2010). The main difference
#' lies in it implementing only two out of three tests on the candidates described by the authors. The implemented algorithm has the
#' following steps.
#' \enumerate{
#'   \item Randomly draw \code{nCandidates} candidates from an interval given by the user.
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
#' This function will write a CSV file to the working directory summarising progress. This file is called \code{modelName}_searchStart.csv .
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
#'                                   \itemize{
#'                                     \item \code{nCandidates}: Numeric scalar. Number of candidate sets of parameters to be used at the start. Should be an integer bigger than 1. Default is 100.
#'                                     \item \code{smartStart}: Boolean. If TRUE, candidates are randomly generated with more chances in the directions the Hessian indicates improvement of the LL function. Default is FALSE.
#'                                     \item \code{apolloBetaMin}: Vector. Minimum possible value of parameters when generating candidates. Ignored if smartStart is TRUE. Default is \code{apollo_beta - 0.1}.
#'                                     \item \code{apolloBetaMax}: Vector. Maximum possible value of parameters when generating candidates. Ignored if smartStart is TRUE. Default is \code{apollo_beta + 0.1}.
#'                                     \item \code{maxStages}: Numeric scalar. Maximum number of search stages. The algorithm will stop when there is only one candidate left, or if it reaches this number of stages. Default is 5.
#'                                     \item \code{dTest}: Numeric scalar. Tolerance for test 1. A candidate is discarded if its distance in parameter space to a better one is smaller than \code{dTest}. Default is 1.
#'                                     \item \code{gTest}: Numeric scalar. Tolerance for test 2. A candidate is discarded if the norm of its gradient is smaller than \code{gTest} AND its LL is further than \code{llTest} from a better candidate. Default is \code{10^(-3)}.
#'                                     \item \code{llTest}: Numeric scalar. Tolerance for test 2. A candidate is discarded if the norm of its gradient is smaller than \code{gTest} AND its LL is further than \code{llTest} from a better candidate. Default is 3.
#'                                     \item \code{bfgsIter}: Numeric scalar. Number od BFGS iterations to perform at each stage to each remaining candidate. Default is 20.
#'                                   }
#' @return named vector of model parameters. These are the best values found.
#' @export
#' @importFrom numDeriv hessian
#' @importFrom maxLik maxLik
apollo_searchStart <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, searchStart_settings=NA){
   
   # # # # # # # # # # #
   #### Initialise  ####
   # # # # # # # # # # #
   
   ### Load defaults
   default <- list(nCandidates=100, apolloBetaMin=apollo_beta - 0.1, apolloBetaMax=apollo_beta + 0.1,
                   smartStart=FALSE, maxStages=5, dTest=1, gTest=10^-3, llTest=3, bfgsIter=20)
   if(length(searchStart_settings)==1 && is.na(searchStart_settings)) searchStart_settings <- default
   tmp <- names(default)[!(names(default) %in% names(searchStart_settings))] # options missing in searchStart_settings
   for(i in tmp) searchStart_settings[[i]] <- default[[i]]
   
   ### Extract setings
   nCandidates   = searchStart_settings[["nCandidates"]]
   apolloBetaMin = searchStart_settings[["apolloBetaMin"]]
   apolloBetaMax = searchStart_settings[["apolloBetaMax"]]
   smartStart    = searchStart_settings[["smartStart"]]
   maxStages     = searchStart_settings[["maxStages"]]
   dTest         = searchStart_settings[["dTest"]]
   gTest         = searchStart_settings[["gTest"]]
   llTest        = searchStart_settings[["llTest"]]
   bfgsIter      = searchStart_settings[["bfgsIter"]]
   if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
   
   ### Checks
   if(nCandidates<2) stop("Argument 'nCandidates' should be at least 2.")
   if(maxStages<1) stop("Argument 'maxStages' should be at least 1.")
   if(anyNA(c(apolloBetaMin,apolloBetaMax)) & !smartStart) stop("Invalid 'apolloBetaMin' and/or 'apolloBetaMax' parameters.")
   apollo_print("Testing probability function (apollo_probabilities)")
   apollo_inputs$apollo_control$noDiagnostics <- TRUE
   apollo_probabilitiesVal <- apollo_insertComponentName(apollo_probabilities)
   apollo_probabilitiesVal(apollo_beta, apollo_inputs, functionality="validate")
   rm(apollo_probabilitiesVal)
   
   ### Pre-process likelihood function
   if(!silent) apollo_print("Pre-processing likelihood function...")
   # Create multi-core version of likelihood (if needed)
   apollo_logLike <- apollo_makeLogLike(apollo_beta, apollo_fixed, 
                                        apollo_probabilities, apollo_inputs, 
                                        list(estimationRoutine='BHHH'))
   on.exit(if(exists('apollo_logLike') && apollo_inputs$apollo_control$nCores>1) parallel::stopCluster(environment(apollo_logLike)$cl))
   # Create gradient function if required (and possible)
   if(!is.null(apollo_inputs$apollo_control$analyticGrad) && apollo_inputs$apollo_control$analyticGrad){
      # apollo_makeGrad will create gradient function ONLY IF all components have analytical gradient.
      grad <- apollo_makeGrad(apollo_beta, apollo_fixed, apollo_logLike, validateGrad=TRUE)
      if(is.null(grad)) apollo_inputs$apollo_control$analyticGrad <- FALSE
   } else grad <- NULL
   
   ### Separate beta into variable and fixed part
   beta_var_val <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
   beta_fix_val <- apollo_beta[apollo_fixed]
   K <- length(beta_var_val)
   
   # # # # # # # # # # # # # # # # #
   #### Generate candidate set  ####
   # # # # # # # # # # # # # # # # #
   
   ### Only non-fixed parameters are considered.
   # Candidates AT THE START OF EACH ITERATION are stored in a list of matrices called 'candidates'. 
   # Each element of the list is a matrix where each row is a candidate
   set.seed(2)
   candidates <- list()
   apollo_print(paste("Creating initial set of",nCandidates,"candidate values."))
   if(smartStart){
      # Create neighbours of starting value using Hessian (Bierlaire et al. 2007)
      apollo_print(" (using Hessian, this might take a while).")
      if(!is.null(grad)){
         sumGradLL <- function(theta) colSums( grad(theta) )
         H <- numDeriv::jacobian(func=sumGradLL, x=beta_var_val)
      } else H <- numDeriv::hessian(apollo_logLike, beta_var_val, sumLL=TRUE)
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
   
   ### Calculate LL of all candidates AT THE START OF THE ITERATION
   LL <- rep(0, nCandidates)
   cat("Calculating LL of candidates 0%")
   for(j in 1:nCandidates){
      beta_test        <- as.vector(candidates[[1]][j,])
      names(beta_test) <- names(beta_var_val)
      LL[j]            <- apollo_logLike(beta_test, sumLL=TRUE)
      if(j==ceiling(nCandidates/2)) cat("50%") else if(j%%ceiling(nCandidates/10)==0) cat(".")
   }
   cat("100%\n")
   
   ### Remove candidates with infinite or NA loglike
   if(any(!is.finite(LL))){
      removeRows <- which(!is.finite(LL))
      candidates[[1]] <- candidates[[1]][-removeRows,]
      LL <- LL[-removeRows]
      nCandidates <- nCandidates - length(removeRows)
      apollo_print(paste0(length(removeRows), " candidates removed due to non-finite starting LL. ", nCandidates, " remain."))
      if(nCandidates<1) stop("All initial candidates had non-finite starting LL.")
      if(nCandidates==1) return(candidates[[1]])
   }
   
   ### Write candidates
   fileName <- paste0(apollo_inputs$apollo_control$modelName,"_searchStart.csv")
   utils::write.csv(cbind(candidates[[1]], LL=LL, stage=1), fileName, row.names=FALSE)
   
   # # # # # # # # # # # # # # # #
   #### Loop over candidates  ####
   # # # # # # # # # # # # # # # #
   
   active    <- rep(TRUE, nCandidates)
   converged <- rep(FALSE, nCandidates)
   s <- 1
   while(sum(active)>1 & s<=maxStages){
      activeRows <- which(active)
      
      apollo_print('\n')
      apollo_print(paste0("Stage ", s, ", ", length(activeRows), " active candidates."))
      apollo_print(paste("Estimating", bfgsIter, "BFGS iteration(s) for each active candidate."))
      cat(" Candidate......LLstart.....LLfinish.....GradientNorm...Converged")
      
      # Apply BFGS for each active candidate
      LL               <- cbind(LL, NA)
      colnames(LL)     <- paste0("LL", 1:ncol(LL))
      candidates[[s+1]]<- matrix(NA, nrow=nCandidates, ncol=K)
      gradientNorm     <- rep(NA, nCandidates)
      for(j in activeRows){
         LLstart <- as.character(round(LL[j,s],0))
         cat("\n ",rep(" ",9-nchar(as.character(j))), j,
             " ", rep(" ",12-nchar(LLstart)),LLstart, sep="")
         
         if(!converged[j]){ # If it hasn't converged yet
            candidateParam <- as.vector(candidates[[s]][j,])
            names(candidateParam) <- names(beta_var_val)
            model <- maxLik::maxLik(apollo_logLike, start=candidateParam,
                                    method="bfgs", print.level=0, grad=grad, 
                                    finalHessian=FALSE, iterlim=bfgsIter)
            candidates[[s+1]][j,] <- model$estimate
            LL[j,s+1]             <- model$maximum
            gradientNorm[j]       <- sqrt(sum(model$gradient^2))
            converged[j]          <- ifelse(model$code==0, TRUE, FALSE)
         } else {# if it has converged already
            candidates[[s+1]][j,] <- candidates[[s]][j,]
            LL[j,s+1]             <- LL[j,s]
            gradientNorm[j]       <- 0 # approximation
            converged[j]          <- TRUE
         }
         
         LLfinish <- as.character(round(LL[j,s+1],0))
         gradFin  <- as.character(round(gradientNorm[j],3))
         cat(" ", rep(" ",12-nchar(LLfinish)), LLfinish,
             " ", rep(" ",16-nchar(gradFin)), gradFin,
             " ", rep(" ",10), 1*converged[j], sep="")
      }; cat('\n')
      rm(LLstart, candidateParam, model, LLfinish, gradFin)
      
      # Update active list
      for(j in activeRows){
         candParam    <- as.vector(candidates[[s+1]][j,])
         candLL       <- LL[j,s+1]
         betterLLRows <- activeRows[LL[activeRows,s+1]>=candLL & activeRows!=j]
         
         # Test 0: Converged to a worst solution
         if(converged[j] & length(betterLLRows)>0) active[j] <- FALSE
         
         # Apply tests only if it hasn't converged and is not the best LL of the stage
         if(length(betterLLRows)>0 & !converged[j]){
            # Test 1: Distance in parameter space to a better one
            betterParams <- candidates[[s+1]][betterLLRows,]
            if(length(betterLLRows)==1) betterParams <- matrix(betterParams, nrow=1)
            distParam <- apply(betterParams, MARGIN=1, function(x) sqrt(sum((x-candParam)^2)) )
            failedT1Rows <- betterLLRows[ distParam<dTest ]
            
            # Test 2: small gradient norm and close to another
            distLL <- abs(as.vector(LL[betterLLRows,s+1] - candLL))
            failedT2Rows <- betterLLRows[ gradientNorm[betterLLRows]<gTest & distLL>=llTest ]
            
            if(length(failedT1Rows)>0 | length(failedT2Rows)>0){
               active[j] <- FALSE
               apollo_print(paste0("Candidate ", j, " dropped."))
               if(length(failedT1Rows)>0) apollo_print(paste0("- Failed test 1: Too close to ", 
                                                              paste0(failedT1Rows, collapse=', '), 
                                                              " in parameter space."))
               if(length(failedT2Rows)>0) apollo_print(paste0("- Failed test 2: Converging to a worse solution than", 
                                                              paste0(failedT2Rows, collapse=', ')))
            }
            
            rm(betterParams, distParam, failedT1Rows, failedT2Rows)
         }
         rm(candParam, candLL, betterLLRows)
      }
      
      
      # Print best candidate to screen
      bestCandRow <- which.max(LL[,s+1])
      bestCandParam <- as.vector(candidates[[s+1]][bestCandRow,])
      names(bestCandParam) <- names(beta_var_val)
      bestCandParam <- c(bestCandParam, apollo_beta[apollo_fixed])
      bestCandParam <- bestCandParam[names(apollo_beta)]
      bestCandParam <- matrix(bestCandParam, ncol=1, dimnames = list(names(apollo_beta), 'Value'))
      apollo_print(paste0("Best candidate so far (LL=", round(LL[bestCandRow,s+1],1), ")"))
      print(as.matrix(round(bestCandParam,4)))
      
      # Write current state to file
      tryCatch(utils::write.csv(cbind(candidates[[1]], LL), fileName, row.names=FALSE), 
               error=function(e) apollo_print(paste0("Stage update ",s," could not be written to file", fileName)))
      
      # Next iteration
      s <- s+1
   }
   
   tmp <- as.vector(bestCandParam)
   names(tmp) <- rownames(bestCandParam)
   invisible(return(tmp))
}
