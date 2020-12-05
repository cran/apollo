#' Prints estimation results to console
#' 
#' Prints estimation results to console. Amount of information presented can be adjusted through arguments.
#' 
#' Prints to screen the output of a model previously estimated by apollo_estimate()
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param modelOutput_settings List of options. It can include the following.
#'                             \itemize{
#'                               \item \code{printFixed}: Logical. TRUE for printing fixed parameters among estimated parameter. TRUE by default.
#'                               \item \code{printClassical}: Logical. TRUE for printing classical standard errors. TRUE by default.
#'                               \item \code{printPVal}: Logical or Scalar. TRUE or 1 for printing p-values for one-sided test, 2 for printing p-values for two-sided test, FALSE for not printing p-values. FALSE by default.
#'                               \item \code{printT1}: Logical. If TRUE, t-test for H0: apollo_beta=1 are printed. FALSE by default.
#'                               \item \code{printDataReport}: Logical. TRUE for printing summary of choices in database and other diagnostics. FALSE by default.
#'                               \item \code{printModelStructure}: Logical. TRUE for printing model structure. TRUE by default.
#'                               \item \code{printCovar}: Logical. TRUE for printing parameters covariance matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. FALSE by default.
#'                               \item \code{printCorr}: Logical. TRUE for printing parameters correlation matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. FALSE by default.
#'                               \item \code{printOutliers}: Logical or Scalar. TRUE for printing 20 individuals with worst average fit across observations. FALSE by default. If Scalar is given, this replaces the default of 20.
#'                               \item \code{printChange}: Logical. TRUE for printing difference between starting values and estimates. FALSE by default.
#'                               \item \code{printFunctions}: Logical. TRUE for printing apollo_control, apollo_randCoeff (when available), apollo_lcPars (when available) and apollo_probabilities. FALSE by default.
#'                             }
#' @return A matrix of coefficients, s.d. and t-tests (invisible)
#' @export
#' @importFrom coda geweke.diag
#' @importFrom utils capture.output
apollo_modelOutput=function(model, modelOutput_settings=NA){
  if(length(modelOutput_settings)==1 && is.na(modelOutput_settings)) modelOutput_settings=list()
  if(is.null(modelOutput_settings[["printModelStructure"]])) modelOutput_settings[["printModelStructure"]] = TRUE
  if(is.null(modelOutput_settings[["printDataReport"]])) modelOutput_settings[["printDataReport"]] = FALSE
  if(!is.null(modelOutput_settings[["printDiagnostics"]])){
    modelOutput_settings[["printModelStructure"]] = modelOutput_settings[["printDiagnostics"]]
    modelOutput_settings[["printDataReport"]]     = modelOutput_settings[["printDiagnostics"]]
  } 
  if(is.null(modelOutput_settings[["printClassical"]])) modelOutput_settings[["printClassical"]] = TRUE
  if(is.null(modelOutput_settings[["printPVal"]])) modelOutput_settings[["printPVal"]] = FALSE
  if(is.null(modelOutput_settings[["printT1"]])) modelOutput_settings[["printT1"]] = FALSE
  if(is.null(modelOutput_settings[["printCovar"]])) modelOutput_settings[["printCovar"]] = FALSE
  if(is.null(modelOutput_settings[["printCorr"]])) modelOutput_settings[["printCorr"]] = FALSE
  if(is.null(modelOutput_settings[["printOutliers"]])) modelOutput_settings[["printOutliers"]] = FALSE
  if(is.null(modelOutput_settings[["printChange"]])) modelOutput_settings[["printChange"]] = FALSE
  if(is.null(modelOutput_settings[["printFunctions"]])) modelOutput_settings[["printFunctions"]] = FALSE
  if(is.null(modelOutput_settings[["printFixed"]])) modelOutput_settings[["printFixed"]] = TRUE
  
  printClassical   = modelOutput_settings[["printClassical"]]
  printPVal        = modelOutput_settings[["printPVal"]]
  printT1          = modelOutput_settings[["printT1"]]
  printDiagnostics = modelOutput_settings[["printDiagnostics"]]
  printCovar       = modelOutput_settings[["printCovar"]]
  printCorr        = modelOutput_settings[["printCorr"]]
  printOutliers    = modelOutput_settings[["printOutliers"]]
  printChange      = modelOutput_settings[["printChange"]]
  printFunctions   = modelOutput_settings[["printFunctions"]]
  printFixed       = modelOutput_settings[["printFixed"]]
  
  if(length(model$scaling)>0 && !anyNA(model$scaling) && !all(model$scaling==1)){
    scaling_used=TRUE
  }else{
    scaling_used=FALSE
  }
  # ####################### #
  #### MODEL DESCRIPTION ####
  # ####################### #
  
  apollo_control <- model$apollo_control
  nParams     <- length(model$apollo_beta)
  nFreeParams <- nParams
  if(!is.null(model$apollo_fixed)) nFreeParams <- nFreeParams - length(model$apollo_fixed)
  
  
  # Printing model information
  if("package:apollo" %in% search()){
    apolloVersion <- tryCatch(utils::packageDescription("apollo", fields = "Version"),
                              warning=function(w) return("alpha"),
                              error=function(e) return("alpha"))
  } else apolloVersion <- "alpha"
  
  cat('Model run using Apollo for R, version', apolloVersion, 'on', Sys.info()['sysname'], 'by', Sys.info()['user'], '\n')
  cat("www.ApolloChoiceModelling.com\n\n")
  cat("Model name                       : ", model$apollo_control$modelName,"\n", sep="")
  cat("Model description                : ", model$apollo_control$modelDescr,"\n", sep="")
  cat("Model run at                     : ", paste(model$startTime),"\n", sep="")
  cat("Estimation method                : ", model$estimationRoutine, "\n", sep="")
  if(!apollo_control$HB) cat("Model diagnosis                  : ",model$message,"\n", sep="")
  cat("Number of individuals            : ", model$nIndivs,"\n", sep="")
  cat("Number of rows in database       : ", model$nObs,"\n", sep="")
  if(!is.null(model$nObsTot)){
    cat('Number of modelled outcomes      : ', sum(model$nObsTot), '\n', sep='')
    if(length(model$nObsTot)>1) for(i in names(model$nObsTot)){
      if(nchar(i)>26) txt <- substr(i, 1, 26) else txt <- i # max 26 characters for the component name
      if(nchar(txt)<26) txt <- paste0(c(rep(' ', 26 - nchar(txt)), txt), collapse='')
      cat('       ', txt, ': ', model$nObsTot[i], '\n', sep='')
    }; cat('\n')
  }
  cat("Number of cores used             : ",model$apollo_control$nCores,"\n")
  if(model$apollo_control$mixing){
    d <- model$apollo_draws
    if(d$interNDraws>0 && length(c(d$interUnifDraws, d$interNormDraws))>0){
      cat("Number of inter-individual draws : ", d$interNDraws, ' (', d$interDrawsType, ')', "\n", sep='')
    }
    if(d$intraNDraws>0 && length(c(d$intraUnifDraws, d$intraNormDraws))>0){
      cat("Number of intra-individual draws     : ", d$intraNDraws, ' (', d$intraDrawsType, ')', "\n", sep='')
    }
    if(!model$apollo_control$panelData & model$apollo_control$mixing & d$interNDraws>0){
      cat("WARNING: Inter-individual draws were used\n")
      cat("         without a panel data structure.\n")
    }
    rm(d)
  } else { if(!apollo_control$HB) cat("Model without mixing\n") }
  cat("\n")
  
  
  # ####################### #
  #### HB OUTPUT         ####
  # ####################### #
  if(apollo_control$HB){
    apollo_HB <- model$apollo_HB
    cat("Estimation carried out using RSGHB\n")
    cat("Burn-in iterations               : ",model$gNCREP,"\n", sep="")
    cat("Post burn-in iterations          : ",model$gNEREP,"\n", sep="")
    cat("LL(start)                        : ",model$LLStart,"\n", sep="")
    cat("LL(0)                            : ",model$LL0[1],"\n", sep="")
    cat("Average post. LL post burn-in    : ",mean(colSums(log(model$cmcLLout))),"\n",sep="")
    cat("Average post. RLH post burn-in   : ",round(mean(colMeans((model$cmcRLHout))),4),"\n",sep="")
    f <- function(t){
      tmpH <- floor(t/60^2)
      tmpM <- floor((t-tmpH*60^2)/60)
      tmpS <- round(t-tmpH*60^2-tmpM*60,2)
      paste(formatC(tmpH,width=2,format='d',flag=0),
            formatC(tmpM,width=2,format='d',flag=0),
            tmpS,sep=':')
    }
    cat("Time taken (hh:mm:ss)            : ",f(model$timeTaken),"\n")
    cat("     pre-estimation              : ",f(model$timePre),"\n")
    cat("     estimation                  : ",f(model$timeEst),"\n")
    cat("     post-estimation             : ",f(model$timePost),"\n")
    cat("\n\n")
    
    cat("Chain convergence report (Geweke test)\n\n")
    if(!is.null(model$F)){
      cat("Fixed (non random) parameters (t-test value for Geweke test)\n")
      #tmp <- coda::geweke.diag(model$F[,2:(ncol(model$F))], frac1=0.1, frac2=0.5)[[1]]
      #names(tmp) <- model$params.fixed
      #print( round(tmp, 4) )
      print( round(model$F_convergence, 4) )
      cat("\n")
    }
    if(!is.null(model$A)){
      cat("Random parameters (t-test value for Geweke test)\n")
      #tmp <- coda::geweke.diag(model$A[,2:(ncol(model$A))], frac1=0.1, frac2=0.5)[[1]]
      #print( round(tmp, 4) )
      print( round(model$A_convergence, 4) )
      cat("\n")
    }
    if(!is.null(model$D)){
      cat("Covariances of random parameters (t-test value for Geweke test)\n")
      # This assumes the matrix is square
      #tmp <- c()
      #for(i in 1:dim(model$D)[1]) for(j in 1:i){
      #  if(i==1 & j==1) Dmatrix <- as.matrix(model$D[i,j,]) else Dmatrix <- cbind(Dmatrix, as.vector(model$D[i,j,]))
      #  tmp <- c(tmp, paste(colnames(model$A)[i+1],colnames(model$A)[j+1], sep="_"))
      #}
      #colnames(Dmatrix) <- tmp
      #tmp <- coda::geweke.diag(Dmatrix, frac1=0.1, frac2=0.5)[[1]]
      #print( round(tmp, 4) )
      print( round(model$D_convergence, 4) )
    }
    cat("\n\n")
    
    cat("Summary of parameter chains\n\n")
    ans <- list()
    
    if(length(apollo_HB$gVarNamesFixed)>0 | length(model$apollo_fixed)>0){
      cat("Non-random coefficients","\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n\n")
      print(round(model$chain_non_random,4))
      cat("\n")
      ans[["non_random"]] <- model$chain_non_random
    }
    
    apollo_HB$gVarNamesFixed <- model$params.fixed
    apollo_HB$gVarNamesNormal <- model$params.vary
    if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
      cat("Upper level model results for mean parameters for underlying Normals","\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      print(round(model$chain_random_mean,4))
      cat("\n")
      
      cat("Upper level model results for covariance matrix for underlying Normals (means across iterations)","\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      print(round(model$chain_random_cov_mean,4))
      cat("\n")
      
      cat("Upper level model results for covariance matrix for underlying Normals (SD across iterations)","\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      print(round(model$chain_random_cov_sd,4))
      cat("\n")
      
      cat("Summary of distributions of random coeffients (after distributional transforms)","\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      print(round(model$random_coeff_summary,4))
      cat("\n")
      
      if(length(apollo_HB$gVarNamesNormal)>1){
        cat("Covariance matrix of random coeffients (after distributional transforms)","\n")
        if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
        print(round(model$random_coeff_covar,4))
        cat("\n")
        
        cat("Correlation matrix of random coeffients (after distributional transforms)","\n")
        if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
        print(round(model$random_coeff_corr,4))
        cat("\n")
      }
      
      
      cat("Results for posterior means for random coefficients","\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      print(round(model$posterior_mean,4))
      cat("\n")
      
      ans[["random_mean"]]   <- model$chain_random_mean
      ans[["random_cov_mean"]] <- model$chain_random_cov_mean
      ans[["random_cov_sd"]] <- model$chain_random_cov_sd
      ans[["random_coeff_summary"]]     <- model$random_coeff_summary
      ans[["posterior"]]     <- model$posterior_mean
      ans[["random_coeff_covar"]]     <- model$random_coeff_covar
      ans[["random_coeff_corr"]]     <- model$random_coeff_corr
    }    
    return(invisible(ans))
  }
  
  
  # ####################### #
  #### CLASSICAL OUTPUT  ####
  # ####################### #
  
  ### change 7 August (next line, and then replacing 2 by multiplier in several lines belwo)
  if(printPVal==2) pMult <- 2 else pMult <- 1
  output=cbind(model$estimate,
               model$se,
               model$estimate/model$se,
               pMult*(1-stats::pnorm(abs(model$estimate/model$se))),
               (model$estimate-1)/model$se,
               pMult*(1-stats::pnorm(abs((model$estimate-1)/model$se))),
               model$robse,
               model$estimate/model$robse,
               pMult*(1-stats::pnorm(abs(model$estimate/model$robse))),
               (model$estimate-1)/model$robse,
               pMult*(1-stats::pnorm(abs((model$estimate-1)/model$robse))))
  #output <- signif(output,4)
  ### change 7
  if(pMult==2){
    colnames(output) <- c('Estimate',
                          's.e.', 't.rat.(0)', 'p(2-sided)','t.rat(1)','p(2-sided)', 
                          'Rob.s.e.','Rob.t.rat.(0)','p(2-sided)','Rob.t.rat.(1)','p(2-sided)')
  } else {
    colnames(output) <- c('Estimate',
                          's.e.', 't.rat.(0)', 'p(1-sided)','t.rat(1)','p(1-sided)', 
                          'Rob.s.e.','Rob.t.rat.(0)','p(1-sided)','Rob.t.rat.(1)','p(1-sided)')
  }
  rownames(output) <- names(model$estimate)
  # If there is a bootstrap covariance matrix
  if(!is.null(model$bootvarcov)){
    tmp <- model$bootse
    tmp <- cbind(`Bootstrap.s.e.`   = tmp,
                 `Bootstrap.t.rat.(0)` = model$estimate/tmp,
                 `p(2-sided)`       = 2*(1-stats::pnorm(abs(model$estimate/tmp))),
                 `p(1-sided)`       = 1*(1-stats::pnorm(abs(model$estimate/tmp))),
                 `Bootstrap.t.rat.(1)` = (model$estimate-1)/tmp,
                 `p(2-sided)`       = 2*(1-stats::pnorm(abs((model$estimate-1)/tmp))),
                 `p(1-sided)`       = 1*(1-stats::pnorm(abs((model$estimate-1)/tmp))) )
    if(!printPVal) tmp <- tmp[,-grep('p(', colnames(tmp), fixed=TRUE)]
    if(printPVal & pMult==1) tmp <- tmp[,-grep('p(2', colnames(tmp), fixed=TRUE)]
    if(printPVal & pMult==2) tmp <- tmp[,-grep('p(1', colnames(tmp), fixed=TRUE)]
    if(!printT1) tmp <- tmp[,-grep('t.rat.(1', colnames(tmp), fixed=TRUE)]
    output <- cbind(output, tmp)
  }
  
  dropcolumns = NULL
  if(printClassical==FALSE) dropcolumns = c(dropcolumns,2,3,4,5,6)
  if(printT1==FALSE) dropcolumns = c(dropcolumns,5,6,10,11,15,16)
  if(printPVal==FALSE) dropcolumns = c(dropcolumns,4,6,9,11,14,16)
  dropcolumns = unique(dropcolumns)
  if(length(dropcolumns)>0) output = output[,-dropcolumns, drop=FALSE]
  
  cat("LL(start)                        : ",model$LLStart,"\n", sep="")
  if(length(model$LLout)==1) cat("LL(0)                            : ",ifelse(!anyNA(model$LL0[1]),model$LL0[1],"Not applicable"),"\n", sep="")
  if(length(model$LLout)>1)  cat("LL(0, whole model)               : ",ifelse(!anyNA(model$LL0[1]),model$LL0[1],"Not applicable"),"\n", sep="")
  if(length(model$LLout)==1) cat("LL(final)                        : ",model$maximum,"\n",sep="")
  if(length(model$LLout)>1)  cat("LL(final, whole model)           : ",model$maximum,"\n",sep="")
  test <- !is.null(model$modelTypeList) && all(tolower(model$modelTypeList) %in% c("mnl", "nl", "cnl", "el", "dft", "lc"))
  test <- test && !anyNA(model$LL0[1])
  if(test){
    cat("Rho-square (0)                   : ",round(1-(model$maximum/model$LL0[1]),4),"\n")
    cat("Adj.Rho-square (0)               : ",round(1-((model$maximum-nFreeParams)/model$LL0[1]),4),"\n")
  } else {
    cat("Rho-square (0)                   : Not applicable\n")
    cat("Adj.Rho-square (0)               : Not applicable\n")
  }
  cat("AIC                              : ",round(-2*model$maximum + 2*nFreeParams,2),"\n")
  if(!is.null(model$nObsTot)) tmp <- sum(model$nObsTot) else tmp <- model$nObs
  cat("BIC                              : ",round(-2*model$maximum + nFreeParams*log(tmp),2),"\n")
  
  cat("\n")
  if(length(model$LLout)>1){
    for(j in 2:length(model$LLout)){
      nam <- names(model$LLout)[j]
      sp1 <- ifelse(6+nchar(nam)<33, paste0(rep(" ",33-nchar(nam)-6), collapse=""), "")
      sp2 <- ifelse(9+nchar(nam)<33, paste0(rep(" ",33-nchar(nam)-10), collapse=""), "")
      cat("LL(0,",nam,")",sp1,": ",ifelse(is.finite(model$LL0[j]),model$LL0[j],"Not applicable"),"\n", sep="")
      cat("LL(final,",nam,")",sp2,": ",model$LLout[j],"\n", sep="")
    }; rm(nam, sp1, sp2)
  }
  
  cat("\n")
  cat("Estimated parameters             :  ", nFreeParams,"\n", sep="")
  #cat("Norm of the gradient at optimum  : ",round( sqrt(sum(model$gradient^2)),2), "\n\n")
  f <- function(t){
    tmpH <- floor(t/60^2)
    tmpM <- floor((t-tmpH*60^2)/60)
    tmpS <- round(t-tmpH*60^2-tmpM*60,2)
    paste(formatC(tmpH,width=2,format='d',flag=0),
          formatC(tmpM,width=2,format='d',flag=0),
          tmpS,sep=':')
  }
  cat("Time taken (hh:mm:ss)            : ",f(model$timeTaken),"\n")
  cat("     pre-estimation              : ",f(model$timePre),"\n")
  cat("     estimation                  : ",f(model$timeEst),"\n")
  cat("     post-estimation             : ",f(model$timePost),"\n")
  if(length(grep('successful', model$message))==0) tmp <- paste0("(",model$message,")") else tmp <- ""
  cat("Iterations                       : ",model$nIter, tmp, "\n")
  if(model$bootstrapSE>0) cat("Number of bootstrap repetitions  : ", model$bootstrapSE, "\n")
  if(!is.null(model$eigValue) && !anyNA(model$eigValue)){
    cat("Min abs eigenvalue of Hessian    : ",round(min(abs(model$eigValue)),6),"\n")
    if(any(model$eigValue>0)) apollo_print("Some eigenvalues of Hessian are positive, indicating potential problems!") 
  }
  cat("\n")
  
  if(!printClassical & anyNA(model$se[!(names(model$estimate) %in% model$apollo_fixed)]) ){
    apollo_print('WARNING: Classical standard errors could not be calculated for some parameters. This could point to an identification or estimation problem.')
  }
  
  if(scaling_used) apollo_print("These outputs have had the scaling used in estimation applied to them.")
  cat("Estimates:\n")
  if(nrow(output)>options("max.print")) options(max.print=nrow(output)+100)
  if(!printFixed && length(model$apollo_fixed)>0) output <- output[!(rownames(output) %in% model$apollo_fixed),]
  apollo_print(output) #print(output, digits=4)
  cat('\n')
  
  ### Print diagnostics
  # Model structure
  if(!is.null(model$componentReport)) for(r in model$componentReport){
    test <- !is.null(r$param) && length(r$param)>0 && modelOutput_settings$printModelStructure
    if(test) for(j in r$param) cat(j, '\n', sep='')
    if(test) cat('\n')
  }
  if(!is.null(model$componentReport)) for(r in model$componentReport){
    test <- !is.null(r$data ) && length(r$data )>0 && modelOutput_settings$printDataReport
    if(test) for(j in r$data ) cat(j, '\n', sep='')
    if(test) cat('\n')
  }
  if(exists('r')) rm(r)
  if(exists('j')) rm(j)
  
  ### Fill shorter param names with spaces
  longNames <- names( model$estimate )
  if(length(model$apollo_fixed)>0) longNames <- longNames[-which(longNames %in% model$apollo_fixed)]
  maxLen    <- max(nchar(longNames))
  for(i in 1:length(longNames)) longNames[i] <- paste0(paste0(rep(" ", maxLen-nchar(longNames[i])), collapse=""), longNames[i])
  
  if(printCovar){
    if(printClassical==TRUE){
      cat("\n")
      cat("Classical covariance matrix:\n")
      tmp <- model$varcov
      colnames(tmp) <- longNames
      apollo_print(tmp) #print(tmp, digits=4)
    }
    cat("\n")
    cat("Robust covariance matrix:\n")
    tmp <- model$robvarcov
    colnames(tmp) <- longNames
    apollo_print(tmp) #print(tmp, digits=4)
    if(model$bootstrapSE>0){
      cat("\n")
      cat("Bootstrap covariance matrix:\n")
      tmp <- model$bootvarcov
      colnames(tmp) <- longNames
      apollo_print(tmp) #print(tmp, digits=4)
    }
  }
  
  if(printCorr){
    if(printClassical==TRUE){
      cat("\n")
      cat("Classical correlation matrix:\n")
      tmp <- model$corrmat
      colnames(tmp) <- longNames
      apollo_print(tmp) #print(tmp, digits=4)
    }
    cat("\n")
    cat("Robust correlation matrix:\n")
    tmp <- model$robcorrmat
    colnames(tmp) <- longNames
    apollo_print(tmp) #print(tmp, digits=4)
    if(model$bootstrapSE>0){
      cat("\n")
      cat("Bootstrap correlation matrix:\n")
      tmp <- model$bootcorrmat
      colnames(tmp) <- longNames
      apollo_print(tmp) #print(tmp, digits=4)
    }
  }
  
  if(printOutliers>0){
    outliers <- data.frame(ID=names(model$avgCP), avgChoiceProb=model$avgCP)
    colnames(outliers) <- c("ID","Avg prob per choice")
    if(!model$apollo_control$panelData) colnames(outliers) <- c("row", "Avg prob per choice")
    outliers <- outliers[order(outliers[,2]),]
    if(printOutliers==TRUE) printOutliers=20
    printOutliers=floor(min(printOutliers,nrow(outliers)))
    cat("\n",printOutliers,"worst outliers in terms of lowest average per choice prediction:\n")
    print(outliers[(1:printOutliers),], row.names=FALSE)
  }
  
  if(printChange){
    cat("\nChanges in parameter estimates from starting values:\n")
    tmp <- cbind(model$apollo_beta, model$estimate[names(model$apollo_beta)],
                 model$estimate[names(model$apollo_beta)]-model$apollo_beta)
    colnames(tmp) <- c("Initial", "Estimate", "Difference")
    apollo_print(tmp)
  }
  
  if(printFunctions){
    cat("\nSettings and functions used in model definition:\n")
    cat("\napollo_control")
    cat("\n--------------\n")
    ### change 27 July
    ##txt=t(data.frame(model$apollo_control))
    tmp=model$apollo_control
    tmp$cpp=NULL
    ##tmp$analyticGrad=NULL
    tmp$matrixMult=NULL
    tmp$subMaxV=NULL
    txt=t(data.frame(tmp))
    ### end change 27 July
    colnames(txt)="Value"
    print(txt)
    cat("\nHessian routines attempted")
    cat("\n--------------\n")
    cat(model$hessianMethodsAttempted)
    cat("\n")
    if(!is.null(model$scaling) && !all(model$scaling==1)){
      cat("\nScaling in estimation")
      cat("\n--------------\n")
      txt=as.matrix(model$scaling)
      colnames(txt)="Value"
      print(txt)}
    if(!is.null(model$hessianScaling) && !all(model$hessianScaling==1)){
      cat("\nScaling used in computing Hessian")
      cat("\n--------------\n")
      txt=as.matrix(model$hessianScaling)
      colnames(txt)="Value"
      print(txt)}
    if(is.function(model$apollo_randCoeff)){
      cat("\n\napollo_randCoeff")
      cat("\n----------------\n")
      txt=capture.output(print(model$apollo_randCoeff))
      cat(txt[1:(length(txt)-1)],sep="\n")}
    if(is.function(model$apollo_lcPars)){
      cat("\n\napollo_lcPars")
      cat("\n-------------\n")
      txt=capture.output(print(model$apollo_lcPars))
      cat(txt[1:(length(txt)-1)],sep="\n")}
    cat("\n\napollo_probabilities")
    cat("\n--------------------\n")
    txt=capture.output(print(model$apollo_probabilities))
    cat(txt[1:(length(txt)-1)],sep="\n")
  }
  
  invisible(output)
}
