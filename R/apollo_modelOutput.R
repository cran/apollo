#' Prints estimation results to console
#' 
#' Prints estimation results to console. Amount of information presented can be adjusted through arguments.
#' 
#' Prints to screen the output of a model previously estimated by apollo_estimate()
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param modelOutput_settings List of options. It can include the following.
#'                             \itemize{
#'                               \item printClassical: Boolean. TRUE for printing classical standard errors. TRUE by default.
#'                               \item printPVal: Boolean. TRUE for printing p-values. FALSE by default.
#'                               \item printT1: Boolean. If TRUE, t-test for H0: apollo_beta=1 are printed. FALSE by default.
#'                               \item printDiagnostics: Boolean. TRUE for printing summary of choices in database and other diagnostics. TRUE by default.
#'                               \item printCovar: Boolean. TRUE for printing parameters covariance matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. FALSE by default.
#'                               \item printCorr: Boolean. TRUE for printing parameters correlation matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. FALSE by default.
#'                               \item printOutliers: Boolean. TRUE for printing 20 individuals with worst average fit across observations. FALSE by default.
#'                               \item printChange: Boolean. TRUE for printing difference between starting values and estimates. FALSE by default.
#'                             }
#' @return A matrix of coefficients, s.d. and t-tests (invisible)
#' @export
#' @importFrom coda geweke.diag
apollo_modelOutput=function(model, modelOutput_settings=NA){
  if(length(modelOutput_settings)==1 && is.na(modelOutput_settings)) modelOutput_settings=list()
  if(is.null(modelOutput_settings[["printClassical"]])) modelOutput_settings[["printClassical"]] = TRUE
  if(is.null(modelOutput_settings[["printPVal"]])) modelOutput_settings[["printPVal"]] = FALSE
  if(is.null(modelOutput_settings[["printT1"]])) modelOutput_settings[["printT1"]] = FALSE
  if(is.null(modelOutput_settings[["printDiagnostics"]])) modelOutput_settings[["printDiagnostics"]] = TRUE
  if(is.null(modelOutput_settings[["printCovar"]])) modelOutput_settings[["printCovar"]] = FALSE
  if(is.null(modelOutput_settings[["printCorr"]])) modelOutput_settings[["printCorr"]] = FALSE
  if(is.null(modelOutput_settings[["printOutliers"]])) modelOutput_settings[["printOutliers"]] = FALSE
  if(is.null(modelOutput_settings[["printChange"]])) modelOutput_settings[["printChange"]] = FALSE
  
  printClassical   = modelOutput_settings[["printClassical"]]
  printPVal        = modelOutput_settings[["printPVal"]]
  printT1          = modelOutput_settings[["printT1"]]
  printDiagnostics = modelOutput_settings[["printDiagnostics"]]
  printCovar       = modelOutput_settings[["printCovar"]]
  printCorr        = modelOutput_settings[["printCorr"]]
  printOutliers    = modelOutput_settings[["printOutliers"]]
  printChange      = modelOutput_settings[["printChange"]]
  
  apollo_control <- model$apollo_control
  nParams     <- length(model$apollo_beta)
  nFreeParams <- nParams
  if(!is.null(model$apollo_fixed)) nFreeParams <- nFreeParams - length(model$apollo_fixed)
  
  
  apolloVersion <- tryCatch(utils::packageDescription("apollo", fields = "Version"),
                            warning=function(w) return("alpha"),
                            error=function(e) return("alpha"))
  cat("Model run using Apollo for R, version", apolloVersion,"\n")
  cat("www.apollochoicemodelling.com\n\n")
  cat("Model name                       : ", model$apollo_control$modelName,"\n", sep="")
  cat("Model description                : ", model$apollo_control$modelDescr,"\n", sep="")
  cat("Model run at                     : ", paste(model$startTime),"\n", sep="")
  cat("Estimation method                : ", model$estimationRoutine, "\n", sep="")
  if(!apollo_control$HB) cat("Model diagnosis                  : ",model$message,"\n", sep="")
  cat("Number of individuals            : ", model$nIndivs,"\n", sep="")
  cat("Number of observations           : ", model$nObs,"\n", sep="")
  cat("Estimated parameters             : ", nFreeParams,"\n", sep="")
  if(model$apollo_control$mixing){
    d <- model$apollo_draws
    if(d$interNDraws>0 && length(c(d$interUnifDraws, d$interNormDraws))>0){
      cat("Number of inter-person draws     : ", d$interNDraws, ' (', d$interDrawsType, ')', "\n", sep='')
    }
    if(d$intraNDraws>0 && length(c(d$intraUnifDraws, d$intraNormDraws))>0){
      cat("Number of intra-person draws     : ", d$intraNDraws, ' (', d$intraDrawsType, ')', "\n", sep='')
    }
    if(!model$apollo_control$panelData & model$apollo_control$mixing & d$interNDraws>0){
      cat("WARNING: Inter-person draws were used\n")
      cat("         without a panel data structure.\n")
    }
    rm(d)
  } else { if(!apollo_control$HB) cat("Model with no mixing\n") }
  cat("\n")
  
  
  if(apollo_control$HB){
    apollo_HB <- model$apollo_HB
    cat("Average posterior log-likelihood post burn-in",mean(colSums(log(model$cmcLLout))),"\n")
    cat("Average posterior RLH post burn-in",mean(colMeans((model$cmcRLHout))),"\n")
    cat("\n\n")
    
    cat("Chain convergence report:\n")
    if(!is.null(model$F)){
      cat("Fixed (non random) parameters:\n")
      tmp <- coda::geweke.diag(model$F[,2:(ncol(model$F))], frac1=0.1, frac2=0.5)[[1]]
      names(tmp) <- model$params.fixed
      print( round(tmp, 4) )
      cat("\n")
    }
    if(!is.null(model$A)){
      cat("Random parameters:\n")
      tmp <- coda::geweke.diag(model$A[,2:(ncol(model$A))], frac1=0.1, frac2=0.5)[[1]]
      print( round(tmp, 4) )
      cat("\n")
    }
    if(!is.null(model$D)){
      cat("Covariances of random parameters:\n")
      tmp <- c()
      for(i in 1:dim(model$D)[1]) for(j in 1:i){
        if(i==1 & j==1) Dmatrix <- as.vector(model$D[i,j,]) else Dmatrix <- cbind(Dmatrix, as.vector(model$D[i,j,]))
        tmp <- c(tmp, paste(colnames(model$A)[i+1],colnames(model$A)[j+1], sep="_"))
      }
      colnames(Dmatrix) <- tmp
      tmp <- coda::geweke.diag(Dmatrix, frac1=0.1, frac2=0.5)[[1]]
      print( round(tmp, 4) )
    }
    cat("\n\n")
    
    
    cat("Summary of parameters chains:\n")
    ans <- list()
    
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
      cat("Non-random coefficients","\n")
      colnames(non_random)=c("Mean","SD")
      originalOrder <- names(model$apollo_beta)[names(model$apollo_beta) %in% rownames(non_random)]
      print(round(non_random[originalOrder,],4))
      cat("\n")
      
      ans[["non_random"]] <- non_random
    }
    
    apollo_HB$gVarNamesFixed <- model$params.fixed
    apollo_HB$gVarNamesNormal <- model$params.vary
    if(any(!is.null(apollo_HB$gVarNamesNormal)) && length(apollo_HB$gVarNamesNormal)>0){
      random_mean     = matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
      random_mean[,1] = colMeans(model$A)[2:ncol(model$A)]
      random_mean[,2] = apply(model$A,FUN=stats::sd,2)[2:ncol(model$A)]
      rownames(random_mean)=apollo_HB$gVarNamesNormal
      cat("Upper level model results for mean parameters for underlying Normals","\n")
      colnames(random_mean)=c("Mean","SD")
      print(round(random_mean,4))
      cat("\n")
      
      random_cov_mean           = apply(model$D,FUN=mean,c(1,2))
      random_cov_sd             = apply(model$D,FUN=stats::sd,c(1,2))
      rownames(random_cov_mean) = apollo_HB$gVarNamesNormal
      colnames(random_cov_mean) = apollo_HB$gVarNamesNormal
      cat("Upper level model results for covariance matrix for underlying Normals (means across iterations)","\n")
      print(round(random_cov_mean,4))
      cat("\n")
      
      rownames(random_cov_sd) = apollo_HB$gVarNamesNormal
      colnames(random_cov_sd) = apollo_HB$gVarNamesNormal
      cat("Upper level model results for covariance matrix for underlying Normals (SD across iterations)","\n")
      print(round(random_cov_sd,4))
      cat("\n")
      
      posterior=matrix(0,nrow=length(apollo_HB$gVarNamesNormal),2)
      posterior[,1]=colMeans(model$C)[3:ncol(model$C)]
      posterior[,2]=apply(model$C,FUN=stats::sd,2)[3:ncol(model$C)]
      rownames(posterior)=apollo_HB$gVarNamesNormal
      cat("Results for posterior means for random coefficients","\n")
      colnames(posterior)=c("Mean","SD")
      print(round(posterior,4))
      cat("\n")
      
      ans[["random_mean"]]   <- random_mean
      ans[["random_cov_mean"]] <- random_cov_mean
      ans[["random_cov_sd"]] <- random_cov_sd
      ans[["posterior"]]     <- posterior
    }
    
    return(invisible(ans))
  }
  
  
  output=cbind(round(model$estimate,4),
               round(model$se,4),
               round(model$estimate/model$se,2),
               round(2*(1-stats::pnorm(abs(model$estimate/model$se))),3),
               round((model$estimate-1)/model$se,2),
               round(2*(1-stats::pnorm(abs((model$estimate-1)/model$se))),3),
               round(model$robse,4),
               round(model$estimate/model$robse,2),
               round(2*(1-stats::pnorm(abs(model$estimate/model$robse))),3),
               round((model$estimate-1)/model$robse,2),
               round(2*(1-stats::pnorm(abs((model$estimate-1)/model$robse))),3))
  colnames(output) <- c('Estimate','Std.err.','t.ratio(0)','p-val(0)','t.ratio(1)','p-val(1)','Rob.std.err.','Rob.t.ratio(0)','Rob.p-val(0)','Rob.t.ratio(1)','Rob.p-val(1)')
  rownames(output) <- names(model$estimate)
  
  dropcolumns=NULL
  if(printClassical==FALSE) dropcolumns = c(dropcolumns,2,3,4,5,6)
  if(printT1==FALSE) dropcolumns = c(dropcolumns,5,6,10,11)
  if(printPVal==FALSE) dropcolumns = c(dropcolumns,4,6,9,11)
  dropcolumns = unique(dropcolumns)
  if(length(dropcolumns)>0) output = output[,-dropcolumns, drop=FALSE]
  
  cat("LL(start)                        : ",model$LLStart,"\n", sep="")
  if(!anyNA(model$LL0)) cat("LL(0)                            : ",model$LL0,"\n",sep="")
  if(length(model$LLout)==1) cat("LL(final)                        : ",model$maximum,"\n",sep="")
  if(length(model$LLout)>1) cat("LL(final, whole model)           : ",model$maximum,"\n",sep="")
  if(length(model$LLout)>1){
    j=2
    nameList <- names(model$LLout)
    while(j<=length(model$LLout)){
      spaces <- 33-(6+nchar(nameList[j]))
      if(spaces>0) spaces <- paste(rep(" ",spaces),sep='') else spaces <- ""
      cat("  LL(",nameList[j],")",spaces,": ",model$LLout[j],"\n",sep="")
      j=j+1
    }
  }
  if(!anyNA(model$LL0)){
    cat("Rho-square (0)                   : ",round(1-(model$maximum/model$LL0),4),"\n")
    cat("Adj.Rho-square (0)               : ",round(1-((model$maximum-nFreeParams)/model$LL0),4),"\n")}
  
  
  cat("AIC                              : ",round(-2*model$maximum + 2*nFreeParams,2),"\n")
  cat("BIC                              : ",round(-2*model$maximum + nFreeParams*log(model$nObs),2),"\n")
  tmpH <- floor(model$timeTaken/60^2)
  tmpM <- floor((model$timeTaken-tmpH*60^2)/60)
  tmpS <- round(model$timeTaken-tmpH*60^2-tmpM*60,2)
  timeTaken <- paste(formatC(tmpH,width=2,format='d',flag=0),
                     formatC(tmpM,width=2,format='d',flag=0),
                     tmpS,sep=':')
  cat("Time taken (hh:mm:ss)            : ",timeTaken,"\n")
  cat("Iterations                       : ",model$nIter,"\n")
  cat("Number of cores used             : ",model$apollo_control$nCores,"\n")
  cat("\n")
  
  if(!printClassical & anyNA(model$se[!(names(model$estimate) %in% model$apollo_fixed)]) ){
    cat("\nWARNING: Some parameters classical standard errors could not be calculated.")
    cat("\n         This could point to an identification or estimation problem.")
  }
  
  cat("Estimates:\n")
  if(nrow(output)>options("max.print")) options(max.print=nrow(output)+100)
  print(output)
  cat('\n')
    cat("\n")
    
    fileName <- paste(model$apollo_control$modelName, "_tempOutput.txt", sep="")
    fileName <- file.path(tempdir(),fileName)
    if(file.exists(fileName)){
      fileConn <- tryCatch(file(fileName, open="rt"), error=function(e) NA)
      if(!is.na(fileConn)){
        txt <- readLines(fileConn)
        for(l in txt) cat(l ,"\n")
        close(fileConn)
      } else cat('Could not read additional output from temp file.\n')
    }
  
  if(printCovar){
    if(printClassical==TRUE){
      cat("\n")
      cat("Classical covariance matrix:\n")
      print(round(model$varcov,6))
    }
    cat("\n")
    cat("Robust covariance matrix:\n")
    print(round(model$robvarcov,6))
  }
  
  if(printCorr){
    if(printClassical==TRUE){
      cat("\n")
      cat("Classical correlation matrix:\n")
      print(round(model$corrmat,6))
    }
    cat("\n")
    cat("Robust correlation matrix:\n")
    print(round(model$robcorrmat,6))
  }
  
  if(printOutliers){
    outliers <- data.frame(ID=names(model$avgCP), avgChoiceProb=model$avgCP)
    colnames(outliers) <- c("ID","Avg prob per choice")
    outliers <- outliers[order(outliers[,2]),]
    
    cat("\n20 worst outliers in terms of lowest average per choice prediction:\n")
    print(outliers[(1:20),], row.names=FALSE)
  }
  
  if(printChange){
    cat("\nChanges in parameter estimates from starting values:\n")
    x <- cbind(model$apollo_beta, model$estimate[names(model$apollo_beta)],
               model$estimate[names(model$apollo_beta)]-model$apollo_beta)
    x <- round(x,4)
    colnames(x) <- c("Initial", "Estimate", "Difference")
    print(x)
  }
  
  invisible(output)
}
