#' Prints estimation results to console
#'
#' Prints estimation results to console. Amount of information presented can be adjusted through arguments.
#'
#' Prints to screen the output of a model previously estimated by apollo_estimate()
#' @param model Model object. An estimated model as returned by \code{apollo_estimate}
#' @param printClassical Boolean. TRUE (default) for printing classical standard errors.
#' @param printPVal Boolean. TRUE for printing p-values. FALSE by default.
#' @param printT1 Boolean. TRUE for printing t-test for H0: apollo_beta=1 in addition to H0: apollo_beta=0. FALSE by default.
#' @param printDiagnostics Boolean. TRUE (default) for printing summary of choices in database and other diagnostics.
#' @param printCovar Boolean. TRUE (default) for printing parameter covariance matrix.
#'                   If \code{printClassical=TRUE}, both the classical and robust covariance matrices are printed.
#' @param printCorr Boolean. TRUE (default) for printing parameter correlation matrix.
#'                  If \code{printClassical=TRUE}, both the classical and robust correlation matrices are printed.
#' @param printOutliers Boolean. TRUE (default) for printing 20 individuals with worst average fit across observations.
#' @param printChange Boolean. TRUE (default) for printing difference between starting values and estimates.
#' @return A matrix of coefficients, s.d. and t-tests (invisible)
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
apollo_modeloutput=function(model, printClassical=TRUE, printPVal=FALSE, printT1=FALSE, printDiagnostics=TRUE, printCovar=FALSE, printCorr=FALSE,
                         printOutliers=FALSE, printChange=FALSE)
{
  # ####################### #
  #### MODEL DESCRIPTION ####
  # ####################### #

  apollo_control <- model$apollo_control
  nParams     <- length(model$apollo_beta)
  nFreeParams <- nParams
  if(!is.null(model$apollo_beta_fixed)) nFreeParams <- nFreeParams - length(model$apollo_beta_fixed)


  # Printing model information
  cmcRcodeVersion <- tryCatch(utils::packageDescription("cmcRcode", fields = "Version"),
                              warning=function(w) return("alpha"),
                              error=function(e) return("alpha"))
  cat("Model run using CMC choice modelling code for R, version", cmcRcodeVersion,"\n")
  cat("www.cmc.leeds.ac.uk\n\n")
  cat("Model name                       : ", model$apollo_control$modelName,"\n", sep="")
  cat("Model description                : ", model$apollo_control$modelDescr,"\n", sep="")
  cat("Model run at                     : ", paste(model$startTime),"\n", sep="")
  cat("Estimation method                : ", model$estimation_routine, "\n", sep="")
  if(!apollo_control$HB) cat("Model diagnosis                  : ",model$message,"\n", sep="")
  cat("Number of decision makers        : ", model$nIndivs,"\n", sep="")
  cat("Number of observations           : ", model$nObs,"\n", sep="")
  cat("Estimated parameters             : ", nFreeParams,"\n", sep="")
  if(model$apollo_control$mixing){
    cat("Number of inter-person draws     : ", model$apollo_draws$inter_nDraws, ' (', model$apollo_draws$inter_drawsType, ')', "\n", sep='')
    cat("Number of intra-person draws     : ", model$apollo_draws$intra_nDraws, ' (', model$apollo_draws$intra_drawsType, ')', "\n", sep='')
    if(!model$apollo_control$panelData & model$apollo_control$mixing & model$apollo_draws$inter_nDraws>0){
      cat("WARNING: Inter-person draws were used\n")
      cat("         without a panel structure.\n")
    }
  } else { if(!apollo_control$HB) cat("Model with no mixing\n") }
  cat("\n")


  # ####################### #
  #### HB OUTPUT         ####
  # ####################### #
  if(apollo_control$HB){
    cat("Average posterior log-likelihood post burn-in",mean(colSums(log(model$cmcLLout))),"\n")
    cat("Average posterior RLH post burn-in",mean(colMeans((model$cmcRLHout))),"\n")
    cat("\n\n")

    cat("MCMC convergence report:\n")
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
      # This assumes the matrix is square
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
    HB_control <- model$HB_control
    if(length(HB_control$gVarNamesFixed)>0 | length(model$apollo_beta_fixed)>0){
      if(length(HB_control$gVarNamesFixed)>0){
        non_random=matrix(0,nrow=length(HB_control$gVarNamesFixed),2)
        non_random[,1]=colMeans(model$F)[2:ncol(model$F)]
        non_random[,2]=apply(model$F,FUN=stats::sd,2)[2:ncol(model$F)]
        rownames(non_random)=HB_control$gVarNamesFixed}
      if(length(model$apollo_beta_fixed)>0){
        if(length(HB_control$gVarNamesFixed)>0){
          non_random=rbind(non_random,cbind(matrix(model$apollo_beta[model$apollo_beta_fixed]),NA))
          rownames(non_random)[(length(HB_control$gVarNamesFixed)+1):nrow(non_random)]=model$apollo_beta_fixed
        } else{
          non_random=cbind(matrix(model$apollo_beta[model$apollo_beta_fixed]),NA)
          rownames(non_random)=model$apollo_beta_fixed
        }
      }
      cat("Non-random coefficients","\n")
      colnames(non_random)=c("Mean","SD")
      originalOrder <- names(model$apollo_beta)[names(model$apollo_beta) %in% rownames(non_random)]
      print(round(non_random[originalOrder,],4))
      cat("\n")

      ans[["non_random"]] <- non_random
    }

    if(any(!is.null(HB_control$gVarNamesNormal)) && length(HB_control$gVarNamesNormal)>0){
      random_mean=matrix(0,nrow=length(HB_control$gVarNamesNormal),2)
      random_mean[,1]=colMeans(model$A)[2:ncol(model$A)]
      random_mean[,2]=apply(model$A,FUN=stats::sd,2)[2:ncol(model$A)]
      rownames(random_mean)=HB_control$gVarNamesNormal
      cat("Upper level model results for mean parameters for underlying Normals","\n")
      colnames(random_mean)=c("Mean","SD")
      print(round(random_mean,4))
      cat("\n")

      random_cov_mean=apply(model$D,FUN=mean,c(1,2))
      random_cov_sd=apply(model$D,FUN=stats::sd,c(1,2))
      rownames(random_cov_mean)=HB_control$gVarNamesNormal
      colnames(random_cov_mean)=HB_control$gVarNamesNormal
      cat("Upper level model results for covariance matrix for underlying Normals (means across iterations)","\n")
      print(round(random_cov_mean,4))
      cat("\n")

      rownames(random_cov_sd)=HB_control$gVarNamesNormal
      colnames(random_cov_sd)=HB_control$gVarNamesNormal
      cat("Upper level model results for covariance matrix for underlying Normals (SD across iterations)","\n")
      print(round(random_cov_sd,4))
      cat("\n")

      posterior=matrix(0,nrow=length(HB_control$gVarNamesNormal),2)
      posterior[,1]=colMeans(model$C)[3:ncol(model$C)]
      posterior[,2]=apply(model$C,FUN=stats::sd,2)[3:ncol(model$C)]
      rownames(posterior)=HB_control$gVarNamesNormal
      cat("Results for posterior means for random coefficients","\n")
      colnames(posterior)=c("Mean","SD")
      print(round(posterior,4))
      cat("\n")

      ans[["random_mean"]]   <- random_mean
      ans[["random_cov_sd"]] <- random_cov_sd
      ans[["posterior"]]     <- posterior
    }

    return(invisible(ans))
  }


  # ####################### #
  #### CLASSICAL OUTPUT  ####
  # ####################### #
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
  colnames(output) <- c('Estimate','Std.err.','t.ratio(0)','p-value(0)','t.ratio(1)','p-value(1)','Rob.std.err.','Rob.t.ratio(0)','Rob.p-value(0)','Rob.t.ratio(1)','Rob.p-value(1)')
  rownames(output) <- names(model$estimate)

  dropcolumns=NULL
  if(printClassical==FALSE) dropcolumns = c(dropcolumns,2,3,4,5,6)
  if(printT1==FALSE) dropcolumns = c(dropcolumns,5,6,10,11)
  if(printPVal==FALSE) dropcolumns = c(dropcolumns,4,6,9,11)
  dropcolumns = unique(dropcolumns)
  if(length(dropcolumns)>0) output = output[,-dropcolumns, drop=FALSE]

  cat("LL(start)                        : ",model$LLStart,"\n", sep="")
  if(!is.na(model$LL0)) cat("LL(0)                            : ",model$LL0,"\n",sep="")
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
  if(!is.na(model$LL0)){
    cat("Rho-square (0)                   : ",round(1-(model$maximum/model$LL0),4),"\n")
    cat("Adj.Rho-square (0)               : ",round(1-((model$maximum-nFreeParams)/model$LL0),4),"\n")}


  cat("AIC                              : ",round(-2*model$maximum + 2*nFreeParams,2),"\n")
  cat("BIC                              : ",round(-2*model$maximum + nFreeParams*log(model$nObs),2),"\n\n")
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

  if(!printClassical & anyNA(model$se[!(names(model$estimate) %in% model$apollo_beta_fixed)]) ){
    cat("\nWARNING: Some parameters classical standard errors could not be calculated.")
    cat("\          This could point to an identification or estimation problem.")
  }

  cat("Estimates:\n")
  if(nrow(output)>options("max.print")) options(max.print=nrow(output)+100)
  print(output)
  cat('\n')
  if(length(model$apollo_beta_fixed)>0) cat("The following parameters were fixed (they have no std.err.):\n",
                                         paste(model$apollo_beta_fixed, collapse=", ", sep=""),
                                         "\n", sep="")

  if(printDiagnostics==TRUE){
    cat("\n")

    fileName <- paste(model$apollo_control$modelName, "_tempOutput.txt", sep="")
    fileName <- file.path(tempdir(),fileName)
    if(file.exists(fileName)){
      openSuccesfuly <- TRUE
      tryCatch(fileConn <- file(fileName, open="rt"), error=function(e) openSuccesfuly <- FALSE)
      if(openSuccesfuly){
        nTree <- 1
        txt <- readLines(fileConn, n=1)
        while(length(txt)>0){
          if( grepl("Tree structure", txt) ){
            cat('[Tree structure ',nTree,"]\n", sep="")
            nTree <- nTree + 1
          } else cat(txt,'\n',sep="")
          txt <- readLines(fileConn, n=1)
        }
        close(fileConn)
      } else cat('Could not read tree structure(s) from temp file.\n')

    }
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
