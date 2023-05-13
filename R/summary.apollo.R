#' Prints summary of Apollo model
#' 
#' Receives an estimated model object and prints a summary using the generic summary function.
#' @param object Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param ... further arguments passed to or from other methods.
#' @param pTwoSided Logical. Should two-sided p-values be printed instead of one-sided p-values. FALSE by default.
#' #' @return nothing.
#' @export
#' @importFrom stats printCoefmat
summary.apollo <- function(object, ... ,pTwoSided=FALSE){
  model = object
  cat("\nApollo model summary\n\n")
  
  cat("Model name                   : ", model$apollo_control$modelName,"\n", sep="")
  cat("Model description            : ", model$apollo_control$modelDescr,"\n", sep="")
  cat("Estimation method            : ", model$estimationRoutine, "\n", sep="")
  if(!is.null(model$nObsTot)) cat("Modelled outcomes            : ",sum(model$nObsTot), "\n", sep="")
  
  if(!model$apollo_control$HB){
    cat("\n")
    cat("LL(final)                    : ",round(model$maximum,2),"\n",sep="")
    cat("Estimated parameters         : ", model$nFreeParams,"\n", sep="")
  }else{
    cat("\n")
    cat("LL(final)                    : ",round(model$maximum,2),"\n",sep="")
    cat("Equiv. estimated parameters  : ", model$nFreeParams,"\n", sep="")
    # if(model$HB_n_nonrandom>0){
    #   cat(" (non-random parameters      : ", model$HB_n_nonrandom,")\n", sep="")
    # }
    # if(model$HB_n_random_means>0){
    #   cat(" (means of random parameters : ", model$HB_n_random_means,")\n", sep="")
    # }
    # if(model$HB_n_random_covar>0){
    #   cat(" (covariance matrix terms    : ", model$HB_n_random_covar,")\n", sep="")
    # }
    cat("\nClassical model fit statistics were calculated at parameter values obtained using averaging across the post burn-in iterations.\n")
  }  
  
  if(!model$apollo_control$HB){
    if(pTwoSided) pMult <- 2 else pMult <- 1
    if(!is.null(model$bootse)){
      se=model$bootse
      if(pMult==1){
        errors="(bootstrap covariance matrix, 1-sided p-values)"
      }else {
        errors="(bootstrap covariance matrix, 2-sided p-values)"
      }
    } else if(!is.null(model$robse)){
      se=model$robse
      if(pMult==1){
        errors="(robust covariance matrix, 1-sided p-values)"
      }else {
        errors="(robust covariance matrix, 2-sided p-values)"
      }
    } else {
      se=model$se
      if(pMult==1){
        errors="(classical covariance matrix, 1-sided p-values)"
      }else {
        errors="(classical covariance matrix, 2-sided p-values)"
      }
    }
    
    if(all(is.na(se))){
      cat("\n\nEstimates (no standard errors computed):\n")
      signif(model$estimate[!names(model$estimate)%in%model$apollo_fixed],4)
    }else{
      cat(paste0("\n\nEstimates ",errors,":\n\n"))
      output=cbind(model$estimate,
                   se,
                   model$estimate/se,
                   pMult*(1-stats::pnorm(abs(model$estimate/se))))
      tmp=c("estimate","std. error","t-ratio")
      #tmp=c(tmp,"p-value")
      if(pMult==1){
        tmp=c(tmp,"p (1-sided)")
      } else{
        tmp=c(tmp,"p (2-sided)")
      }
      colnames(output)=tmp
      output=output[!(rownames(output)%in%model$apollo_fixed),]
      #printCoefmat(output,digits= max(3L, getOption("digits") - 3L),P.values=TRUE,has.Pvalue=TRUE)
      stats::printCoefmat(output,digits= 2,P.values=TRUE,has.Pvalue=TRUE)
    }
  }else{
    cat("\nSummary of parameter chains\n")
    
    if(length(model$apollo_HB$gVarNamesFixed)>0 | length(model$apollo_fixed)>0){
      cat("\nNon-random coefficients","\n")
      print(round(model$HB_chains_non_random,4))
    }
    
    if(any(!is.null(model$apollo_HB$gVarNamesNormal)) && length(model$apollo_HB$gVarNamesNormal)>0){
      cat("\nResults for posterior means for random coefficients","\n")
      print(round(model$HB_posterior_means_summary,4))
    }
  }
  
  cat("\n")
  apollo_print("For more detailed output, use apollo_modelOutput")
  
}
