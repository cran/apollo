#' Prints brief summary of Apollo model
#' 
#' Receives an estimated model object and prints a brief summary using the generic print function.
#' @param x Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param ... further arguments passed to or from other methods.
#' @return nothing.
#' @export
print.apollo <- function(x, ...){
  model=x
  cat("\nApollo model summary\n\n")
  
  cat("Model name                   : ", model$apollo_control$modelName,"\n", sep="")
  cat("Model description            : ", model$apollo_control$modelDescr,"\n", sep="")
  cat("Estimation method            : ", model$estimationRoutine, "\n", sep="")
  
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
    cat("\nEstimates:\n")
    print(signif(model$estimate[!names(model$estimate)%in%model$apollo_fixed],4))
  }else{
    cat("\nSummary of parameter chains\n")
    
    if(length(model$apollo_HB$gVarNamesFixed)>0 | length(model$apollo_fixed)>0){
      cat("\nNon-random coefficients","\n")
      print(round(model$HB_chains_non_random[,1],4))
    }
    
    if(any(!is.null(model$apollo_HB$gVarNamesNormal)) && length(model$apollo_HB$gVarNamesNormal)>0){
      cat("\nResults for posterior means for random coefficients","\n")
      print(round(model$HB_posterior_means_summary[,1],4))
    }
  }
  cat("\n")
  apollo_print("For more detailed output, use summary, or apollo_modelOutput for full outputs")
}