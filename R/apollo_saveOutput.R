#' Saves estimation results to files.
#' 
#' Writes files in the working directory with the estimation results.
#' 
#' Estimation results are printed to different files in the working directory:
#' \describe{
#'   \item{(modelName)_output.txt}{Text file with the output produced by function \code{apollo_modelOutput}.}
#'   \item{(modelName)_estimates.csv}{CSV file with the estimated parameter values, their standars errors, and t-ratios.}
#'   \item{(modelName)_covar.csv}{CSV file with the estimated classical covariance matrix. Only when bayesian estimation was not used.}
#'   \item{(modelName)_robcovar.csv}{CSV file with the estimated robust covariance matrix. Only when bayesian estimation was not used.}
#'   \item{(modelName)_corr.csv}{CSV file with the estimated classical correlation matrix. Only when bayesian estimation was not used.}
#'   \item{(modelName)_robcorr.csv}{CSV file with the estimated robust correlation matrix. Only when bayesian estimation was not used.}
#'   \item{(modelName).F12}{F12 file with model results. Compatible with ALOGIT.}
#' }
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param saveOutput_settings List of options. Valid options are the following.
#'                            \itemize{
#'                               \item printClassical: Boolean. TRUE for printing classical standard errors. TRUE by default.
#'                               \item printPVal: Boolean. TRUE for printing p-values. FALSE by default.
#'                               \item printT1: Boolean. If TRUE, t-test for H0: apollo_beta=1 are printed. FALSE by default.
#'                               \item printDiagnostics: Boolean. TRUE for printing summary of choices in database and other diagnostics. TRUE by default.
#'                               \item printCovar: Boolean. TRUE for printing parameters covariance matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. TRUE by default.
#'                               \item printCorr: Boolean. TRUE for printing parameters correlation matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. TRUE by default.
#'                               \item printOutliers: Boolean or Scalar. TRUE for printing 20 individuals with worst average fit across observations. FALSE by default. If Scalar is given, this replaces the default of 20.
#'                               \item printChange: Boolean. TRUE for printing difference between starting values and estimates. TRUE by default.
#'                               \item saveEst: Boolean. TRUE for saving estimated parameters and standard errors to a CSV file. TRUE by default.
#'                               \item saveCov: Boolean. TRUE for saving estimated correlation matrix to a CSV file. TRUE by default.
#'                               \item saveCorr: Boolean. TRUE for saving estimated correlation matrix to a CSV file. TRUE by default.
#'                               \item saveModelObject: Boolean. TRUE to save the R model object to a file (use \link{apollo_loadModel} to load it to memory). TRUE by default.
#'                               \item writeF12: Boolean. TRUE for writing results into an F12 file (ALOGIT format). FALSE by default.
#'                            }
#' @return nothing
#' @export
#' @importFrom RSGHB writeModel
apollo_saveOutput=function(model, saveOutput_settings=NA){
  if(length(saveOutput_settings)==1 && is.na(saveOutput_settings)) saveOutput_settings=list()
  default <- list(printClassical   = TRUE,
                  printPVal        = FALSE,
                  printT1          = FALSE,
                  printDiagnostics = TRUE,
                  printCovar       = TRUE,
                  printCorr        = TRUE,
                  printOutliers    = TRUE,
                  printChange      = TRUE,
                  saveEst          = TRUE,
                  saveCov          = TRUE,
                  saveCorr         = TRUE,
                  saveModelObject  = TRUE,
                  writeF12         = FALSE)
  tmp <- names(default)[ !(names(default) %in% names(saveOutput_settings)) ]
  for(i in tmp)  saveOutput_settings[[i]] <- default[[i]]
  rm(tmp, default)
  printClassical  = saveOutput_settings[["printClassical"]]
  printPVal       = saveOutput_settings[["printPVal"]]
  printT1         = saveOutput_settings[["printT1"]]
  printDiagnostics= saveOutput_settings[["printDiagnostics"]]
  printCovar      = saveOutput_settings[["printCovar"]]
  printCorr       = saveOutput_settings[["printCorr"]]
  printOutliers   = saveOutput_settings[["printOutliers"]]
  printChange     = saveOutput_settings[["printChange"]]
  saveEst         = saveOutput_settings[["saveEst"]]
  saveCov         = saveOutput_settings[["saveCov"]]
  saveCorr        = saveOutput_settings[["saveCorr"]]
  saveModelObject = saveOutput_settings[["saveModelObject"]]
  writeF12        = saveOutput_settings[["writeF12"]]
  
  if(length(model$scaling)>0 && !is.na(model$scaling)){
    scaling_used=TRUE
  }else{
    scaling_used=FALSE
  }
  
  sink(paste(model$apollo_control$modelName,"_output.txt",sep=""))
  tmp <- apollo_modelOutput(model,saveOutput_settings)
  sink()
  
  # ################################## #
  #### HB only components           ####
  # ################################## #
  
  if(model$apollo_control$HB){
    if(saveEst){
      RSGHB::writeModel(model, writeDraws = FALSE, path = getwd())
      cat("RSGHB output saved in following files\n")
      cat("\nOutputs at iteration level (post burn-in chains)\n")
      if(!is.null(tmp$non_random)) cat("Non-random parameters:",paste(model$apollo_control$modelName, "_F"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      cat("\n")
      if(!is.null(tmp$random_mean)) cat("Means for underlying normals:",paste(model$apollo_control$modelName, "_A"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      cat("\n")
      cat("\nPosteriors\n")
      if(!is.null(tmp$random_mean)) cat("Mean individual-level draws for underlying normals:",paste(model$apollo_control$modelName, "_B"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      cat("\n")
      if(!is.null(tmp$random_mean)) cat("SD of individual-level draws for underlying normals:",paste(model$apollo_control$modelName, "_Bsd"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      cat("\n")
      if(!is.null(tmp$random_mean)) cat("Mean individual-level draws after transformations to underlying normals:",paste(model$apollo_control$modelName, "_C"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      cat("\n")
      if(!is.null(tmp$random_mean)) cat("SD of individual-level draws after transformations to underlying normals:",paste(model$apollo_control$modelName, "_Csd"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      cat("\n")
      if(!is.null(tmp$random_mean)) cat("Sample variance-covariance matrix for underlying normals:",paste(model$apollo_control$modelName, "_D"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      cat("\n")
      cat("\nRSGHB log file saved to",paste(model$apollo_control$modelName, ".log", sep=""),"\n")
      cat("\n")
      if(!is.null(tmp$non_random)){
      cat("\nAdditional output files:\n")
      utils::write.csv(tmp$non_random   , paste(model$apollo_control$modelName, "_param_non_random"   ,".csv", sep=""))
      cat("Summary of chains for non-random parameters:",paste(model$apollo_control$modelName, "_param_non_random"   ,".csv", sep=""),"\n")}
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      if(!is.null(tmp$random_mean)){utils::write.csv(tmp$random_mean  , paste(model$apollo_control$modelName, "_param_random_mean"  ,".csv", sep=""))
      cat("\n")
      cat("Summary of chains for means of normals:",paste(model$apollo_control$modelName, "_param_random_mean"   ,".csv", sep=""),"\n")}
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      if(!is.null(tmp$random_cov_mean)){utils::write.csv(tmp$random_cov_mean, paste(model$apollo_control$modelName, "_param_random_cov_mean",".csv", sep=""))
      cat("\n")
      cat("Means of chains for covariance of normals:",paste(model$apollo_control$modelName, "_param_random_cov_mean"   ,".csv", sep=""),"\n")}
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      if(!is.null(tmp$random_cov_sd)){utils::write.csv(tmp$random_cov_sd, paste(model$apollo_control$modelName, "_param_random_cov_sd",".csv", sep=""))
      cat("\n")
      cat("SDs of chains for covariance of normals:",paste(model$apollo_control$modelName, "_param_random_cov_sd"   ,".csv", sep=""),"\n")}
      if(scaling_used) cat("These outputs have NOT had the scaling used in estimation applied to them\n")
      if(!is.null(tmp$posterior)){utils::write.csv(tmp$posterior    , paste(model$apollo_control$modelName, "_param_posterior"    ,".csv", sep=""))
      cat("\n")
      cat("Summary of posteriors for random parameters:",paste(model$apollo_control$modelName, "_param_posterior"   ,".csv", sep=""),"\n")}
      if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")
      cat("\n")
     }
    if(saveModelObject){
      tryCatch( {
        saveRDS(model, file=paste0(model$apollo_control$modelName,"_model.rds"))
        cat("\nModel object saved to",paste(model$apollo_control$modelName, ".rds", sep=""),"\n")
        cat("\n")
        }, error=function(e) cat("Model object could not be written to file."))
    }
    return(invisible(TRUE))
  }
  
  # ################################## #
  #### for all                      ####
  # ################################## #
  
  if(saveEst){
    # Build matrix with results table
    output <- matrix(model$estimate, nrow=length(model$estimate), ncol=1, dimnames=list(names(model$estimate)))
    colnames(output) <- c("Estimate")
    if(printClassical){
      output <- cbind(output, Std.err.=model$se, `t-ratio(0)`=model$estimate/model$se)
      if(printPVal) output <- cbind(output, `p-val(0)`=2*(1-stats::pnorm(abs(model$estimate/model$se))) )
      if(printT1){
        output <- cbind(output, `t-ratio(1)`=(model$estimate-1)/model$se)
        if(printPVal) output <- cbind(output, `p-val(1)`=2*(1-stats::pnorm(abs((model$estimate-1)/model$se))) )
      }
    }
    output <- cbind(output, Rob.std.err.=model$robse, `Rob.t-ratio(0)`=model$estimate/model$robse)
    if(printPVal) output <- cbind(output, `Rob.p-val(0)`=2*(1-stats::pnorm(abs(model$estimate/model$robse))) )
    if(printT1){
      output <- cbind(output, `Rob.t-ratio(1)`=(model$estimate-1)/model$robse)
      if(printPVal) output <- cbind(output, `Rob.p-val(1)`=2*(1-stats::pnorm(abs((model$estimate-1)/model$robse))) )
    }
    
    # Write to file
    utils::write.csv(output,paste(model$apollo_control$modelName,"_estimates.csv",sep=""))
    cat("Estimates saved to",paste(model$apollo_control$modelName, "_estimates.csv"   , sep=""),"\n")
  }
  if(saveCov){
    if(printClassical==TRUE){
      utils::write.csv(model$varcov,paste(model$apollo_control$modelName,"_covar.csv",sep=""))
      cat("Classical covariance matrix saved to",paste(model$apollo_control$modelName, "_covar.csv"   , sep=""),"\n")
      }
    utils::write.csv(model$robvarcov,paste(model$apollo_control$modelName,"_robcovar.csv",sep=""))
    cat("Robust covariance matrix saved to",paste(model$apollo_control$modelName, "_robcovar.csv"   , sep=""),"\n")
    if(!is.null(model$bootstrapSE) && model$bootstrapSE>0){
      utils::write.csv(model$bootvarcov,paste(model$apollo_control$modelName,"_bootcovar.csv",sep=""))
      cat("Bootstrap covariance matrix saved to",paste(model$apollo_control$modelName, "_bootcovar.csv"   , sep=""),"\n")
    }
  }
  if(saveCorr){
    if(printClassical==TRUE){ 
      utils::write.csv(model$corrmat,paste(model$apollo_control$modelName,"_corr.csv",sep=""))
      cat("Classical correlation matrix saved to",paste(model$apollo_control$modelName, "_covar.csv"   , sep=""),"\n")
    }
    utils::write.csv(model$robcorrmat,paste(model$apollo_control$modelName,"_robcorr.csv",sep=""))
    cat("Robust correlation matrix saved to",paste(model$apollo_control$modelName, "_robcorr.csv"   , sep=""),"\n")
    if(!is.null(model$bootstrapSE) && model$bootstrapSE>0){
      utils::write.csv(model$bootcorrmat, paste(model$apollo_control$modelName,"_bootcorr.csv",sep=""))
      cat("Bootstrap correlation matrix saved to",paste(model$apollo_control$modelName, "_bootcorr.csv"   , sep=""),"\n")
    }
  }
  
  if(saveModelObject){
    tryCatch( {
      saveRDS(model, file=paste0(model$apollo_control$modelName,"_model.rds"))
      cat("Model object saved to",paste(model$apollo_control$modelName, ".rds", sep=""),"\n")
      }, error=function(e) cat("Model object could not be written to file."))
  }
  
  if(writeF12) apollo_writeF12(model)
}
