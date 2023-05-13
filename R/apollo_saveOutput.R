#' Saves estimation results to files.
#' 
#' Writes files in the working/output directory with the estimation results.
#' 
#' Estimation results are saved different files in the working/output directory:
#' \itemize{
#'   \item \strong{\code{(modelName)_corr.csv}} CSV file with the estimated classical correlation matrix. Only when bayesian estimation was not used.
#'   \item \strong{\code{(modelName)_covar.csv}} CSV file with the estimated classical covariance matrix. Only when bayesian estimation was not used.
#'   \item \strong{\code{(modelName)_estimates.csv}} CSV file with the estimated parameter values, their standars errors, and t-ratios.
#'   \item \strong{\code{(modelName).F12}} F12 file with model results. Compatible with ALOGIT.
#'   \item \strong{\code{(modelName)_output.txt}} Text file with the output produced by function \code{apollo_modelOutput}.
#'   \item \strong{\code{(modelName)_robcorr.csv}} CSV file with the estimated robust correlation matrix. Only when bayesian estimation was not used.
#'   \item \strong{\code{(modelName)_robcovar.csv}} CSV file with the estimated robust covariance matrix. Only when bayesian estimation was not used.
#' }
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param saveOutput_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                            \itemize{
#'                               \item \strong{\code{printChange}}: Boolean. TRUE for printing difference between starting values and estimates. TRUE by default.
#'                               \item \strong{\code{printClassical}}: Boolean. TRUE for printing classical standard errors. TRUE by default.
#'                               \item \strong{\code{printCorr}}: Boolean. TRUE for printing parameters correlation matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. For Bayesian estimation, this setting is used for the covariane of random parameters. TRUE by default.
#'                               \item \strong{\code{printCovar}}: Boolean. TRUE for printing parameters covariance matrix. If \code{printClassical=TRUE}, both classical and robust matrices are printed. For Bayesian estimation, this setting is used for the correlation of random parameters. TRUE by default.
#'                               \item \strong{\code{printDataReport}}: Boolean. TRUE for printing summary of choices in database and other diagnostics. FALSE by default.
#'                               \item \strong{\code{printFixed}}: Logical. TRUE for printing fixed parameters among estimated parameter. TRUE by default.
#'                               \item \strong{\code{printFunctions}}: Boolean. TRUE for printing apollo_control, apollo_randCoeff (when available), apollo_lcPars (when available) and apollo_probabilities. TRUE by default.                               
#'                               \item \strong{\code{printHBconvergence}}: Boolean. TRUE for printing Geweke convergence tests. TRUE by default.                               
#'                               \item \strong{\code{printHBiterations}}: Boolean. TRUE for printing an iterations report for HB estimation. TRUE by default.                               
#'                               \item \strong{\code{printModelStructure}}: Boolean. TRUE for printing model structure. TRUE by default.
#'                               \item \strong{\code{printOutliers}}: Boolean or Scalar. TRUE for printing 20 individuals with worst average fit across observations. FALSE by default. If Scalar is given, this replaces the default of 20.
#'                               \item \strong{\code{printPVal}}: Boolean or Scalar. TRUE or 1 for printing p-values for one-sided test, 2 for printing p-values for two-sided test, FALSE for not printing p-values. FALSE by default.
#'                               \item \strong{\code{printT1}}: Boolean. If TRUE, t-test for H0: apollo_beta=1 are printed. FALSE by default.
#'                               \item \strong{\code{saveEst}}: Boolean. TRUE for saving estimated parameters and standard errors to a CSV file. TRUE by default.
#'                               \item \strong{\code{saveCorr}}: Boolean. TRUE for saving estimated correlation matrix to a CSV file. FALSE by default.
#'                               \item \strong{\code{saveCov}}: Boolean. TRUE for saving estimated covariance matrix to a CSV file. FALSE by default.
#'                               \item \strong{\code{saveHBiterations}}: Boolean. TRUE for including HB iterations in the saved model object. FALSE by default.
#'                               \item \strong{\code{saveModelObject}}: Boolean. TRUE to save the R model object to a file (use \link{apollo_loadModel} to load it to memory). TRUE by default.
#'                               \item \strong{\code{writeF12}}: Boolean. TRUE for writing results into an F12 file (ALOGIT format). FALSE by default.
#'                            }
#' @return nothing
#' @export
#' @importFrom RSGHB writeModel
#' @importFrom utils capture.output
apollo_saveOutput=function(model, saveOutput_settings=NA){
  if(length(saveOutput_settings)==1 && is.na(saveOutput_settings)) saveOutput_settings=list()
  default <- list(printClassical   = TRUE,
                  printChange      = TRUE,
                  printCorr        = TRUE,
                  printCovar       = TRUE,
                  printDataReport  = TRUE,
                  printFixed       = TRUE,
                  printFunctions   = TRUE,
                  printHBconvergence = TRUE,
                  printHBiterations = TRUE,
                  printModelStructure = TRUE, 
                  printOutliers    = TRUE,
                  printPVal        = FALSE,
                  printT1          = FALSE,
                  saveEst          = TRUE,
                  saveCov          = FALSE, #TRUE,
                  saveCorr         = FALSE, #TRUE,
                  saveHBiterations = FALSE,
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
  printFunctions  = saveOutput_settings[["printFunctions"]]
  saveEst         = saveOutput_settings[["saveEst"]]
  saveCov         = saveOutput_settings[["saveCov"]]
  saveCorr        = saveOutput_settings[["saveCorr"]]
  saveHBiterations= saveOutput_settings[["saveHBiterations"]]
  saveModelObject = saveOutput_settings[["saveModelObject"]]
  writeF12        = saveOutput_settings[["writeF12"]]
  scaling_used    = length(model$scaling)>0 && !anyNA(model$scaling)
  
  # Check if files exists, and if they do, rename them as OLD
  modName <- paste0(model$apollo_control$outputDirectory,model$apollo_control$modelName)
  if(file.exists( paste0(modName, "_output.txt") )){
    # Figure out corresponding OLD version
    n <- 1
    while( file.exists( paste0(modName, "_OLD", n, "_output.txt") ) ) n <- n + 1
    modNameOld <- paste0(modName, "_OLD", n)
    # Rename files
    outFiles <- c("_output.txt", "_estimates.csv", 
                  "_covar.csv", "_robcovar.csv", "_bootcovar.csv", 
                  "_corr.csv", "_robcorr.csv", "_bootcorr.csv",
                  "_model.rds", ".F12", 
                  "_F.csv", "_A.csv", "_B.csv", "_Bsd.csv", "_C.csv", "_Csd.csv", "_D.csv", 
                  ".log", "_params_non_random.csv", "_params_random_mean.csv", "_params_random_cov_mean.csv",
                  "_params_random_cov_sd.csv", "_params_posterior.csv", "_HB_random_params_covar.csv", "_HB_random_params_corr.csv",
                  "_params_non_random.csv",
                  "_posterior_means_summary.csv",
                  "_posterior_means.csv",
                  "_posterior_sd.csv",
                  "_params_random.csv",
                  "_params_random_covar.csv",
                  "_params_random_corr.csv"
                  )
    for(i in outFiles) if(file.exists(paste0(modName, i))){
      file.rename(from=paste0(modName, i), to=paste0(modNameOld, i))
      cat("\nOld result file \"", paste0(modName, i), "\" \n renamed to: \"", paste0(modNameOld, i), "\"", sep="")
    }
    cat("\n")
  }
  
  
  #sink( paste0(modName, "_output.txt") )
  #tmp <- apollo_modelOutput(model,saveOutput_settings)
  capture.output(tmp <- apollo_modelOutput(model,saveOutput_settings),
                 file=paste0(modName, "_output.txt"))
  cat("Model output saved to",paste0(modName, "_output.txt"),"\n")
  #sink()
  
  # ################################## #
  #### HB only components           ####
  # ################################## #
  
  if(model$apollo_control$HB){
    
    if(!saveHBiterations){
      model$HB_iterations_means=NULL
      model$HB_iterations_covar=NULL
      model$HB_iterations_non_random=NULL
    }

    if(saveEst){
      if(!is.null(model$HB_chains_non_random)){
        utils::write.csv(model$HB_chains_non_random, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_non_random"   ,".csv", sep=""))
        cat("Summary of chains for non-random parameters saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_non_random"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}

      if(!is.null(model$HB_posterior_means_summary)){
        utils::write.csv(model$HB_posterior_means_summary, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_means_summary"   ,".csv", sep=""))
        cat("Summary of distribution of posterior means for random parameters saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_means_summary"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
      
      if(!is.null(model$HB_posterior_means)){
        utils::write.csv(model$HB_posterior_means, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_means"   ,".csv", sep=""))
        cat("Individual-level posterior means for random parameters saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_means"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
      
      if(!is.null(model$HB_posterior_sd)){
        utils::write.csv(model$HB_posterior_sd, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_sd"   ,".csv", sep=""))
        cat("Individual-level posterior standard deviation for random parameters saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_posterior_means"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}

      if(!is.null(model$HB_random_params_mean_sd)){
      utils::write.csv(model$HB_random_params_mean_sd, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random"   ,".csv", sep=""))
      cat("Summary for randomly distributed parameters (with distributional transforms) saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random"   ,".csv", sep=""),"\n")
      if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    }
    
    if(saveCov){
      if(!is.null(model$HB_random_params_covar)){      
        utils::write.csv(model$HB_random_params_covar    , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_covar"    ,".csv", sep=""))
        cat("Covariance matrix of random coeffients (after distributional transforms) saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_covar"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    }
    
    
    if(saveCorr){
      if(!is.null(model$HB_random_params_corr)){      
        utils::write.csv(model$HB_random_params_corr    , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_corr"    ,".csv", sep=""))
        cat("Correlation matrix of random coeffients (after distributional transforms) saved to:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_corr"   ,".csv", sep=""),"\n")
        if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    }
    
    # if(saveEst){
    #   
    #   
    #   
    #   currentWD <- getwd()
    #   if(dir.exists(model$apollo_control$outputDirectory)) setwd(model$apollo_control$outputDirectory)
    #   RSGHB::writeModel(model, writeDraws = FALSE, path = getwd())
    #   setwd(currentWD)
    #   cat("\n\nRSGHB output saved in following files\n")
    #   cat("\nOutputs at iteration level (post burn-in chains)\n")
    #   if(!is.null(tmp$non_random)){ cat("Non-random parameters:",paste(model$apollo_control$modelName, "_F"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("These outputs have had the scaling used in estimation applied to them\n")}
    # 
    #   if(!is.null(tmp$random_mean)){ cat("Means for underlying normals:",paste(model$apollo_control$modelName, "_A"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    # 
    #   cat("\nPosteriors\n")
    #   if(!is.null(tmp$random_mean)){ cat("Mean individual-level draws for underlying normals:",paste(model$apollo_control$modelName, "_B"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   
    #   if(!is.null(tmp$random_mean)){ cat("SD of individual-level draws for underlying normals:",paste(model$apollo_control$modelName, "_Bsd"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   
    #   if(!is.null(tmp$random_mean)){ cat("Mean individual-level draws after transformations to underlying normals:",paste(model$apollo_control$modelName, "_C"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   
    #   if(!is.null(tmp$random_mean)){ cat("SD of individual-level draws after transformations to underlying normals:",paste(model$apollo_control$modelName, "_Csd"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   
    #   if(!is.null(tmp$random_mean)){ cat("Sample variance-covariance matrix for underlying normals:",paste(model$apollo_control$modelName, "_D"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   
    #   cat("\nRSGHB log file saved to",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, ".log", sep=""),"\n")
    #   
    #   cat("\nAdditional output files:\n")
    #   if(!is.null(tmp$non_random)){
    #   utils::write.csv(tmp$non_random   , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_non_random"   ,".csv", sep=""))
    #   cat("Summary of chains for non-random parameters:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_non_random"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   if(!is.null(tmp$random_mean)){utils::write.csv(tmp$random_mean  , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_mean"  ,".csv", sep=""))
    #   
    #   cat("Summary of chains for means of normals:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_mean"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   if(!is.null(tmp$random_cov_mean)){utils::write.csv(tmp$random_cov_mean, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_cov_mean",".csv", sep=""))
    #   
    #   cat("Means of chains for covariance of normals:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_cov_mean"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   if(!is.null(tmp$random_cov_sd)){utils::write.csv(tmp$random_cov_sd, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_cov_sd",".csv", sep=""))
    #   
    #   cat("SDs of chains for covariance of normals:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_random_cov_sd"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have NOT had the scaling used in estimation applied to them\n")}
    #   if(!is.null(tmp$posterior)){utils::write.csv(tmp$posterior    , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_posterior"    ,".csv", sep=""))
    #   
    #   cat("Summary of posteriors for random parameters:",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_params_posterior"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   
    # 
    #   if(!is.null(tmp$HB_random_params_covar)){utils::write.csv(tmp$HB_random_params_covar, paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_HB_random_params_covar",".csv", sep=""))
    #   
    #   cat("Covariance matrix of random coeffients (after distributional transforms):",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_HB_random_params_covar"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   if(!is.null(tmp$HB_random_params_corr)){utils::write.csv(tmp$HB_random_params_corr    , paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_HB_random_params_corr"    ,".csv", sep=""))
    #   
    #   cat("Correlation matrix of random coeffients (after distributional transforms):",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, "_HB_random_params_corr"   ,".csv", sep=""),"\n")
    #   if(scaling_used) cat("   These outputs have had the scaling used in estimation applied to them\n")}
    #   
    #  }
    if(saveModelObject){
      tryCatch( {
        #model$cmcLLout=NULL
        #model$cmcRLHout=NULL
        saveRDS(model, file=paste0(model$apollo_control$outputDirectory,model$apollo_control$modelName,"_model.rds"))
        cat("\nModel object saved to",paste(model$apollo_control$outputDirectory,model$apollo_control$modelName, ".rds", sep=""),"\n")
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
      if(printPVal==1) output <- cbind( output, `p(1-sided)`=printPVal*(1-stats::pnorm(abs(model$estimate/model$se))) )
      if(printPVal==2) output <- cbind( output, `p(2-sided)`=printPVal*(1-stats::pnorm(abs(model$estimate/model$se))) )
      if(printT1){
        output <- cbind(output, `t-ratio(1)`=(model$estimate-1)/model$se)
        if(printPVal==1) output <- cbind( output, `p(1-sided)`=printPVal*(1-stats::pnorm(abs((model$estimate-1)/model$se))) )
        if(printPVal==2) output <- cbind( output, `p(2-sided)`=printPVal*(1-stats::pnorm(abs((model$estimate-1)/model$se))) )
      }
    }
    output <- cbind(output, Rob.std.err.=model$robse, `Rob.t-ratio(0)`=model$estimate/model$robse)
    if(printPVal) output <- cbind(output, `Rob.p-val(0)`=2*(1-stats::pnorm(abs(model$estimate/model$robse))) )
    if(printT1){
      output <- cbind(output, `Rob.t-ratio(1)`=(model$estimate-1)/model$robse)
      if(printPVal) output <- cbind(output, `Rob.p-val(1)`=2*(1-stats::pnorm(abs((model$estimate-1)/model$robse))) )
    }
    
    # Write to file
    utils::write.csv(output,paste0(modName,"_estimates.csv"))
    cat("Estimates saved to",paste0(modName, "_estimates.csv"),"\n")
  }
  if(saveCov){
    if(printClassical==TRUE){
      utils::write.csv(model$varcov,paste0(modName,"_covar.csv"))
      cat("Classical covariance matrix saved to",paste0(modName, "_covar.csv"),"\n")
      }
    utils::write.csv(model$robvarcov,paste0(modName,"_robcovar.csv"))
    cat("Robust covariance matrix saved to",paste0(modName, "_robcovar.csv"),"\n")
    if(!is.null(model$bootstrapSE) && model$bootstrapSE>0){
      utils::write.csv(model$bootvarcov,paste0(modName,"_bootcovar.csv"))
      cat("Bootstrap covariance matrix saved to",paste(modName, "_bootcovar.csv"   , sep=""),"\n")
    }
  }
  if(saveCorr){
    if(printClassical==TRUE){ 
      utils::write.csv(model$corrmat,paste0(modName,"_corr.csv"))
      cat("Classical correlation matrix saved to",paste(modName, "_covar.csv"   , sep=""),"\n")
    }
    utils::write.csv(model$robcorrmat,paste0(modName,"_robcorr.csv"))
    cat("Robust correlation matrix saved to",paste0(modName, "_robcorr.csv"),"\n")
    if(!is.null(model$bootstrapSE) && model$bootstrapSE>0){
      utils::write.csv(model$bootcorrmat, paste0(modName,"_bootcorr.csv"))
      cat("Bootstrap correlation matrix saved to",paste(modName, "_bootcorr.csv"   , sep=""),"\n")
    }
  }
  
  if(saveModelObject){
    tryCatch( {
      saveRDS(model, file=paste0(modName,"_model.rds"))
      cat("Model object saved to",paste0(modName, ".rds"),"\n")
      }, error=function(e) cat("Model object could not be written to file."))
  }
  
  if(writeF12) apollo_writeF12(model)
}
