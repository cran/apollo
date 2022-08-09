#' Checks and modifies Apollo user-defined functions
#'
#' Checks and enhances user defined functions apollo_probabilities, apollo_randCoeff and apollo_lcPars.
#'
#' Internal use only. Called by \code{apollo_estimate} before estimation.
#' Checks include: no re-definition of variables, no (direct) calls to database, 
#' calling of apollo_weighting if weights are defined.
#' 
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_probabilities Function. Returns probabilities of the model to 
#'                             be estimated. Must receive three arguments:
#'                             \itemize{
#'                              \item \strong{\code{apollo_beta}}: Named numeric vector. 
#'                                                                 Names and values of 
#'                                                                 model parameters.
#'                              \item \strong{\code{apollo_inputs}}: List containing 
#'                                                                   options of the model. 
#'                                                                   See \link{apollo_validateInputs}.
#'                              \item \strong{\code{functionality}}: Character. Can be 
#'                                                                   either \code{"components"}, 
#'                                                                   \code{"conditionals"}, 
#'                                                                   \code{"estimate"} (default), 
#'                                                                   \code{"gradient"}, 
#'                                                                   \code{"output"}, 
#'                                                                   \code{"prediction"}, 
#'                                                                   \code{"preprocess"}, 
#'                                                                   \code{"raw"}, \code{"report"}, 
#'                                                                   \code{"shares_LL"}, 
#'                                                                   \code{"validate"} or 
#'                                                                   \code{"zero_LL"}.
#'                             }
#' @param apollo_fixed Character vector. Names of parameters inside apollo_beta
#'                     whose values should be kept constant throughout estimation.
#' @param apollo_inputs List grouping most common inputs. Created by function 
#'                      \link{apollo_validateInputs}.
#' @param validate Logical. If TRUE, the original and modified 
#'                 \code{apollo_probabilities} functions are estimated. If their 
#'                 results do not match, then the original functions are 
#'                 returned, and \code{success} is set to \code{FALSE} inside 
#'                 the returned list.
#' @return List with four elements: apollo_probabilities, apollo_randCoeff, 
#'         apollo_lcPars and a dummy called success (TRUE if modification was
#'         successful, FALSE if not. FALSE will be only be returnes if
#'         the modifications are validated).
#' @export
apollo_modifyUserDefFunc <- function(apollo_beta, apollo_fixed, 
                                     apollo_probabilities, apollo_inputs, 
                                     validate=TRUE){
  silent <- apollo_inputs$silent
  debug  <- apollo_inputs$apollo_control$debug
  
  # # # #  # # # # 
  #### Checks ####
  # # # #  # # # # 
  
  ### Check for apollo_prepareProb() and return() in apollo_probabilities
  tmp  <- deparse(apollo_probabilities)
  if(!any(grepl("apollo_prepareProb", tmp))) stop("Your 'apollo_probabilities' function should include a call to 'apollo_prepareProb'!")
  if(!any(grepl("return", tmp))) stop("Your 'apollo_probabilities' function should include a 'return' statement at the end, usually 'return(P)'!")
  
  ### Check that names of params in apollo_beta, database, apollo_randCoeff & apollo_lcPars are not re-defined
  tmp <- as.character(body(apollo_probabilities))
  tmp <- gsub("(", "", tmp, fixed=TRUE)
  tmp <- gsub(")", "", tmp, fixed=TRUE)
  tmp <- gsub(" ", "", tmp, fixed=TRUE)
  # check for apollo_beta
  for(i in names(apollo_beta)){
    test <- grep(paste0("^",i,"="), tmp)
    test <- c(test, grep(paste0("^",i,"<-"), tmp))
    if(length(test)>0) stop("Parameter ", i, " from apollo_beta was re-defined ",
                            "inside apollo_probabilities. This is not allowed.")
  }
  for(i in names(apollo_inputs$database)){
    test <- grep(paste0("^",i,"(=|<-)"), tmp)
    if(length(test)>0) stop("Variable ", i, " from database is re-defined ", 
                            "inside apollo_probabilities. This is not allowed")
  }; rm(i, test)
  #check for apollo_randCoeff
  if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
    env <- c(apollo_inputs$database, apollo_inputs$draws, as.list(apollo_beta))
    env <- list2env(env, hash=TRUE, parent=parent.frame())
    rnd <- apollo_inputs$apollo_randCoeff; environment(rnd) <- env
    rnd <- rnd(apollo_beta, apollo_inputs)
    for(i in names(rnd)){
      test <- grep(paste0('^', i, '=|^', i, '<-'), tmp)
      if(length(test)>0) stop("Parameter ", i, " from apollo_randCoeff was re-defined ",
                              "inside apollo_probabilities. This is not allowed.")
    }; rm(env, i, test)
  }
  #check for apollo_lcPars
  if(is.function(apollo_inputs$apollo_lcPars)){
    env <- c(apollo_inputs$database, as.list(apollo_beta))
    if(exists('rnd', inherits=FALSE)) env <- c(env, apollo_inputs$draws, rnd)
    env <- list2env(env, hash=TRUE, parent=parent.frame())
    lcp <- apollo_inputs$apollo_lcPars; environment(lcp) <- env
    lcp <- names(lcp(apollo_beta, apollo_inputs))
    for(i in lcp){
      test <- grep(paste0('^', i, '=|^', i, '<-'), tmp)
      if(length(test)>0) stop("Parameter ", i, " from apollo_lcPars was re-defined ",
                              "inside apollo_probabilities. This is not allowed.")
    }; rm(env, lcp, i, test)
  }; if(exists('rnd')) rm(rnd)
  rm(tmp)
  
  ### Check there are no references to database inside apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- gsub("apollo_inputs$database", " ", tmp, fixed=TRUE)
    tmp <- grep("database", tmp, fixed=TRUE)
    if(length(tmp)>0) stop("The database object is 'attached' and elements should thus be called",
                           " directly in apollo_probabilities without the 'database$' prefix.")
    rm(tmp)
  }
  
  ### Check apollo_weighting is called if apollo_control$weights are defined (unless apollo_inputs$EM is TRUE)
  w <- apollo_inputs$apollo_control[['weights']]
  test <- is.null(apollo_inputs$EM) || (is.logical(apollo_inputs$EM) && !apollo_inputs$EM)
  test <- test && !is.null(w) && !is.null(apollo_inputs$database) && (w %in% names(apollo_inputs$database))
  if(test){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- grep('apollo_weighting', tmp, fixed=TRUE)
    if(length(tmp)==0) stop('When using weights, apollo_weighting should be called inside apollo_probabilities.')
    rm(tmp)
  }; rm(w)
  
  
  # # # # # # # # # # # # # # # # # # #
  #### Validate and prepare scaling ####
  # # # # # # # # # # # # # # # # # # #
  scaling <- setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                      names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  test <- is.vector(apollo_inputs$apollo_scaling) && is.numeric(apollo_inputs$apollo_scaling)
  test <- test && !is.null(names(apollo_inputs$apollo_scaling)) && !anyNA(apollo_inputs$apollo_scaling)
  if(test){
    # If the user did not provide any scaling
    # (apollo_inputs$apollo_scaling is NULL, NA, or doesn't have names)
    apollo_inputs$apollo_scaling <- scaling
    apollo_inputs$manualScaling  <- FALSE
  } else {
    # If the user did provide scaling
    # Check that all scaling correspond to existing parameters
    if(!all(names(apollo_inputs$apollo_scaling) %in% names(apollo_beta))){
      txt <- names(apollo_inputs$apollo_scaling)[!(names(apollo_inputs$apollo_scaling) %in% names(apollo_beta))]
      stop(paste0("Some parameters included in 'scaling' (", paste0(txt, collapse=", "), 
                  ") are not included in 'apollo_beta'."))
    }
    # Check for duplicates
    if(anyDuplicated(names(apollo_inputs$apollo_scaling))){
      txt <- names(apollo_inputs$apollo_scaling)[duplicated(names(apollo_inputs$apollo_scaling))]
      stop(paste0("The \"scaling\" setting contains duplicate elements (", paste0(txt, collapse=", "), ")."))
    }
    # Copy user provided scales into "scaling"
    scaling[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling
    # Check no fixed params are scaled
    if(any(names(scaling) %in% apollo_fixed)) stop("Parameters in 'apollo_fixed' should not be included in 'scaling'")
    # Check there are no negative scaling values. If there are, take their abs value.
    if(any(scaling<0)){
      scaling <- abs(scaling)
      txt <- 'Some negative values in "scaling" were replaced by their absolute value'
      if(!silent) apollo_print(paste0('WARNING: ', txt, '.')) else warning(txt)
    }; if(any(scaling<=0)) stop('All terms in "scaling" should be strictly positive!')
    txt <- "During estimation, parameters will be scaled using the values in estimate_settings$scaling"
    if(!all(scaling==1)){ if(!silent) apollo_print(txt) else warning(txt)}
    rm(txt)
    apollo_inputs$apollo_scaling <- scaling
    apollo_inputs$manualScaling  <- TRUE
  }
  rm(scaling)
  
  
  # # # # # # # # # # # #
  #### Modifications ####
  # # # # # # # # # # # #
  
  if(!silent) apollo_print("Preparing user-defined functions.")
  
  ### Evaluate apollo_probabilities before changes, return immediately if it doesn't work
  if(validate){
    apollo_beta_shifted <- apollo_beta + 0.0001
    test1 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), 
                      error=function(e) if(grepl("not found",e$message)) stop(e) else NULL )
    if(anyNA(test1)) test1 <- NULL
    if(is.null(test1)) return(list(apollo_probabilities = apollo_probabilities, 
                                   apollo_randCoeff     = apollo_inputs$apollo_randCoeff, 
                                   apollo_lcPars        = apollo_inputs$apollo_lcPars,
                                   apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                                                                   names(apollo_inputs$apollo_scaling)), 
                                   success              = FALSE))
  }
  
  ### Insert componentName if missing
  if(debug) apollo_print("- Inserting component name in apollo_probabilities")
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  
  ### Store unaltered version of functions in case modification fails
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_inputs$apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_inputs$apollo_lcPars
  
  ### Expand loop
  if(is.function(apollo_inputs$apollo_randCoeff) && apollo_inputs$apollo_control$mixing){
    if(debug) apollo_print("- Expanding loops in apollo_randCoeff.")
    tmp <- tryCatch( apollo_expandLoop(apollo_inputs$apollo_randCoeff, apollo_inputs), error=function(e) NULL)
    if(is.function(tmp)) apollo_inputs$apollo_randCoeff <- tmp
    rm(tmp)
  }
  if(is.function(apollo_inputs$apollo_lcPars)){
    if(debug) apollo_print("- Expanding loops in apollo_lcPars.")
    tmp <- tryCatch( apollo_expandLoop(apollo_inputs$apollo_lcPars, apollo_inputs), error=function(e) NULL)
    if(is.function(tmp)) apollo_inputs$apollo_lcPars <- tmp
    rm(tmp)
  }
  if(debug) apollo_print("- Expanding loops in apollo_probabilities.")
  tmp <- tryCatch( apollo_expandLoop(apollo_probabilities, apollo_inputs), error=function(e) NULL)
  if(is.function(tmp)) apollo_probabilities <- tmp
  rm(tmp)
  
  ### Insert scaling (only if no apollo_inputs$apollo_scaling is found inside)
  test <- as.character(body(apollo_probabilities))
  test <- grepl('apollo_inputs$apollo_scaling', test, fixed=TRUE)
  if(any(test)){
    if(debug) apollo_print("- Scaling not inserted because it is already present in apollo_probabilities.")
  } else {
    if(debug) apollo_print("- Inserting scaling in apollo_probabilities")
    apollo_probabilities <- apollo_insertScaling(apollo_probabilities, apollo_inputs$apollo_scaling)
    if(apollo_inputs$apollo_control$mixing && is.function(apollo_inputs$apollo_randCoeff)){
      if(debug) apollo_print("- Inserting scaling in apollo_randCoeff")
      apollo_inputs$apollo_randCoeff <- apollo_insertScaling(apollo_inputs$apollo_randCoeff, 
                                                             apollo_inputs$apollo_scaling)
    }
    if(is.function(apollo_inputs$apollo_lcPars)){
      if(debug) apollo_print("- Inserting scaling in apollo_lcPars")
      apollo_inputs$apollo_lcPars <- apollo_insertScaling(apollo_inputs$apollo_lcPars, apollo_inputs$apollo_scaling)
    }
  }
  
  ### Introduce quotes into apollo_rrm
  if(debug) apollo_print("- Inserting quotes in settings for apollo_rrm (if present)")
  apollo_probabilities <- apollo_insertRRMQuotes(apollo_probabilities)
  
  ### Introduce 'function()' at the beginning of definitions (only if using analytic gradients)
  if(debug) apollo_print('- Inserting function() in user-defined functions')
  if(apollo_inputs$apollo_control$analyticGrad){
    tmp <- apollo_insertFunc(apollo_probabilities, like=TRUE)
    if(is.function(tmp)) apollo_probabilities <- tmp
    if(is.function(apollo_inputs$apollo_randCoeff)){
      tmp <- apollo_insertFunc(apollo_inputs$apollo_randCoeff, randCoeff=TRUE)
      if(is.function(tmp))apollo_inputs$apollo_randCoeff <- tmp
    } 
    if(is.function(apollo_inputs$apollo_lcPars)){
      tmp <- apollo_insertFunc(apollo_inputs$apollo_lcPars, lcPars=TRUE)
      if(is.function(tmp)) apollo_inputs$apollo_lcPars <- tmp
    }
  }
  
  ### Evaluate apollo_probabilities after changes and compare to result before them
  if(validate){
    apollo_beta_shifted[names(apollo_inputs$apollo_scaling)] <- apollo_beta_shifted[names(apollo_inputs$apollo_scaling)]/apollo_inputs$apollo_scaling
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The output of 'apollo_probabilities' has changed", 
                   " after the optimisation of the code carried out by Apollo.", 
                   " This indicates a problem, and the unoptimised and unscaled", 
                   " version will be used instead. Please contact the", 
                   " developers for assistance!"), highlight=TRUE)
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                                                  names(apollo_inputs$apollo_scaling)), 
                  success              = FALSE))
    }
  }
  
  # If functions were modified and in debug mode, then write them to file
  if(debug){
    test <- !is.null(apollo_inputs) && is.list(apollo_inputs) && !is.null(apollo_inputs$apollo_control)
    test <- test && is.list(apollo_inputs$apollo_control) && !is.null(apollo_inputs$apollo_control$outputDirectory)
    test <- test && is.character(apollo_inputs$apollo_control$outputDirectory)
    test <- test && length(apollo_inputs$apollo_control$outputDirectory)==1
    if(test) outputDirectory <- apollo_inputs$apollo_control$outputDirectory else outputDirectory=getwd()
    if(substr(outputDirectory, nchar(outputDirectory), nchar(outputDirectory))!="/") outputDirectory <- paste0(outputDirectory, "/")
    txt <- utils::capture.output(print(apollo_probabilities))
    fileConn <- file(paste0(outputDirectory, apollo_inputs$apollo_control$modelName, "_apollo_probabilities_modified.txt"))
    writeLines(txt, fileConn)
    close(fileConn)
    if(is.function(apollo_inputs$apollo_randCoeff)){
      txt <- utils::capture.output(print(apollo_inputs$apollo_randCoeff))
      fileConn <- file(paste0(outputDirectory, apollo_inputs$apollo_control$modelName, "_apollo_randCoeff_modified.txt"))
      writeLines(txt, fileConn)
      close(fileConn)
    }
    if(is.function(apollo_inputs$apollo_lcPars)){
      txt <- utils::capture.output(print(apollo_inputs$apollo_lcPars))
      fileConn <- file(paste0(outputDirectory, apollo_inputs$apollo_control$modelName, "_apollo_lcPars_modified.txt"))
      writeLines(txt, fileConn)
      close(fileConn)
    }
    rm(test, txt)
  }
  
  ### Return
  return(list(apollo_probabilities = apollo_probabilities, 
              apollo_randCoeff     = apollo_inputs$apollo_randCoeff, 
              apollo_lcPars        = apollo_inputs$apollo_lcPars, 
              apollo_scaling       = apollo_inputs$apollo_scaling, 
              success              = TRUE))
}