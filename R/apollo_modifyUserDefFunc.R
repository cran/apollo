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
#' @param noModification Logical. If TRUE, loop expansion etc are skipped.
#' @return List with four elements: apollo_probabilities, apollo_randCoeff, 
#'         apollo_lcPars and a dummy called success (TRUE if modification was
#'         successful, FALSE if not. FALSE will be only be returnes if
#'         the modifications are validated).
#' @export
apollo_modifyUserDefFunc <- function(apollo_beta, apollo_fixed, 
                                     apollo_probabilities, apollo_inputs, 
                                     validate=TRUE, noModification=FALSE){
  silent <- apollo_inputs$silent
  debug  <- apollo_inputs$apollo_control$debug
  
  # # # #  # # # # 
  #### Checks ####
  # # # #  # # # # 
  
  ### Check for apollo_attach in apollo_probabilities
  tmp  <- deparse(apollo_probabilities)
  usesAttach <- FALSE
  if(any(grepl("apollo_attach", tmp))) usesAttach <- TRUE
  if(!usesAttach && !silent) apollo_print("You are not using apollo_attach, this may affect performance and capabilities", type="w", pause=5)
  
  ### Check for apollo_prepareProb() and return() in apollo_probabilities
  if(!any(grepl("apollo_prepareProb", tmp))) stop("SYNTAX ISSUE - The 'apollo_probabilities' function should include a call to 'apollo_prepareProb'!")
  if(!any(grepl("return", tmp))) stop("SYNTAX ISSUE - The 'apollo_probabilities' function should include a 'return' statement at the end, usually 'return(P)'!")
  
  ### Check that apollo_beta is not called, unless not using attach
  tmp <- as.character(body(apollo_probabilities))
  if(usesAttach && any(grepl("apollo_beta[", tmp, fixed=TRUE))){
    stop("SYNTAX ISSUE - The apollo_beta object is 'attached' and elements should thus be called",
         " directly in apollo_probabilities without the 'apollo_beta[...]' syntax.")
  }
  
  ### Check that names of params in apollo_beta, database, apollo_randCoeff & apollo_lcPars are not re-defined
  tmp <- gsub("(", "", tmp, fixed=TRUE)
  tmp <- gsub(")", "", tmp, fixed=TRUE)
  tmp <- gsub(" ", "", tmp, fixed=TRUE)
  # check for apollo_beta
  for(i in names(apollo_beta)){
    test <- grep(paste0("^",i,"="), tmp)
    test <- c(test, grep(paste0("^",i,"<-"), tmp))
    if(length(test)>0) stop("SYNTAX ISSUE - Parameter ", i, " from apollo_beta was re-defined ",
                            "inside apollo_probabilities. This is not allowed.")
  }
  for(i in names(apollo_inputs$database)){
    test <- grep(paste0("^",i,"(=|<-)"), tmp)
    if(length(test)>0) stop("SYNTAX ISSUE - Variable ", i, " from database is re-defined ", 
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
      if(length(test)>0) stop("SYNTAX ISSUE - Parameter ", i, " from apollo_randCoeff was re-defined ",
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
      if(length(test)>0) stop("SYNTAX ISSUE - Parameter ", i, " from apollo_lcPars was re-defined ",
                              "inside apollo_probabilities. This is not allowed.")
    }; rm(env, lcp, i, test)
  }; if(exists('rnd')) rm(rnd)
  rm(tmp)
  
  ### Check there are no references to database inside apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- gsub("apollo_inputs$database", " ", tmp, fixed=TRUE)
    tmp <- grep("database", tmp, fixed=TRUE)
    if(length(tmp)>0) stop("SYNTAX ISSUE - The database object is 'attached' and elements should thus be called",
                           " directly in apollo_probabilities without the 'database$' prefix. If there is a specific",
                           " reason for doing so, the  'apollo_inputs$database$' prefix has to be used.")
    rm(tmp)
  }
  
  ### Check apollo_weighting is called if apollo_control$weights are defined (unless apollo_inputs$EM is TRUE)
  w <- apollo_inputs$apollo_control[['weights']]
  test <- is.null(apollo_inputs$EM) || (is.logical(apollo_inputs$EM) && !apollo_inputs$EM)
  test <- test && !is.null(w) && !is.null(apollo_inputs$database) && (w %in% names(apollo_inputs$database))
  if(test){
    tmp <- as.character(body(apollo_probabilities))
    tmp <- grep('apollo_weighting', tmp, fixed=TRUE)
    if(length(tmp)==0) stop('SYNTAX ISSUE - When using weights, apollo_weighting should be called inside apollo_probabilities.')
    rm(tmp)
  }; rm(w)
  
  
  # # # # # # # # # # # # # # # # # # #
  #### Validate and prepare scaling ####
  # # # # # # # # # # # # # # # # # # #
  scaling <- setNames(rep(1, length(apollo_beta)-length(apollo_fixed)), 
                      names(apollo_beta)[!(names(apollo_beta) %in% apollo_fixed)])
  test <- is.vector(apollo_inputs$apollo_scaling) && is.numeric(apollo_inputs$apollo_scaling)
  test <- test && !is.null(names(apollo_inputs$apollo_scaling)) && !anyNA(apollo_inputs$apollo_scaling)
  if(!test){
    # If the user did not provide any scaling
    # (apollo_inputs$apollo_scaling is NULL, NA, or doesn't have names)
    apollo_inputs$apollo_scaling <- scaling
    apollo_inputs$manualScaling  <- FALSE
  } else {
    # If the user did provide scaling
    # Check that all scaling correspond to existing parameters
    if(!all(names(apollo_inputs$apollo_scaling) %in% names(apollo_beta))){
      txt <- names(apollo_inputs$apollo_scaling)[!(names(apollo_inputs$apollo_scaling) %in% names(apollo_beta))]
      stop(paste0("SYNTAX ISSUE - Some parameters included in 'scaling' (", paste0(txt, collapse=", "), 
                  ") are not included in 'apollo_beta'."))
    }
    # Check for duplicates
    if(anyDuplicated(names(apollo_inputs$apollo_scaling))){
      txt <- names(apollo_inputs$apollo_scaling)[duplicated(names(apollo_inputs$apollo_scaling))]
      stop(paste0("SYNTAX ISSUE - The \"scaling\" setting contains duplicate elements (", paste0(txt, collapse=", "), ")."))
    }
    # Copy user provided scales into "scaling"
    scaling[names(apollo_inputs$apollo_scaling)] <- apollo_inputs$apollo_scaling
    # Check no fixed params are scaled
    if(any(names(scaling) %in% apollo_fixed)) stop("SYNTAX ISSUE - Parameters in 'apollo_fixed' should not be included in 'scaling'")
    # Check there are no negative scaling values. If there are, take their abs value.
    if(any(scaling<0)){
      scaling <- abs(scaling)
      txt <- 'Some negative values in "scaling" were replaced by their absolute value'
      if(!silent) apollo_print(paste0(txt, '.'), type="w") else warning(txt)
    }; if(any(scaling<=0)) stop('SYNTAX ISSUE - All terms in "scaling" should be strictly positive!')
    txt <- "During estimation, parameters will be scaled using the values in estimate_settings$scaling"
    if(!all(scaling==1)){ if(!silent) apollo_print(paste0(txt, '.'), type="i") else warning(txt)}
    rm(txt)
    apollo_inputs$apollo_scaling <- scaling
    apollo_inputs$manualScaling  <- TRUE
  }
  rm(scaling)
  
  # # # # # # # # # # # # # # # # # # # # #
  #### Create testing values for beta ####
  # # # # # # # # # # # # # # # # # # # # #
  
  if(!apollo_inputs$apollo_control$HB){
    apollo_beta_shifted <- apollo_beta + 0.0001  
  }else{
    apollo_HB=apollo_inputs$apollo_HB
    apollo_test_beta=apollo_beta
    if(!is.null(apollo_HB$gVarNamesFixed)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesFixed )
      r <- names(apollo_beta)[r]
      apollo_test_beta[r] <- apollo_HB$FC[r]
    }
    if(!is.null(apollo_HB$gVarNamesNormal)){
      r <- ( names(apollo_beta) %in% apollo_HB$gVarNamesNormal )
      r <- names(apollo_beta)[r]
      dists_normal=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==1])
      dists_lnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==2])
      dists_lnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==3])
      dists_cnp=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==4])
      dists_cnn=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==5])
      dists_sb=names(apollo_HB$gDIST[r][apollo_HB$gDIST[r]==6])
      if(length(dists_normal)>0) apollo_test_beta[dists_normal] <- apollo_HB$svN[dists_normal]
      if(length(dists_lnp)>0) apollo_test_beta[dists_lnp] <- exp(apollo_HB$svN[dists_lnp])
      if(length(dists_lnn)>0) apollo_test_beta[dists_lnn] <- -exp(apollo_HB$svN[dists_lnn])
      if(length(dists_cnp)>0) apollo_test_beta[dists_cnp] <- apollo_HB$svN[dists_cnp]*(apollo_HB$svN[dists_cnp]>0)
      if(length(dists_cnn)>0) apollo_test_beta[dists_cnn] <- apollo_HB$svN[dists_cnn]*(apollo_HB$svN[dists_cnn]<0)
      if(length(dists_sb)>0){
        names(apollo_HB$gMINCOEF)=names(apollo_HB$svN)
        names(apollo_HB$gMAXCOEF)=names(apollo_HB$svN)
        apollo_test_beta[dists_sb] <- apollo_HB$gMINCOEF[dists_sb]+(apollo_HB$gMAXCOEF[dists_sb]-apollo_HB$gMINCOEF[dists_sb])/(1+exp(-apollo_HB$svN[dists_sb]))
      }
      rm(dists_normal, dists_lnp, dists_lnn, dists_cnp, dists_cnn, dists_sb)
    }
    apollo_beta_shifted <- apollo_test_beta + 0.0001  
    rm(apollo_test_beta)
  }

  # # # # # # # # # # # #
  #### Modifications ####
  # # # # # # # # # # # #
  
  if(!silent) apollo_print("Preparing user-defined functions.")
  
  ### Store unaltered version of functions in case modification fails
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_inputs$apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_inputs$apollo_lcPars
  
  ### Insert componentName if missing
  if(debug) apollo_print("- Inserting component name in apollo_probabilities")
  apollo_probabilities <- apollo_insertComponentName(apollo_probabilities)
  
  ### Evaluate apollo_probabilities before changes, return immediately if it doesn't work
  # NOTE: We should simplify this by calling apollo_validate BEFORE modifying the functions.
  #       But to do that, we would need to either not check for component names
  #       or insert the component names beforehand. Or call apollo_probabilities(..., "validate")
  #       right here.
  if(validate || noModification){
    #apollo_beta_shifted <- apollo_beta + 0.0001
    test1 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), 
                      error=function(e) if(grepl("not found|ISSUE|INCORRECT",e,fixed=FALSE)){
                        stop(e)
                      }else{
                        NULL
                      } )
    if(anyNA(test1)) test1 <- NULL
    test <- is.null(test1) || noModification
    if(test){
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed in initial testing.",
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities, # returns version with names inserted
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                               names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
                  success              = FALSE))
    } 
  }
  ### Update ORIG version to keep changes so far
  apollo_probabilities_ORIG <- apollo_probabilities
  
  ### Change c to list if using OL
  if(debug) apollo_print("- Replacing tau=c(...) by tau=list(...) in calls to apollo_ol.")
  test <- any(grepl("apollo_ol", as.character(body(apollo_probabilities))))
  if(test) apollo_probabilities <- apollo_insertOLList(apollo_probabilities)
  test <- is.function(apollo_inputs$apollo_randCoeff) && 
    apollo_inputs$apollo_control$mixing && 
    any(grepl("apollo_ol", as.character(body(apollo_inputs$apollo_randCoeff))))
  if(test) apollo_inputs$apollo_lcPars <- apollo_insertOLList(apollo_inputs$apollo_lcPars)
  test <- is.function(apollo_inputs$apollo_lcPars) &&
    any(grepl("apollo_ol", as.character(body(apollo_inputs$apollo_lcPars))))
  if(test) apollo_inputs$apollo_lcPars <- apollo_insertOLList(apollo_inputs$apollo_lcPars)
  ### Evaluate apollo_probabilities after current changes, return immediately if it doesn't work
  if(validate){
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed during syntax checking.", 
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                                names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
                  success              = FALSE))
    }
  }
  ### Update ORIG version to keep changes so far
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
  
  ### Evaluate apollo_probabilities after current changes, return immediately if it doesn't work
  if(validate){
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed during loop expansion.", 
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                                names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
                  success              = FALSE))
    }
  }
  ### Update ORIG version to keep changes so far
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_inputs$apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_inputs$apollo_lcPars
  
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
  
  ### Evaluate apollo_probabilities after current changes, return immediately if it doesn't work
  if(validate){
    apollo_beta_shifted[names(apollo_inputs$apollo_scaling)] <- apollo_beta_shifted[names(apollo_inputs$apollo_scaling)]/apollo_inputs$apollo_scaling
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed after inserting parameter scaling.", 
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                                names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
                  success              = FALSE))
    }
  }
  ### Update ORIG version to keep changes so far
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_inputs$apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_inputs$apollo_lcPars
  
  ### Introduce quotes into apollo_rrm
  if(debug) apollo_print("- Inserting quotes in settings for apollo_rrm (if present)")
  apollo_probabilities <- apollo_insertRRMQuotes(apollo_probabilities)
  
  ### Evaluate apollo_probabilities after current changes, return immediately if it doesn't work
  if(validate){
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed after additional syntax processing.", 
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                                names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
                  success              = FALSE))
    }
  }
  ### Update ORIG version to keep changes so far
  apollo_probabilities_ORIG <- apollo_probabilities
  apollo_randCoeff_ORIG     <- apollo_inputs$apollo_randCoeff
  apollo_lcPars_ORIG        <- apollo_inputs$apollo_lcPars
  
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
    test2 <- tryCatch(apollo_probabilities(apollo_beta_shifted, apollo_inputs), error=function(e) NULL)
    test <- !is.null(test1) && !is.null(test2) && is.numeric(test1) && is.numeric(test2)
    test <- test && !any(is.nan(test1)) && !any(is.nan(test2)) && abs(sum(test2)/sum(test1) - 1) < 0.001
    if(!test){
      # If they are different or evaluation of test2 failed, then undo changes
      apollo_print(paste0("The pre-processing of 'apollo_probabilities' failed after inserting functions.", 
                          " Your model may still run, but this indicates a potential problem. Please contact the", 
                          " developers for assistance!"),  pause=5, type="w")
      return(list(apollo_probabilities = apollo_probabilities_ORIG, 
                  apollo_randCoeff     = apollo_randCoeff_ORIG, 
                  apollo_lcPars        = apollo_lcPars_ORIG,
                  #apollo_scaling       = setNames(rep(1, length(apollo_inputs$apollo_scaling)),
                  #                                names(apollo_inputs$apollo_scaling)), 
                  apollo_scaling       = apollo_inputs$apollo_scaling,
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
              manualScaling        = apollo_inputs$manualScaling, 
              success              = TRUE))
}