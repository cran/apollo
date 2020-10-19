#' Prepares input for \code{apollo_estimate}
#'
#' Searches the user work space (.GlobalEnv) for all necessary input to run \code{apollo_estimate}, and packs it in a single list.
#'
#' All arguments to this function are optional. If the function is called without arguments, then it it will look in
#' the user workspace (i.e. the global environment) for variables with the same name as its ommited arguments.
#' We strongly recommend users to visit www.ApolloChoiceModelling.com for examples on how to use Apollo.
#' In the website, users will also find a detailed manual and a user-group for help and further reference.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                       \itemize{
#'                         \item \code{modelName}: Character. Name of the model. Used when saving the output to files. Avoid characters not allowed in file names, such as \code{\\}, \code{*}, \code{:}, etc.
#'                         \item \code{modelDescr}: Character. Description of the model. Used in output files.
#'                         \item \code{indivID}: Character. Name of column in the database with each decision maker's ID.
#'                         \item \code{mixing}: Boolean. TRUE for models that include random parameters.
#'                         \item \code{nCores}: Numeric>0. Number of threads (processors) to use in estimation of the model.
#'                         \item \code{workInLogs}: Boolean. TRUE for higher numeric stability at the expense of computational time.
#'                                                  Useful for panel models only. Default is FALSE.
#'                         \item \code{seed}: Numeric. Seed for random number generation.
#'                         \item \code{HB}: Boolean. TRUE if using RSGHB for Bayesian estimation of model.
#'                         \item \code{noValidation}: Boolean. TRUE if user does not wish model input to be validated before estimation - FALSE by default.
#'                         \item \code{noDiagnostics}: Boolean. TRUE if user does not wish model diagnostics to be printed - FALSE by default.
#'                         \item \code{panelData}: Boolean. TRUE if using panelData data (created automatically by \code{apollo_validateControl}).
#'                         \item \code{weights}: Character. Name of column in database containing weights for estimation.
#'                       }
#' @param apollo_HB List. Contains options for bayesian estimation. See \code{?RSGHB::doHB} for details.
#'                   Parameters \code{modelname}, \code{gVarNamesFixed}, \code{gVarNamesNormal},
#'                   \code{gDIST}, \code{svN} and \code{FC} are automatically set based on the
#'                   other arguments of this function.
#'                   Other settings to include are the following.
#'                   \itemize{
#'                     \item \strong{\code{hbDist}} \emph{Mandatory} setting. A named character vector determining
#'                                                  the distribution of each parameter to be estimated. Possible 
#'                                                  values are as follows.
#'                                                  \itemize{
#'                                                    \item \code{"DNE"}: Parameter kept at its starting value (not estimated).
#'                                                    \item \code{"F"}: Fixed (as in non-random) parameter.
#'                                                    \item \code{"N"}: Normal.
#'                                                    \item \code{"LN+"}: Positive log-normal.
#'                                                    \item \code{"LN-"}: Negative log-normal.
#'                                                    \item \code{"CN+"}: Positive censored normal.
#'                                                    \item \code{"CN-"}: Negative censored normal.
#'                                                    \item \code{"JSB"}: Johnson SB.
#'                                                  }
#'                     \item \strong{\code{constraintNorm}} Character vector. Constraints for \emph{random} coefficients 
#'                                                          in bayesian estimation. Constraints can be written as 
#'                                                          "b1>b2", "b1<b2", "b1>0", or "b1<0".
#'                     \item \strong{\code{fixedA}} Named numeric vector. Contains the names and fixed mean values of 
#'                                                  random parameters. For example, c(b1=0) fixes the mean of b1 to zero.
#'                     \item \strong{\code{fixedD}} Named numeric vector. Contains the names and fixed variance of 
#'                                                  random parameters. For example, c(b1=1) fixes the variance of b1 to zero.
#'                   }
#' @param apollo_draws List of arguments describing the inter and intra individual draws. Required only if \code{apollo_control$mixing = TRUE}. Unused elements can be ommited.
#'                  \itemize{
#'                    \item \code{interDrawsType}: Character. Type of inter-individual draws ('halton','mlhs','pmc','sobol','sobolOwen',
#'                                                 'sobolFaureTezuka', 'sobolOwenFaureTezuka' or the name of an object loaded in memory,
#'                                                 see manual in www.ApolloChoiceModelling.com for details).
#'                    \item \code{interNDraws}: Numeric scalar (>=0). Number of inter-individual draws per individual. Should be set to 0 if not using them.
#'                    \item \code{interUnifDraws}: Character vector. Names of uniform-distributed inter-individual draws.
#'                    \item \code{interNormDraws}: Character vector. Names of normaly distributed inter-individual draws.
#'                    \item \code{intraDrawsType}: Character. Type of intra-individual draws ('halton','mlhs','pmc','sobol','sobolOwen','sobolFaureTezuka', 'sobolOwenFaureTezuka' or the name of an object loaded in memory).
#'                    \item \code{intraNDraws}: Numeric scalar (>=0). Number of intra-individual draws per individual. Should be set to 0 if not using them.
#'                    \item \code{intraUnifDraws}: Character vector. Names of uniform-distributed intra-individual draws.
#'                    \item \code{intraNormDraws}: Character vector. Names of normaly distributed intra-individual draws.
#'                  }
#' @param apollo_randCoeff Function. Used with mixing models. Constructs the random parameters of a mixing model. Receives two arguments:
#'                      \itemize{
#'                        \item \code{apollo_beta}: Named numeric vector. Names and values of model parameters. 
#'                        \item \code{apollo_inputs}: The output of this function (\code{apollo_validateInputs}).
#'                      }
#' @param apollo_lcPars Function. Used with latent class models. Constructs a list of parameters for each latent class. Receives two arguments:
#'                      \itemize{
#'                        \item \code{apollo_beta}: Named numeric vector. Names and values of model parameters. 
#'                        \item \code{apollo_inputs}: The output of this function (\code{apollo_validateInputs}).
#'                      }
#' @param recycle Logical. If TRUE, an older version of apollo_inputs is looked for in the calling environment (parent frame), and any
#'                element in that old version created by the user is copied into the new apollo_inputs returned by this function.
#'                For \code{recycle=TRUE} to work, the old version of apollo_inputs \strong{must} be named "apollo_inputs".
#'                If FALSE, nothing is copied from any older version of apollo_inputs. FALSE is the default.
#' @param silent Logical. TRUE to keep the function from printing to the console. Default is FALSE.
#' @return List grouping several required input for model estimation.
#' @export
apollo_validateInputs <- function(apollo_beta=NA, apollo_fixed=NA, database=NA,
                                  apollo_control=NA, 
                                  apollo_HB=NA, apollo_draws=NA,
                                  apollo_randCoeff=NA, apollo_lcPars=NA,
                                  recycle=FALSE, silent=FALSE){
  
  ### Try to recover mandatory variables from global environment if not provided
  tmp <- c("database", paste0("apollo_", c("beta", "fixed", "control"))) 
  for(i in tmp){
    x <- get(i, envir=environment(), inherits=FALSE)
    if(length(x)==1 && is.na(x)) x <- tryCatch( get(i, envir=globalenv()), error=function(e) NA )
    if(length(x)==1 && is.na(x)) stop("No variable called ", i, " found in user workspace (i.e. global environment).") else assign(i, x, envir=environment())
  }; rm(tmp, x, i)
  
  ### Validate apollo_beta & apollo_fixed
  if(!is.numeric(apollo_beta) | !is.vector(apollo_beta) | is.null(names(apollo_beta))) stop("The \"apollo_beta\" argument needs to be a named vector")
  if(length(apollo_fixed)>0 && !is.character(apollo_fixed)) stop("'apollo_fixed' is not an empty vector nor a vector of names.")
  if(anyDuplicated(names(apollo_beta))){
    txt <- paste0(names(apollo_beta)[duplicated(names(apollo_beta))], collapse=", ")
    txt <- paste0("The \"apollo_beta\" argument contains duplicate elements (", txt, ").")
    stop(txt)
  }
  if(anyDuplicated(apollo_fixed)){
    txt <- paste0(apollo_fixed[duplicated(apollo_fixed)], collapse=", ")
    txt <- paste0("The \"apollo_fixed\" argument contains duplicate elements (", txt, ").")
    stop(txt)
  }
  if(!all(apollo_fixed %in% names(apollo_beta))){
    txt <- apollo_fixed[!(apollo_fixed %in% names(apollo_beta))]
    txt <- paste0(txt, collapse=", ")
    txt <- paste0("Some parameters included in 'apollo_fixed' (", txt, ") are not included in 'apollo_beta'.")
    stop(txt)
  }
  
  ### Validate apollo_control, database
  apollo_control <- apollo_validateControl(database, apollo_control, silent=silent)
  database       <- apollo_validateData(database, apollo_control, silent=silent)
  
  ### Try to recover apollo_HB if appropiate, and sets the default value for the missing parts
  if(!apollo_control$HB) apollo_HB <- NA else{
    if(length(apollo_HB)==1 && is.na(apollo_HB)) apollo_HB <- tryCatch( get("apollo_HB", envir=globalenv()), error=function(e) NA )
    if(length(apollo_HB)==1 && is.na(apollo_HB)) stop("No variable called apollo_HB found in user workspace (i.e. global environment)!")
    apollo_HB <- apollo_validateHBControl(apollo_HB, apollo_beta, apollo_fixed, apollo_control, silent)
    # Check that database is sorted by indiv
    ordID <- sort(database[,apollo_control$indivID])
    if(any(ordID!=database[,apollo_control$indivID])) stop('For estimation using HB the database needs to be sorted by ID!')
  }
  if(apollo_control$HB && anyNA(apollo_HB)) stop("Argument apollo_HB must be provided when using Bayesian estimation!")
  
  ### Try to recover apollo_draws and apollo_randCoeff if appropiate, and sets the default value for the missing parts
  if(!apollo_control$mixing){
    if(!is.function(apollo_randCoeff)) apollo_randCoeff <- tryCatch( get("apollo_randCoeff", envir=globalenv()), error=function(e) NA )
    if(is.function(apollo_randCoeff)) apollo_print("Function called apollo_randCoeff found in user workspace will be ignored as model not using mixing.")
    if(length(apollo_draws)==1 && is.na(apollo_draws)) apollo_draws <- tryCatch( get("apollo_draws", envir=globalenv()), error=function(e) NA )
    if(length(apollo_draws)==1 && !is.na(apollo_draws)) apollo_print("Variable called apollo_draws found in user workspace will be ignored as model not using mixing.")
    apollo_draws     <- NA
    draws            <- NA
    apollo_randCoeff <- NA
  } else{
    if(length(apollo_draws)==1 && is.na(apollo_draws)) apollo_draws <- tryCatch( get("apollo_draws", envir=globalenv()), error=function(e) NA )
    if(length(apollo_draws)==1 && is.na(apollo_draws)) stop("Mixing set to TRUE in apollo_control, but no variable called apollo_draws found in user workspace (i.e. global environment).")
    default <- list(interDrawsType="halton", interNDraws=0, interUnifDraws=c(), interNormDraws=c(), 
                    intraDrawsType='halton', intraNDraws=0, intraUnifDraws=c(), intraNormDraws=c())
    for(i in names(default)) if(!(i %in% names(apollo_draws))) apollo_draws[[i]] <- default[[i]]
    
    if(!is.function(apollo_randCoeff)) apollo_randCoeff <- tryCatch( get("apollo_randCoeff", envir=globalenv()), error=function(e) NA )
    if(!is.function(apollo_randCoeff)) stop("Mixing set to TRUE in apollo_control, but no function called apollo_randCoeff found in user workspace (i.e. global environment).")
  }
  
  ### Try to recover apollo_lcPars if not provided
  if(length(apollo_lcPars)==1 && is.na(apollo_lcPars)) apollo_lcPars <- tryCatch( get("apollo_lcPars", envir=globalenv()), error=function(e) NA )
  
  ### Create apolloLog
  #apolloLog      <- new.env(parent=emptyenv())
  #apolloLog$cppV <- list()
  
  ### Pack everything into a single list
  apollo_inputs <- list(database=database, apollo_control=apollo_control, 
                        apollo_HB=apollo_HB, apollo_draws=apollo_draws, apollo_randCoeff=apollo_randCoeff,
                        apollo_lcPars=apollo_lcPars, draws=NA, #apolloLog=apolloLog, 
                        apollo_fixed=apollo_fixed, silent=silent, class_specific=0,
                        apollo_beta_names=names(apollo_beta))
  
  ### Make draws
  if(apollo_control$mixing) apollo_inputs$draws <- apollo_makeDraws(apollo_inputs, silent=silent)
  
  ### If recycle is TRUE, fetch "old" apollo_inputs from calling enironment
  ### and copy all elements not in the new version, from the old version.
  old_inputs <- tryCatch(get("apollo_inputs", envir=parent.frame(), inherits=TRUE),
                         error=function(e) return(NULL))
  if(recycle && !is.null(old_inputs)){
    toCopy <- names(old_inputs)[!(names(old_inputs) %in% names(apollo_inputs))]
    if(length(toCopy)>0) for(i in toCopy) apollo_inputs[[i]] <- old_inputs[[i]]
  }
  
  ### Check for repeated names across apollo_beta, database, apollo_draws, apollo_randCoeff, apollo_lcPars
  namesList <- list(apollo_beta = names(apollo_beta), # store apollo_beta names
                    database    = names(database))    # store database names
  if(apollo_control$mixing){ # run apollo_randCoeff and store its names
    namesList$draws <- names(apollo_inputs$draws)
    env <- c(apollo_inputs$database, as.list(apollo_beta), apollo_inputs$draws)
    f   <- apollo_randCoeff
    environment(f) <- list2env(env, hash=TRUE, parent=parent.frame())
    randCoeff <- f(apollo_beta, apollo_inputs)
    namesList$apollo_randCoeff <- names(randCoeff)
    rm(env, f)
  }
  if(is.function(apollo_inputs$apollo_lcPars)){ # run apollo_lcPars and store its names
    env <- c(apollo_inputs$database, as.list(apollo_beta))
    if(exists('randCoeff', inherits=FALSE)) env <- c(env, apollo_inputs$draws, randCoeff)
    f   <- apollo_lcPars
    environment(f) <- list2env(env, hash=TRUE, parent=parent.frame())
    lc_pars <- f(apollo_beta, apollo_inputs)
    namesList$apollo_lcPars <- names(lc_pars)
    rm(env, f, lc_pars)
  }; if(exists('randCoeff', inherits=FALSE)) rm(randCoeff)
  for(i in 1:(length(namesList)-1)) for(j in (i+1):length(namesList)){ # check for name overlaps
    tmp <- namesList[[i]][which(namesList[[i]] %in% namesList[[j]])]
    if(length(tmp)>0){
      one <- length(tmp)==1
      if(!one) tmp <- paste0(tmp, collapse='", "')
      tmp <- paste0('The variable name', ifelse(one,'','s'), ' "', tmp, '" ', ifelse(one, 'is', 'are'),
                    ' used in both ', names(namesList)[i], ' and ', names(namesList)[j], '. ',
                    'Variable names must be unique.')
      stop(tmp)
    } 
  }
  
  #### Check for reserved names
  #reserved <- c("alpha", "alternatives", "altStart", "apollo_beta", 
  #              "apolloBetaMax", "apolloBetaMin", "apollo_fixed", 
  #              "apollo_lcPars", "apollo_probabilities", "apollo_randCoeff", 
  #              "attrScaling", "attrValues", "attrWeights", 
  #              "avail", "baseModel", "bfgsIter", 
  #              "bootstrapSE", "bootstrapSeed", "budget", 
  #              "choiceVar", "choiceVars", "classProb", 
  #              "cnlNests", "cnlStructure", "coding", 
  #              "componentName", "constraints", "constraintsNorm", 
  #              "continuousChoice", "cost", "database", 
  #              "dTest", "estimateDigits", "estimationRoutine", 
  #              "explanators", "functionality", "gamma", 
  #              "generalModel", "gFULLCV", "gINFOSKIP", 
  #              "gNCREP", "gNEREP", "gTest", "HB", "hbDist", 
  #              "hessianRoutine", "inClassProb", "indivID", 
  #              "inputModelName", "interDrawsType", "interNDraws", 
  #              "interNormDraws", "interUnifDraws", "intraDrawsType", 
  #              "intraNDraws", "intraNormDraws", "intraUnifDraws", 
  #              "llTest", "maxIterations", "maxStages", "mdcnevNests", 
  #              "mdcnevStructure", "minConsumption", "mixing", 
  #              "modelComponent", "modelDescr", "modelName", 
  #              "modelNames", "mu", "multiPar1", "multiPar2", 
  #              "nCandidates", "nCores", "nCoresTry", "nDrawsTry", 
  #              "nlNests", "nlStructure", "noDiagnostics", 
  #              "noValidation", "nRep", "numDeriv_settings", "operation", 
  #              "outcomeNormal", "outcomeOrdered", "outside", "overwriteFixed", 
  #              "P", "panelData", "parName1", "parName2", "pDigits", 
  #              "printChange", "printClassical", "printCorr", "printCovar", 
  #              "printDiagnostics", "printLevel", "printOutliers", "printPVal", 
  #              "printT1", "procPars", "rows", "samples", "saveCorr", "saveCov", 
  #              "saveEst", "saveModelObject", "scales", "scaling", "seed", 
  #              "sigma", "silent", "smartStart", "sortByDate", "subsamples", 
  #              "tau", "tDigits", "V", "validationSize", "weights", "workInLogs", 
  #              "writeF12", "writeIter", "xNormal", "apollo_control", 
  #              "apollo_draws", "apollo_HB", "apollo_inputs", "bootstrap_settings", 
  #              "choiceAnalysis_settings", "cnl_settings", "combineResults_settings", 
  #              "deltaMethod_settings", "dft_settings", "el_settings", 
  #              "estimate_settings", "fitsTest_settings", "lc_settings", 
  #              "mdcev_settings", "mdcnev_settings", "mnl_settings", "model", 
  #              "modelOutput_settings", "nl_settings", "normalDensity_settings", 
  #              "ol_settings", "op_settings", "outOfSample_settings", 
  #              "prediction_settings", "saveOutput_settings", "searchStart_settings", 
  #              "sharesTest_settings", "speedTest_settings", "apollo_attach", 
  #              "apollo_avgInterDraws", "apollo_avgIntraDraws", "apollo_bootstrap", 
  #              "apollo_choiceAnalysis", "apollo_cnl", "apollo_combineModels", 
  #              "apollo_combineResults", "apollo_conditionals", "apollo_deltaMethod", 
  #              "apollo_detach", "apollo_dft", "apollo_el", "apollo_estimate", 
  #              "apollo_firstRow", "apollo_fitsTest", "apollo_initialise", 
  #              "apollo_lcConditionals", "apollo_lc", "apollo_lcUnconditionals", 
  #              "apollo_llCalc", "apollo_loadModel", "apollo_lrTest", 
  #              "apollo_mdcev", "apollo_mdcnev", "apollo_mnl", "apollo_modelOutput", 
  #              "apollo_nl", "apollo_normalDensity", "apollo_ol", "apollo_op", 
  #              "apollo_outOfSample", "apollo_panelProd", "apollo_prediction", 
  #              "apollo_prepareProb", "apollo_readBeta", "apollo_saveOutput", 
  #              "apollo_searchStart", "apollo_sharesTest", "apollo_speedTest", 
  #              "apollo_unconditionals", "apollo_validateInputs", "apollo_weighting")
  
  return(apollo_inputs)
}
