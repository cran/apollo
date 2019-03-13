#' Prepares input for \code{apollo_estimate}
#'
#' Searches the user work space for all necessary input to run \code{apollo_estimate}, and packs it in a single list.
#'
#' All arguments to this function are optional. If the function is called without arguments, then it it will look in
#' the user workspace (i.e. the global environment) for variables with the same name as its ommited arguments.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                       \itemize{
#'                         \item modelName: Character. Name of the model. Used when saving the output to files.
#'                         \item modelDescr: Character. Description of the model. Used in output files.
#'                         \item indivID: Character. Name of column in the database with each decision maker's ID.
#'                         \item mixing: Boolean. TRUE for models that include random parameters.
#'                         \item nCores: Numeric>0. Number of threads (processors) to use in estimation of the model.
#'                         \item workInLogs: Boolean. TRUE for higher numeric stability at the expense of computational time.
#'                                           Useful for panel models only. Default is FALSE.
#'                         \item seed: Numeric. Seed for random number generation.
#'                         \item HB: Boolean. TRUE if using RSGHB for Bayesian estimation of model.
#'                         \item noValidation: Boolean. TRUE if user does not wish model input to be validated before estimation - FALSE by default.
#'                         \item noDiagnostics: Boolean. TRUE if user does not wish model diagnostics to be printed - FALSE by default.
#'                         \item panelData: Boolean. TRUE if using panelData data (created automatically by \code{apollo_validateControl}).
#'                       }
#' @param apollo_HB List. Contains options for bayesian estimation. See \code{?RSGHB::doHB} for details.
#'                   Parameters \code{modelname}, \code{gVarNamesFixed}, \code{gVarNamesNormal},
#'                   \code{gDIST}, \code{svN} and \code{FC} are automatically set based on the
#'                   other arguments of this function.
#'                   It should also include a named character vector called \code{hbDist} identifying 
#'                   the distribution of each parameter to be estimated. Possible values are as follows.
#'                   \itemize{
#'                     \item "DNE": Parameter kept at its starting value (not estimated).
#'                     \item "F": Fixed (as in non-random) parameter.
#'                     \item "N": Normal.
#'                     \item "LN+": Positive log-normal.
#'                     \item "LN-": Negative log-normal.
#'                     \item "CN+": Positive censored normal.
#'                     \item "CN-": Negative censored normal.
#'                     \item "JSB": Johnson SB.
#'                   }
#' @param apollo_draws List of arguments describing the inter and intra individual draws.
#'                  \itemize{
#'                    \item interDrawsType: Character. Type of inter-individual draws ('halton','mlhs','pmc','sobol','sobolOwen','sobolFaureTezuka', 'sobolOwenFaureTezuka' or the name of an object loaded in memory).
#'                    \item interNDraws: Numeric scalar (>=0). Number of inter-individual draws per individual. Should be set to 0 if not using them.
#'                    \item interUnifDraws: Character vector. Names of uniform-distributed inter-individual draws.
#'                    \item interNormDraws: Character vector. Names of normaly distributed inter-individual draws.
#'                    \item intraDrawsType: Character. Type of intra-individual draws ('halton','mlhs','pmc','sobol','sobolOwen','sobolFaureTezuka', 'sobolOwenFaureTezuka' or the name of an object loaded in memory).
#'                    \item intraNDraws: Numeric scalar (>=0). Number of intra-individual draws per individual. Should be set to 0 if not using them.
#'                    \item intraUnifDraws: Character vector. Names of uniform-distributed intra-individual draws.
#'                    \item intraNormDraws: Character vector. Names of normaly distributed intra-individual draws.
#'                  }
#' @param apollo_randCoeff Function. Used with mixing models. Constructs the random parameters of a mixing model. Receives two arguments:
#'                      \itemize{
#'                        \item apollo_beta: Named numeric vector. Names and values of model parameters. 
#'                        \item apollo_inputs: The output of this function (\code{apollo_validateInputs}).
#'                      }
#' @param apollo_lcPars Function. Used with latent class models. Constructs a list of parameters for each latent class. Receives two arguments:
#'                      \itemize{
#'                        \item apollo_beta: Named numeric vector. Names and values of model parameters. 
#'                        \item apollo_inputs: The output of this function (\code{apollo_validateInputs}).
#'                      }
#' @param silent Boolean. TRUE to keep the function from printing to the console. Default is FALSE.
#' @return List grouping several required input for model estimation.
#' @export
apollo_validateInputs <- function(apollo_beta=NA, apollo_fixed=NA, database=NA,
                                 apollo_control=NA, 
                                 apollo_HB=NA, apollo_draws=NA,
                                 apollo_randCoeff=NA, apollo_lcPars=NA,
                                 silent=FALSE){
  
  tmp <- c("database", paste0("apollo_", c("beta", "fixed", "control"))) 
  for(i in tmp){
    x <- get(i, envir=environment(), inherits=FALSE)
    if(length(x)==1 && is.na(x)) x <- tryCatch( get(i, envir=globalenv()), error=function(e) NA )
    if(length(x)==1 && is.na(x)) stop("No variable called ", i, " found in user workspace (i.e. global environment).") else assign(i, x, envir=environment())
  }; rm(tmp, x, i)

  apollo_control <- apollo_validateControl(database, apollo_control, silent=silent)
  database       <- apollo_validateData(database, apollo_control, silent=silent)

  if(!apollo_control$HB) apollo_HB <- NA else{
    if(length(apollo_HB)==1 && is.na(apollo_HB)) apollo_HB <- tryCatch( get("apollo_HB", envir=globalenv()), error=function(e) NA )
    if(length(apollo_HB)==1 && is.na(apollo_HB)) stop("No variable called 'apollo_HB' found in user workspace (i.e. global environment).")
    apollo_HB <- apollo_validateHBControl(apollo_HB, apollo_beta, apollo_fixed, apollo_control)
  }
  if(apollo_control$HB && anyNA(apollo_HB)) stop("Argument 'apollo_HB' must be provided when using Bayesian estimation.")

  if(!apollo_control$mixing){
    apollo_draws <- NA
    draws <- NA
    apollo_randCoeff <- NA
  } else{
    if(length(apollo_draws)==1 && is.na(apollo_draws)) apollo_draws <- tryCatch( get("apollo_draws", envir=globalenv()), error=function(e) NA )
    if(length(apollo_draws)==1 && is.na(apollo_draws)) stop("No variable called 'apollo_draws' found in user workspace (i.e. global environment).")
    if(!is.function(apollo_randCoeff)) apollo_randCoeff <- tryCatch( get("apollo_randCoeff", envir=globalenv()), error=function(e) NA )
    if(!is.function(apollo_randCoeff)) stop("No variable called 'apollo_randCoeff' found in user workspace (i.e. global environment).")
  }

  if(length(apollo_lcPars)==1 && is.na(apollo_lcPars)) apollo_lcPars <- tryCatch( get("apollo_lcPars", envir=globalenv()), error=function(e) NA )

  apollo_inputs <- list(database=database, apollo_control=apollo_control, 
                       apollo_HB=apollo_HB, apollo_draws=apollo_draws, apollo_randCoeff=apollo_randCoeff,
                       apollo_lcPars=apollo_lcPars, draws=NA)

  if(apollo_control$mixing) apollo_inputs$draws <- apollo_makeDraws(apollo_inputs, silent=silent)

  if(apollo_control$mixing){
    if(anyNA(apollo_inputs$draws)) stop("Argument 'draws' must be provided when estimating mixture models. Use apollo_makeDraws.")
    if(!is.function(apollo_inputs$apollo_randCoeff)) stop("Argument 'apollo_randCoeff' must be provided when estimating mixture models.")
    if(!apollo_inputs$apollo_control$panel & dim(apollo_inputs$draws[[1]])[2]>1) warning('Inter-person draws are used without a panel structure. This is unusual.')
  }

  return(apollo_inputs)
}
