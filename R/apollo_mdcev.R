#' Calculates MDCEV likelihoods.
#'
#' Calculates the likelihood of a Multiple Discrete Continuous Extreme Value (MDCEV) model.
#'
#' @param mdcev_settings List of settings for the MDCEV model. It must include the following.
#'                       \itemize{
#'                         \item V: Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#'                         \item alternatives: Character vector. Names of alternatives, elements must match the names in list 'V'.
#'                         \item alpha: Named list. Alpha parameters for each alternative, including for any outside good. As many elements as alternatives.
#'                         \item gamma: Named list. Gamma parameters for each alternative, excluding any outside good. As many elements as inside good alternatives.
#'                         \item sigma: Numeric scalar. Scale parameter of the model extreme value type I error.
#'                         \item cost: Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item avail: Named list. Availabilities of alternatives, one element per alternative. Names of elements must match those in argument 'alternatives'. Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation. If all alternatives are always available, then user can just omit this argument.
#'                         \item continuousChoice: Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item budget: Numeric vector. Budget for each observation.
#'                         \item minConsumption: Named list of scalars or numeric vectors. Minimum consumption of the alternatives, if consumed. As many elements as alternatives. Names must match those in \code{alternatives}.
#'                         \item outside: Character. Optional name of the outside good.
#'                         \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                         #'                       }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate" Used for model estimation.
#'                        \item "prediction" Used for model predictions.
#'                        \item "validate" Used for validating input.
#'                        \item "zero_LL" Used for calculating null likelihood.
#'                        \item "conditionals" Used for calculating conditionals.
#'                        \item "output" Used for preparing output after model estimation.
#'                        \item "raw" Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": A matrix with one row per observation, and means and s.d. of predicted consumptions.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": Not applicable.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @export
apollo_mdcev <- function(mdcev_settings,functionality){
  if(is.null(mdcev_settings[["alternatives"]])) stop("The mdcev_settings list needs to include an object called \"alternatives\"!")
  if(is.null(mdcev_settings[["avail"]])) stop("The mdcev_settings list needs to include an object called \"avail\"!")
  if(is.null(mdcev_settings[["continuousChoice"]])) stop("The mdcev_settings list needs to include an object called \"continuousChoice\"!")
  if(is.null(mdcev_settings[["V"]])) stop("The mdcev_settings list needs to include an object called \"V\"!")
  if(is.null(mdcev_settings[["alpha"]])) stop("The mdcev_settings list needs to include an object called \"alpha\"!")
  if(is.null(mdcev_settings[["gamma"]])) stop("The mdcev_settings list needs to include an object called \"gamma\"!")
  if(is.null(mdcev_settings[["sigma"]])) stop("The mdcev_settings list needs to include an object called \"sigma\"!")
  if(is.null(mdcev_settings[["cost"]])) stop("The mdcev_settings list needs to include an object called \"cost\"!")
  if(is.null(mdcev_settings[["budget"]])) stop("The mdcev_settings list needs to include an object called \"budget\"!")
  if(is.null(mdcev_settings[["minConsumption"]])) mdcev_settings[["minConsumption"]]=NA
  if(is.null(mdcev_settings[["rows"]])) mdcev_settings[["rows"]]="all"
  if(is.null(mdcev_settings[["outside"]])) mdcev_settings[["outside"]]=NA
  
  alternatives     = mdcev_settings[["alternatives"]]
  avail            = mdcev_settings[["avail"]]
  continuousChoice = mdcev_settings[["continuousChoice"]]
  V                = mdcev_settings[["V"]]
  alpha            = mdcev_settings[["alpha"]]
  gamma            = mdcev_settings[["gamma"]]
  sigma            = mdcev_settings[["sigma"]]
  cost             = mdcev_settings[["cost"]]
  budget           = mdcev_settings[["budget"]]
  minConsumption   = mdcev_settings[["minConsumption"]]
  rows             = mdcev_settings[["rows"]]
  outsideName      = mdcev_settings[["outside"]]
  if(anyNA(outsideName) && "outside" %in% alternatives) outsideName = "outside"
  if(!anyNA(outsideName) && !(outsideName %in% alternatives)) stop("Name provided for outside good does not correspond to an alternative!")
  
  apollo_inputs <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE ), error=function(e) return(NA))
  
  if(any(budget<=0)) stop("Budget for some rows in data is less than or equal to zero!")
  
  
  if(!anyNA(outsideName)){
    if(is.null(gamma[[outsideName]])) gamma[[outsideName]]=1
    # Reorder 'alternatives'
    alternatives=c(outsideName,alternatives[-which(alternatives==outsideName)])
    # Reorder V and avail to match 'alternatives' order, if necessary
    if(any(alternatives != names(V))) V <- V[alternatives]
    if(any(alternatives != names(avail))) avail <- avail[alternatives]
    if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
    if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
    if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
    if(any(alternatives != names(cost))) cost <- cost[alternatives]
    if(!anyNA(minConsumption)) if(any(alternatives != names(minConsumption))) minConsumption <- minConsumption[alternatives]
    # Rename first element
    alternatives = c("outside", alternatives[2:length(alternatives)])
    names(V    ) = alternatives
    names(avail) = alternatives
    names(alpha) = alternatives
    names(gamma) = alternatives
    names(cost ) = alternatives
    names(continuousChoice) = alternatives
    if(!anyNA(minConsumption)) names(minConsumption) = alternatives
    # Call mdcev outside
    ans <- apollo_mdcevOutside(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice, budget, functionality, minConsumption, outsideName, rows)
  } else {
    ans <- apollo_mdcevInside(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice, budget, functionality, minConsumption, rows)
  }

  return(ans)
}
