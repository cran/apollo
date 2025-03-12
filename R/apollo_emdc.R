#' MDC model with exogenous budget
#' 
#' Calculates the likelihood function of the MDC model with exogenous budget. Can also predict and validate inputs.
#' 
#' This model extends the traditional multiple discrete-continuous (MDC) framework by (i) making the 
#' marginal utility of the outside good deterministic, and (ii) including complementarity and 
#' substitution in the model formulation. See the following working paper for more details:
#' 
#' Palma, D. & Hess, S. (2022) Extending the Multiple Discrete Continuous (MDC) modelling 
#' framework to consider complementarity, substitution, and an unobserved budget. Transportation 
#' Reserarch 161B, 13 - 35. https://doi.org/10.1016/j.trb.2022.04.005
#' 
#' @param emdc_settings List of settings for the model. It includes the following.
#'                        \itemize{
#'                          \item \strong{\code{avail}}: Named list of numeric vectors. Availability of each product. Can also be called "A".
#'                          \item \strong{\code{budget}}: Optional numeric vector. Budget. Must be bigger that the expenditure on all inside goods. Can also be called "B".
#'                          \item \strong{\code{cost}}: Named list of numeric vectors. Price of each product.
#'                          \item \strong{\code{delta}}: Lower triangular numeric matrix, or list of lists. Complementarity/substitution parameter.
#'                          \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount consumed of each inside good. Outside good must not be included. Can also be called "X".
#'                          \item \strong{\code{gamma}}: Named list of numeric vectors. Satiation parameter of each product.
#'                          \item \strong{\code{nRep}}: Scalar positive integer. Number of repetitions used when prediction
#'                          \item \strong{\code{sigma}}: Numeric vector or scalar. Standard deviation of the error term. Default is one.
#'                          \item \strong{\code{timeLimit}}: Positive scalar. Maximum amount of seconds the optimiser can spend calculating a prediction before setting it to NA.
#'                          \item \strong{\code{tol}}: Positive scalar. Tolerance of the prediction algorithm.
#'                          \item \strong{\code{utilities}}: Named list of numeric vectors (or matrices or arrays). Base utility of each product. Can also be called "V".
#'                          \item \strong{\code{utilityOutside}}: Numeric vector (or matrix or array). Shadow price of the budget. Must be normalised to 0 for at least one individual. Default is 0 for every observation. Can also be called "V0".
#'                        }
#' @param functionality Character. Either "validate", "zero_LL", "estimate", "conditionals", "raw", "output" or "prediction"
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @export
apollo_emdc <- function(emdc_settings, functionality="estimate"){
  
  if(!all(sort(names(emdc_settings$continuousChoice))==sort(names(emdc_settings$utilities)))) stop("Names of alternatives do not match between continuousChoice and utilities!")
  if(!all(sort(names(emdc_settings$continuousChoice))==sort(names(emdc_settings$gamma)))) stop("Names of alternatives do not match between continuousChoice and gamma!")
  if(!is.null(emdc_settings$avail) && (!all(sort(names(emdc_settings$continuousChoice))==sort(names(emdc_settings$avail))))) stop("Names of alternatives do not match between continuousChoice and avail!")

    
  # Rename input if necessary
  map <- c(X = "continuousChoice", B = "budget", A = "avail", 
           V0= "utilityOutside",   V = "utilities")
  for(i in 1:length(map)) if(!is.null(emdc_settings[[map[i]]])){
    emdc_settings[[names(map)[i]]] <- emdc_settings[[map[i]]]
    emdc_settings[[map[i]]]        <- NULL
  }; rm(i, map)
  
  if("B" %in% names(emdc_settings)) return(apollo_emdc1(emdc_settings, functionality))
  return(apollo_emdc2(emdc_settings, functionality))
}