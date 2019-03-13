#' Attaches predefined variables.
#'
#' Attaches parameters and data to allow users to refer to individual variables by name without reference to the object they are contained in.
#'
#' This function should be called at the beginning of \code{apollo_probabilities}
#' to make writing the log-likelihood more user-friendly. If used, then \link{apollo_detach} should
#' be called at the end \code{apollo_probabilities}, or more conveniently, using
#' \link{on.exit}.
#' \code{apollo_attach} attaches \code{apollo_beta}, \code{database}, \code{draws},
#' and the output of \code{apollo_randCoeff} and \code{apollo_lcPars}, if they are
#' defined by the user.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return Nothing.
#' @examples
#' apollo_beta  <- c(b1=0.3, b2=-0.5)
#' apollo_fixed <- c()
#' apollo_control <- list(indivID="id", mixing = FALSE, panelData = FALSE)
#' database <- data.frame(id=1:100, x1=stats::runif(100), x2=stats::runif(100))
#' apollo_inputs <- apollo_validateInputs()
#' apollo_attach(apollo_beta, apollo_inputs)
#' V = b1*x1 + b2*x2
#' apollo_detach(apollo_beta, apollo_inputs)
#' @export
apollo_attach=function(apollo_beta, apollo_inputs){
  apollo_control   = apollo_inputs[["apollo_control"]]
  database         = apollo_inputs[["database"]]
  draws            = apollo_inputs[["draws"]]
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars    = apollo_inputs[["apollo_lcPars"]]

  attach(as.list(apollo_beta))
  attach(database)

  if(apollo_control$HB==FALSE && apollo_control$mixing){
    if(anyNA(draws)) stop("Random draws have not been specified despite setting apollo_control$mixing==TRUE!")
    if(!is.function(apollo_randCoeff)) stop("apollo_randCoeff function has not been defined despite setting apollo_control$mixing==TRUE!")
    if("draws" %in% search()) detach("draws")
    attach(draws)
    randcoeff = apollo_randCoeff(apollo_beta, apollo_inputs)
    if("randcoeff" %in% search()) detach("randcoeff")
    attach(randcoeff)
  }

  if(is.function(apollo_lcPars)){
    lcpars = apollo_lcPars(apollo_beta, apollo_inputs)
    if("lcpars" %in% search()) detach("lcpars")
    attach(lcpars)
  }

}
