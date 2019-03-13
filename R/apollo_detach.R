#' Detaches parameters and the database.
#'
#' Detaches variables attached by \link{apollo_attach}.
#'
#' This function detaches the variables attached by \link{apollo_attach}. 
#' It should be called at the end of \code{apollo_probabilities}, only if 
#' \link{apollo_attach} was called and the beginning. This can be achieved 
#' by adding the line \code{on.exit(apollo_detach(apollo_beta, apollo_inputs))} 
#' right after calling \link{apollo_attach}.
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
apollo_detach=function(apollo_beta, apollo_inputs){
  apollo_control=apollo_inputs[["apollo_control"]]
  
  if(apollo_control$mixing){
    detach(randcoeff)
    detach(draws)
  }
  detach("as.list(apollo_beta)")
  detach(database)
  if("lcpars" %in% search()){
    detach("lcpars")
  }
}
