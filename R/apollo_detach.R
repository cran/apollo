#' Detaches parameters and the database.
#'
#' Detaches variables attached by \link{apollo_attach}.
#'
#' This function detaches the variables attached by \link{apollo_attach}. 
#' It should be called at the end of \code{apollo_probabilities}, only if 
#' \link{apollo_attach} was called and the beginning. This can be achieved 
#' by adding the line \code{on.exit(apollo_detach(apollo_beta, apollo_inputs))} 
#' right after calling \link{apollo_attach}.
#' This function can also be called without any arguments, i.e. \code{apollo_detach()}.
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
apollo_detach=function(apollo_beta=NA, apollo_inputs=NA){
  ### Detach things if necessary
  tmp <- c("database", "as.list(apollo_beta)", "randcoeff", "draws", "lcpars")
  tmp <- tmp[tmp %in% search()]
  if(length(tmp)>0) for(i in tmp) detach(i, character.only=TRUE)
}
