#' Prepares environment
#' 
#' Prepares environment (the global environment if called by the user) for model definition and estimation.
#' 
#' This function detaches variables and makes sure that output is directed to console. It does not delete variables from the working environment.
#' 
#' @return Nothing.
#' @export
apollo_initialise <- function()
{
  doDetach <- !grepl("^(.GlobalEnv|package:|tools:|Autoloads|CheckExEnv|TempEnv)", search())
  doDetach <- (search())[which(doDetach)]
  if(length(doDetach)>0) for(i in 1:length(doDetach)) detach(pos=(which(doDetach[i]==search()))[1])
  if(sink.number()>0) sink()
  cat("Apollo ignition sequence completed\n")
}
