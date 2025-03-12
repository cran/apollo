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
  
  ### WARNING if more than six months old
  releaseDate <- as.POSIXct('2025-03-15', format='%Y-%m-%d')
  isOld <- as.POSIXct(Sys.Date()) > (releaseDate + 6*30.4*24*60*60)
  #if(isOld) txt2 <- paste0(txt, 
  #if(isOld) txt2 <- paste0('\n\nYour version of Apollo is more than six months old.',
  if(isOld) txt2 <- paste0('\nYour version of Apollo is more than six months old.',
                           '\nUsing the latest version will ensure you have all',
                           '\n current functionality and bug fixes.', 
                           '\nYou can update to the latest version by typing:', 
                           '\n install.packages("apollo")')
  
  if(isOld){#cat("\n")
    apollo_print("\n")
    apollo_print(txt2, pause=0, type="i")
  } 
  
  cat("Apollo ignition sequence completed\n")
}
