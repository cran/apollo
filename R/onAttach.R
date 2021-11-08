#' Prints package startup message
#' 
#' This function is only called by R when attaching the package.
#' 
#' @param libname Name of library.
#' @param pkgname Name of package.
#' @return Nothing
.onAttach <- function(libname, pkgname) {
  
  ### Welcome message
  apolloVersion <- tryCatch(utils::packageDescription("apollo", fields = "Version"),
                            warning=function(w) return("alpha"),
                            error=function(e) return("alpha"))
  txt <- paste0("Apollo ", apolloVersion,
  "\nwww.ApolloChoiceModelling.com",
  "\nSee url for a detailed manual, examples and a help forum.",
  "\nSign up to our mailing list for updates on new releases.")
  
  ### Warning if more than six months old
  releaseDate <- as.POSIXct('2021-11-01', format='%Y-%m-%d')
  isOld <- as.POSIXct(Sys.Date()) > (releaseDate + 6*30.4*24*60*60)
  if(isOld) txt <- paste0(txt, 
                          '\n\nYour version of Apollo is more than six months old.',
                          '\nUsing the latest version will ensure you have all',
                          '\n current functionality and bug fixes.', 
                          '\nYou can update to the latest version by typing:', 
                          '\n install.packages("apollo")')
  
  ### Print message
  packageStartupMessage(txt)
  if(isOld) Sys.sleep(5)
}