#' Prints package startup message
#' 
#' This function is only called by R when attaching the package.
#' 
#' @param libname Name of library.
#' @param pkgname Name of package.
#' @return Nothing
.onAttach <- function(libname, pkgname) {
  apolloVersion <- tryCatch(utils::packageDescription("apollo", fields = "Version"),
                            warning=function(w) return("alpha"),
                            error=function(e) return("alpha"))
  txt <- paste0("Apollo ", apolloVersion,
  "\nwww.ApolloChoiceModelling.com",
  "\nSee url for a detailed manual, examples and a help forum.",
  "\nSign up to our mailing list for updates on new releases.")
  packageStartupMessage(txt)
}