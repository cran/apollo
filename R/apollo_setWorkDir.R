#' Automatically sets working directory to active file directory
#' 
#' This function only works in Rstudio. If called outside RStudio
#' will just print a message to screen saying it could not set the 
#' working directory.
#' @return (invisibly) TRUE if it manages to set the working directory, 
#'         FALSE if not.
#' @export
#' @importFrom rstudioapi getActiveDocumentContext
apollo_setWorkDir <- function(){
  path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                   error=function(e) "")
  if(path==""){
    apollo_print("Could not set working directory, please 
                            do it manually.", type="w", pause=0)
    return(invisible(FALSE))
  } 
  if(path!="") setwd(dirname(path))
  return(invisible(TRUE))
}