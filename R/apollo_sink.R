#' Writes or stops writing output to a text file
#'
#' Starts or stops writing the output shown in the console to a file named "modelName_additional_output.txt".
#'
#' After the first time this function is called, all output shown in the console will also be written to a
#' text file called "modelName_additional_output.txt", where "modelName" is the modelName set inside
#' apollo_control. 
#' The second time this function is called, it stops writing the console output to the file. The user 
#' should always call this function an even number of times to close the output file and prevents data loss.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#'                      If not provided, it will be looked for in the global environment.
#' @return Nothing.
#' @export
apollo_sink <- function(apollo_inputs=NULL){
  
  # Checks if sinks are open, if so closes them and exists
  if(sink.number()>0){
    sink()
    apollo_print(paste0('Output is no longer being written to file.'))
    return(invisible(NULL))
  } 
  
  # Try to fetch apollo_inputs if it is missing
  if(is.null(apollo_inputs)) apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                                                       error=function(e) NULL)
  if(is.null(apollo_inputs)) stop('apollo_inputs could not be retrieved')
  
  # Check apollo_inputs has everything needed
  test <- is.list(apollo_inputs) && !is.null(apollo_inputs$apollo_control)
  test <- test && !is.null(apollo_inputs$apollo_control$outputDirectory) && is.character(apollo_inputs$apollo_control$outputDirectory)
  test <- test && !is.null(apollo_inputs$apollo_control$modelName) && is.character(apollo_inputs$apollo_control$modelName)
  if(!test) stop('apollo_inputs is corrupted')
  
  # Do sink()
  f <- paste0(apollo_inputs$apollo_control$outputDirectory, 
              apollo_inputs$apollo_control$modelName,
              "_additional_output.txt")
  if(!apollo_inputs$silent) apollo_print(paste0('Writing output to file ', f, '. ', 
                                                'Please run "apollo_sink()" again after finishing writing results.'))
  sink(f, split=TRUE)
}