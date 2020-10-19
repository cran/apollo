#' Reports market share for subsamples
#'
#' Compares market shares across subsamples in dataset, and writes results to a file.
#'
#' Saves the output to a csv file in the working directory.
#' @param choiceAnalysis_settings List containing settings for this function. The settings must be:
#'                                \itemize{
#'                                  \item \strong{alternatives}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                                  \item \strong{avail}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                                  \item \strong{choiceVar}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                                  \item \strong{explanators}: data.frame. Variables determining subsamples of the database. Values in each column must describe a group or groups of individuals (e.g. socio-demographics). Most usually a subset of columns from database.
#'                                  \item \strong{rows}: Boolean vector. Consideration of rows to include, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                }
#' @param apollo_control List. Options controlling the running of the code. See \link{apollo_validateInputs}.
#' @param database data.frame. Data used by model.
#' @return Silently returns a matrix containg the mean ehen chosen and un chose for each explanator, 
#'         as well as the t-test comparing those means (H0: equivalence).
#'         The table is also writen to a file called \code{modelName_choiceAnalysis.csv}.
#' @export
apollo_choiceAnalysis=function(choiceAnalysis_settings, apollo_control, database){
  if(is.null(choiceAnalysis_settings[["rows"]])) choiceAnalysis_settings[["rows"]]="all"
  tmp <- c("alternatives", "avail", "choiceVar", "explanators")
  for(i in tmp) if(is.null(choiceAnalysis_settings[[i]])) stop("The choiceAnalysis_settings list needs to include an object called \"",i,"\"!")
  
  ### Validate control & data as this function can be run before validateInputs
  apollo_control <- apollo_validateControl(database, apollo_control, silent=TRUE)
  database       <- apollo_validateData(database, apollo_control, silent=TRUE)
  
  ### Extract useful values
  modelName   = apollo_control$modelName
  alternatives= choiceAnalysis_settings[["alternatives"]]
  avail       = choiceAnalysis_settings[["avail"]]
  choiceVar   = choiceAnalysis_settings[["choiceVar"]]
  explanators = choiceAnalysis_settings[["explanators"]]
  if(is.vector(explanators)) explanators = data.frame(explanators)
  rows         = choiceAnalysis_settings[["rows"]]
  
  ### Filter rows
  if(!is.vector(rows)) stop("Argument 'rows', when provided, must be a boolean vector or the words 'all'.")
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nrow(database))
  if(any(!rows)){
    database <- database[rows,]
    avail <- lapply(avail,function(x) x[rows])
    choiceVar  <- choiceVar[rows] 
    explanators <- explanators[rows,] 
  }
  
  
  ### Create and matrix of outputs and sets its col and row names
  output           = matrix(0,nrow=length(alternatives),ncol=ncol(explanators)*3)
  rownames(output) = names(alternatives)
  outputnames=c(rep(0,ncol(output)))
  for(s in 1:ncol(explanators)){
    outputnames[(s-1)*3+1] = paste("Mean for",colnames(explanators)[s],"if chosen")
    outputnames[(s-1)*3+2] = paste("Mean for",colnames(explanators)[s],"if not chosen")
    outputnames[(s-1)*3+3] = paste("t-test for difference")
  }
  colnames(output)=outputnames
  
  ### Expand availability if needed
  if(length(avail)==1 && avail==1){
    avail <- rep(1, length(alternatives))
    names(avail) <- names(alternatives)
    avail <- as.list(avail)
  }
  
  ### Populate output matrix
  for(j in 1:length(alternatives)){
    r               = avail[[j]]==1
    database_sub    = subset(database, r)
    explanators_sub = subset(explanators, r)
    chosen          = subset(choiceVar, r)==alternatives[j]
    for(s in 1:ncol(explanators)){
      x = tapply(explanators_sub[,s],chosen, mean,na.rm=TRUE)
      output[j,((s-1)*3+1)] = x[2]
      output[j,((s-1)*3+2)] = x[1]
      if(x[1]==x[2]){
        output[j,((s-1)*3+3)] = NA
      }else{
        if(length(subset(explanators_sub[,s], !chosen))>1 && length(subset(explanators_sub[,s], chosen))>1){
          output[j,((s-1)*3+3)] = stats::t.test(subset(explanators_sub[,s], !chosen),subset(explanators_sub[,s], chosen))[["statistic"]]
        } else {
          output[j,((s-1)*3+3)] = NA
        }
        
      }
    }
    }
  
  filename = paste(modelName,"_choiceAnalysis.csv",sep="")
  
  output_new=c()
  for(s in 1:ncol(explanators)){
    output_new=cbind(output_new,round(output[,((s-1)*3+1)],4))
    output_new=cbind(output_new,round(output[,((s-1)*3+2)],4))
    output_new=cbind(output_new,round(output[,((s-1)*3+3)],2))
  }
  
  colnames(output_new)=colnames(output)
  utils::write.csv(output_new, filename)
  cat("Ouputs of apollo_choiceAnalysis saved to ",filename,"\n",sep="")
  invisible(output_new)
}


