#' Reports market share for subsamples
#'
#' Compares market shares across subsamples in dataset, and conducts statistical tests.
#'
#' Saves the output to a csv file in the working/output directory.
#' @param choiceAnalysis_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                                \itemize{
#'                                  \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}. Note that these need not necessarily be the alternatives as defined in the model, but could e.g. relate to cheapest/most expensive.
#'                                  \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                                  \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                                  \item \strong{\code{explanators}}: data.frame. Variables determining subsamples of the database. Values in each column must describe a group or groups of individuals (e.g. socio-demographics). Most usually a subset of columns from the database.
#'                                  \item \strong{\code{printToScreen}}: Logical. TRUE for returning output to screen as well as file. TRUE by default.
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                }
#' @param apollo_control List. Options controlling the running of the code. See \link{apollo_validateInputs}.
#' @param database data.frame. Data used by model.
#' @return Silently returns a matrix containing the mean value for each explanator for those cases where an alternative is chosen and where it is not chosen, 
#'         as well as the t-test comparing those means (H0: equivalence).
#'         The table is also written to a file called \code{modelName_choiceAnalysis.csv} and printed to screen.
#' @export
apollo_choiceAnalysis=function(choiceAnalysis_settings, apollo_control, database){
  if(is.null(choiceAnalysis_settings[["printToScreen"]])) choiceAnalysis_settings[["printToScreen"]]=TRUE
  if(is.null(choiceAnalysis_settings[["rows"]])) choiceAnalysis_settings[["rows"]]="all"
  tmp <- c("alternatives", "choiceVar", "explanators")
  for(i in tmp) if(is.null(choiceAnalysis_settings[[i]])) stop("SYNTAX ISSUE - The choiceAnalysis_settings list needs to include an object called \"",i,"\"!")
  
  ### Validate control & data as this function can be run before validateInputs
  apollo_control <- apollo_validateControl(database, apollo_control, silent=TRUE)
  database       <- apollo_validateData(database, apollo_control, silent=TRUE)
  
  ### Extract useful values
  modelName   = apollo_control$modelName
  alternatives= choiceAnalysis_settings[["alternatives"]]
  if(is.null(choiceAnalysis_settings[["avail"]])){
    choiceAnalysis_settings[["avail"]]=1
    apollo_print("Availabilities not provided for 'apollo_choiceAnalysis', so full availability is assumed.")
  } 
  avail       = choiceAnalysis_settings[["avail"]]
  choiceVar   = choiceAnalysis_settings[["choiceVar"]]
  explanators = choiceAnalysis_settings[["explanators"]]
  if(is.vector(explanators)) explanators = data.frame(explanators)
  rows         = choiceAnalysis_settings[["rows"]]
  
  ### Filter rows
  if(!is.vector(rows)) stop("SYNTAX ISSUE - Argument 'rows', when provided, must be a boolean vector or the words 'all'.")
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nrow(database))
  if(any(!rows)){
    database <- database[rows,]
    #avail <- lapply(avail,function(x) x[rows])
    #if(is.list(avail)) for(j in alternatives) if(length(avail[[j]])!=1) avail[[j]]=avail[[j]][rows]
    if(is.list(avail)) for(j in 1:length(alternatives)) if(length(avail[[j]])!=1) avail[[j]]=avail[[j]][rows]
    choiceVar  <- choiceVar[rows] 
    explanators <- explanators[rows,] 
  }
  
  ### Make sure there are no strange values in choiceVar
  if(!all(choiceVar %in% alternatives)) stop('SYNTAX ISSUE - Some values in "choiceVar" are not defined in "alternatives".')
  
  ### Check for unchosen and unavailable choices
  #for(j in alternatives){
  for(j in 1:length(alternatives)){
    if(sum(choiceVar==alternatives[j])==0) apollo_print(paste0("WARNING: Alternative ",names(alternatives)[j]," never chosen in your data."))
    if(is.list(avail)&&(length(avail[[j]])!=1)){
     if(any((choiceVar==alternatives[j])&(avail[[j]]==0))) stop(paste0("INPUT ISSUE - Your data contains rows where alternative ",names(alternatives)[j]," is chosen when not available!"))
    }
  }
  
  ### Create and matrix of outputs and sets its col and row names
  output           = matrix(0,nrow=length(alternatives),ncol=ncol(explanators)*3)
  rownames(output) = names(alternatives)
  outputnames=c(rep(0,ncol(output)))
  for(s in 1:ncol(explanators)){
    outputnames[(s-1)*3+1] = paste0("Explanator ",s," (",colnames(explanators)[s],"), mean when alt is chosen:")
    outputnames[(s-1)*3+2] = paste0("Explanator ",s," (",colnames(explanators)[s],"), mean when alt is not chosen:")
    outputnames[(s-1)*3+3] = paste0("Explanator ",s," (",colnames(explanators)[s],"), t-test (mean if chosen - mean if not chosen)")
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
      if(length(x)==1){
       if(names(x)=="TRUE"){
        x=c(nFALSE=0,nTRUE=x)
        }else{
        x=c(nFALSE=x,nTRUE=0)
      }}
       output[j,((s-1)*3+1)] = x[2]
      output[j,((s-1)*3+2)] = x[1]
      if(x[1]==x[2]){
        output[j,((s-1)*3+3)] = NA
      }else{
        if(length(subset(explanators_sub[,s], !chosen))>1 && length(subset(explanators_sub[,s], chosen))>1){
          output[j,((s-1)*3+3)] = stats::t.test(subset(explanators_sub[,s], chosen),subset(explanators_sub[,s], !chosen))[["statistic"]]
        } else {
          output[j,((s-1)*3+3)] = NA
        }
      }
    }
  }
  
  ### Determine name of output file
  filename = paste(modelName,"_choiceAnalysis.csv",sep="")
  test <- !is.null(apollo_control$outputDirectory) && is.character(apollo_control$outputDirectory)
  test <- test && apollo_control$outputDirectory!=''
  if(test){
    n <- nchar(apollo_control$outputDirectory)
    tmp <- ''
    if(!(substr(apollo_control$outputDirectory, n, n) %in% c('/','\\'))) tmp <- '/'
    filename <- paste0(apollo_control$outputDirectory, tmp, filename)
    rm(n, tmp)
  }; rm(test)
  
  output_new=c()
  for(s in 1:ncol(explanators)){
    output_new=cbind(output_new,round(output[,((s-1)*3+1)],4))
    output_new=cbind(output_new,round(output[,((s-1)*3+2)],4))
    output_new=cbind(output_new,round(output[,((s-1)*3+3)],2))
  }
  
  ### Print output to screen, if requested
  colnames(output_new)=colnames(output)
  if(choiceAnalysis_settings$printToScreen){
    for(s in 1:ncol(explanators)){
      print(t(output_new[,(1+(s-1)*3):(s*3)]))
      cat("\n")
    } 
  } 
  
  ### Write file
  utils::write.csv(t(output_new), filename)
  cat("Ouputs of apollo_choiceAnalysis saved to ",filename,"\n",sep="")
  invisible(t(output_new))
  
}


