#' Reports market share for subsamples
#'
#' Compares market shares across subsamples in dataset, and writes results to a file.
#'
#' Saves the output to a csv file in the working directory.
#' @param choiceAnalysis_settings List containing settings for this function. The settings must be:
#'                                \itemize{
#'                                  \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                                  \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                                  \item choiceVar: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                                  \item explanators: data.frame. Variables determining subsamples of the database. Values in each column must describe a group or groups of individuals (e.g. socio-demographics). Most usually a subset of columns from database.
#'                                }
#' @param apollo_control List. Options controlling the running of the code. See \link{apollo_validateInputs}.
#' @param database data.frame. Data used by model.
#' @return nothing, but prints the output to screen and writes a csv file to the working directory.
#' @export
apollo_choiceAnalysis=function(choiceAnalysis_settings, apollo_control, database){
  tmp <- c("alternatives", "avail", "choiceVar", "explanators")
  for(i in tmp) if(is.null(choiceAnalysis_settings[[i]])) stop("The choiceAnalysis_settings list needs to include an object called \"",i,"\"!")

  ### added these to run checks as this function can now run before validateInputs
  apollo_control <- apollo_validateControl(database, apollo_control, silent=TRUE)
  database       <- apollo_validateData(database, apollo_control, silent=TRUE)
  
  ### changed next line to not use apollo_inputs
  modelName   = apollo_control$modelName
  alternatives= choiceAnalysis_settings[["alternatives"]]
  avail       = choiceAnalysis_settings[["avail"]]
  choiceVar   = choiceAnalysis_settings[["choiceVar"]]
  explanators = choiceAnalysis_settings[["explanators"]]
  if(is.vector(explanators)) explanators <- data.frame(explanators)

  ### dropped this
  #database=apollo_inputs[["database"]]

  output=matrix(0,nrow=length(alternatives),ncol=ncol(explanators)*3)
  rownames(output)=names(alternatives)
  
  
  outputnames=c(rep(0,ncol(output)))
  s=1
  while(s<=ncol(explanators)){
    outputnames[(s-1)*3+1]=paste("Mean for",colnames(explanators)[s],"if chosen")
    outputnames[(s-1)*3+2]=paste("Mean for",colnames(explanators)[s],"if not chosen")
    outputnames[(s-1)*3+3]=paste("p-val for difference")
    s=s+1
  }
  
  colnames(output)=outputnames
  
  j=1
  s=1

  if(length(avail)==1 && avail==1){
    avail <- rep(1, length(alternatives))
    names(avail) <- names(alternatives)
    avail <- as.list(avail)
  }

  while(j<=length(alternatives)){
    database_sub=subset(database,avail[[j]]==1)
    explanators_sub=subset(explanators,avail[[j]]==1)
    chosen=subset(choiceVar,avail[[j]]==1)==alternatives[j]
    s=1
    while(s<=ncol(explanators)){
      x=tapply(explanators_sub[,s],chosen, mean,na.rm=TRUE)
      output[j,((s-1)*3+1)]=x[2]
      output[j,((s-1)*3+2)]=x[1]
      if(x[1]==x[2]){
        output[j,((s-1)*3+3)]=NA
      }else{
        output[j,((s-1)*3+3)]=stats::t.test(subset(explanators_sub[,s],chosen==0),subset(explanators_sub[,s],chosen==1))[["p.value"]]
      }
      s=s+1}
    j=j+1}

  filename=paste(modelName,"_choiceAnalysis.csv",sep="")

  utils::write.csv(round(output,2), filename)
  cat("Ouputs of apollo_choiceAnalysis saved to ",filename,"\n",sep="")
}


