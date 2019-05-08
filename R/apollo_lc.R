#' Calculates the likelihood of a latent class model
#' 
#' Using the conditional likelihoods of each latent class, as well as teir classification probabilities, calculate the weighted likelihood of the whole model.
#' 
#' @param lc_settings List of arguments used by \code{apollo_lc}. It must include the following.
#'                  \itemize{
#'                    \item inClassProb: List of probabilities. Conditional likelihood for each class. One element per class, in the same order as \code{classProb}.
#'                    \item classProb: List of probabilities. Allocation probability for each class. One element per class, in the same order as \code{inClassProb}.
#'                  }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate" Used for model estimation.
#'                        \item "prediction" Used for model predictions.
#'                        \item "validate" Used for validating input.
#'                        \item "zero_LL" Used for calculating null likelihood.
#'                        \item "conditionals" Used for calculating conditionals.
#'                        \item "output" Used for preparing output after model estimation.
#'                        \item "raw" Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @examples
#' data(apollo_modeChoiceData)
#' database <- apollo_modeChoiceData
#' rm(apollo_modeChoiceData)
#' N <- nrow(database)
#' lc_settings <- list(inClassProb=list(rnorm(N), rnorm(N)),
#'                     classProb=list(stats::runif(N), stats::runif(N)))
#' apollo_control <- list(indivID="ID")
#' x <- apollo_lc(lc_settings, apollo_control, functionality="estimate")
#' summary(x)
#' @export
apollo_lc <- function(lc_settings, apollo_inputs, functionality){
  if(is.null(lc_settings[["inClassProb"]])) stop("The lc_settings list needs to include an object called \"inClassProb\"!")
  if(is.null(lc_settings[["classProb"]])) stop("The lc_settings list needs to include an object called \"classProb\"!")
  
  inClassProb=lc_settings[["inClassProb"]]
  classProb=lc_settings[["classProb"]]
  
  apollo_control=apollo_inputs[["apollo_control"]]
  
  # Count number of rows in classProb
  nRowsClassProb <- 0
  if(is.vector(classProb[[1]])) nRowsClassProb <- length(classProb[[1]]) else {
    nRowsClassProb <- dim(classProb[[1]])[1]
  }
  
  if(functionality=="validate"){
    if(length(inClassProb)!=length(classProb)) stop("Arguments 'inClassProb' and 'classProb' must have the same length.")
    class_summary=matrix(0,nrow=length(classProb),ncol=1)
    rownames(class_summary)=paste0("class_", 1:length(classProb))
    colnames(class_summary)=c("Mean prob.")
    c=1
    while(c<=length(classProb)){
      class_summary[c,1]=mean(classProb[[c]])
      c=c+1
    }
    #cat("\nAll checks passed for latent class structure.")
    #cat("\nClass allocation summary for LC before estimation:\n")
    #print(round(class_summary,2))
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Class allocation summary for LC before estimation:", round(class_summary,2), apolloLog)
    return(TRUE)
  } 
  
  if(functionality %in% c("estimate", "conditionals", "output", "zero_LL")){
    
    # Count number of rows in inClassProb
    nRowsInClassProb <- 0
    if(is.vector(inClassProb[[1]])) nRowsInClassProb <- length(inClassProb[[1]]) else {
      nRowsInClassProb <- dim(inClassProb[[1]])[1]
    }
    
    # Match the number of rows in each list
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- as.vector(table(indivID))
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        classProb <- lapply(classProb, rowsum, group=indivID)
        classProb <- lapply(classProb, '/', nObsPerIndiv)
        warning('Class probability averaged across observations of each individual.')
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        classProb <- lapply(classProb, rep, nObsPerIndiv)
        ### warning('Class probability repeated for all observations of the same individual.')
      }
    }
    
    # Returns inClassProb*classProb
    Pout <- mapply(function(p,pi) pi*p, inClassProb, classProb, SIMPLIFY=FALSE)
    Pout <- Reduce('+', Pout)
    if(functionality=="output"){
      #if(!apollo_control$noDiagnostics){
      class_summary=matrix(0,nrow=length(classProb),ncol=1)
      rownames(class_summary)=paste0("class_", 1:length(classProb))
      colnames(class_summary)=c("Mean prob.")
      c=1
      while(c<=length(classProb)){
        class_summary[c,1]=mean(classProb[[c]])
        c=c+1
      }
      
      ## Write class allocation summary to temp file
      #on.exit(sink())
      #fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
      #fileName <- file.path(tempdir(),fileName)
      #fileConn <- tryCatch( file(fileName, open="at"),
      #                      error=function(e){
      #                        cat('apollo_lc could not write diagnostics to temporary file. No diagnostics in output.\n')
      #                        return(NA)
      #                      })
      #if(!anyNA(fileConn)){
      #  sink(fileConn)
      #  on.exit({if(sink.number()>0) sink(); close(fileConn)})
      #  cat('Class allocation summary for LC after estimation:\n')
      #  print(round(class_summary,2))
      #  cat("\n")
      #}
      
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Class allocation summary for LC after estimation:", round(class_summary,2), apolloLog)
    }
    return( Pout )
  }
  
  if((functionality=="prediction")|(functionality=="raw")){
    
    # Count number of rows in inClassProb
    nRowsInClassProb <- 0
    if(is.vector(inClassProb[[1]][[1]])) nRowsInClassProb <- length(inClassProb[[1]][[1]]) else {
      nRowsInClassProb <- dim(inClassProb[[1]][[1]])[1]
    }
    
    # Match the number of rows in each list
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb  & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- as.vector(table(indivID))
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        stop('Class-probability variable has more elements than in-class-probability.')
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        S=length(classProb)
        s=1
        while(s<=S){
          isMat <- is.matrix(classProb[[s]])
          isCub <- is.array(classProb[[s]]) && !isMat && length(dim(classProb[[s]]))==3
          if(isCub){
            # first average out over intra
            classProb[[s]]=colSums(aperm(classProb[[s]], perm=c(3,1,2)))/dim(classProb[[s]])[3]
          } 
          # now check what's left
          isVec <- is.vector(classProb[[s]])
          isMat <- is.matrix(classProb[[s]])
          if(isVec) classProb[[s]]=rep(classProb[[s]],nObsPerIndiv)
          if(isMat){
            tmp <- matrix(0, nrow=sum(nObsPerIndiv), ncol=ncol(classProb[[s]]))
            for(n in 1:length(nObsPerIndiv)){
              a <- ifelse(n==1, 1, sum(nObsPerIndiv[1:(n-1)]) + 1)
              b <- a + nObsPerIndiv[n] - 1 
              tmp[a:b,] <- rep(classProb[[s]][n,], each=nObsPerIndiv[n])
            }
            classProb[[s]] <- tmp
          } 
          s=s+1
        }
        
      }
    }
    
    nClass<- length(classProb)
    nAlts <- length(inClassProb[[1]])
    Pout <- vector(mode="list", length=nAlts)
    for(i in 1:nAlts){
      Pout[[i]] <- 0*inClassProb[[1]][[1]]
      for(k in 1:nClass) Pout[[i]] <- Pout[[i]] + inClassProb[[k]][[i]]*classProb[[k]]
    }
    return(Pout)
    
  }
}