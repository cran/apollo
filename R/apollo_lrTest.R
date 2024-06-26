#' Likelihood ratio test
#' 
#' Calculates the likelihood ratio test value between two models and reports the corresponding p-value. 
#' 
#' The two models need to have been estimated on the same data, and one model needs to be nested within the other model.
#' 
#' @param model1 Either a character variable with the name of a previously estimated model, or an estimated 
#'               model in memory, as returned by \link{apollo_estimate}.
#' @param model2 Either a character variable with the name of a previously estimated model, or an estimated 
#'               model in memory, as returned by \link{apollo_estimate}.
#' @return LR-test p-value (invisibly)
#' @export
apollo_lrTest = function(model1, model2){
  
  modelNames = list()
  LL         = list()
  obs        = list()
  k          = list()
  inputs     = list(model1, model2)
  
  for(i in 1:2){
    modeluse=inputs[[i]]
    if(is.character(modeluse)){
      ### INPUT IS NAME OF MODEL ####
      
      ### Search for output directory
      outputDirectory <- ''
      tmp  <- get('apollo_control', envir=parent.frame(), inherits=FALSE)
      test <- is.list(tmp) && !is.null(tmp$outputDirectory) && is.character(tmp$outputDirectory)
      if(test) outputDirectory <- tmp$outputDirectory else {
        tmp <- get('apollo_inputs', envir=parent.frame(), inherits=FALSE)
        test <- is.list(tmp) && !is.null(tmp$apollo_control) && !is.null(tmp$apollo_control$outputDirectory)
        test <- test && is.character(tmp$apollo_control$outputDirectory)
        if(test) outputDirectory <- tmp$apollo_control$outputDirectory
      }
      test <- outputDirectory!=''
      test <- test && !(substr(outputDirectory, nchar(outputDirectory), nchar(outputDirectory)) %in% c('/','\\'))
      if(test) outputDirectory <- paste0(outputDirectory, '/')
      
      ### Try to read file
      filename = paste0(outputDirectory, modeluse, "_output.txt", collapse='')
      if(!file.exists(filename)) filename = paste0(modeluse,"_output.txt", collapse='')
      if(!file.exists(filename)){
        if(outputDirectory=='') stop('INPUT ISSUE - File ', filename, ' not found in working directory.') else 
          stop('INPUT ISSUE - File ', filename, ' not found in working directory, nor in ', outputDirectory, '.') 
      } 
      lines = tryCatch(readLines(filename), 
                       warning=function(w) x=FALSE,
                       error=function(e) x=FALSE)
      if(is.logical(lines) && lines==FALSE) stop("INPUT ISSUE - Could not open file ",filename) 
      
      ### Read model name
      id <- grepl(paste0("Model name"), lines)
      value=lines[which(id)] 
      position=gregexpr(pattern=":",value)[[1]][1]
      modelNames[[i]]=(substr(value,position+2,nchar(value)))
      
      id <- grepl(paste0("Number of modelled outcomes"), lines)
      if(any(id)){
        value=lines[which(id)] 
      } else {
        stop("INPUT ISSUE - Number of observations not found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      obs[[i]]=as.double(substr(value,position+1,nchar(value)))

      id1 <- grepl("LL\\(final, whole model)", lines)
      id2 <- grepl("LL\\(final)", lines) 
      if(any(id1)){
        value=lines[which(id1)] 
      } else if(any(id2)){
        value=lines[which(id2)] 
      } else {
        stop("INPUT ISSUE - No final LL found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      LL[[i]]=as.double(substr(value,position+1,nchar(value)))

      id <- grepl("Estimated parameters", lines)
      if(any(id)){
        value=lines[which(id)] 
      } else {
        stop("INPUT ISSUE - Number of estimated parameters not found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      k[[i]]=as.double(substr(value,position+1,nchar(value)))

    } else {
      ### INPUT IS MODEL OBJECT ####
      
      modelNames[[i]]=modeluse$apollo_control$modelName
      if(is.null(modeluse$maximum)){
        stop("INPUT ISSUE - No LL found in ",paste0(modeluse))
      } else {
        LL[[i]]=modeluse$maximum  
      }
      #if(is.null(modeluse$nObs) || anyNA(modeluse$nObs[1])){
      #  stop("INPUT ISSUE - Number of observations not found in ",paste0(modeluse))
      #} else {
      #  obs[[i]]=modeluse$nObs[1]
      #} 
      if(is.null(modeluse$nObsTot) || anyNA(modeluse$nObsTot)){
        stop("INPUT ISSUE - Number of observations not found in ",paste0(modeluse))
      } else {
        obs[[i]]=sum(modeluse$nObsTot)
      }  
      nParams     <- length(modeluse$apollo_beta)
      nFreeParams <- nParams
      if(!is.null(modeluse$apollo_fixed)) nFreeParams <- nFreeParams - length(modeluse$apollo_fixed)
      k[[i]]=nFreeParams
    }
  }

  if(obs[[1]]!=obs[[2]]) stop("INCORRECT FUNCTION/SETTING USE - The two models to be compared were not estimated on the same number of observations. A likelihood ratio test cannot be used!")
  if(k[[1]]==k[[2]]) stop("INCORRECT FUNCTION/SETTING USE - The two models to be compared have the same number of parameters. A likelihood ratio test cannot be used!")
  if(((k[[1]]-k[[2]])*(LL[[1]]-LL[[2]]))<0) stop("INCORRECT FUNCTION/SETTING USE - The model with more parameters does not have a better log-likelihood. A likelihood ratio test cannot be used!")

  if(LL[[2]]<LL[[1]]){
    apollo_print(paste0("The order of your two models will be reversed in the output as model 1 has better fit than model 2."))
    apollo_print("\n")
    LL=LL[c(2,1)]
    k=k[c(2,1)]
    inputs=inputs[c(2,1)]
    modelNames=modelNames[c(2,1)]
  }
  
  output=matrix(0,nrow=3,ncol=2)
  output[1:2,1]=round(unlist(LL),2)
  output[1:2,2]=round(unlist(k),2)
  output[3,]=output[2,]-output[1,]
  colnames(output)=c("LL","par")
  rownames(output)=c(unlist(modelNames),"Difference")
  
  LR_test_value     = 2*output[3,1]
  df                = output[3,2]
  p=stats::pchisq(LR_test_value, df, lower.tail=FALSE)

  print(output)
  cat("\nLikelihood ratio test-value:   ",round(LR_test_value,2),"\n")
  cat("Degrees of freedom:            ",df,"\n")
  cat("Likelihood ratio test p-value: ",formatC(p),"\n")
  apollo_print("\nThe p-value from the test is returned invisibly as an output from this function. Calling the function via result=apollo_lrTest(...) will save this output in an object called result (or otherwise named object).", type="i")
  return(invisible(p))
}
