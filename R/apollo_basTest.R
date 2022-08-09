#' Ben-Akiva & Swait test
#' 
#' Carries out the Ben-Akiva & Swait test for non-nested models and reports the corresponding p-value. 
#' 
#' The two models need to both be discrete choice, and need to have been estimated on the same data.
#' 
#' @param model1 Either a character variable with the name (and possibly path) of a previously estimated model, or an estimated 
#'               model in memory, as returned by \link{apollo_estimate}.
#' @param model2 Either a character variable with the name (and possibly path) of a previously estimated model, or an estimated 
#'               model in memory, as returned by \link{apollo_estimate}.
#' @return Ben-Akiva & Swait test p-value (invisibly)
#' @export
apollo_basTest = function(model1,model2){
  
  modelNames = list()
  LL         = list()
  LL0        = list()
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
        if(outputDirectory=='') stop('File ', filename, ' not found in working directory.') else 
          stop('File ', filename, ' not found in working directory, nor in ', outputDirectory, '.') 
      } 
      lines = tryCatch(readLines(filename), 
                       warning=function(w) x=FALSE,
                       error=function(e) x=FALSE)
      if(is.logical(lines) && lines==FALSE) stop("Could not open file ",filename) 
      
      ### Read file
      id <- grepl(paste0("Model name"), lines)
      value=lines[which(id)] 
      position=gregexpr(pattern=":",value)[[1]][1]
      modelNames[[i]]=(substr(value,position+2,nchar(value)))
      
      #id <- grepl("Adj.Rho-square (0)", lines, fixed=TRUE)
      id <- grepl("Adj.Rho-squared vs equal shares", lines, fixed=TRUE)
      if(any(id)){
        value=lines[which(id)] 
      } else {
        stop("No adjusted rho2 found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      if(((substr(value,position+1,nchar(value))))==" Not applicable") stop("No adjusted rho2 found in ",filename)

      id1 <- grepl("LL\\(final, whole model)", lines)
      id2 <- grepl("LL\\(final)", lines) 
      if(any(id1)){
        value=lines[which(id1)] 
      } else if(any(id2)){
        value=lines[which(id2)] 
      } else {
        stop("No final LL found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      LL[[i]]=as.double(substr(value,position+1,nchar(value)))

      id1 <- grepl("LL\\(0, whole model)", lines)
      id2 <- grepl("LL\\(0)", lines) 
      if(any(id1)){
        value=lines[which(id1)] 
      } else if(any(id2)){
        value=lines[which(id2)] 
      } else {
        stop("No LL(0) found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      LL0[[i]]=as.double(substr(value,position+1,nchar(value)))

      id <- grepl("Estimated parameters", lines)
      if(any(id)){
        value=lines[which(id)] 
      } else {
        stop("Number of estimated parameters not found in ",filename)
      }
      position=gregexpr(pattern=":",value)[[1]][1]
      k[[i]]=as.double(substr(value,position+1,nchar(value)))

    } else {
      modelNames[[i]]=modeluse$apollo_control$modelName
      test <- !is.null(modeluse$modelTypeList) && all(tolower(modeluse$modelTypeList) %in% c("mnl", "nl", "cnl", "el", "dft", "lc", "rrm"))
      if(!test) stop("Adjusted rho2 calculation not possible for model ",paste0(modeluse$apollo_control$modelName))
      if(is.null(modeluse$maximum)){
        stop("No LL found in ",paste0(modeluse$apollo_control$modelName))
      } else {
        LL[[i]]=modeluse$maximum  
      }
      if(is.null(modeluse$LL0) || anyNA(modeluse$LL0[1])){
        stop("No LL(0) found in ",paste0(modeluse$apollo_control$modelName))
      } else {
        LL0[[i]]=modeluse$LL0[1]
      }      
      nParams     <- length(modeluse$apollo_beta)
      nFreeParams <- nParams
      if(!is.null(modeluse$apollo_fixed)) nFreeParams <- nFreeParams - length(modeluse$apollo_fixed)
      k[[i]]=nFreeParams
    }
  }

  if(round(LL0[[1]],2)!=round(LL0[[2]],2)) stop("The two models to be compared do not have the same LL(0). This suggests they were not estimated on the same data!")
  
  adjrho=list(1-(LL[[1]]-k[[1]])/LL0[[1]],
              1-(LL[[2]]-k[[2]])/LL0[[2]])
  
  flag_reverse=FALSE
  if(adjrho[[2]]<adjrho[[1]]){
    flag_reverse=TRUE
    adjrho=adjrho[c(2,1)]
    LL=LL[c(2,1)]
    k=k[c(2,1)]
    inputs=inputs[c(2,1)]
    modelNames=modelNames[c(2,1)]
  }
  
  output=matrix(0,nrow=3,ncol=4)
  output[1:2,1]=round(unlist(LL0),2)
  output[1:2,2]=round(unlist(LL),2)
  output[1:2,3]=round(unlist(k),2)
  output[1,4]=round(1-(LL[[1]]-k[[1]])/LL0[[1]],4)
  output[2,4]=round(1-(LL[[2]]-k[[2]])/LL0[[2]],4)
  output[3,]=output[2,]-output[1,]
  colnames(output)=c("LL0","LL","par","adj.rho2")
  rownames(output)=c(unlist(modelNames),"Difference")
  
  if(output[3,3]<0) stop("Test not run as better fitting model has fewer parameters!")
  p=pnorm(-sqrt(-2*output[3,4]*output[1,1]+output[3,3]))
  
  if(flag_reverse){
    apollo_print(paste0("The order of your two models will be reversed in the output as model 1 has better fit than model 2."))
    apollo_print("\n")}

  print(output)
  apollo_print("\n")
  apollo_print(paste0("p-value for Ben-Akiva & Swait test: ",formatC(p)))
  
  return(invisible(p))
}
