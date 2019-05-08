#' Write model results to file
#' 
#' Writes results from various models to a single CSV file.
#' 
#' @param combineResults_settings List of options. It can include the following.
#'                                \itemize{
#'                                  \item modelNames: Character vector. List of names of models to combine. Use an empty vector to combine results from all models in the directory.
#'                                  \item printClassical: Boolean. TRUE for printing classical standard errors. TRUE by default.
#'                                  \item printPVal: Boolean. TRUE for printing p-values. FALSE by default.
#'                                  \item printT1: Boolean. If TRUE, t-test for H0: apollo_beta=1 are printed. FALSE by default.
#'                                  \item estimateDigits: Numeric scalar. Number of decimal places to print for estimates. Default is 4.
#'                                  \item tDigits: Numeric scalar. Number of decimal places to print for t-ratios values. Default is 2.
#'                                  \item pDigits: Numeric scalar. Number of decimal places to print for p-values. Default is 2.
#'                                }
#' @return Nothing, but writes a file called 'model_comparison|_[date].csv' in the working directory.
#' @export
apollo_combineResults = function(combineResults_settings=NULL){
   if(is.null(combineResults_settings[["modelNames"]])){
     cat("The combineResults_settings does not include an object called \"modelNames\" so all files will be used!\n")
     combineResults_settings[["modelNames"]]=list.files(pattern="*estimates.csv")
     j=1
     while(j<=length(combineResults_settings[["modelNames"]])){
       l=nchar(combineResults_settings[["modelNames"]][j])
       l=l-14  
       combineResults_settings[["modelNames"]][j]=substr(combineResults_settings[["modelNames"]][j],1,l)
       j=j+1
     } }
   if(is.null(combineResults_settings[["printClassical"]])) combineResults_settings[["printClassical"]]=FALSE
   if(is.null(combineResults_settings[["printPVal"]])) combineResults_settings[["printPVal"]]=FALSE
   if(is.null(combineResults_settings[["printT1"]])) combineResults_settings[["printT1"]]=FALSE
   if(is.null(combineResults_settings[["estimateDigits"]])) combineResults_settings[["estimateDigits"]]=4
   if(is.null(combineResults_settings[["tDigits"]])) combineResults_settings[["tDigits"]]=2
   if(is.null(combineResults_settings[["pDigits"]])) combineResults_settings[["pDigits"]]=2

   modelNames=combineResults_settings[["modelNames"]]
   printClassical=combineResults_settings[["printClassical"]]
   printPVal=combineResults_settings[["printPVal"]]
   printT1=combineResults_settings[["printT1"]]
   estimateDigits=combineResults_settings[["estimateDigits"]]
   tDigits=combineResults_settings[["tDigits"]]
   pDigits=combineResults_settings[["pDigits"]]
### END NEW LINES

  outputfile_check=paste(modelNames,"_output.txt",sep="")
  estimatefile_check=paste(modelNames,"_estimates.csv",sep="")
  files=list.files()
  if(!(all(outputfile_check%in%files))) stop("Model output files missing for some models!")
  if(!(all(estimatefile_check%in%files))) stop("Model estimate files missing for some models!")     

  estimateDigits = max(1,estimateDigits)
  tDigits        = max(1,tDigits)
  pDigits        = max(1,pDigits)
  
  if(!is.character(modelNames)) stop("Argument 'modelNames' must be a character vector.")
  
  estimates    = list()
  otheroutputs = data.frame(matrix(0,nrow=8,ncol=length(modelNames)))
  rownames(otheroutputs) = c("Model name","Number of individuals","Number of observations","Estimated parameters","LL(final)","Adj.Rho-square (0)","AIC","BIC")
  values = 1 + ( 1 + printClassical ) * ( 1 + printT1) * ( 1 + printPVal )
  
  j=1
  while(j<=length(modelNames)){
    inputs=utils::read.csv(paste(modelNames[[j]],"_estimates.csv",sep=""))
    
    if(printClassical==TRUE){
      if(!("t.ratio.0." %in% colnames(inputs))){
        warning("Classical t.ratios not available in ",paste(modelNames[[j]],"_estimates.csv",sep=""))
        inputs[,"t.ratio.0."]=NA
      }}
    if(printT1==TRUE){
      if(!("Rob.t.ratio.1." %in% colnames(inputs))){
        warning("t.ratios against 1 not available in ",paste(modelNames[[j]],"_estimates.csv",sep=""))
        inputs[,"Rob.t.ratio.1."]=NA
        if(printClassical==TRUE) inputs[,"t.ratio.1."]=NA
      }}
    if(printPVal==TRUE){
      if(!("Rob.p.val.0." %in% colnames(inputs))){
        warning("p-values not available in ",paste(modelNames[[j]],"_estimates.csv",sep=""))
        inputs[,"Rob.p.val.0."]=NA
        if(printT1==TRUE) inputs[,"Rob.p.val.1."]=NA
        if(printClassical==TRUE){
          inputs[,"p.val.0."]=NA
          if(printT1==TRUE) inputs[,"p.val.1."]=NA
        }
      }}
    
    estimates[[j]]=round(inputs[,"Estimate"],estimateDigits)
    
    if(printClassical==TRUE){
      estimates[[j]]=cbind(estimates[[j]],round(inputs[,"t.ratio.0."],tDigits))
      if(printPVal==TRUE) estimates[[j]]=cbind(estimates[[j]],round(inputs[,"p.val.0."],pDigits))
      if(printT1==TRUE){ 
        estimates[[j]]=cbind(estimates[[j]],round(inputs[,"t.ratio.1."],tDigits))
        if(printPVal==TRUE) estimates[[j]]=cbind(estimates[[j]],round(inputs[,"p.val.1."],pDigits))
      }
    }
    estimates[[j]]=cbind(estimates[[j]],round(inputs[,"Rob.t.ratio.0."],tDigits))
    if(printPVal==TRUE) estimates[[j]]=cbind(estimates[[j]],round(inputs[,"Rob.p.val.0."],pDigits))
    if(printT1==TRUE){ 
      estimates[[j]]=cbind(estimates[[j]],round(inputs[,"Rob.t.ratio.1."],tDigits))
      if(printPVal==TRUE) estimates[[j]]=cbind(estimates[[j]],round(inputs[,"Rob.p.val.1."],pDigits))
    }
    
    rownames(estimates[[j]])=inputs[,1]
    if(j==1){
      combined_names=rownames(estimates[[j]])
    }else{
      combined_names=c(combined_names,rownames(estimates[[j]]))
    }
    j=j+1
  }
  
  combined_names=unique(combined_names)
  
  combined_outputs=data.frame(matrix(0,nrow=length(combined_names)+1,ncol=values*length(modelNames)))
  row.names(combined_outputs)=c("",combined_names)
  
  outputnames=c("estimate")
  if(printClassical==TRUE){
    outputnames=c(outputnames,"t-ratio(0)")
    if(printPVal==TRUE) outputnames=c(outputnames,"p-val(0)")
    if(printT1==TRUE){
      outputnames=c(outputnames,"t-ratio(1)")
      if(printPVal==TRUE) outputnames=c(outputnames,"p-val(1)")
    }
  }
  outputnames=c(outputnames,"Rob.t-ratio(0)")
  if(printPVal==TRUE) outputnames=c(outputnames,"Rob.p-val(0)")
  if(printT1==TRUE){
    outputnames=c(outputnames,"Rob.t-ratio(1)")
    if(printPVal==TRUE) outputnames=c(outputnames,"Rob.p-val(1)")
  }
  
  combined_outputs[1,]=rep(outputnames,length(modelNames))
  
  j=1
  while(j<=length(modelNames)){
    applicable_names=combined_names[combined_names%in%row.names(estimates[[j]])]  
    combined_outputs[applicable_names,((j-1)*values+1):(j*values)]=estimates[[j]][applicable_names,]
    inapplicable_names=combined_names[!(combined_names%in%row.names(estimates[[j]]))] 
    combined_outputs[inapplicable_names,((j-1)*values+1):(j*values)]=NA
    j=j+1
  }
  
  j=1
  while(j<=length(modelNames)){
    
    lines <- readLines(paste(modelNames[[j]],"_output.txt",sep="")) 
    
    inputvar = grep("Model name", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position = gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[1,j]=substr(inputvar,position+1,nchar(inputvar))
    } else otheroutputs[1,j] = NA
    
    k=2
    inputvar = grep("Number of individuals", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("Number of observations", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("Estimated parameters", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("LL\\(final\\)", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("Adj.Rho-square \\(0\\)", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("AIC", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    k=k+1
    inputvar = grep("BIC", lines) 
    if(length(inputvar)!=0){
      inputvar = lines[inputvar]
      position=gregexpr(pattern=":",inputvar)[[1]][1]
      otheroutputs[k,j]=as.double(substr(inputvar,position+1,nchar(inputvar)))
    } else otheroutputs[k,j] = NA
    
    
    j=j+1
  }
  
  
  otherouputs_new=data.frame(matrix("",nrow=8,ncol=values*length(modelNames)))
  rownames(otherouputs_new)=rownames(otheroutputs)
  j=1
  while(j<=length(modelNames)){
    otherouputs_new[,((j-1)*values+1)]=otheroutputs[,j]  
    colnames(otherouputs_new)[((j-1)*values+1)]=colnames(otheroutputs)[j]
    colnames(otherouputs_new)[((j-1)*values+2)]=""
    if(printClassical==TRUE){
      colnames(otherouputs_new)[((j-1)*values+3)]=""  
    }
    j=j+1
  }
  
  filename=paste("model_comparison_",gsub("[: -]", "" , Sys.time(), perl=TRUE),".csv",sep="")
  
  utils::write.table(otherouputs_new, filename, sep = ",", col.names = F, append = T)
  utils::write.table("",filename, sep = ",", col.names = F, append = T)
  utils::write.table(combined_outputs, filename, sep = ",", col.names = F, append = T)
  
  cat("Outputs of apollo_combineResults saved to",filename,"\n")
  return(invisible(TRUE))
}