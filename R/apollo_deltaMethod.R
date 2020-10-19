#' Delta method
#'
#' Applies the delta method to calculate the standard errors of transformations of parameters.
#' If the bootstrap covariance matrix is available, it is used. If not, the robust covariance matrix is used.
#'
#' \code{apollo_deltaMethod} supports the following five operations.
#' \describe{
#'   \item{sum}{Calculates the s.e. of \code{parName1} + \code{parName2}}
#'   \item{diff}{Calculates the s.e. of \code{parName1} - \code{parName2} and \code{parName2} - \code{parName1}}
#'   \item{ratio}{Calculates the s.e. of \code{parName1}/\code{parName2} and \code{parName2}/\code{parName1}}
#'   \item{exp}{Calculates the s.e. of \code{exp(parName1)}}
#'   \item{logistic}{If only \code{parName1} is provided, it calculates the s.e. of
#'         \code{exp(parName1)/(1+exp(parName1))} and \code{1/(1+exp(parName1))}.
#'         If \code{parName1} and \code{parName2} are provided, it calculates
#'         \code{exp(par_i)/(1+exp(parName1)+exp(parName2))} for i=1, 2, and 3 (par_3 = 1).}
#'   \item{lognormal}{Calculates the mean and s.d. of a lognormal distribution based on the mean (\code{parName1}) and s.d. (\code{parName2}) of the underlying normal.}
#' }
#'
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param deltaMethod_settings List of arguments. It must contain the following.
#'                             \itemize{
#'                               \item operation: Character. Function to calculate the delta method for. See details.
#'                               \item parName1: Character. Name of the first parameter.
#'                               \item parName2: Character. Name of the second parameter. Optional depending on \code{operation}.
#'                               \item multPar1: Numeric scalar. A value to scale \code{parName1}.
#'                               \item multPar2: Numeric scalar. A value to scale \code{parName2}.
#'                             }
#' @return Matrix contating calue, s.e. and t-ratio resulting from the operation. This is also printed to screen.
#' @export
apollo_deltaMethod=function(model, deltaMethod_settings){
  
  # Check the model was not estimated using HB
  apollo_control = tryCatch(get("apollo_inputs", envir=parent.frame())$apollo_control,
                            error=function(e){
                              cat("WARNING: Could not retrieve apollo_control.\n         Assuming this model was not estimated using HB.")
                              return(list(HB=FALSE))
                            })
  if(is.null(apollo_control$HB)) apollo_control$HB=FALSE
  if(apollo_control$HB) stop("The function \'apollo_deltaMethod\' is not applicables for models estimated using HB!") 
  
  if(is.null(model$estimate)) stop("Estimates missing from model object!")
  if(is.null(model$robvarcov) && is.null(model$bootvarcov)) stop("Covariance matrix missing from model object!")
  
  if(is.null(deltaMethod_settings[["operation"]])) stop("The deltaMethod_settings list needs to include an object called \"operation\"!")
  if(is.null(deltaMethod_settings[["parName1"]])) stop("The deltaMethod_settings list needs to include an object called \"parName1\"!")
  if(is.null(deltaMethod_settings[["parName2"]])) deltaMethod_settings[["parName2"]]=NA
  if(is.null(deltaMethod_settings[["multPar1"]])) deltaMethod_settings[["multPar1"]]=1
  if(is.null(deltaMethod_settings[["multPar2"]])) deltaMethod_settings[["multPar2"]]=1
  
  operation = deltaMethod_settings[["operation"]]
  parName1  = deltaMethod_settings[["parName1"]]
  parName2  = deltaMethod_settings[["parName2"]]
  multPar1  = deltaMethod_settings[["multPar1"]]
  multPar2  = deltaMethod_settings[["multPar2"]]
  
  operation <- tolower(operation)
  
  if( !(operation %in% c("sum","diff","ratio","exp","logistic","lognormal","prod")) ) stop("Invalid value of 'operation' parameter. See ?apollo_deltaMethod.")
  if(is.na(parName2) & !(operation %in% c("logistic","exp"))) stop("Need two parameters if using operation: ",operation)
  if(!(parName1 %in% names(model$estimate))) stop("parName1=", parName1, " not found among model estimates.")
  if(!is.na(parName2) && !(parName2 %in% names(model$estimate))) stop("parName2=", parName2, " not found among model estimates.")
  if(!is.na(parName2) & (operation=="exp")) stop("Should only have one parameters if using operation: ",operation)
  
  est <- model$estimate
  if(!is.null(model$est)) est=model$est
  if(!is.null(model$bootvarcov)) robvarcov = model$bootvarcov else robvarcov=model$robvarcov
  est[parName1]=multPar1*est[parName1]
  robvarcov[parName1,]=multPar1*robvarcov[parName1,]
  robvarcov[,parName1]=multPar1*robvarcov[,parName1]
  if(!is.na(parName2)){
    robvarcov[parName2,]=multPar2*robvarcov[parName2,]
    robvarcov[,parName2]=multPar2*robvarcov[,parName2]
    est[parName2]=multPar2*est[parName2]
  }
  
  if(multPar1!=1){
    parName1name=paste(parName1," (multiplied by ",multPar1,")",sep = "")
  } else {
    parName1name=parName1
  }
  if(multPar2!=1){
    parName2name=paste(parName2," (multiplied by ",multPar2,")",sep = "")
  } else {
    parName2name=parName2
  }
  
  if(operation=="sum"){
    v=est[parName1]+est[parName2]
    se=sqrt(robvarcov[parName1,parName1]+robvarcov[parName2,parName2]+2*robvarcov[parName1,parName2])
    operation_name=paste("Sum of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="diff"){
    v=est[parName1]-est[parName2]
    se=sqrt(robvarcov[parName1,parName1]+robvarcov[parName2,parName2]-2*robvarcov[parName1,parName2])
    operation_name=paste("Difference between ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="ratio"){
    v=est[parName1]/est[parName2]
    se=sqrt(v^2*(robvarcov[parName1,parName1]/(est[parName1]^2)+robvarcov[parName2,parName2]/(est[parName2]^2)-2*robvarcov[parName1,parName2]/(est[parName1]*est[parName2])))
    operation_name=paste("Ratio of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="prod"){
    v=est[parName1]*est[parName2]
    se=sqrt((est[parName2]^2)*robvarcov[parName1,parName1]+(est[parName1]^2)*robvarcov[parName2,parName2]+2*est[parName1]*est[parName2]*robvarcov[parName1,parName2])
    operation_name=paste("Product of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="exp"){
    v=exp(est[parName1])
    se=sqrt(robvarcov[parName1,parName1])*v
    operation_name=paste("Exponential of ",parName1name,": ",sep="")}
  
  if((operation=="logistic")&is.na(parName2)){
    v1=1/(1+exp(-est[parName1]))
    v2=1-v1
    se1=sqrt(robvarcov[parName1,parName1])*exp(-est[parName1])/(1+exp(-est[parName1]))^2
    se2=se1
    operation_name1=paste("exp(",parName1name,")/(1+exp(",parName1name,")): ",sep="")
    operation_name2=paste("1/(1+exp(",parName1name,")): ",sep="")
    v=rbind(v1,v2)
    se=rbind(se1,se2)
    operation_name=rbind(operation_name1,operation_name2)}
  
  if((operation=="logistic")&!is.na(parName2)){
    v1=exp(est[parName1])/(1+exp(est[parName1])+exp(est[parName2]))
    v2=exp(est[parName2])/(1+exp(est[parName1])+exp(est[parName2]))
    v3=1-v1-v2
    phi1=(exp(est[parName1])*(1+exp(est[parName2])))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    phi2=(-exp(est[parName1])*exp(est[parName2]))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    se1=sqrt(robvarcov[parName1,parName1]*phi1^2+robvarcov[parName2,parName2]*phi2^2+2*robvarcov[parName1,parName2]*phi1*phi2)
    
    phi1=(-exp(est[parName1])*exp(est[parName2]))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    phi2=(exp(est[parName2])*(1+exp(est[parName1])))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    se2=sqrt(robvarcov[parName1,parName1]*phi1^2+robvarcov[parName2,parName2]*phi2^2+2*robvarcov[parName1,parName2]*phi1*phi2)
    
    phi1=-exp(est[parName1])/((1+exp(est[parName1])+exp(est[parName2]))^2)
    phi2=-exp(est[parName2])/((1+exp(est[parName1])+exp(est[parName2]))^2)
    se3=sqrt(robvarcov[parName1,parName1]*phi1^2+robvarcov[parName2,parName2]*phi2^2+2*robvarcov[parName1,parName2]*phi1*phi2)
    
    operation_name1=paste("exp(",parName1name,")/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    operation_name2=paste("exp(",parName2name,")/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    operation_name3=paste("1/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    v=rbind(v1,v2,v3)
    se=rbind(se1,se2,se3)
    operation_name=rbind(operation_name1,operation_name2,operation_name3)}
  
  if((operation=="lognormal")){
    v1  = exp( est[parName1] + est[parName2]^2/2 )
    v2  = v1*sqrt( exp(est[parName2]^2) - 1 )
    se1 = sqrt( v1^2*(robvarcov[parName1,parName1] + est[parName2]^2*robvarcov[parName2,parName2] + 2*robvarcov[parName1,parName2]*est[parName2]) ) 
    se2 = sqrt( v2^2*robvarcov[parName1,parName1] + 
                  2*est[parName2]*(exp(2*est[parName1] + 2*est[parName2]^2) + v2^2)*robvarcov[parName1,parName2] + 
                  est[parName2]^2*((exp(2*est[parName1] + 2*est[parName2]^2))/v2 + v2)^2*robvarcov[parName2,parName2] ) 
    operation_name1 = paste("Mean for exp(N(",parName1name,",",parName2name,")")
    operation_name2 = paste("SD for exp(N(",parName1name,",",parName2name,")")
    v               = rbind(v1,v2)
    se              = rbind(se1,se2)
    operation_name  = rbind(operation_name1,operation_name2)
  }
  
  
  t = v/se
  
  delta_output=cbind(v,se,t)
  colnames(delta_output)=c("Value","Robust s.e.","Rob t-ratio (0)")
  rownames(delta_output)=operation_name
  
  cat("\nRunning Delta method computations\n")
  print(delta_output, digits=4)
  cat("\n")
  return(invisible(delta_output))
}
