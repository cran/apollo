#' Delta method for Apollo models
#'
#' Applies the Delta method to calculate the standard errors of transformations of parameters. 
#' 
#' \code{apollo_deltaMethod} can be used in two ways. The first and recommended way is to provide an 
#' element called \code{expression} inside its argument \code{deltaMethod_settings}. \code{expression} 
#' should contain the expression or expressions for which the standard error is/are to be calculated, as text. For 
#' example, to calculate the ratio between parameters b1 and b2, \code{expression=c(vtt="b1/b2")} should be used.
#' 
#' The second method is to provide the name of a specific operation inside \code{deltaMethod_settings}.
#' The following five operations are supported.
#' \itemize{
#'   \item \strong{\code{sum}}: {Calculates the s.e. of \code{parName1} + \code{parName2}}
#'   \item \strong{\code{diff}}: {Calculates the s.e. of \code{parName1} - \code{parName2} and \code{parName2} - \code{parName1}}
#'   \item \strong{\code{prod}}: {Calculates the s.e. of \code{parName1}*\code{parName2}}
#'   \item \strong{\code{ratio}}: {Calculates the s.e. of \code{parName1}/\code{parName2} and \code{parName2}/\code{parName1}}
#'   \item \strong{\code{exp}}: {Calculates the s.e. of \code{exp(parName1)}}
#'   \item \strong{\code{logistic}}: {If only \code{parName1} is provided, it calculates the s.e. of
#'                   \code{exp(parName1)/(1+exp(parName1))} and \code{1/(1+exp(parName1))}.
#'                   If \code{parName1} and \code{parName2} are provided, it calculates
#'                   \code{exp(par_i)/(1+exp(parName1)+exp(parName2))} for i=1, 2, and 3 (par_3 = 1).}
#'   \item \strong{\code{lognormal}}: {Calculates the mean and s.d. of a lognormal distribution based on the mean 
#'                    (\code{parName1}) and s.d. (\code{parName2}) of the underlying normal.}
#' }
#' 
#' By default, \code{apollo_deltaMethod} uses the robust covariance matrix. However, the user can change this through the \code{varcov} setting.
#' 
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param deltaMethod_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                             \itemize{
#'                               \item \strong{\code{expression}}: Character vector. A character vector with a single or multiple arbitrary functions of the estimated parameters, as text. 
#'                                                 For example: \code{c(VTT="b1/b2*60")}. Each expression can only contain model parameters (estimated or fixed), 
#'                                                 numeric values, and operands. At least one of the parameters used needs to not have been fixed in estimation. Variables in the database
#'                                                 cannot be included. If the user does not provide a name for an expression, then the expression itself is used in the output. If this setting is provided, then \code{operation}, 
#'                                                 \code{parName1}, \code{parName2}, \code{multPar1} and \code{multPar2} are
#'                                                 ignored.
#'                               \item \strong{\code{varcov}}: Character. Type of variance-covariance matrix to use in calculations. 
#'                                                        It can take values \code{"classical"}, \code{"robust"} and 
#'                                                        \code{"bootstrap"}. Default is \code{"robust"}.
#'                               \item \strong{\code{operation}}: Character. Function to calculate the delta method for. See details. Not used if \code{expression} is provided.
#'                               \item \strong{\code{parName1}}: Character. Name of the first parameter if \code{operation} is used. See details. Not used if \code{expression} is provided.
#'                               \item \strong{\code{parName2}}: Character. Name of the second parameter if \code{operation} is used. See details. Not used if \code{expression} is provided.. Optional depending on \code{operation}.
#'                               \item \strong{\code{multPar1}}: Numeric scalar. An optional value to scale \code{parName1}. Not used if \code{expression} is provided.
#'                               \item \strong{\code{multPar2}}: Numeric scalar. An optional value to scale \code{parName2}. Not used if \code{expression} is provided.
#'                             }
#' @return Matrix containing value, s.e. and t-ratio resulting from the requested expression or operation. This is also printed to screen.
#' @export
#' @importFrom stringr str_replace
apollo_deltaMethod <- function(model, deltaMethod_settings){
  
  #### Validation ####
  
  ### Validate model
  test <- is.list(model) && !is.null(names(model)) && !is.null(model$apollo_control)
  test <- test && any(c("varcov", "robvarcov", "bootvarcov") %in% names(model))
  test <- test && all("estimate" %in% names(model))
  if(!test) stop("Input \"model\" is not a valid Apollo model object.")
  test <- !is.null(model$apollo_control$HB) && !model$apollo_control$HB
  if(!test) stop("The function \'apollo_deltaMethod\' is not applicables for models estimated using HB!")
  
  ### Validate settings
  test <- is.list(deltaMethod_settings) && !is.null(names(deltaMethod_settings))
  if(!test) stop("Input \"deltaMethod_settings\" should be a named list. Type ?apollo_deltaMethod for help.")
  if(is.null(deltaMethod_settings[["varcov"]])) deltaMethod_settings[["varcov"]] <- "robust"
  test <- deltaMethod_settings[["varcov"]] %in% c("classical", "robust", "bootstrap")
  if(!test) stop("Setting \"varcov\" can only take values \"classical\", \"robust\" or \"bootstrap\".")
  if(is.null(deltaMethod_settings[["expression"]])){
    # Checks if "expression" is not provided
    nm <- names(deltaMethod_settings)
    mandatory <- c("operation", "parName1")
    for(i in mandatory) if(!(i %in% nm)) stop("The deltaMethod_settings list needs to include an element called \"", 
                                              i, "\" if \"expression\" is not provided.")
    nonMand <- c(parName1=NA, multPar1=1, multPar2=1)
    for(i in names(nonMand)) if(!(i %in% nm)) deltaMethod_settings[[i]] <- nonMand[i]
    rm(nm, mandatory, nonMand, i)
  } else {
    # Checks if "expression" is provided
    e <- deltaMethod_settings[["expression"]]
    test <- is.vector(e) && is.character(e)
    if(!test) stop("Element \"expression\" must be text.")
    tmp_e = e
    for(s in 1:length(tmp_e)){
      e <- tmp_e[s]
      e_name = ifelse((is.null(names(e))||names(e)==""),e,names(e))
      e <- tryCatch(as.expression(e), error=function(e) NULL)
      if(is.null(e)) stop("The expression ",e_name," could not be parsed into a valid R expression.")
      test <- all.vars(str2lang(deltaMethod_settings[["expression"]][s]))
      ### then expanded the below line
      test1 <- test[!(test %in% names(model$estimate))]
      if(length(test1)>0) stop("The expression ",e_name," includes variables that are not parameters from the model: ",
                               paste(test1, collapse=", "))
      test2 <- test[!(test %in% names(model$estimate)[!(names(model$estimate) %in% model$apollo_fixed)])]
      if(length(test2)==length(test)) stop("The expression ",e_name," needs to include at least one parameter that was not fixed in estimation!")
      if(length(test2)>0){
        for(tmp in test2) tmp_e[s]=stringr::str_replace(tmp_e[s],tmp,as.character(model$estimate[tmp]))
        apollo_print(paste0("The expression ",e_name," includes parameters that were fixed in estimation: ",
                            paste(test2, collapse=", ")))
        apollo_print(paste0("These have been replaced by their fixed values, giving:\n",
                            tmp_e[s]))
        cat("\n")
      } 
    }
    
    deltaMethod_settings[["expression"]]=tmp_e

    rm(e, test, tmp_e, e_name)
  }
  
  ### Check for availability of requested covariance matrix availability, and store it in vcov
  vcov <- deltaMethod_settings[["varcov"]]
  if(vcov=="classical") vcov <- "varcov"
  if(vcov=="robust"   ) vcov <- "robvarcov"
  if(vcov=="bootstrap") vcov <- "bootvarcov"
  test <- vcov %in% names(model)
  if(!test) stop("The requested variance covariance matrix cannot be found inside \"model\".")
  vcov <- model[[vcov]]
  rm(test)
  
  #### Processing using expression ####
  if(!is.null(deltaMethod_settings[["expression"]])){
    # Useful functions
    is.val <- function(e) if(is.symbol(e) || is.numeric(e) || is.character(e) || is.logical(e) ) return(TRUE) else return(FALSE)
    insertAB <- function(e, apollo_beta_names){
      if(is.character(e)){
        e <- str2lang(paste0("function(apollo_beta) ", e))
        e <- eval(e)
      } 
      isFunction <- is.function(e)
      if(isFunction){f <- e; e <- body(e)}
      # Case 1: x
      test1 <- is.symbol(e) && (as.character(e) %in% apollo_beta_names)
      if(test1) e <- str2lang(paste0("apollo_beta[\"", as.character(e), "\"]"))
      # Case 2: expression
      if( !test1 && (is.call(e) || is.expression(e)) ){
        for(i in 1:length(e)) if(!is.null(e[[i]])) e[[i]] <- insertAB(e[[i]], apollo_beta_names)
      } 
      # Return
      if(isFunction){body(f) <- e; e <- f}
      return(e)
    }
    delta_from_expression=function(f){
      b <- model$estimate[!(names(model$estimate) %in% model$apollo_fixed)]
      tmp=list()
      for(k in 1:length(names(b))) tmp[[k]]=paste0(names(b)[k],"=",Deriv::Deriv(f, names(b)[k]))
      tmp=paste0("c(",toString(tmp),")")
      g=tmp
      ###
      g <- insertAB(g, names(b))
      g <- setNames(g(b), names(b))
      if(length(g)!=dim(vcov)[1]) stop("Dimensions of \"model$estimate\" and the covariance matrix do not match.")
      # Prepare output
      return(data.frame(`Expression`     = ifelse((is.null(names(f))||names(f)==""),f,names(f)),
                       `Value`           = round(eval(str2lang(f), envir=list2env(as.list(b))),4),
                       `Robust s.e.`     = round(sqrt(t(g)%*%vcov%*%g),4),
                       `Rob t-ratio (0)` = round(eval(str2lang(f), envir=list2env(as.list(b)))/sqrt(t(g)%*%vcov%*%g),2),
                       check.names       = FALSE))
    }
    
    cat("Running Delta method computation for user-defined function:\n\n")
    out = delta_from_expression(deltaMethod_settings[["expression"]][1])
    if(length(deltaMethod_settings[["expression"]])>1){
      for(s in 2:length(deltaMethod_settings[["expression"]])) out=rbind(out,delta_from_expression(deltaMethod_settings[["expression"]][s]))
    }
    print(out,row.names=FALSE)
    return(invisible(out))
  }
  
  
  #### Processing using operation ####
  ### SH bug fix 22 Dec - these 3 lines were missing
  if(is.null(deltaMethod_settings[["parName2"]])) deltaMethod_settings[["parName2"]]=NA
  if(is.null(deltaMethod_settings[["multPar1"]])) deltaMethod_settings[["multPar1"]]=1
  if(is.null(deltaMethod_settings[["multPar2"]])) deltaMethod_settings[["multPar2"]]=1
  
  operation = deltaMethod_settings[["operation"]]
  parName1  = deltaMethod_settings[["parName1"]]
  parName2  = deltaMethod_settings[["parName2"]]
  multPar1  = deltaMethod_settings[["multPar1"]]
  multPar2  = deltaMethod_settings[["multPar2"]]
  
  if(parName1%in%model$apollo_fixed) stop(paste0("Parameter ",parName1," cannot be used with apollo_deltaMethod as it was fixed in estimation!"))
  if(parName2%in%model$apollo_fixed) stop(paste0("Parameter ",parName2," cannot be used with apollo_deltaMethod as it was fixed in estimation!"))
  operation <- tolower(operation)
  
  if( !(operation %in% c("sum","diff","ratio","exp","logistic","lognormal","prod")) ) stop("Invalid value of 'operation' parameter. See ?apollo_deltaMethod.")
  if(is.na(parName2) & !(operation %in% c("logistic","exp"))) stop("Need two parameters if using operation: ",operation)
  if(!(parName1 %in% names(model$estimate))) stop("parName1=", parName1, " not found among model estimates.")
  if(!is.na(parName2) && !(parName2 %in% names(model$estimate))) stop("parName2=", parName2, " not found among model estimates.")
  if(!is.na(parName2) & (operation=="exp")) stop("Should only have one parameters if using operation: ",operation)
  
  est <- model$estimate
  if(!is.null(model$est)) est=model$est
  if(!is.null(model$bootvarcov)) vcov = model$bootvarcov else vcov=model$robvarcov
  est[parName1]=multPar1*est[parName1]
  vcov[parName1,]=multPar1*vcov[parName1,]
  vcov[,parName1]=multPar1*vcov[,parName1]
  if(!is.na(parName2)){
    vcov[parName2,]=multPar2*vcov[parName2,]
    vcov[,parName2]=multPar2*vcov[,parName2]
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
    se=sqrt(vcov[parName1,parName1]+vcov[parName2,parName2]+2*vcov[parName1,parName2])
    operation_name=paste("Sum of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="diff"){
    v=est[parName1]-est[parName2]
    se=sqrt(vcov[parName1,parName1]+vcov[parName2,parName2]-2*vcov[parName1,parName2])
    operation_name=paste("Difference between ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="ratio"){
    v=est[parName1]/est[parName2]
    se=sqrt(v^2*(vcov[parName1,parName1]/(est[parName1]^2)+vcov[parName2,parName2]/(est[parName2]^2)-2*vcov[parName1,parName2]/(est[parName1]*est[parName2])))
    operation_name=paste("Ratio of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="prod"){
    v=est[parName1]*est[parName2]
    se=sqrt((est[parName2]^2)*vcov[parName1,parName1]+(est[parName1]^2)*vcov[parName2,parName2]+2*est[parName1]*est[parName2]*vcov[parName1,parName2])
    operation_name=paste("Product of ",parName1name," and ",parName2name,": ",sep="")}
  
  if(operation=="exp"){
    v=exp(est[parName1])
    se=sqrt(vcov[parName1,parName1])*v
    operation_name=paste("Exponential of ",parName1name,": ",sep="")}
  
  if((operation=="logistic")&is.na(parName2)){
    v1=1/(1+exp(-est[parName1]))
    v2=1-v1
    se1=sqrt(vcov[parName1,parName1])*exp(-est[parName1])/(1+exp(-est[parName1]))^2
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
    se1=sqrt(vcov[parName1,parName1]*phi1^2+vcov[parName2,parName2]*phi2^2+2*vcov[parName1,parName2]*phi1*phi2)
    
    phi1=(-exp(est[parName1])*exp(est[parName2]))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    phi2=(exp(est[parName2])*(1+exp(est[parName1])))/((1+exp(est[parName1])+exp(est[parName2]))^2)
    se2=sqrt(vcov[parName1,parName1]*phi1^2+vcov[parName2,parName2]*phi2^2+2*vcov[parName1,parName2]*phi1*phi2)
    
    phi1=-exp(est[parName1])/((1+exp(est[parName1])+exp(est[parName2]))^2)
    phi2=-exp(est[parName2])/((1+exp(est[parName1])+exp(est[parName2]))^2)
    se3=sqrt(vcov[parName1,parName1]*phi1^2+vcov[parName2,parName2]*phi2^2+2*vcov[parName1,parName2]*phi1*phi2)
    
    operation_name1=paste("exp(",parName1name,")/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    operation_name2=paste("exp(",parName2name,")/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    operation_name3=paste("1/(1+exp(",parName1name,")+exp(",parName2name,")): ",sep="")
    v=rbind(v1,v2,v3)
    se=rbind(se1,se2,se3)
    operation_name=rbind(operation_name1,operation_name2,operation_name3)}
  
  if((operation=="lognormal")){
    v1  = exp( est[parName1] + est[parName2]^2/2 )
    v2  = v1*sqrt( exp(est[parName2]^2) - 1 )
    se1 = sqrt( v1^2*(vcov[parName1,parName1] + est[parName2]^2*vcov[parName2,parName2] + 2*vcov[parName1,parName2]*est[parName2]) ) 
    se2 = sqrt( v2^2*vcov[parName1,parName1] + 
                  2*est[parName2]*(exp(2*est[parName1] + 2*est[parName2]^2) + v2^2)*vcov[parName1,parName2] + 
                  est[parName2]^2*((exp(2*est[parName1] + 2*est[parName2]^2))/v2 + v2)^2*vcov[parName2,parName2] ) 
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
