#' Calculates product of panel observations.
#' 
#' Multiplies likelihood of observations from the same individual, or adds the log of them.
#' 
#' This function should be called inside apollo_probabilities only if the data has a panel structure.
#' It should be called after apollo_avgIntraDraws if intra-individual draws are used.
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"} For model estimation, returns probabilities of chosen alternatives.
#'                        \item \code{"prediction"} For model predictions, returns probabilities of all alternatives.
#'                        \item \code{"validate"} Validates input.
#'                        \item \code{"zero_LL"} Return probabilities with all parameters at zero.
#'                        \item \code{"conditionals"} For conditionals, returns probabilities of chosen alternatives.
#'                        \item \code{"output"} Checks that the model is well defined.
#'                        \item \code{"raw"} For debugging, returns probabilities of all alternatives.
#'                      }
#' @return Probabilities at the individual level.
#' @export
apollo_panelProd <- function(P, apollo_inputs, functionality){
  apollo_control = apollo_inputs[["apollo_control"]]
  workInLogs     = apollo_control$workInLogs
  
  inputIsList <- is.list(P)
  if(inputIsList && is.null(P[["model"]]) && functionality!="prediction") stop("The list \"P\" should contain an element called \"model\"!")
  
  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #
  
  if(apollo_control$HB==TRUE) stop('Function apollo_panelProd should not be used when apollo_control$HB==TRUE!')
  
  # ############################### #
  #### pre-checks                ####
  # ############################### #
  
  if(!apollo_control$panelData){
    stop('Panel data setting not used, so multiplying over choices not applicable!')
    return(P)
  } 
  
  # ########################################## #
  #### functionality="gradient"             ####
  # ########################################## #
  if(functionality=="gradient"){
    if(!apollo_control$panelData) return(P)
    if(is.list(P)) P <- P[[1]]
    g <- function(b, db){
      #like  <- get("apollo_probabilities", envir=parent.frame(), inherits=TRUE)
      input <- get("apollo_inputs", envir=parent.frame(), inherits=TRUE)
      indiv <- db[,input$apollo_control$indivID]
      #input$apollo_control$panelData <- FALSE
      #dLL <- P(b, db)
      #LL  <- like(b, input, functionality="estimate")
      #if(!input$apollo_control$workInLogs) LL <- log(LL)
      #dLL <- rowsum(dLL/LL, group=indiv)
      #LL  <- stats::aggregate(as.data.frame(LL), by=list(ind=indiv), FUN=prod)[,2]
      #return( LL*dLL )
      dLL <- P(b, db)
      dLL <- rowsum(dLL, group=indiv)
      return(dLL)
    }
    environment(g) <- new.env(parent=baseenv())
    assign("P", P, envir=environment(g))
    return(g)
  }
  
  # ############################################# #
  #### functionality="prediction/validate/raw" ####
  # ############################################# #
  
  ### NEW LINES SH 28/2
  if(functionality %in% c("prediction","raw")) return(P)
  #if(functionality %in% c("prediction","raw","validate")) return(P)
  
  # ########################################################## #
  #### functionality="estimate/conditionals/zero_LL/output" ####
  # ########################################################## #
  
  indivID <- apollo_inputs$database[, apollo_control$indivID]
  
  #if(functionality=="zero_LL"){
  #  if(inputIsList) P <- P[["model"]]
  #  Pout = rowsum(log(P), group=indivID)
  #  if(!workInLogs) Pout <- exp(Pout)
  #  if(is.matrix(Pout) && ncol(Pout)==1) Pout=as.vector(Pout)
  #  if(inputIsList) Pout <- list(model=Pout)
  #  return(Pout)
  #}
  
  ### NEW LINES SH 28/2
  #if(functionality %in% c("estimate", "conditionals")){
  if(functionality %in% c("estimate", "conditionals","validate")){
    if(inputIsList) P <- P[["model"]]
    if(is.array(P) && length(dim(P))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
    if(is.vector(P) || (is.matrix(P) && !workInLogs) ){
      Pout <- rowsum(log(P), group=indivID)
      if(!workInLogs) Pout <- exp(Pout)
    }
    if(apollo_control$panelData && is.matrix(P) && workInLogs){
      # approach to use if working in logs with mixing
      B    <- rowsum(log(P), group=indivID) # nIndiv x nDraws
      Bbar <- rowMeans(B) # nIndiv x 1
      Pout <- Bbar + log( rowMeans(exp(B-Bbar)) ) # nIndiv x 1
    }
    if(is.matrix(Pout) && ncol(Pout)==1) Pout <- as.vector(Pout) ### Added on 16 April due to Fangqing spotted bug
    if(inputIsList) Pout <- list(model=Pout)
    return(Pout)
  }
  
  if(functionality%in%c("output","zero_LL")){
    if(!inputIsList) P <- list(model=P)
    j=1
    Pout=P
    while(j<= length(P)){
      if(is.array(P[[j]]) && length(dim(P[[j]]))==3) stop('Need to average over intra-individual draws first before multiplying over choices!')
      if(is.vector(P[[j]]) || (is.matrix(P[[j]]) && !workInLogs) ){
        Pout[[j]] <- rowsum(log(P[[j]]), group=indivID)
        if(!workInLogs) Pout[[j]] <- exp(Pout[[j]])
      }
      if(apollo_control$panelData && is.matrix(P[[j]]) && workInLogs){
        # approach to use if working in logs with mixing
        B    <- rowsum(log(P[[j]]), group=indivID) # nIndiv x nDraws
        Bbar <- rowMeans(B) # nIndiv x 1
        Pout[[j]] <- Bbar + log( rowMeans(exp(B-Bbar)) ) # nIndiv x 1
      }
      if(is.matrix(Pout[[j]]) && ncol(Pout[[j]])==1) Pout[[j]] <- as.vector(Pout[[j]])
      j=j+1
    }
    if(!inputIsList) Pout <- Pout[[1]]
    return(Pout)
  }
  
  
}