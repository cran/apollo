#' Averages intra-individual draws
#'
#' Averages observation-specific likelihood across intra-individual draws.
#'
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Description of the desired output from \code{apollo_probabilities}. Can take the values: "estimate", "prediction", "validate", "zero_LL", "conditionals", "output", "raw".
#' @return Likelihood averaged over intra-individual draws (shape depends on argument \code{functionality}).
#'         \itemize{
#'           \item"estimate": Returns the likelihood of the model averaged across intra-individual draws.
#'           \item"prediction": Returns the likelihood of all alternatives and all model components across intra-individual draws.
#'           \item"validate": Returns P without changes.
#'           \item"zero_LL": Returns P without changes.
#'           \item"conditionals": Returns P without changes.
#'           \item"output": Returns the same than "estimate", but also prints a summary of estimation data.
#'           \item"raw": Returns P without changes.
#'         }
#' @export
apollo_avgIntraDraws <- function(P, apollo_inputs, functionality){
  apollo_control=apollo_inputs[["apollo_control"]]
  
  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #
  
  if(apollo_control$HB==TRUE) return(P)
  
  # ############################### #
  #### pre-checks                ####
  # ############################### #
  
  if(!apollo_control$mixing) stop('No mixing used in model!')
  
  inputIsList <- is.list(P)
  
  if(inputIsList && functionality!="prediction" && is.null(P[["model"]])) stop('Element called "model" is missing in list P!')
  
  # ########################################## #
  #### functionality="zero_LL/raw/validate" ####
  # ########################################## #
  
  if(functionality %in% c("zero_LL","raw","validate")) return(P)
  
  # ########################################## #
  #### functionality="estimate/conditionals" ####
  # ########################################## #
  
  if(functionality=="estimate" | functionality=="conditionals"){
    if(!is.list(P)){
      if(is.array(P)){
        if(length(dim(P))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          return( colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3] )
        } else stop('No intra-individual draws present to average over!')
      } else stop('No draws present to average over!')
    } else {
      if(is.array(P[["model"]])){
        if(length(dim(P[["model"]]))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          return( list(model=colSums(aperm(P[["model"]], perm=c(3,1,2)))/dim(P[["model"]])[3] ))
        } else stop('No intra-individual draws present to average over!')
      } else stop('No draws present to average over!')
    }
  }
  
  # ########################################## #
  #### functionality="prediction"           ####
  # ########################################## #
  
  if(functionality=="prediction"){
    if(!is.list(P)){
      if(is.array(P) && length(dim(P))==3){
        # Returns a matrix of dimensions nObs x nDrawsInter
        output=colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3]
      }
      return(output)
    } else {
      output_list=P
      j=1
      while(j<= length(P)){
        k=1
        while(k<= length(P[[j]])){
          if(is.array(P[[j]][[k]]) && length(dim(P[[j]][[k]]))==3){
            # Returns a matrix of dimensions nObs x nDrawsInter
            output_list[[j]][[k]]=colSums(aperm(P[[j]][[k]], perm=c(3,1,2)))/dim(P[[j]][[k]])[3]
          }
          k=k+1}
        j=j+1}
      return(output_list)
    }
    
  }
  
  # ########################################## #
  #### functionality="output"               ####
  # ########################################## #
  
  if(functionality=="output"){
    if(!is.list(P)){
      if(is.array(P) && length(dim(P))==3){
        # Returns a matrix of dimensions nObs x nDrawsInter
        output=colSums(aperm(P, perm=c(3,1,2)))/dim(P)[3]
      }
      return(output)
    } else {
      output_list=P
      j=1
      while(j<= length(P)){
        if(is.array(P[[j]]) && length(dim(P[[j]]))==3){
          # Returns a matrix of dimensions nObs x nDrawsInter
          output_list[[j]]=colSums(aperm(P[[j]], perm=c(3,1,2)))/dim(P[[j]])[3]}
        j=j+1}
      return(output_list)
    }
  }
  
}
