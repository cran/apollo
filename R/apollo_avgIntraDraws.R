#' Averages intra-individual draws
#'
#' Averages observation-specific likelihood across intra-individual draws.
#'
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Description of the desired output from \code{apollo_probabilities}. Can take the values: "estimate", "prediction", "validate", "zero_LL", "conditionals", "output", "raw".
#' @return Likelihood averaged over intra-individual draws (shape depends on argument \code{functionality}).
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: Returns the likelihood of the model averaged across intra-individual draws. Drops all components but \code{"model"}.
#'           \item \strong{\code{"prediction"}}: Returns the likelihood of all alternatives and all model components averaged across intra-individual draws.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"zero_LL"}}: Returns \code{P} without changes.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"output"}}: Returns the likelihood of all model components averaged across intra-individual draws.
#'           \item \strong{\code{"raw"}}: Returns \code{P} without changes.
#'         }
#' @export
apollo_avgIntraDraws <- function(P, apollo_inputs, functionality){
  apollo_control=apollo_inputs[["apollo_control"]]
  
  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #
  
  if(apollo_control$HB==TRUE) stop('Function apollo_avgIntraDraws should not be used when apollo_control$HB==TRUE!')
  
  # ########################################## #
  #### functionality="gradient"             ####
  # ########################################## #
  if(functionality=="gradient") stop("Analytical gradient only implemented for models without mixing.")
  
  # ############################### #
  #### pre-checks                ####
  # ############################### #
  
  if(!apollo_control$mixing) stop('No mixing used in model!')
  
  inputIsList <- is.list(P)
  
  if(inputIsList && functionality!="prediction" && is.null(P[["model"]])) stop('Element called "model" is missing in list P!')
  
  # ########################################## #
  #### functionality="zero_LL/raw"          ####
  # ########################################## #
  
  if(functionality %in% c("zero_LL","raw")) return(P)
  
  # #################################################### #
  #### functionality="estimate/conditionals/validate" ####
  # #################################################### #
  
  if(functionality %in% c("estimate", "conditionals", "validate")){
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
