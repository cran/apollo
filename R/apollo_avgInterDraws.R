#' Averages inter-individual draws
#'
#' Averages individual-specific likelihood across inter-individual draws.
#'
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Description of the desired output from \code{apollo_probabilities}. Can take the values: "estimate", "prediction", "validate", "zero_LL", "conditionals", "output", "raw".
#' @return Likelihood averaged over inter-individual draws (shape depends on argument \code{functionality}).
#'         \itemize{
#'           \item"estimate": Returns the likelihood of the model averaged across inter-individual draws.
#'           \item"prediction": Returns the likelihood of all alternatives and all model components across inter-individual draws.
#'           \item"validate": Returns P without changes.
#'           \item"zero_LL": Returns P without changes.
#'           \item"conditionals": Returns P without changes.
#'           \item"output": Returns the same as "estimate", but also prints a summary of the estimation data.
#'           \item"raw": Returns P without changes.
#'         }
#' @export
apollo_avgInterDraws <- function(P, apollo_inputs, functionality){
  apollo_control = apollo_inputs[["apollo_control"]]
  workInLogs    <- apollo_control$workInLogs

  # ############################### #
  #### ignored for HB estimation ####
  # ############################### #

  if(apollo_control$HB==TRUE) return(P)

  # ############################### #
  #### pre-checks                ####
  # ############################### #

  if(!apollo_control$mixing) stop('No mixing used in model!')

  if(apollo_control$panelData){
    indivID <- apollo_inputs$database[, apollo_control$indivID]
    nIndiv <- length(unique(indivID))
    if(is.list(P)){
      if(is.array(P[[1]])) pRows <- dim(P[[1]])[1] else pRows <- length(P[[1]])
    } else {
      if(is.array(P)) pRows <- dim(P)[1] else pRows <- length(P)
    }
    if(nIndiv!=pRows & !(functionality %in% c("zero_LL", "raw", "validate", "prediction"))) stop("Observations from the same individual must be aggregated (e.g. multiplied) before averaging over inter-individual draws.")
  }

  inputIsList <- is.list(P)

  if(inputIsList && functionality!="prediction" && functionality!="raw" && is.null(P[["model"]])) stop('Element called "model" is missing in list P!')

  # ########################################## #
  #### functionality="zero_LL/raw/validate" ####
  # ########################################## #

  if(functionality %in% c("zero_LL","raw","validate")) return(P)

  # ########################################## #
  #### functionality="estimate" ####
  # ########################################## #

  if(functionality=="estimate"){
    if(inputIsList) P <- P[["model"]]
    if(is.vector(P) && !workInLogs ) stop('No Inter-individuals draws to average over!')
    if(is.array(P) && length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(is.matrix(P)) P <- rowMeans(P)
    if(inputIsList) P <- list(model=P)
    return(P)
  }

  # ########################################## #
  #### functionality="prediction"           ####
  # ########################################## #

  if(functionality=="prediction"){
    if(!inputIsList){
      if(is.matrix(P)) output=rowMeans(P)
      if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
      return(output)
    }

    if(inputIsList){
      output_list=P
      j=1
      while(j <= length(P)){
        if(is.list(P[[j]])){
        k=1
        while(k<= length(P[[j]])){
          if(is.array(P[[j]][[k]]) && length(dim(P[[j]][[k]]))==3) stop('Intra-individual draws still present to average over!')
          if(is.matrix(P[[j]][[k]])) output_list[[j]][[k]]=rowMeans(P[[j]][[k]])
          k=k+1
        }}else{
          if(is.array(P[[j]]) && length(dim(P[[j]]))==3) stop('Intra-individual draws still present to average over!')
          if(is.matrix(P[[j]])) output_list[[j]]=rowMeans(P[[j]])
        }
        j=j+1
      }
      return(output_list)
    }
  }

  # ########################################## #
  #### functionality="conditionals"         ####
  # ########################################## #

  if(functionality=="conditionals"){
    if(inputIsList) P <- P[["model"]]
    if(!is.array(P)) stop('No draws present to average over!')
    if(is.array(P) & length(dim(P))==3) stop('Intra-individual draws still present to average over!')
    if(inputIsList) P <- list(model=P)
    return(P)
  }

  # ########################################## #
  #### functionality="output"               ####
  # ########################################## #

  if(functionality=="output"){
    if(!inputIsList) P <- list(model=P)
    output_list=P
    j=1
    while(j<= length(P)){
      if(is.array(P[[j]]) & length(dim(P[[j]]))==3) stop('Intra-individual draws still present to average over!')
      if(is.matrix(P[[j]])) output_list[[j]]=rowMeans(P[[j]])
      j=j+1
    }
    if(!inputIsList) output_list <- output_list[[1]]
    return(output_list)
  }

}
