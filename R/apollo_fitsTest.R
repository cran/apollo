#' Compares log-likelihood of model across categories
#' 
#' Given the estimates of a model, it compares the log-likelihood at the observation level across categories of observations.
#' 
#' Prints a table comparing the average log-likelihood at the observation level for each category.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item apollo_beta: Named numeric vector. Names and values of model parameters.
#'                            \item apollo_inputs: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item functionality: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param fitsTest_settings List of arguments. It must contain the following elements.
#'                         \itemize{
#'                              \item \strong{subsamples}: Named list of boolean vectors. Each element of the list defines whether a given observation belongs to a given subsample (e.g. by sociodemographics).
#'                              \item \strong{modelComponent}: Name of model component. Set to model by default.
#'                         }
#' @return Matrix with average log-likelihood at observation level per category (invisibly).
#' @export
apollo_fitsTest=function(model,apollo_probabilities,apollo_inputs,fitsTest_settings){
  
  if(is.null(fitsTest_settings[["modelComponent"]])) fitsTest_settings$modelComponent="model"
  if(is.null(fitsTest_settings[["subsamples"]])) fitsTest_settings[["subsamples"]]=NULL
  
  fits=model$avgLL
  
  
  if(is.null(fitsTest_settings[["subsamples"]])){
    iterations=0
  } else {
    categories = fitsTest_settings[["subsamples"]]
    iterations = length(categories)
  }
  
  output=matrix(0,nrow=6,ncol=(iterations+1))
  output[1,1]=min(fits)
  output[2,1]=mean(fits)
  output[3,1]=stats::median(fits)
  output[4,1]=max(fits)
  output[5,1]=stats::sd(fits)
  output[6,1]=NA
  
  for(j in 1:iterations){
    categories[[j]]=apollo_firstRow(categories[[j]],apollo_inputs)
    
    tmp = fits[categories[[j]]]
    output[1,j+1]=min(tmp)
    output[2,j+1]=mean(tmp)
    output[3,j+1]=stats::median(tmp)
    output[4,j+1]=max(tmp)
    output[5,j+1]=stats::sd(tmp)
    output[6,j+1]=mean(tmp-mean(fits[!categories[[j]]]))
    
  }
  if(iterations>0){
    colnames(output)=c("All data",names(categories))
  } else {
    colnames(output)="All data"
  }
  rownames(output)=c("Min LL per obs","Mean LL per obs","Median LL per obs","Max LL per obs","SD LL per obs","mean vs mean of all others")
  if(iterations==0) output=output[1:5,]
  print(round(output,2))
  invisible(output)
}
