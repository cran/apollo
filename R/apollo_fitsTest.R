#' Compares fit of model across categories
#' 
#' Given the predictions of a model, it compares the fit across categories of observations.
#' 
#' Prints a table comparing the average fit for each category.
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
#'                              \item subsamples: Named list of boolean vectors. Each element of the list defines whether a given observation belongs to a given subsample (e.g. by sociodemographics).
#'                              \item modelComponent: Name of model component. Set to model by default.
#'                         }
#' @return Matrix with average fit per category (invisibly).
#' @export
apollo_fitsTest=function(model,apollo_probabilities,apollo_inputs,fitsTest_settings){
  
  if(is.null(fitsTest_settings[["modelComponent"]])) fitsTest_settings$modelComponent="model"
  if(is.null(fitsTest_settings[["subsamples"]])) fitsTest_settings[["subsamples"]]=NULL

   fits       = apollo_prediction(model, 
                                  apollo_probabilities, apollo_inputs,fitsTest_settings$modelComponent)
   
   fits=fits[,ncol(fits)]

  if(is.null(fitsTest_settings[["subsamples"]])){
    iterations=0
  } else {
    categories = fitsTest_settings[["subsamples"]]
    iterations=length(categories)
  }
  
  output=matrix(0,nrow=6,ncol=(iterations+1))
  output[1,1]=min(fits)
  output[2,1]=mean(fits)
  output[3,1]=stats::median(fits)
  output[4,1]=max(fits)
  output[5,1]=stats::sd(fits)
  output[6,1]=NA
  
  j=1
  while(j<=iterations){
    temp1=as.data.frame(cbind(fits,categories[[j]]))
    output[1,j+1]=min(temp1[temp1$V2==1,1])
    output[2,j+1]=mean(temp1[temp1$V2==1,1])
    output[3,j+1]=stats::median(temp1[temp1$V2==1,1])
    output[4,j+1]=max(temp1[temp1$V2==1,1])
    output[5,j+1]=stats::sd(temp1[temp1$V2==1,1])
    output[6,j+1]=mean(temp1[temp1$V2==1,1])-mean(temp1[temp1$V2!=1,1])
    j=j+1
  }
  if(iterations>0){
    colnames(output)=c("All data",names(categories))
  } else {
    colnames(output)="All data"
  }
  rownames(output)=c("Min P(chosen)","Mean P(chosen)","Median P(chosen)","Max P(chosen)","SD P(chosen)","mean vs mean of all other")
  if(iterations==0) output=output[1:5,]
  print(round(output,2))
  invisible(output)
}
