#' Likelihood ratio test
#' 
#' Calculates the likelihood ratio test and prints result.
#' 
#' @param baseModel Character. Name of a previously estimated model whose results were written to disk by \link{apollo_saveOutput}.
#'                   This is the restricted model, i.e. the one with fewer parameters.
#' @param generalModel Either a character variable with the name of a previously estimated model, or an estimated model in memory, as
#'                      returned by \link{apollo_estimate}. This model nests \code{baseModel}, and it should have more parameters than it.
#' @return LL ratio test statistic (invisibly)
#' @export
apollo_lrTest=function(baseModel,generalModel){
lines <- readLines(paste(baseModel,"_output.txt",sep="")) 
id <- grep("LL\\(final", lines) 
fit=lines[id]
nchar(fit)
position=gregexpr(pattern=":",fit)[[1]][1]
base_LL=as.double(substr(fit,position+1,nchar(fit)))

id <- grep("Estimated parameters", lines) 
npar=lines[id]
position=gregexpr(pattern=":",npar)[[1]][1]
dfbase=as.double(substr(npar,position+1,nchar(npar)))

if(is.character(generalModel)){
  lines <- readLines(paste(generalModel,"_output.txt",sep="")) 
  id <- grep("LL\\(final", lines) 
  fit=lines[id]
  nchar(fit)
  position=gregexpr(pattern=":",fit)[[1]][1]
  general_LL=as.double(substr(fit,position+1,nchar(fit)))
  
  id <- grep("Estimated parameters", lines) 
  npar=lines[id]
  position=gregexpr(pattern=":",npar)[[1]][1]
  dfgeneral=as.double(substr(npar,position+1,nchar(npar)))
} else {
  general_LL=generalModel$maximum
  nParams     <- length(generalModel$apollo_beta)
  nFreeParams <- nParams
  if(!is.null(generalModel$apollo_fixed)) nFreeParams <- nFreeParams - length(generalModel$apollo_fixed)
  dfgeneral=nFreeParams
}

if(general_LL<base_LL) stop('General model LL should be superior to base model LL!')
df=dfgeneral-dfbase
if(df<1) stop('No degrees of freedom for LR test')

LR_test_value     = -2*( base_LL - general_LL )

cat("\nLikelihood ratio test-value:  ",round(LR_test_value,2))
cat("\nDegrees of freedom:           ",df)
cat("\nLikelihood ratio test p-value:",round(stats::pchisq(LR_test_value, df, lower.tail=FALSE),4),"\n")
return(invisible(LR_test_value))
}
