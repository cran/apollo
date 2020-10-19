#' Prints message to terminal
#'
#' Prints message to terminal if apollo_inputs$silent is FALSE
#'
#' @param txt Character, what to print.
#' @param nSignifD Optional numeric integer. Minimum number of significant digits when printing numeric matrices. Default is 4.
#' @param widthLim Optional numeric integer. Minimum width (in characters) of each column when printing numeric matrices. Default is 11
#' @param highlight Optional logical. If TRUE, the message will be highlighted and will pause excecution for 5 seconds.
#' @return Nothing
#' @export
apollo_print <- function(txt, nSignifD=4, widthLim=11, highlight=FALSE){
  # If highlighted
  if(highlight){ cat(paste0(rep('#', getOption("width")-1), collapse='')); cat('\n')}
  
  # If character
  if(is.character(txt)) writeLines(strwrap(txt, exdent=2))
  
  #if(is.vector(txt)) writeLines(strwrap(paste0(txt, collapse=", ")))
  
  # If matrix
  if(is.matrix(txt)){
    if(is.numeric(txt)){ # If numeric matrix
      # Initialise
      M  <- txt
      M2 <- matrix('', nrow=nrow(M), ncol=ncol(M), dimnames=list(rownames(M), colnames(M)))
      # Go over each column
      for(j in 1:ncol(M)){
        m <- formatC(M[,j], digits=nSignifD, width=widthLim, format='g', flag='#')
        scie <- regexpr('e', m, fixed=TRUE)>0
        nDec <- regexpr('.', m, fixed=TRUE)
        nDec <- ifelse(nDec>0, nchar(m)-nDec, 0)
        # Convert to scientific notation if it is too long
        test <- !scie & nDec>6
        if(any(test)) m[test] <- formatC(M[test,j], format='e', digits=nSignifD, width=widthLim)
        if(any(test)) scie <- scie | test
        # Increase decimals if there are too few
        if(any(!scie)) objDec <- min(max(nDec[!scie]),6) else objDec <- 6
        test <- !scie & nDec<objDec
        if(any(test)) m[test] <- formatC(M[test,j], digits=objDec, width=widthLim, format='f')
        # Add empty spaces at the beginning if lacking
        test <- nchar(m)<max(nchar(m))
        if(any(test)) for(i in which(test)) m[i] <- paste0(rep(' ', max(nchar(m))-nchar(m[i])), m[i], collapse='')
        M2[,j] <- m
      }
      print(as.data.frame(M2), max=(nrow(M2)+1)*(ncol(M2)+1))
    } else print(txt) # If not numeric matrix
  }
  
  # If highlighted
  if(highlight){
    cat(paste0(rep('#', getOption("width")-1-5), collapse=''))
    for(i in 1:5){ Sys.sleep(1); cat('#')}
    cat('\n')
  } 
  return(invisible(NULL))
}