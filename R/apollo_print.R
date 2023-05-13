#' Prints message to terminal
#'
#' Prints message to terminal if \code{apollo_inputs$silent} is FALSE
#'
#' @param txt Character, what to print.
#' @param nSignifD Optional numeric integer. Minimum number of significant digits when printing numeric matrices. Default is 4.
#' @param widthLim Optional numeric integer. Minimum width (in characters) of each column when printing numeric matrices. Default is 11
#' @param pause Scalar integer. Number of seconds the execution will pause after printing the message. Default is 0.
#' @param type Character. "t" for regular text (default), "w" for warning, "i" for information.
#' @return Nothing
#' @export
apollo_print <- function(txt, nSignifD=4, widthLim=11, pause=0, type="t"){
  
  ### Validate input
  test <- is.numeric(pause) && is.vector(pause) && length(pause)==1
  if(!test) stop("INTERNAL ISSUE - Argument 'pause' must be a scalar integer.")
  pause <- as.integer(round(pause,0))
  test <- is.character(type) && is.vector(type) && length(type)==1 && type %in% c("t", "w", "i")
  if(!test) stop("INTERNAL ISSUE - Argument 'type' must be characters 't', 'w', or 'i'.")
  test <- is.character(txt) || is.matrix(txt)
  if(!test) stop("INTERNAL ISSUE - Argument 'txt' must be a matrix or a character")
  
  # If txt is character
  if(is.character(txt)){
    if(type=="t") writeLines(strwrap(txt, exdent=2))
    if(type=="w"){
      # Attach WARNING in red text and yellow background
      ###txt <- strwrap(paste0("WARNING: ", txt), exdent=2)
      txt <- strwrap(txt, exdent=2)
      #if(nchar(txt[1])>8) txt[1] <- substr(txt[1], 9, nchar(txt[1]))
      txt <- paste(txt, collapse="\n")
      # Print to screen
      #cat(cli::make_ansi_style("#ddcc77", bg=TRUE)( cli::make_ansi_style("#882255", bg=FALSE)("WARNING:") ),
      #    txt, "\n", sep="")
      cat(cli::make_ansi_style("#ddcc77", bg=TRUE)( cli::make_ansi_style("#882255", bg=FALSE)("WARNING:") ),
          txt, "\n", sep=" ")
    }
    if(type=="i"){
      # Attach INFORMATION in red text and yellow background
      ###txt <- strwrap(paste0("INFORMATION: ", txt), exdent=2)
      txt <- strwrap(txt, exdent=2)
      #if(nchar(txt[1])>12) txt[1] <- substr(txt[1], 13, nchar(txt[1]))
      txt <- paste(txt, collapse="\n")
      # Print to screen
      cat(cli::bg_br_black(cli::col_br_white("INFORMATION:")),
          txt, "\n", sep=" ")
    }
  }
  
  # If matrix
  if(is.matrix(txt)){
    if(type=="w") cat(cli::make_ansi_style("#ddcc77", bg=TRUE)( cli::make_ansi_style("#882255", bg=FALSE)("WARNING:") ), "\n")
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
  if(pause>0){
    cat("\n")
    txt <- paste0("Current process will resume in ", pause, " seconds unless interrupted by the user")
    cat(paste(strwrap(txt, indent=2, exdent=2), collapse="\n"))
    for(i in 1:pause){ Sys.sleep(1); cat('.')}
    cat("\n\n")
  }
  return(invisible(NULL))
}