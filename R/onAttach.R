#' Prints package startup message
#' 
#' This function is only called by R when attaching the package.
#' 
#' @param libname Name of library.
#' @param pkgname Name of package.
#' @return Nothing
.onAttach <- function(libname, pkgname) {
  
  ### Welcome message
  apolloVersion <- tryCatch(utils::packageDescription("apollo", fields = "Version"),
                            warning=function(w) return("alpha"),
                            error=function(e) return("alpha"))
  
  txt0 <- paste0("\n",  
                 "\n             . ,,                                                            ",
                 "\n            ,      ,,                                                        ",
                 "\n ,,,,,,    ,         ,,                                                      ",
                 "\n,     ,,  ,            ,,,,.                                                 ",
                 "\n,,     , ,,   ,,,,,,    ,,,                                 //  //           ",
                 "\n  ,     ,,,.   ,,,,,.   ,,      ////                        //  //           ",
                 "\n,,     ,,,,,.           ,,     // //     //////    /////    //  //    /////  ",
                 "\n,,,        ,,           ,      //  //    /    //  //   //   //  //   //   // ",
                 "\n              ,,       ,      ////////   /    //  //   //   //  //   //   // ",
                 "\n                ,,   ,,      //     //   /   ///  //   //   //  //   //   // ",
                 "\n                   ,         //      //  /////      ///      //  //    ///   ",
                 "\n                                         //                                  ",
                 "\n                                         //                                  ",
                 "\n")                                       
    
  txt1 <- paste0("\nApollo ", apolloVersion,
  "\nhttp://www.ApolloChoiceModelling.com",
  "\nSee url for a detailed manual, examples and a user forum.",
  "\nSign up to the user forum to receive updates on new releases.",
  "\n",
  "\nPlease cite Apollo in all written material you produce:",
  "\nHess S, Palma D (2019). \"Apollo: a flexible, powerful and customisable",
  "\nfreeware package for choice model estimation and application.\" Journal",
  "\nof Choice Modelling, 32. doi.org/10.1016/j.jocm.2019.100170",
  "\n",
  "\nThe developers of Apollo acknowledge the substantial support provided by",
  "\nthe European Research Council (ERC) through the consolidator grant DECISIONS,",
  "\nthe proof of concept grant APOLLO, and the advanced grant SYNERGY."
  )
  
  ### Print message
  txt=paste0(txt0,txt1)
  
  txt=cli::make_ansi_style("#0097ba",bg=FALSE)(txt)
  #cat(cli::make_ansi_style("#0097ba",bg=FALSE)(txt))
  packageStartupMessage(txt)

  #if(isOld) txt=paste0(txt,txt2)
  #packageStartupMessage(txt)
  #if(isOld) Sys.sleep(5)
}