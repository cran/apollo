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
  "\nSign up to the user forum to receive updates on new releases.")
  
  ### WARNING if more than six months old
  releaseDate <- as.POSIXct('2023-01-19', format='%Y-%m-%d')
  isOld <- as.POSIXct(Sys.Date()) > (releaseDate + 6*30.4*24*60*60)
  #if(isOld) txt2 <- paste0(txt, 
  #if(isOld) txt2 <- paste0('\n\nYour version of Apollo is more than six months old.',
  if(isOld) txt2 <- paste0('\nYour version of Apollo is more than six months old.',
                          '\nUsing the latest version will ensure you have all',
                          '\n current functionality and bug fixes.', 
                          '\nYou can update to the latest version by typing:', 
                          '\n install.packages("apollo")')
  
  ### Print message
  txt=paste0(txt0,txt1)
  
  txt=cli::make_ansi_style("#0097ba",bg=FALSE)(txt)
  #cat(cli::make_ansi_style("#0097ba",bg=FALSE)(txt))
  packageStartupMessage(txt)
  if(isOld){#cat("\n")
    apollo_print("\n")
    apollo_print(txt2, pause=5, type="w")
  } 
  
  #if(isOld) txt=paste0(txt,txt2)
  #packageStartupMessage(txt)
  #if(isOld) Sys.sleep(5)
}