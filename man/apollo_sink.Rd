% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_sink.R
\name{apollo_sink}
\alias{apollo_sink}
\title{Starts or stops writing output to a text file.}
\usage{
apollo_sink(apollo_inputs = NULL)
}
\arguments{
\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.
If not provided, it will be looked for in the global environment.}
}
\value{
Nothing.
}
\description{
Starts or stops writing the output shown in the console to a file named "modelName_additional_output.txt".
}
\details{
After the first time this function is called, all output shown in the console will also be written to a
text file called "modelName_additional_output.txt", where "modelName" is the modelName set inside
apollo_control. 
The second time this function is called, it stops writing the console output to the file. The user 
should always call this function an even number of times to close the output file and prevents data loss.
}
