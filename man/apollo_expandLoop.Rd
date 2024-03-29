% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_expandLoop.R
\name{apollo_expandLoop}
\alias{apollo_expandLoop}
\title{Expands loops in a function or expression}
\usage{
apollo_expandLoop(f, apollo_inputs, validate = TRUE)
}
\arguments{
\item{f}{function (usually \code{apollo_probabilities}) inside which the name of the components are inserted.}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{validate}{Logical. If TRUE, the new function will be validated before being returned}
}
\value{
A function or an expression (same type as input \code{f})
}
\description{
Expands loops replacing the index by its value. It also evaluates \code{paste} and \code{paste0}, and removes \code{get}.
}
\details{
For example, the expression
\code{for(j in 1:3) V[[paste0('alt',j)]] = b1*get(paste0('x',j)) + b2*X[,j]}

would be expanded into:

\code{
V[[alt1]] = b1*x1 + b2*X[,1]
V[[alt2]] = b1*x2 + b2*X[,2]
V[[alt3]] = b1*x3 + b2*X[,3]
}
}
