% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apollo_basTest.R
\name{apollo_basTest}
\alias{apollo_basTest}
\title{Ben-Akiva & Swait test}
\usage{
apollo_basTest(model1, model2)
}
\arguments{
\item{model1}{Either a character variable with the name of a previously estimated model, or an estimated 
model in memory, as returned by \link{apollo_estimate}.}

\item{model2}{Either a character variable with the name of a previously estimated model, or an estimated 
model in memory, as returned by \link{apollo_estimate}.}
}
\value{
Ben-Akiva & Swait test p-value (invisibly)
}
\description{
Calculates the p-value for the Ben-Akiva & Swait test for non-nested models. 
The two models need to both be discrete choice, and estimated on the same data.
}
