#' Dataset of route choice.
#'
#' A Stated Preference dataset containing 3,492 route choices among two alternatives.
#'
#' This dataset is to be used for discrete choice modelling.
#' Data comes from 388 individuals who participated in a Stated
#' Choice (SC) survey, providing a total of 3,492 observations.
#' Each choice scenario includes two alternatives described in
#' terms of travel time, cost, headway and interchanges. Additional
#' information on respondents is available.
#' This dataset comes from the following publication.
#' Vrtic, M. & Axhausen, K.W. (2003), The impact of tilting trains in 
#' Switzerland: A route choice model of regional and long
#' distance public transport trips. 82nd annual meeting of 
#' the transportation research board, Washington, DC.
#' @format A data frame with 3,492 rows and 16 variables:
#' \describe{
#'   \item{ID}{Numeric. Identification number of the individual.}
#'   \item{choice}{Numeric. Choice indicator, 1 for alternative 1, and 2 for alternative 2.}
#'   \item{tt1}{Numeric. Travel time (in minutes) for alternative 1.}
#'   \item{tc1}{Numeric. Travel cost (in CHF) for alternative 1.}
#'   \item{hw1}{Numeric. Headway time (in minutes) for alternative 1.}
#'   \item{ch1}{Numeric. Number of interchanges for alternative 1.}
#'   \item{tt2}{Numeric. Travel time (in minutes) for alternative 2.}
#'   \item{tc2}{Numeric. Travel cost (in CHF) for alternative 2.}
#'   \item{hw2}{Numeric. Headway time (in minutes) for alternative 2.}
#'   \item{ch2}{Numeric. Number of interchanges for alternative 2.}
#'   \item{hh_inc_abs}{Numeric. Household income (in CHF per annum).}
#'   \item{car_availability}{Numeric. 1 if respondent has a car available, 0 otherwise.}
#'   \item{commute}{Numeric. 1 if the purpose of the trip is commuting. 0 otherwise.}
#'   \item{shopping}{Numeric. 1 if the purpose of the trip is shopping. 0 otherwise.}
#'   \item{business}{Numeric. 1 if the purpose of the trip is business. 0 otherwise.}
#'   \item{leisure}{Numeric. 1 if the purpose of the trip is leisure. 0 otherwise.}
#' }
#' @source \url{http://www.apollochoicemodelling.com/}
"apollo_swissRouteChoiceData"
