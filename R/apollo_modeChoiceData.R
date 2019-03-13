#' Simulated dataset of mode choice.
#'
#' A simulated dataset containing 8000 mode choices among four alternatives.
#'
#' This dataset is to be used for discrete choice modelling.
#' Data comes from 500 individuals, each with one revealed
#' preferences (RP) observation, and 15 stated preferences
#' (SP) observations. There are 8000 choices in total.
#' Data is simulated.
#' Each observation contains attributes of the alternatives,
#' availability of alternatives, and characteristics of the
#' individuals.
#' @format A data frame with 8000 rows and 25 variables:
#' \describe{
#'   \item{ID}{Numeric. Identification number of the individual.}
#'   \item{RP}{Numeric. 1 if the row corresponds to a revealed preference (RP) observation. 0 otherwise.}
#'   \item{RP_journey}{Numeric. Consecutive ID of RP observation.}
#'   \item{SP}{Numeric. 1 if the row corresponds to a stated preference (SP) observation. 0 otherwise.}
#'   \item{SP_task}{Numeric. Consecutive ID of SP choice task.}
#'   \item{access_air}{Numeric. Access time (in minutes) of mode air.}
#'   \item{access_bus}{Numeric. Access time (in minutes) of mode bus.}
#'   \item{access_rail}{Numeric. Access time (in minutes) of mode rail.}
#'   \item{av_air}{Numeric. 1 if the mode air (plane) is available. 0 otherwise.}
#'   \item{av_bus}{Numeric. 1 if the mode bus is available. 0 otherwise.}
#'   \item{av_car}{Numeric. 1 if the mode car is available. 0 otherwise.}
#'   \item{av_rail}{Numeric. 1 if the mode rail (train) is available. 0 otherwise.}
#'   \item{business}{Numeric. Purpose of the trip. 1 for business, 0 for other.}
#'   \item{choice}{Numeric. Choice indicator, 1=car, 2=bus, 3=air, 4=rail.}
#'   \item{cost_air}{Numeric. Cost (in GBP) of mode air.}
#'   \item{cost_bus}{Numeric. Cost (in GBP) of mode bus.}
#'   \item{cost_car}{Numeric. Cost (in GBP) of mode car.}
#'   \item{cost_rail}{Numeric. Cost (in GBP) of mode rail.}
#'   \item{female}{Numeric. Sex of individual. 1 for female, 0 for male.}
#'   \item{income}{Numeric. Income (in GBP per annum) of the individual.}
#'   \item{service_air}{Numeric. Additional services in the air mode. 0 for none, 1 for a meal, 2 for wifi, 3 for meal and wifi.}
#'   \item{service_rail}{Numeric. Additional services in the rail mode. 0 for none, 1 for a meal, 2 for wifi, 3 for meal and wifi.}
#'   \item{time_air}{Numeric. Travel time (in minutes) of mode air.}
#'   \item{time_bus}{Numeric. Travel time (in minutes) of mode bus.}
#'   \item{time_car}{Numeric. Travel time (in minutes) of mode car.}
#'   \item{time_rail}{Numeric. Travel time (in minutes) of mode rail.}
#' }
#' @source \url{http://www.apollochoicemodelling.com/}
"apollo_modeChoiceData"
