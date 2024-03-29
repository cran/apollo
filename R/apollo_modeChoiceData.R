#' Simulated dataset of mode choice.
#'
#' A simulated dataset containing 8,000 mode choices among four alternatives.
#'
#' This dataset is to be used for discrete choice modelling.
#' Data comes from 500 individuals, each with two revealed
#' preferences (RP) observation, and 14 stated stated
#' (SC) observations. There are 8,000 choices in total.
#' Data is simulated.
#' Each observation contains attributes for the alternatives,
#' availability of alternatives, and characteristics of the
#' individuals.
#' @format A data.frame with 8,000 rows and 25 variables:
#' \describe{
#'   \item{ID}{Numeric. Identification number of the individual.}
#'   \item{RP}{Numeric. 1 if the row corresponds to a revealed preference (RP) observation. 0 otherwise.}
#'   \item{RP_journey}{Numeric. Consecutive ID of RP observations. 0 if SP observation.}
#'   \item{SP}{Numeric. 1 if the row corresponds to a stated preference (SP) observation. 0 otherwise.}
#'   \item{SP_task}{Numeric. Consecutive ID of SP choice tasks. 0 if RP observation.}
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
#'   \item{service_air}{Numeric. Additional services for the air alternative. 1 for no-frills, 2 for wifi, 3 for food. This is not used in the RP data, where it is set to 0.}
#'   \item{service_rail}{Numeric. Additional services for the rail alternative. 1 for no-frills, 2 for wifi, 3 for food. This is not used in the RP data, where it is set to 0.}
#'   \item{time_air}{Numeric. Travel time (in minutes) of mode air.}
#'   \item{time_bus}{Numeric. Travel time (in minutes) of mode bus.}
#'   \item{time_car}{Numeric. Travel time (in minutes) of mode car.}
#'   \item{time_rail}{Numeric. Travel time (in minutes) of mode rail.}
#' }
#' @source \url{http://www.apollochoicemodelling.com/}
"apollo_modeChoiceData"
