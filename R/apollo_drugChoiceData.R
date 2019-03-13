#' Simulated dataset of medication choice.
#'
#' A simulated dataset containing 10,000 stated medication choices among four alternatives.
#'
#' This dataset is to be used for discrete choice modelling.
#' Data comes from 1,000 individuals, each with ten stated 
#' preferences (SP) observations among headache medication. 
#' There are 10,000 choices in total. Data is simulated.
#' Each observation contains attributes of the alternatives,
#' characteristics of the respondent, and their answers to
#' four attitudinal questions. All four alternatives are
#' always available for all individuals. Alternatives 1 and
#' 2 are branded, while alternatives 3 and 4 are generic.
#' Respondents provide a full ranking of alternatives for
#' each choice task (i.e. observation).
#' @format A data frame with 10000 rows and 33 variables:
#' \describe{
#'   \item{ID}{Numeric. Identification number of the individual.}
#'   \item{task}{Numeric. 1 if the row corresponds to a revealed preference (RP) observation. 0 otherwise.}
#'   \item{best}{Numeric. Consecutive ID of RP observation.}
#'   \item{second_pref}{Numeric. 1 if the row corresponds to a stated preference (SP) observation. 0 otherwise.}
#'   \item{third_pref}{Numeric. Consecutive ID of SP choice task.}
#'   \item{worst}{Numeric. Access time (in minutes) of mode air.}

#'   \item{brand_1}{Character. Name of alternative's brand.}
#'   \item{country_1}{Character. Name of alternative's country of origin.}
#'   \item{char_1}{Character. Characteristics of the alternative (standard, fast acting, or double strength).}
#'   \item{side_effects_1}{Numeric. Chance of suffering negative side effects if this alternative is consumed.}
#'   \item{price_1}{Numeric. Cost of this alternative in Pounds sterling (GBP).}

#'   \item{brand_2}{Character. Name of alternative's brand.}
#'   \item{country_2}{Character. Name of alternative's country of origin.}
#'   \item{char_2}{Character. Characteristics of the alternative (standard, fast acting, or double strength).}
#'   \item{side_effects_2}{Numeric. Chance of suffering negative side effects if this alternative is consumed.}
#'   \item{price_2}{Numeric. Cost of this alternative in Pounds sterling (GBP).}

#'   \item{brand_3}{Character. Name of alternative's brand.}
#'   \item{country_3}{Character. Name of alternative's country of origin.}
#'   \item{char_3}{Character. Characteristics of the alternative (standard, fast acting, or double strength).}
#'   \item{side_effects_3}{Numeric. Chance of suffering negative side effects if this alternative is consumed.}
#'   \item{price_3}{Numeric. Cost of this alternative in Pounds sterling (GBP).}

#'   \item{brand_4}{Character. Name of alternative's brand.}
#'   \item{country_4}{Character. Name of alternative's country of origin.}
#'   \item{char_4}{Character. Characteristics of the alternative (standard, fast acting, or double strength).}
#'   \item{side_effects_4}{Numeric. Chance of suffering negative side effects if this alternative is consumed.}
#'   \item{price_4}{Numeric. Cost of this alternative in Pounds sterling (GBP).}

#'   \item{regular_user}{Numeric. 1 if the respondent is a regular user of headache medicine, 0 otherwise.}
#'   \item{university_educated}{Numeric. 1 if the respondent holds a university degree, 0 otherwise.}
#'   \item{over_50}{Numeric. 1 if the respondent is 50 years old or older, 0 otherwise.}
#'   \item{attitude_quality}{Numeric. Level of agreement from 1 (strongly disagree) to 5 (strongly agree) with the phrase 'I am concerned about the quality of drugs developed by unknown companies'.}
#'   \item{attitude_ingredients}{Numeric. Level of agreement from 1 (strongly disagree) to 5 (strongly agree) with the phrase 'I believe that ingredients are the same no matter what brand'.}
#'   \item{attitude_patent}{Numeric. Level of agreement from 1 (strongly disagree) to 5 (strongly agree) with the phrase 'The original patent holders have valuable experience with their medicines'.}
#'   \item{attitude_dominance}{Numeric. Level of agreement from 1 (strongly disagree) to 5 (strongly agree) with the phrase 'I believe the dominance of big pharmaceutical companies is unhelpful'.}
#' }
#' @source \url{http://www.apollochoicemodelling.com/}
"apollo_drugChoiceData"
