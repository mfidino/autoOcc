#' Virginia opossum detection / non-detection data from Chicago, IL
#'
#' A dataset containing weekly detection / non-detection data
#' of Virginia opossum (Didelphis virginiana) at 102 unique
#' camera trapping locations during four primary sampling periods
#' in 2019 throughout the Chicago greater metropolitan area.
#' Camera traps were deployed for roughly 4 weeks over
#' deployment in January, April, July, and October. This is a
#' subset of the data used in a Urban Wildlife Information
#' Network (UWIN) publication (see source for link to
#' article). However, the data have been slightly formatted
#' from the original data so that there are detection weeks (
#' the UWIN publication used a Binomial distribution for the
#' detection model not a Bernoulli like \code{\link{auto_occ}})
#' does and so some formatting was needed).
#'
#' @format A data frame with 408 rows and 6 variables:
#' \describe{
#'   \item{Site}{Abbreviation for the location sampled}
#'   \item{Season}{Codes for the different primary sampling periods. JA19 = January 2019,
#'   AP19 = April 2019, JU19 = July 2019, and OC19 = October 2019}
#'   \item{Week_1}{For a given Site and Season whether the camera was active detected opossum (1),
#'   active and did not detect opossum (0), or was not active (NA)
#'   }
#'   \item{Week_2}{For a given Site and Season whether the camera was active detected opossum (1),
#'   active and did not detect opossum (0), or was not active (NA)
#'   }
#'   \item{Week_3}{For a given Site and Season whether the camera was active detected opossum (1),
#'   active and did not detect opossum (0), or was not active (NA)
#'   }
#'   \item{Week_4}{For a given Site and Season whether the camera was active detected opossum (1),
#'   active and did not detect opossum (0), or was not active (NA)
#'   }
#'   ...
#' }
#' @source \url{https://doi.org/10.1111/gcb.15800}
"opossum_det_hist"



#' Spatial environmental covariates throughout Chicago, IL
#'
#'
#' This is a companion dataset to go along with \code{\link{opossum_det_hist}},
#' which is detection / non-detection data for Virginia opossum (Didelphis virginiana).
#' This contains spatial covariates collected within a 1 km
#' buffer of each camera trapping site (i.e., no temporal variation).
#'
#' @format A data frame with 102 rows and 7 variables in Site order:
#' \describe{
#'   \item{Site}{Abbreviation for the location sampled}
#'   \item{Building_age}{The median building age within 1 km of a site, in years. Data comes
#'   from the 2014-2018 American Community Survey.
#'   }
#'   \item{Impervious}{Percent impervious cover (0 - 100 range) within 1 km of a site. Data came
#'   from the 2016 National Landcover Database Developed Imperviousness product.
#'   }
#'   \item{Income}{Per capita income within 1 km of a site in US dollars. Data comes from
#'   the 2014-2018 American Community Survey.
#'   }
#'   \item{Population_density}{people per km^2 within 1 km of a site. Data came from
#'   the Silvis Lab's block-level housing density change database for U.S. cities.
#'   }
#'   \item{Vacancy}{The density of vacant buildings within 1 km of a site (units per km^2). Data
#'   came from the 2014-2018 American Community Survey.}
#'   ...
#' }
#' @source \url{https://doi.org/10.1111/gcb.15800}
"opossum_covariates"
