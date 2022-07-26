% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{opossum_covariates}
\alias{opossum_covariates}
\title{Spatial environmental covariates throughout Chicago, IL}
\format{
A data frame with 96 rows and 6 variables in Site order:
\describe{
  \item{Site}{Abbreviation for the location sampled}
  \item{Building_age}{The median building age within 1 km of a site, in years. Data comes
  from the 2014-2018 American Community Survey.
  }
  \item{Impervious}{Percent impervious cover (0 - 100 range) within 1 km of a site. Data came
  from the 2016 National Landcover Database Developed Imperviousness product.
  }
  \item{Income}{Per capita income within 1 km of a site in US dollars. Data comes from
  the 2014-2018 American Community Survey.
  }
  \item{Population_density}{people per km^2 within 1 km of a site. Data came from
  the Silvis Lab's block-level housing density change database for U.S. cities.
  }
  \item{Vacancy}{The density of vacant buildings within 1 km of a site (units per km^2). Data
  came from the 2014-2018 American Community Survey.}
  ...
}
}
\source{
\url{https://doi.org/10.1111/gcb.15800}
}
\usage{
opossum_covariates
}
\description{
This is a companion dataset to go along with \code{\link{opossum_det_hist}},
which is detection / non-detection data for Virginia opossum (Didelphis virginiana).
This contains spatial covariates collected within a 1 km
buffer of each camera trapping site (i.e., no temporal variation).
}
\keyword{datasets}
