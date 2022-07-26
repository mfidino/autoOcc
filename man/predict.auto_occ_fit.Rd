% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\docType{methods}
\name{predict.auto_occ_fit}
\alias{predict.auto_occ_fit}
\title{Makes predictions from an autologistic occupancy model}
\usage{
\method{predict}{auto_occ_fit}(
  object,
  type,
  newdata = NULL,
  level = 0.95,
  nsim = 3000,
  seed = NULL,
  ...
)
}
\arguments{
\item{object}{Object of class inheriting from \code{"auto_occ_fit"}.}

\item{type}{Either \code{"psi"} for occupancy or \code{"rho"} for detection.
See details for how expected occupancy is derived from this model object.}

\item{newdata}{Optionally, a data.frame of covariates to make predictions from.
If omitted, the fitted linear predictors are used.}

\item{level}{Tolerance / confidence level for predictions. Defaults to \code{0.95}.}

\item{nsim}{The number of parmater simulations to be made via the
models estimated parameters, variance covariance matrix,
and \code{\link[mvtnorm]{rmvnorm}}. Defaults to 3000.}

\item{seed}{The random seed to set for simulations, defaults to \code{NULL}.}

\item{...}{additional arguments. Not used yet.}
}
\value{
If \code{newdata} is supplied, this function will return a \code{data.frame}
with a number of rows equal to \code{newdata} and three columns.
\itemize{
  \item{estimate}{Median estimate from the model}
  \item{lower}{Lower confidence interval, based on confidence level supplied to \code{level}}
  \item{upper}{Upper confidence interval, based on confidence level supplied to \code{level}}
}

If \code{newdata = NULL}, the same data.frame will be output except the covariates
used to fit the model will be added on as additional columns. If there are temporally
varying covariates, the outputted data.frame will be in long format.
}
\description{
Predicted values based on a \code{auto_occ_fit} model object.
}
\details{
This function takes the output of an \code{\link[autoOcc]{auto_occ}}
model and will generate predictions for either occupancy (\code{type = "psi"})
or detection (\code{type = "rho"}). If \code{newdata} is supplied, all
covariates supplied to that level of the model must be present except for
the case of the autologistic occupancy term (that is handled). To
approximate uncertainty in model estimates, this function first uses
\code{\link[mvtnorm]{rmvnorm}} to generate parameter estimates from
multivariate distribution with a mean vector equal to the model parameters
and a variance covariance matrix supplied by the \code{\link[autoOcc]{auto_occ}}
model. This method was chosen as it creates functionally similar uncertainty
estimates as a Bayesian autologistic occupancy model.

For occupancy, each set of simulated parameters is used to
generate two logit-linear predictions for the supplied covariates:

\deqn{\Large \mathrm{logit}(\psi_a) = \beta_0 + \beta_1 \times x_1 + \dots + \beta_n \times x_n}

and

\deqn{\Large \mathrm{logit}(\psi_b) = \beta_0 + \beta_1 \times x_1 + \dots + \beta_n \times x_n + \theta}

where \eqn{\theta} is the estimated autologistic term. Following this, the expected occupancy of an autologisitic occupancy model is

\deqn{\huge
 \frac{
   \mathrm{ilogit}(\psi_a)
 }{
 \mathrm{ilogit}(\psi_a) + (1 - \mathrm{ilogit}(\psi_b))
 }
}
which is similar to the expected occupancy of a dynamic occupancy model (\eqn{\gamma \div (\gamma + \epsilon)})
where \eqn{\gamma} is colonization and \eqn{\epsilon} is extinction. Following this calculation for all simuated
parameter estimates and covariate values, the median estimate and confidence intervals are collected across
simulations for each covariate value.

Detection predictions are more straight-forward given there is no need to
derive the expected value. Simulations are still carried out to create to
generate a vector of values for each parameter.
}
\examples{
data("opossum_det_hist")
data("opossum_covariates")
odh <- opossum_det_hist
oc <- opossum_covariates

# function to generate detection history
opossum_y <- autoOcc::format_y(
  x = odh,
  site_column = "Site",
  time_column = "Season",
  history_columns = "^Week" # start with Week
)

# scale the covariates (with base R)
oc_scaled <- as.data.frame(
  lapply(
    oc,
    function(x){
      if(is.numeric(x)){
        scale(x)
      }else{
        x
      }
    }
  )
)
# dropping site column from oc_scaled
oc_scaled <- oc_scaled[,-1]
# suppress warnings because opossum_y has sites
#  with no data.
m1 <- auto_occ(
  ~Impervious + Income  ~ Impervious + Income,
  y = opossum_y,
  det_covs = oc_scaled,
  occ_covs = oc_scaled
)

# first make the prediction data.frame with a realistic
#   range based on the actual data and not the scaled data.
#   The range(oc$Impervious) is about 18 to 81, so choose 20
#   to 80. We do this so that we have nice numbers for plotting.
#   Likewise, we scaled all of the other data, so we leave Income
#   at it's mean (i.e., 0) for predictions.
imperv_real <- data.frame(
  Impervious = seq(20,80,0.5),
  Income = 0
)

# We will use imperv_real for plotting purposes, but to make predictions
#  we need to scale imperv_real$Impervious in the exact same way we did
#  with the fitted model. Thus, we subtract the mean of the actual data
#  and divide by the standard deviation.
imperv_scaled <- imperv_real
imperv_scaled$Impervious <- (
  imperv_scaled$Impervious - mean(oc$Impervious)
) / sd(oc$Impervious)

# the model prediction across a gradient of Impervious cover
opo_imperv <- predict(
  object = m1,
  type = "psi",
  newdata = imperv_scaled
)

# do the same thing with income
income_real <- data.frame(
  Impervious = 0,
  Income = seq(40000, 160000, by = 500)
)

income_scaled <- income_real
income_scaled$Income <- (income_scaled$Income - mean(oc$Income)) / sd(oc$Income)

opo_income <- predict(
  object = m1,
  type = "psi",
  newdata = income_scaled
)

# plot them out
par(mfrow = c(1,2))
plot(
  opo_imperv$estimate ~ imperv_real$Impervious,
  bty = "l",
  type = "l",
  las = 1,
  ylab = "Occupancy",
  xlab= "Impervious Cover (\%)",
  ylim = c(0,1),
  lwd = 3
)
lines(opo_imperv$lower ~ imperv_real$Impervious, lwd = 2, lty = 2)
lines(opo_imperv$upper ~ imperv_real$Impervious, lwd = 2, lty = 2)
plot(
  opo_income$estimate ~ income_real$Income,
  bty = "l",
  type = "l",
  las = 1,
  ylab = "Occupancy",
  xlab= "Per Capita Income (US Dollar)",
  ylim = c(0,1),
  lwd = 3
)
lines(opo_income$lower ~ income_real$Income, lwd = 2, lty = 2)
lines(opo_income$upper ~ income_real$Income, lwd = 2, lty = 2)
}
