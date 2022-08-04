#' @title Fit an autologistic occupancy model
#'
#' @description Given repeat sampling across multiple primary sampling periods,
#' this model accounts for first-order temporal autocorrelation in species
#' occupancy patterns at a site with an additional parameter added to the
#' latent state of the model (\code{theta}). See Details for additional information.
#'
#'
#' @param formula Double right-hand side formula describing covariates of
#' detection (rho) and occupancy (psi) in that order.
#'
#' @param y A three-dimensional array of species detections. The first dimension
#' is sites, the second dimension is primary sampling periods, and the third
#' dimension is the secondary sampling periods within each primary sampling
#' period. If the species was detected on a given survey, that element would receive
#' a 1, otherwise it is 0. If sampling did not occur for a given sampling period, those
#' elements should be NA.
#'
#' @param det_covs A named list of detection covariates. If a given covariate
#' does not vary temporally, then that list element can be a vector of length nsite (i.e., \code{dim(y)[1]}. If a given
#' covariate varies across primary sampling periods, then that list element
#' should be a site by primary sampling period data.frame or matrix (i.e.,
#' dimensons \code{dim(y)[1]} by \code{dim(y)[2]}. If a given
#' covariate varies across primary and secondary sampling periods, then that
#' list element should be a matrix or data.frame that is nsite rows with a number
#' of columns equal to the product of the number of primary and secondary sampling
#' periods (i.e., \code{dim(y)[1]} by \code{prod(dim(y)[2:3])}. Thus, each
#' primary sampling period has a number of columns equal to the number of
#' secondary sampling periods. See details for additional information.
#'
#' @param occ_covs Either a data.frame or named list of occupancy covariates. If
#' there are temporally varying covariates, use a list. If there are no temporally
#' varying covariates, use a data.frame. If a given covariate
#' does not vary temporally, then that list element can be a vector of length nsite (i.e., \code{dim(y)[1]}. If a given
#' covariate varies across primary sampling periods, then that list element
#' should be a site by primary sampling period data.frame or matrix (i.e.,
#' dimensons \code{dim(y)[1]} by \code{dim(y)[2]}.
#'
#' @param method The optimization method used by \code{\link[stats]{optim}}. Defaults
#' to \code{BFGS}.
#'
#' @param level The confidence interval size to be calculated for model parameters.
#'   Defaults to 0.95.
#'
#' @param ... additional arguments. Not used yet.
#'
#' @details
#'
#' The autologistic occupancy model can be seen as a very simplified
#' version of a dynamic occupancy model that makes inference on species
#' patterns of occupancy from one time period to the next such that species
#' presence at a site during one time step may modify the probability a species
#' occupies that same site in the following time step. This is done through the
#' inclusion of a single logit-scale autoregressive term into the model: \eqn{\theta}.
#' For \eqn{i} in \eqn{1 \dots I} sites and \eqn{t} in \eqn{1 \dots T} primary sampling periods, let
#' \eqn{z_{i,t}} be a species latent occupancy state. For T=1, we do not
#' know if a species occupied a site in a previous time step
#' (as sampling has not occurred). THus, given a vector of parameters (\eqn{\boldsymbol{\beta}};
#'  i.e., the occupancy intercept and slope terms) and a matrix of
#'  covariates whose leading column is a vector of 1's (\eqn{\boldsymbol{X}}),
#'  the probability of occupancy during the first time step is the logit-linear predictor
#' at the first time step is:
#'
#' \deqn{\Large\mathrm{logit}(\psi_{i,t=1}) = \boldsymbol{\beta}\boldsymbol{x}_{i}}
#'
#'where
#'
#' \deqn{\Large z_{i,t=1}\sim \mathrm{Bernoulli}(\psi_{i,t=1})}
#'
#' Following the first time-step, we introduce \eqn{\theta} to help account for changes in occupancy
#' at time t given their presence at the same site at time t-1. Thus, for the remaining
#' time periods the logit-linear predictor is
#'
#' \deqn{\Large \mathrm{logit}(\psi{i,t}) = \boldsymbol{\beta}\boldsymbol{x}_{i} + \theta \times z_{i,t-1} }
#'
#'
#' where
#'
#'\deqn{\Large z_{i,t}\sim \mathrm{Bernoulli}(\psi_{i,t}), t>1}
#'
#' For the data model, let \eqn{\boldsymbol{Y}} be a site by primary sampling period by observation period
#' array (i.e., a three dimensional array). Thus, for \eqn{i} in \eqn{1 \dots I} sites, \eqn{t} in \eqn{1 \dots T} primary sampling periods and
#' \eqn{j} in \eqn{1 \dots J} observation events within a given time period (i.e., secondary samples within a primary sampling period), the scalar
#' \eqn{y_{i,t,j}} can the value of 1 if the species was detected, 0 if it was not detected but sampling occurred, and NA
#' if sampling did not occur. Therefore, you can consider the vector \eqn{y_{i,t,1:J}} the detection history for a species at site \eqn{i} and time
#' \eqn{t}. When creating this array for \code{auto_occ}, it is important that \eqn{J} equal the max number of
#' observation events that happened across all time periods. Thus, if there is variation in the amount of sampling
#' events across time periods you can pad the time periods with less data with NA values.
#'
#' Given a vector of parameters (\eqn{\boldsymbol{\alpha}};
#'  i.e., the detection intercept and slope terms) and a matrix of
#'  covariates whose leading column is a vector of 1's (\eqn{\boldsymbol{W}}),
#'  the logit-linear predictor for the data model is
#'
#' \deqn{\Large\mathrm{logit}(\rho_{i,t,j}) = \boldsymbol{\alpha}\boldsymbol{w}_{i}}
#'
#' where
#'
#' \deqn{\Large y_{i,t,j}\sim \mathrm{Bernoulli}(\rho_{i,t,j} \times z_{i,t})}
#'
#' Note that in this example detection probability only varied by sites, but this
#' level of the model can incorporate covariates that vary be site, time period, or
#' observation.
#'
#' Just like with a dynamic occupancy model,  site occupancy estimates
#' are derived parameters, and are most easily handled with using
#' \code{predict()}. Nevertheless, calculating expected occupancy
#' can be done with the following equation:
#'
#'\deqn{\Large
#'  \frac{
#'    \mathrm{ilogit}(\boldsymbol{\beta}\boldsymbol{x}_{i})
#'  }{
#'  \mathrm{ilogit}(\boldsymbol{\beta}\boldsymbol{x}_{i}) + (1 - \mathrm{ilogit}(\boldsymbol{\beta}\boldsymbol{x}_{i} + \theta))
#'  }
#'}
#'
#' Keeping in mind that your covariate structure may vary (here there
#' is simply variation by sites). Again, \code{predict()} can easily
#' handle this for you.
#'
#' @export
#' @importFrom stats as.formula
#' @importFrom stats optim
#' @importFrom stats qnorm
#' @importFrom methods new
#'
#' @examples
#' \dontrun{
#' # fitting an intercept only model
#' data("opossum_det_hist")
#'
#' # create y array
#' opossum_y <- autoOcc::format_y(
#'   x = opossum_det_hist,
#'   site_column = "Site",
#'   time_column = "Season",
#'   history_columns = "^Week", # regex for starts with Week
#'   report = FALSE # defaults to TRUE, turned off for example
#' )
#'
#' # fit the model
#' m1 <- autoOcc::auto_occ(
#'   formula = ~1 ~1,
#'   y = opossum_y
#' )
#'
#' # Now bring in some covariate data
#'
#' data("opossum_covariates")
#' oc <- opossum_covariates
#'
#' # scale the covariates (with base R)
#' oc_scaled <- as.data.frame(
#'   lapply(
#'     oc,
#'     function(x){
#'       if(is.numeric(x)){
#'         scale(x)
#'       }else{
#'         x
#'       }
#'     }
#'   )
#' )
#'
#' # And fit a single model with covariates that vary spatially
#'
#' m2 <- autoOcc::auto_occ(
#'   formula = ~Impervious ~Impervious,
#'   y = opossum_y,
#'   det_covs = oc_scaled,
#'   occ_covs = oc_scaled
#' )
#' }
#'

auto_occ <- function(formula, y, det_covs = NULL, occ_covs = NULL,
                     method = "BFGS", level = 0.95,...){
  # parse the formulas real quick
  if(length(formula)!= 3){
    stop("Double right-hand side formula required")
  }
  char <- lapply(formula, function(x) {
    paste(deparse(x), collapse = "")
  })
  rho_formula <- as.formula(char[[2]])
  psi_formula <- as.formula(paste("~", char[[3]]))

  # check to see if any sites need to get dropped
  na_count <- apply(
    is.na(y),
    1, sum
  )
  to_go <- which(
    na_count == prod(dim(y)[2:3])
  )
  if(length(to_go) > 0){
    cat(
      "\nSome sites have no data.\nRemoved sites at rows: ",
      ifelse(
        length(to_go) == 1,
          to_go,
          paste0(to_go, collapse = ", ")
      )
    )
    y <- y[-to_go,,]
  }

  # get sites, seasons, and reps from y array
  nsite <- dim(y)[1]
  nseason <- dim(y)[2]
  nrep <- dim(y)[3]

  # get the design matrices together
  if(is.null(occ_covs)){
    occ_covs <- data.frame(intercept = rep(1, nsite + length(to_go)))
  }
  occ_dm <- get_dm(
    x = occ_covs,
    my_formula = psi_formula,
    type = "psi",
    y=y,
    to_drop = to_go
  )

  if(is.null(det_covs)){
    det_covs <- list(
      intercept = data.frame(
        matrix(
          1,
          ncol = nseason,
          nrow = nsite + length(to_go)
        )
      )
    )
  }
  rho_dm <- get_dm(
    det_covs,
    rho_formula,
    "rho",
    y=y,
    to_drop = to_go
  )

  # and the number of parameters

  # Occupancy, +1 for theta term
  nocc_parms <- ncol(occ_dm[[1]]) + 1

  # Detection
  nrho_parms <- ncol(rho_dm[[1]][[1]])

  # total parameters
  nparms <- nocc_parms + nrho_parms

  my_covs <- list(
    psi = occ_dm,
    rho = rho_dm
  )

  initial_parms <- rep(0, nparms)
  fit <- optim(
    par = initial_parms,
    negloglik,
    y = y,
    cov_list = my_covs,
    method = method,
    hessian = TRUE
  )
  mle_table <- data.frame(
    parameter = c(
      paste0("psi - ",colnames(my_covs$psi[[1]])),
      "psi - theta",
      paste0("rho - ",colnames(my_covs$rho[[1]][[1]]))
    ),
    Est = fit$par,
    SE = sqrt(diag(solve(fit$hessian)))
  )
  mle_table$lower <- qnorm(
    (1-level)/2,
    mle_table$Est,
    mle_table$SE
  )
  mle_table$upper <- qnorm(
    1 - (1-level)/2,
    mle_table$Est,
    mle_table$SE
  )
  aic <- 2 * fit$value + 2 * nparms

  if(length(to_go) == 0){
    to_go <- NULL
  }
  to_return <- new(
    "auto_occ_fit",
    fitType = "auto_occ_fit",
    call = match.call(),
    formula = formula,
    y = y,
    estimates = mle_table,
    AIC = aic,
    opt = fit,
    negLogLike = fit$value,
    nllFun = negloglik,
    detcovs = det_covs,
    occcovs = occ_covs,
    sites_removed = to_go
  )
  return(to_return)

}
