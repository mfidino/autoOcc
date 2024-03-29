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
#' @returns
#' An S4 of class \code{auto_occ_fit}. See \code{\linkS4class{auto_occ_fit}}
#' for additional details.
#'
#'
#' @examples
#' ###################################
#' # Example 1 : Intercept only model
#' ###################################
#'
#'  data("opossum_det_hist")
#'
#'  # reducing sample size a bit to speed up example run
#'  opossum_det_hist <- split(
#'    opossum_det_hist,
#'    factor(
#'      opossum_det_hist$Season,
#'      levels = unique(opossum_det_hist$Season)
#'     )
#'   )
#'  opossum_det_hist <- lapply(
#'    opossum_det_hist,
#'    function(x) head(x, 25)
#'  )
#'  opossum_det_hist <- do.call(
#'    "rbind.data.frame",
#'    opossum_det_hist
#'  )
#'  # create y array
#'  opossum_y <- format_y(
#'    x = opossum_det_hist,
#'    site_column = "Site",
#'    time_column = "Season",
#'    history_columns = "^Week", # regex for starts with Week
#'    report = FALSE # defaults to TRUE, turned off for example
#'  )
#'
#'
#'  # fit the model, s
#'  m1 <- auto_occ(
#'    formula = ~1 ~1,
#'    y = opossum_y
#'  )
#'
#'
#' ########################################
#' # Example 2 : Adding spatial covariates,
#' #               no temporal variation
#' ########################################
#'
#'  # Now bring in some covariate data,
#'  #  again reducing sample size here
#'  #  just like we did with opossum_y
#'
#'  data("opossum_covariates")
#'  oc <- opossum_covariates
#'  oc <- head(oc,25)
#'
#'  # scale the covariates (with base R)
#'  oc_scaled <- as.data.frame(
#'    lapply(
#'      oc,
#'      function(x){
#'        if(is.numeric(x)){
#'          scale(x)
#'        }else{
#'          x
#'        }
#'      }
#'    )
#'  )
#'
#'  # And fit a single model with covariates that vary spatially.
#'
#'  m2 <- auto_occ(
#'    formula = ~Impervious ~Impervious,
#'    y = opossum_y,
#'    det_covs = oc_scaled,
#'    occ_covs = oc_scaled
#'  )
#'
#'  summary(m2)
#'
#' ###################################################
#' # Example 3 : Spatial covariate, plus a categorical
#' #               covariate that varies by
#' #               primary sampling period
#' ###################################################
#'
#'  # And a third model with a categorical seasonal covariate.
#'  #  Temporally varying covariates need to be a named list.
#'  #  For this example, the seasonal information is in the
#'  #  opossum detection history (opossum_det_hist).
#'
#'  season_frame <- list(
#'    Season = matrix(
#'      opossum_det_hist$Season,
#'      ncol = dim(opossum_y)[2],
#'      nrow = dim(opossum_y)[1]
#'    ),
#'    Impervious = oc_scaled$Impervious
#'  )
#'  m3 <- auto_occ(
#'    formula =  ~Season + Impervious ~Season + Impervious,
#'    y = opossum_y,
#'    det_covs = season_frame,
#'    occ_covs = season_frame
#'  )
#'
#'  summary(m3)
#'
#' ###################################################
#' # Example 4 :   Detection covariates that vary by
#' #               secondary sampling period.
#' ###################################################
#'
#' # This example is with simulated data.
#'
#' set.seed(3122020)
#'
#' # latent occupancy, just doing intercept only.
#' nsite <- 50
#' nseason <- 6
#' nsecondary <- 4
#'
#' b0 <- -0.7
#' theta <- 1
#'
#' occ_prob <- plogis(b0)
#' occ_prob_theta <- plogis(b0 + theta)
#'
#' z <- matrix(
#'   NA,
#'   ncol = nseason,
#'   nrow = nsite
#' )
#' z[,1] <- rbinom(
#'   nsite,
#'   1,
#'   occ_prob
#' )
#' for(t in 2:nseason){
#'   z[,t] <- rbinom(
#'     nsite,
#'     1,
#'     (1 - z[,t-1]) * occ_prob + z[,t-1] * occ_prob_theta
#'   )
#' }
#'
#' # Let's assume we have different observers that go out and
#' #  sample the sites, and that precipitation at sites can make
#' #  it more difficult to detect the species.
#'
#'
#'
#' # Make the precipitation matrix. Each season needs
#' # nsecondary columns
#'
#' precip <- matrix(
#'   rnorm(nsite * nseason * nsecondary),
#'   nrow = nsite,
#'   ncol = nseason * nsecondary
#' )
#' # give them more helpful names
#' colnames(precip) <- paste0(
#'   "precip_",
#'   rep(1:nseason, each = nsecondary),
#'   "_",
#'   rep(1:nsecondary, nseason)
#' )
#'
#'
#' observers <- matrix(
#'   sample(LETTERS[1:3], nsite * nseason * nsecondary, replace = TRUE),
#'   nrow = nsite,
#'   ncol = nseason * nsecondary
#' )
#' colnames(observers) <- paste0(
#'   "observer_",
#'   rep(1:nseason, each = nsecondary),
#'   "_",
#'   rep(1:nsecondary, nseason)
#' )
#'
#' # simulate detection / non-detection data
#' y <- array(
#'   NA,
#'   dim = c(nsite, nseason, nsecondary)
#' )
#'
#' # create the regression coefficients
#' a0 <- 0.5 #intercept
#' a1 <- -0.6 # precip
#' a2 <- -1 # observer B
#' a3 <- 0.5 # observer C
#'
#' alphas <- c(a0,a1,a2,a3)
#'
#' for(t in 1:nseason){
#'   for(j in 1:nsecondary){
#'     # temporary design matrix
#'     dmat <- data.frame(
#'       precip = precip[,grep(paste0("_",t,"_",j), colnames(precip))],
#'       observer = factor(
#'         observers[,grep(paste0("_",t,"_",j), colnames(observers))],
#'         levels = LETTERS[1:3]
#'       )
#'     )
#'     dmat <- model.matrix(~precip + observer, data = dmat)
#'     logit_prob <- dmat %*% alphas
#'     my_prob <- plogis(logit_prob)
#'     y[,t,j] <- rbinom(
#'       nsite,
#'       1,
#'       my_prob * z[,t]
#'     )
#'   }
#' }
#' # make observers a data.frame
#' observers <- as.data.frame(observers)
#'
#' # make all columns factors with shared levels
#' observers <- as.data.frame(
#'   lapply(
#'     observers,
#'     function(x) factor(x, levels = LETTERS[1:3])
#'   )
#' )
#'
#' season_frame <- list(
#'   precip = precip,
#'   observer = observers
#' )
#'
#'
#' m4 <- auto_occ(
#'   formula = ~precip + observer ~1,
#'   y = y,
#'   det_covs = season_frame
#' )
#'
#' summary(m4)


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

  if(!is.numeric(level)){
    stop("level must be a numeric")
  }
  if(level >= 1 | level <= 0){
    stop("level must be a number between 0 and 1")
  }
  if(!is.array(y) | length(dim(y)) != 3){
    stop("y must be a site by primary sampling period by secondary sampling period array. See ?autoOcc::format_y().")
  }

  # check to see if any sites need to get dropped
  na_count <- apply(
    is.na(y),
    1, sum
  )
  to_go <- which(
    na_count == prod(dim(y)[2:3])
  )
  if(length(to_go) > 0){
    warning(
      paste0(
        "\nSome sites have no data.\nRemoved sites at rows: ",
        ifelse(
          length(to_go) == 1,
            to_go,
            paste0(to_go, collapse = ", ")
        ),"\n"
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
    x = det_covs,
    my_formula = rho_formula,
    type = "rho",
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
  # p value via wald statistic
  mle_table$p <- pnorm(-abs(mle_table$Est)/ mle_table$SE) * 2
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
