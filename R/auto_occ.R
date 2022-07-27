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
#' is sites, the second dimension denotes primary sampling periods, and the third
#' dimension denotes the secondary sampling periods within each primary sampling
#' period. If the species was detected on a given survey, that element would receive
#' a 1, otherwise it is 0. If sampling did not occur for a given sampling period, those
#' elements should be NA.
#'
#' @param det_covs A named list of detection covariates. See details for additional
#'   information on how to format this list.
#'
#' @param occ_covs Either a data.frame or named list of occupancy covariates. If
#' there are temporally varying covariates, use a list. If there are no temporally
#' varying covariates, use a data.frame. See details for additional infomration
#' for how to format this.
#'
#' @param method The optimization method used by \code{\link[stats]{optim}}.
#'
#' @param level The confidence interval size to be calculated for model parameters.
#'   Defaults to 0.95.
#'
#' @param ... additional arguments. Not used yet.
#'
#' @details
#'
#' The autologistic formulation of a standard occupancy model is a simplified
#' version of a dynamic occupancy model that makes inference on species
#' patterns of occupancy from one time period to the next such that species
#' presence at a site during one time step may modify the probability a species
#' occupies that same site in the following time step. This is done through the
#' inclusion of a single logit-scale autoregressive term into the model: \eqn{\theta}.
#' For \eqn{t} in \eqn{1 \dots T} and \eqn{i} in \eqn{1 \dots I} sites, let
#' \eqn{z_{i,t}} be a species latent occupancy state. For T=1, we do not
#' know if a species occupied a site in a previous time step
#' (as sampling has not occurred). Thus, to generate the occupancy probability
#' at the first time step we have two nearly identical logit linear predictors.
#' Assuming an intercept only model this would be:
#'
#' \deqn{logit(\psi_{i,t=1}) = \beta_0}
#' \deqn{logit(\phi_{i,t=1}) = \beta_0 + \theta}
#'
#' where the second equation is just the same as the first but with the added
#' auto-regressive term that is informed by the data from the remaining
#' seasons. From this, we can derive the equilibrium, or expected,
#' occupancy at the first time step such that:
#'
#' \deqn{E(\psi)}
#'
#'  We don't know if
#' .
#' If we assume equilibrium at T = 1, then let \eqn{\psi_{i,t=1}} be the
#'
#' @export
#' @importFrom stats as.formula
#' @importFrom stats optim
#' @importFrom stats qnorm
#' @importFrom methods new

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

  # get sites, seasons, and reps from y array
  nsite <- dim(y)[1]
  nseason <- dim(y)[2]
  nrep <- dim(y)[3]

  # get the design matrices together
  if(is.null(occ_covs)){
    occ_covs <- data.frame(intercept = rep(1, nsite))
  }
  occ_dm <- get_dm(
    occ_covs,
    psi_formula,
    "psi",
    y=y
  )

  if(is.null(det_covs)){
    det_covs <- list(
      intercept = data.frame(
        matrix(
          1,
          ncol = nseason,
          nrow = nsite
        )
      )
    )
  }
  rho_dm <- get_dm(det_covs, rho_formula, "rho", y=y)

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
    occcovs = occ_covs
  )
  return(to_return)

}
