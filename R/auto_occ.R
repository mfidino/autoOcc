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
