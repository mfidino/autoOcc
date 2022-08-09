setClassUnion("dataframe_or_list", c("data.frame", "list"))
setClassUnion("numeric_or_null", c("numeric", "NULL"))

#' @title An S4 class to represent the fitted autologistic occupancy model
#
#' @slot fitType The class of the fitted model object.
#' @slot call The formula call.
#' @slot formula Double right-hand side formula describing covariates of
#' detection (rho) and occupancy (psi) in that order.
#' @slot y A three-dimensional array of species detections. The first dimension
#' is sites, the second dimension denotes primary sampling periods, and the third
#' dimension denotes the secondary sampling periods within each primary sampling
#' period. If the species was detected on a given survey, that element would receive
#' a 1, otherwise it is 0. If sampling did not occur for a given sampling period, those
#' elements should be NA.
#' @slot estimates A data.frame of the estimated model parameters and their uncertainty.
#' @slot AIC the AIC of the fitted model.
#' @slot opt the list returned from the \code{\link[stats]{optim}} call.
#' @slot negLogLike The negative log likelihood of the data given the parameters.
#' @slot nllFun the likelihood function supplied to \code{\link[stats]{optim}}.
#' @slot detcovs A list of the detection covariates supplied to \code{auto_occ}.
#' @slot occcovs A list of the occupancy covariates supplied to \code{auto_occ}.
#' @slot sites_removed An optional vector that denotes which row indices in y had
#' no data across all sampled seasons.
#' @export
#' @rdname auto_occ_fit

#' @importFrom methods new

setClass("auto_occ_fit",
         representation(fitType = "character",
                        call = "call",
                        formula = "formula",
                        y = "array",
                        estimates = "data.frame",
                        AIC = "numeric",
                        opt = "list",
                        negLogLike = "numeric",
                        nllFun = "function",
                        detcovs = "list",
                        occcovs = "dataframe_or_list",
                        sites_removed = "numeric_or_null"))

# constructor for auto_occ_fit objects
auto_occ_fit <- function(fitType, call, formula, y,
                        estimates, AIC, opt, negLogLike, nllFun,
                        detcovs, occcovs, sites_removed){
  aofit <- new(
    "auto_occ_fit",
    fitType = fitType,
    call = call,
    formula = formula,
    y = y,
    estimates = estimates,
    AIC = AIC,
    opt = opt,
    negLogLike = negLogLike,
    nllFun = nllFun,
    detcovs = detcovs,
    occcovs = occcovs,
    sites_removed = sites_removed
  )
  return(aofit)
}


#' @title Summarizing autologistic occupancy models
#'
#' @rdname summary.auto_occ_fit
#'
#' @docType methods
#'
#' @method summary auto_occ_fit
#'
#' @description summary method for class \code{"auto_occ_fit"}
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#'
#' @param ... Additional arguments. Not currently used.
#'
#'
#' @returns
#' \code{summary} returns an object of class \code{"auto_occ_summary"}.
#' See \code{\linkS4class{auto_occ_summary}} for additional
#' details.
#'
#' @examples
#' data("opossum_det_hist")
#' opossum_y <- format_y(
#'   opossum_det_hist,
#'   "Site",
#'   "Season",
#'   "^Week",
#'   report = FALSE
#' )
#' # suppressing warnings as there are sites with 0 data
#' m1 <- suppressWarnings(
#'   auto_occ(~1~1, y = opossum_y)
#' )
#' summary(m1)
#'
#' @export

summary.auto_occ_fit <- function(object, ...)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  summaryOut <- summary(object@estimates)
  cat("\noptim convergence code:", object@opt$convergence)
  cat("\noptim iterations:", object@opt$counts[1], "\n")
  if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
  invisible(summaryOut)
  psis <- object@estimates[grep("^psi - ", object@estimates$parameter),]
  cat("\nOccupancy estimates:\n\n")
  print(psis, digits = 3)
  cat("\nNote: psi - theta is the autologistic term\n")
  rhos <- object@estimates[grep("^rho - ", object@estimates$parameter),]
  cat("\nDetection estimates:\n\n")
  print(rhos, digits = 3)
  cat("\n")
  cat("AIC:", object@AIC,"\n")
  to_return <- new(
    "auto_occ_summary",
    call = object@call,
    optim_convergence_code = object@opt$convergence,
    optim_iterations = object@opt$counts[1],
    psi = psis,
    rho = rhos,
    AIC = object@AIC
  )
  invisible(to_return)
}


#' @importFrom methods show
#' @noRd
setMethod("show", "auto_occ_fit", function(object)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  show(object@estimates)
  cat("\nAIC:", object@AIC,"\n")
  if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
})



