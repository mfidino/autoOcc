#' @title Get variance covariance matrix from autologistic occupancy model.
#'
#' @description Given a \code{"auto_occ_fit"} model. Derive the variance /
#' covariance matrix for either occupancy (\code{"psi"}) or detection (\code{"rho"}).
#' @rdname vcov.auto_occ_fit
#'
#' @method vcov auto_occ_fit
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#'
#' @param type Either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' If type is not supplied, provide entire variance covariance matrix.
#'
#' @param ... Additional arguments. Not currently used.
#' @docType methods
#'
#' @details This function is largely used for making model predictions via
#' \code{\link{predict.auto_occ_fit}}.
#'
#' @returns
#' A square matrix of the estimated covariances between the parameter estimates
#' in the linear predictor of the model.
#'
#' @examples
#'
#'  data("opossum_det_hist")
#'  data("opossum_covariates")
#'  odh <- opossum_det_hist
#'  oc <- opossum_covariates
#'  # only grabbing the first 30 data.points from each season
#'  # for this example
#'  odh <- split(
#'    odh,
#'    factor(odh$Season, levels = unique(odh$Season))
#'  )
#'  odh <- do.call(
#'    "rbind.data.frame",
#'    lapply(
#'      odh,
#'      function(x) head(x,30)
#'    )
#'  )
#'  oc <- head(oc,30)
#'
#'  # function to generate detection history
#'  opossum_y <- autoOcc::format_y(
#'    x = odh,
#'    site_column = "Site",
#'    time_column = "Season",
#'    history_columns = "^Week", # start with Week
#'    report = FALSE
#'  )
#'
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
#'  # drop site names
#'  oc_scaled <- oc_scaled[,-1]
#'  # suppressing warnings because of missing data
#'  m1 <- suppressWarnings(
#'    auto_occ(
#'      ~Impervious~Impervious+Income,
#'      y = opossum_y,
#'      det_covs = oc_scaled,
#'      occ_covs = oc_scaled
#'    )
#'  )
#'  my_vcov <- vcov(m1)
#' @export

   vcov.auto_occ_fit <- function (object, type, ...)
          {
            if (is.null(object@opt$hessian)) {
              stop("Hessian was not computed for this model.")
            }
            v <- solve(object@opt$hessian)
            pnames <- object@estimates$parameter
            pnames <- substr(pnames,1,3)
            tmp <- rle(pnames)
            tmp <- sapply(tmp$lengths, function(x) 1:x)
            tmp <- unlist(tmp)
            pnames <- paste0(
              substr(
                pnames,1,3
              ),"_",tmp
            )
            rownames(v) <- colnames(v) <- pnames
            if (missing(type)) {
              return (v)
            } else {
              inds <- grep(type,colnames(v))
              return (v[inds, inds, drop = FALSE])
            }
          }


#' @title Calculates confidence intervals for model parameters.
#'
#' @rdname confint.auto_occ_fit
#'
#' @method confint auto_occ_fit
#'
#' @docType methods
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#' @param parm a specification of which parameters to be given confidence
#' intervals. See details for additional details.
#' @param level the confidence level required.
#' @param type must be either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' @param ... Additional arguments to be supplied. Not currently used.
#'
#' @details Largely used to fill in the rest of the estimates data.frame
#' for objects of class \code{auto_occ_fit}. However, this could be used
#' to generate confidence intervals of different levels for a given
#' \code{auto_occ_fit}.
#'
#' If \code{parm} is not null, it should be a character vector. While the
#' theta parameter estimate is titled \code{"psi - theta"}, \code{confint}
#' uses the \code{type} argument to specify the first part of a parameter
#' name. Thus, to get only \code{theta} back, \code{type} must equal \code{"psi"}
#' and \code{parm} must be \code{"theta"}. See last example.
#'
#' @returns
#' A data.frame with 5 columns that contains the parameter estimates,
#' standard error, lower, and upper confidence limits for a given
#' logit-linear predictor, which is specified by \code{type}.
#'
#' @examples
#' data("opossum_det_hist")
#' odh <- opossum_det_hist
#'
#' # only grabbing the first 30 data.points from each season
#' # for this example
#' odh <- split(
#'   odh,
#'   factor(odh$Season, levels = unique(odh$Season))
#' )
#' odh <- do.call(
#'   "rbind.data.frame",
#'   lapply(
#'     odh,
#'     function(x) head(x,30)
#'   )
#' )
#'
#'
#' # function to generate detection history
#' opossum_y <- autoOcc::format_y(
#'   x = odh,
#'   site_column = "Site",
#'   time_column = "Season",
#'   history_columns = "^Week", # start with Week
#'   report = FALSE
#' )
#'
#' # suppressing warnings because of missing data
#' m1 <- suppressWarnings(
#'   auto_occ(
#'     ~1~1,
#'     y = opossum_y
#'   )
#' )
#' # get 80% confidence intervals for occupancy
#' psi_80 <- confint(
#'   m1,
#'   level = 0.80,
#'   type = "psi"
#' )
#' # and the same for detection
#' rho_80 <- confint(
#'   m1,
#'   level = 0.80,
#'   type = "rho"
#' )
#'
#' # get just the theta estimate from occupancy,
#' theta_estimate <- confint(
#'   m1,
#'   parm = "theta",
#'   type = "psi"
#' )
#'
#'
#' @docType methods
#' @export
#'
#'
#'

confint.auto_occ_fit <- function(object,parm, level = 0.95,type,...){
  if(!is.numeric(level)){
    stop("level must be a numeric")
  }
  if(level >= 1 | level <=0){
    stop("level must be between 0 and 1")
  }

  if(missing(type)){
    stop("Must specify type as either 'psi' or 'rho'")
  }
  if(missing(parm)){
    my_parm <- grep(type, object@estimates$parameter)
  }else{
    to_grep <- paste0(type," - ", parm, collapse = "|")
    my_parm <- grep(to_grep,object@estimates$parameter)
    if(length(my_parm)!= length(parm)){
      stop(
          "One or more of parm supplied are not in model."
      )
    }
  }

  ests <- object@estimates[my_parm,]
  lwr <- qnorm(
    (1 - level)/2,
    ests$Est,
    ests$SE
  )
  upr <- qnorm(
    (1 - (1 - level)/2),
    ests$Est,
    ests$SE
  )
  return(
    data.frame(
      parameter = ests$parameter,
      Est = ests$Est,
      lower = lwr,
      upper = upr
    )
  )

}

