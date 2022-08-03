#' @title Get variance covariance matrix from autologistic occupancy model.
#'
#' @description Given a \code{"auto_occ_fit"} model. Derive the variance /
#' covariance matrix for either occupancy (\code{"psi}) or detection (\code{"rho}).
#' @rdname vcov-methods
#'
#' @method vcov auto_occ_fit
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#'
#' @param type Either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' If type is not supplied, provide entire variance covariance matrix.
#'
#' @param ... Additional arguments. Not currently used.
#' @aliases vcov,auto_occ_fit-method
#' @docType methods
#'
#' @export


setMethod("vcov", "auto_occ_fit",
          function (object, type, ...)
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
          })


#' @title Calculate confidence intervals for model parameters.
#'
#' @rdname confint-methods
#'
#' @method confint auto_occ_fit
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#' @param parm a specification of which parameters to be given confidence
#' intervals. See details for additional details.
#' @param level the confidence level required.
#' @param type must be either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' @param ... Additional arguments to be supplied.

#' @aliases confint,auto_occ_fit-method
#' @docType methods
#' @export
#'
#'
#'
setMethod("confint","auto_occ_fit",
          function(object,parm, level = 0.95,type,...){

  if(missing(type)){
    stop("Must specify type as either 'psi' or 'rho'")
  }
  if(missing(parm)){
    parm <- grep(type, object@estimates$parameter)
  }else{
    to_grep <- paste0("psi - ", parm, collapse = "|")
    parm <- which(object@estimates$parameter %in% to_grep)
    if(length(parm)<1){
      stop(
        paste0(
          "One or more of parm supplied are not in model.",
        )
      )
    }
  }

  ests <- object@estimates[parm,]
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

})

