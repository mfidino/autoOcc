#' @title Makes predictions from an autologistic occupancy model
#'
#' @rdname predict-methods
#'
#' @docType methods
#'
#' @method predict auto_occ_fit
#'
#' @description Predicted values based on a \code{auto_occ_fit} model object.
#'
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#'
#' @param newdata Optionally, a data.frame of covariates to make predictions from.
#' If omitted, the fitted linear predictors are used.
#'
#' @param type Either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' See details for how expected occupancy is derived from this model object.
#'
#' @param level Tolerance / confidence level for predictions. Defaults to \code{0.95}.
#'
#' @param nsim The number of parmater simulations to be made via the
#' models estimated parameters, variance covariance matrix,
#' and \code{\link[mvtnorm]{rmvnorm}}. Defaults to 3000.
#'
#' @param seed The random seed to set for simulations, defaults to \code{NULL}.
#'
#' @aliases predict,auto_occ_fit-method
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.offset
#' @importFrom stats quantile
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'@details
#'
#' This function takes the output of an \code{\link[autoOcc]{auto_occ}}
#' model and will generate predictions for either occupancy (\code{type = "psi"})
#' or detection (\code{type = "rho"}). If \code{newdata} is supplied, all
#' covariates supplied to that level of the model must be present except for
#' the case of the autologistic occupancy term (that is handled). To
#' approximate uncertainty in model estimates, this function first uses
#' \code{\link[mvtnorm]{rmvnorm}} to generate parameter estimates from
#' multivariate distribution with a mean vector equal to the model parameters
#' and a variance covariance matrix supplied by the \code{\link[autoOcc]{auto_occ}}
#' model. This method was chosen as it creates functionally similar uncertainty
#' estimates as a Bayesian autologistic occupancy model.
#'
#' For occupancy, each set of simulated parameters is used to
#' generate two logit-linear predictions for the supplied covariates:
#'
#' \deqn{\Large \mathrm{logit}(\psi_a) = \beta_0 + \beta_1 \times x_1 + \dots + \beta_n \times x_n}
#'
#'and
#'
#'\deqn{\Large \mathrm{logit}(\psi_b) = \beta_0 + \beta_1 \times x_1 + \dots + \beta_n \times x_n + \theta}
#'
#' where \eqn{\theta} is the estimated autologistic term. Following this, the expected occupancy of an autologisitic occupancy model is
#'
#' \deqn{\huge
#'  \frac{
#'    \mathrm{ilogit}(\psi_a)
#'  }{
#'  \mathrm{ilogit}(\psi_a) + (1 - \mathrm{ilogit}(\psi_b))
#'  }
#'}
#' which is similar to the expected occupancy of a dynamic occupancy model (\eqn{\gamma \div (\gamma + \epsilon)})
#' where \eqn{\gamma} is colonization and \eqn{\epsilon} is extinction. Following this calculation for all simuated
#' parameter estimates and covariate values, the median estimate and confidence intervals are collected across
#' simulations for each covariate value.
#'
#' Detection predictions are more straight-forward given there is no need to
#' derive the expected value. Simulations are still carried out to create to
#' generate a vector of values for each parameter.
#'
#' If \code{newdata} is supplied, then  a \code{data.frame} will be returned that has the
#' same number of rows as \code{newdata} with three columns: \code{estimate},
#' \code{lower}, and \code{upper}. \code{estimate} is the median estimate, \code{lower} is the lower
#' confidence level, \code{upper} is the upper confidence level.
#'
#'  If \code{newdata} is not supplied, then the output will
#' depend on the covariates the model was fitted with:
#'\itemize{
#'  \item{For occupancy with no temporal variation in covariates}{
#'    a data.frame will be returned
#'    that lines up with the covariates provided when the occupancy model was fit.
#'  }
#'  \item{For occupancy with temporal variation across primary sampling periods}{
#'    a list of data frames will
#'    be returned, one for each primary sampling period.
#'  }
#'  \item{For detection with no temporal variation in covariates}{
#'    a data.frame will be returned
#'    that lines up with the covariates provided when the occupancy model was fit.
#'  }
#'  \item{For detection with temporal variation across primary sampling periods}{
#'    a list of data frames will be returned, one for each primary period.
#'  }
#'  \item{For detection with temporal varaition across secondary observation periods}{
#'    a nested list of data.frames will be returned, one for each primary and secondary sampling period. This
#'    make time a substantial amount of time to do these calculations, and therefore it is not
#'    recommended (i.e., use \code{newdata} instead).
#'  }
#'}
#'

setMethod(
  "predict",
  signature(object = "auto_occ_fit"),
  function(object,type, newdata = NULL,level = 0.95, nsim = 3000, seed = NULL){
    if(missing(type)){
      stop("Supply either 'psi' or 'rho' as type for predictions.")
    }
    if(!type %in% c("psi","rho")){
      stop("type must be either' psi' for occupancy or 'rho' for detection." )
    }
    char <- lapply(object@formula, function(x) {
      paste(deparse(x), collapse = "")
    })
    if(type == "psi"){
      my_formula <- as.formula(paste("~", char[[3]]))
    }else{
      my_formula <- as.formula(char[[2]])
    }
    if(
      is.null(newdata)
    ){
     if(type == "psi"){
       newdata <- object@occcovs
     }else{
       newdata <- object@detcovs
     }
    }

#get the old dm
    if(type == "psi"){
      data <- object@occcovs
    } else {
      data <- object@detcovs
    }
    # check if temporally varying psi
    if(all(sapply(sapply(data, nrow),is.null))){
      data <- as.data.frame(data)
    } else{
      stop("Temporally varying covariates not set up for predictions.")
    }

    mf <- get_dm(data, my_formula = my_formula, type = type, y = object@y)
    fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
    xlevs <- lapply(fac_cols, levels)
    xlevs <- xlevs[names(xlevs) %in% names(mf)]
    nmf <- model.frame(my_formula, newdata, na.action=stats::na.pass, xlev = xlevs)
    X <- model.matrix(my_formula, newdata, xlev=xlevs)
    offset <- model.offset(nmf)

    # get variance covariance matrix
    covMat <- vcov(
      object,
      type = type
    )

    est <- object@estimates
    est <- est$Est[grep(type,est$parameter)]
    if (is.null(offset)){
      offset <- rep(0, nrow(X))
    }

    if(!is.null(seed)){
      set.seed(seed = seed)
    }
    mvn_samples <- mvtnorm::rmvnorm(
      nsim,
      mean=est,
      sigma=covMat,
      method="svd"
    )
    # logit-predictions without theta
    if(type == "psi"){
      e1 <- cbind(X,0) %*% t(mvn_samples)
      e1 <- sweep(e1, 1, offset, FUN = "+")
      # logit-predictions with theta
      e2 <- cbind(X,1) %*% t(mvn_samples)
      e2 <- sweep(e2, 1, offset, FUN = "+")
      # calculate expected occupancy
      e3 <- plogis(e1) / (plogis(e1) + (1 - plogis(e2)))
    }
    if(type == "rho"){
      e2 <- X %*% t(mvn_samples)
      e2 <- sweep(e2, 1, offset, FUN = "+")
      e3 <- plogis(e3)
    }
    my_levels <-(1 - level)/2
    predictions <- t(
      apply(
        e3,
        1,
        stats::quantile,
        probs = c(my_levels, 0.5, 1 - my_levels)
      )
    )
    pred_frame <- data.frame(
      estimate = predictions[,2],
      lower = predictions[,1],
      upper = predictions[,3]
    )
    return(pred_frame)

  }
)

