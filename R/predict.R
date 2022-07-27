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
#' @param newdata A data.frame of covariates to make predictions from.
#' See details for more information.
#' @param object Object of class inheriting from \code{"auto_occ_fit"}.
#'
#' @param type Either \code{"psi"} for occupancy or \code{"rho"} for detection.
#' See details for how expected occupancy is derived from this model object.
#'
#' @param backTransform Should predictions be converted back to the probability
#' scale. Defaults to \code{TRUE}.
#'
#' @param level Tolerance / confidence level for predictions. Defaults to \code{0.95}.
#'
#' @aliases predict,auto_occ_fit-method
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.offset
#'
#' @export

setMethod(
  "predict",
  signature(object = "auto_occ_fit"),
  function(object,type, newdata = NULL,backTransform = TRUE,level = 0.95){
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
    data <- ifelse(
      type == "psi",
      object@occcovs,
      object@detcovs
    )
    names(data) <- ifelse(
      type == "psi",
      names(object@occcovs),
      names(object@detcovs)
    )
    # check if temporally varying psi
    if(all(sapply(sapply(data, nrow),is.null))){
      data <- data.frame(data)
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
    # without theta
    e1 <- as.vector(cbind(X,0) %*% est) + offset
    # with theta
    e2 <- as.vector(cbind(X,1) %*% est) + offset
    # calculate expected occupancy
    e3 <- plogis(e1) / (plogis(e1) + (1 - plogis(e2)))
    # and transform back to logit scale
    e <- log(e3/(1-e3))
    # get the SE
    v <- cbind(X,1) %*% covMat %*% t(cbind(X,1))
    pred_se <- sqrt(diag(v))

    my_levels <-(1 - level)/2
    lower_ci <- qnorm(
      my_levels,
      e,
      pred_se
    )
    upper_ci <-qnorm(
      1 - my_levels,
      e,
      pred_se
    )

    if(backTransform){
      pred_frame <- data.frame(
        estimate = plogis(e),
        lower = plogis(lower_ci),
        upper = plogis(upper_ci),
        SE = pred_se
      )
    } else {
      pred_frame <- data.frame(
        estimate = e,
        lower = lower_ci,
        upper = upper_ci,
        SE = pred_se
      )
      return(pred_frame)
    }
  }
)
