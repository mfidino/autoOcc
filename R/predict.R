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
#' @importFrom mvtnorm rmvnorm
#'
#' @export

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
    mvn_samples <- rmvnorm(
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
        quantile,
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

