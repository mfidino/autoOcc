#' @title Makes predictions from an autologistic occupancy model
#'
#' @docType methods
#' @method predict auto_occ_fit

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
#'
#' @param seed The random seed to set for simulations, defaults to \code{NULL}.
#'
#' @param ... additional arguments. Not used yet.
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.offset
#' @importFrom stats quantile
#' @importFrom mvtnorm rmvnorm
#' @rdname predict.auto_occ_fit
#'
#' @export
#' @details
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
#' @returns
#'
#' If \code{newdata} is supplied, this function will return a \code{data.frame}
#' with a number of rows equal to \code{newdata} and three columns.
#' \itemize{
#'   \item{estimate}{Median estimate from the model}
#'   \item{lower}{Lower confidence interval, based on confidence level supplied to \code{level}}
#'   \item{upper}{Upper confidence interval, based on confidence level supplied to \code{level}}
#' }
#'
#' If \code{newdata = NULL}, the same data.frame will be output except the covariates
#' used to fit the model will be added on as additional columns. If there are temporally
#' varying covariates, the outputted data.frame will be in long format.
#'
#' @examples
#' data("opossum_det_hist")
#' data("opossum_covariates")
#' odh <- opossum_det_hist
#' oc <- opossum_covariates
#'
#' # function to generate detection history
#' opossum_y <- autoOcc::format_y(
#'   x = odh,
#'   site_column = "Site",
#'   time_column = "Season",
#'   history_columns = "^Week" # start with Week
#' )
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
#' # dropping site column from oc_scaled
#' oc_scaled <- oc_scaled[,-1]
#' # suppress warnings because opossum_y has sites
#' #  with no data.
#' m1 <- auto_occ(
#'   ~Impervious + Income  ~ Impervious + Income,
#'   y = opossum_y,
#'   det_covs = oc_scaled,
#'   occ_covs = oc_scaled
#' )

#'
#' # first make the prediction data.frame with a realistic
#' #   range based on the actual data and not the scaled data.
#' #   The range(oc$Impervious) is about 18 to 81, so choose 20
#' #   to 80. We do this so that we have nice numbers for plotting.
#' #   Likewise, we scaled all of the other data, so we leave Income
#' #   at it's mean (i.e., 0) for predictions.
#' imperv_real <- data.frame(
#'   Impervious = seq(20,80,0.5),
#'   Income = 0
#' )
#'
#' # We will use imperv_real for plotting purposes, but to make predictions
#' #  we need to scale imperv_real$Impervious in the exact same way we did
#' #  with the fitted model. Thus, we subtract the mean of the actual data
#' #  and divide by the standard deviation.
#' imperv_scaled <- imperv_real
#' imperv_scaled$Impervious <- (
#'   imperv_scaled$Impervious - mean(oc$Impervious)
#' ) / sd(oc$Impervious)
#'
#' # the model prediction across a gradient of Impervious cover
#' opo_imperv <- predict(
#'   object = m1,
#'   type = "psi",
#'   newdata = imperv_scaled
#' )
#'
#' # do the same thing with income
#' income_real <- data.frame(
#'   Impervious = 0,
#'   Income = seq(40000, 160000, by = 500)
#' )
#'
#' income_scaled <- income_real
#' income_scaled$Income <- (income_scaled$Income - mean(oc$Income)) / sd(oc$Income)
#'
#' opo_income <- predict(
#'   object = m1,
#'   type = "psi",
#'   newdata = income_scaled
#' )
#'
#' # plot them out
#' par(mfrow = c(1,2))
#' plot(
#'   opo_imperv$estimate ~ imperv_real$Impervious,
#'   bty = "l",
#'   type = "l",
#'   las = 1,
#'   ylab = "Occupancy",
#'   xlab= "Impervious Cover (%)",
#'   ylim = c(0,1),
#'   lwd = 3
#' )
#' lines(opo_imperv$lower ~ imperv_real$Impervious, lwd = 2, lty = 2)
#' lines(opo_imperv$upper ~ imperv_real$Impervious, lwd = 2, lty = 2)
#' plot(
#'   opo_income$estimate ~ income_real$Income,
#'   bty = "l",
#'   type = "l",
#'   las = 1,
#'   ylab = "Occupancy",
#'   xlab= "Per Capita Income (US Dollar)",
#'   ylim = c(0,1),
#'   lwd = 3
#' )
#' lines(opo_income$lower ~ income_real$Income, lwd = 2, lty = 2)
#' lines(opo_income$upper ~ income_real$Income, lwd = 2, lty = 2)


  predict.auto_occ_fit <- function(object,type, newdata = NULL,level = 0.95, nsim = 3000, seed = NULL,...){
    if(!inherits(object,"auto_occ_fit")){
      stop("model object must be of class auto_occ_fit")
    }
    if(missing(type)){
      stop("Supply either 'psi' or 'rho' as type for predictions.")
    }
    if(!type %in% c("psi","rho")){
      stop("type must be either' psi' for occupancy or 'rho' for detection." )
    }
    if(!is.numeric(level)){
      stop("level must be a numeric")
    }
    if(level >= 1 | level <= 0){
      stop("level must be a number between 0 and 1")
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
      tack_on <- TRUE
     if(type == "psi"){
       newdata <- object@occcovs
     }else{
       newdata <- object@detcovs
     }
      newdata <- lapply(
        newdata,
        as.vector
      )
      newdata <- do.call(
        "cbind.data.frame",
        newdata
      )
      site_vec <- 1:(dim(object@y)[1] + length(object@sites_removed))
      site_vec <- rep(
        site_vec,
        nrow(newdata)/length(site_vec)
      )
      if(length(object@sites_removed)>0){
        newdata <- newdata[-which(site_vec %in% object@sites_removed),]
      }
      newdata <- as.data.frame(newdata)
      newdata <- suppressWarnings(
        factor_df_cols(
          newdata,
          type = type
        )
      )
    } else{
      tack_on <- FALSE
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
    }else{
      data <- lapply(
        data,
        function(k) as.vector(unlist(k))
      )
      data <- do.call(
        "cbind.data.frame",
        data
      )
      data <- as.data.frame(data)
    }
    if(length(object@sites_removed)>0){
      site_vec <- 1:(dim(object@y)[1] + length(object@sites_removed))
      site_vec <- rep(
        site_vec,
        nrow(data)/length(site_vec)
      )
      data <- data[-which(site_vec %in% object@sites_removed),,drop = FALSE]
    }
    data <- suppressWarnings(
      factor_df_cols(
        data,
        type = type
      )
    )

    mf <- suppressWarnings(
      get_dm(
        data,
        my_formula = my_formula,
        type = type,
        y = object@y
      )
    )
    fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
    xlevs <- lapply(fac_cols, levels)
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
      e3 <- plogis(e2)
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
    if(tack_on){
      pred_frame <- cbind.data.frame(
        list(
          pred_frame,
          newdata
        )
      )
    }
    if(identical(paste0(as.character(my_formula),collapse = ""),"~1")){
      pred_frame <- pred_frame[1,c("estimate","lower","upper"),drop = FALSE]
    }
    return(pred_frame)
  }

