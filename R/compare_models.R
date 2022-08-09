#' @title Compare auto_occ models via AIC
#'
#' @description Given a (preferably named) list of models, calculate AIC,
#'  delta AIC, and cumulative model weight.
#'
#'
#' @param model_list A list of models of class auto_off_fit. If the list is not
#' named then the models will be given a numeric identifier based on their location
#' in the list.
#'
#' @param add_formula Whether you want the formula of each model added to the
#' data.frame. Defaults to FALSE.
#'
#' @param digits the number of digits to be reported for each numeric column
#' in the data.frame. Exists as an easy way to generate what could be a 'close
#' to publication ready' AIC results table for a given analysis. Defaults to NULL.
#'
#' @details In order to generate a data.frame of model selection results the
#' response variable (i.e., the y array) must be identical among models. The
#' models added to the \code{model_list} must all be unique as well.
#'
#' @returns
#' A data.frame with a number of rows equal to the number of models in
#' \code{model_list}. If models are not named, they will be named in the
#' order they are input (e.g., \code{m1 == model_list[[1]])}). The default
#' data.frame includes:
#'
#' \itemize{
#'   \item{model}{The model names}
#'   \item{npar}{The number of parameters in the model}
#'   \item{AIC}{AIC score of a fitted model}
#'   \item{delta}{delta AIC}
#'   \item{AICwt}{The AIC weight of each model}
#'   \item{cumltvWt}{The cumulative AIC weight}
#' }
#'
#' If \code{add_formula = TRUE}, an additional formula column is added to the
#' data.frame that includes the fitted model call as a character object.
#' This defaults to \code{FALSE} as formulas can become quite long, making it
#' difficult to view the data.frame easily in R.
#'
#' @examples
#'
#' data("opossum_det_hist")
#' data("opossum_covariates")
#' odh <- opossum_det_hist
#' oc <- opossum_covariates
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
#' oc <- head(oc,30)
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
#' # drop site names
#' oc_scaled <- oc_scaled[,-1]
#' # so many warnings
#' m1 <- suppressWarnings(
#'   auto_occ(
#'     ~1~1,
#'     y = opossum_y
#'   )
#' )
#' m2 <- suppressWarnings(
#'   auto_occ(
#'     ~Impervious~Impervious,
#'     y = opossum_y,
#'     det_covs = oc_scaled,
#'     occ_covs = oc_scaled
#'   )
#' )
#' m3 <- suppressWarnings(
#'   auto_occ(
#'     ~Impervious + Income~Impervious + Income,
#'     y = opossum_y,
#'     det_covs = oc_scaled,
#'     occ_covs = oc_scaled
#'   )
#' )
#' m4 <- suppressWarnings(
#'   auto_occ(
#'     ~Income~Income,
#'     y = opossum_y,
#'     det_covs = oc_scaled,
#'     occ_covs = oc_scaled
#'   )
#' )
#' # Added as an unnamed list
#' my_aic_results <- compare_models(
#'   list(m1, m2,m3,m4)
#' )
#'
#' # or you can do a named list
#' also_my_aic_results <- compare_models(
#'   list(
#'     intercept_only = m1,
#'     impervious = m2,
#'     impervious_income = m3,
#'     income = m4
#'   )
#' )
#'
#' # Add in formulas and reduce down to two significant digits
#' pretty_results <- compare_models(
#'   list(m1, m2, m3, m4),
#'   add_formula = TRUE,
#'   digits = 2
#' )
#'
#'
#' @export

compare_models <- function(
    model_list,
    add_formula = FALSE,
    digits = NULL){

  if(!inherits(add_formula,"logical")){
    stop("add_formula must be a logical statement (i.e., TRUE or FALSE)")
  }
  if(!is.null(digits)){
    if(!inherits(digits,"numeric")){
      stop("digits must be a numeric")
    }
    if(digits <=0){
      stop("digits must be positive")
    }
  }
  if(!inherits(model_list, "list")){
    stop("model_list must be a list of auto_occ models.")
  }
  if(!all(sapply(model_list, function(x) inherits(x,"auto_occ_fit")))){
    stop("All objects in model_list must be of class auto_occ_fit")
  }
  # check if duplicate formulas
  tmp <- lapply(
    model_list,
    function(x) as.character(x@formula)
  )
  tmp <- sapply(
    tmp,
    function(x) paste0(
      c(x[2]," ~", x[3]),
      collapse = ""
    )
  )
  if(any(duplicated(tmp))){
    stop("Models with identical formulas added to model_list")
  }
  # quick check to make sure response variable is the same across
  # models
  my_y <- lapply(
    model_list,
    function(x){
      as.vector(x@y)
    }
  )
  y_length <- length(unique(sapply(my_y, length)))
  if(y_length > 1){
    stop("y array has different dimensions among models")
  }
  my_y <- do.call(
        "cbind",
        my_y
      )
  my_y <- apply(
    my_y,
    1,
    function(x) length(unique(x))
  )
  if(!all(my_y == 1)){
    stop("y array differs among models")
  }


  # Check for convergence
  conv_vals <- sapply(model_list, function(x) x@opt$convergence)

  all_conv <- isTRUE(
    all.equal(
      unname(conv_vals),
      rep(0, length(conv_vals))
    )
  )
  if(!all_conv){
    warning("Removing models in model_list that did not converge.")
    to_go <- which(conv_vals != 0L)
    model_list <- model_list[-to_go]
  }

  # get model names
  mnames <- names(model_list)
  if(is.null(mnames)){
    mnames <- paste0("m", 1:length(model_list))
  }
  model_df <- data.frame(
    model = mnames,
    npar = sapply(model_list, function(x) nrow(x@estimates)),
    AIC = sapply(model_list, function(x) x@AIC)
  )
  if(add_formula){
    model_df$formula <- tmp
  }
  # delta AIC and whatnot
  model_df$delta <- model_df$AIC - min(model_df$AIC)
  model_df <- model_df[order(model_df$AIC),]
  model_df$AICwt <- exp(-0.5 * model_df$delta) /
    sum(exp(-0.5 * model_df$delta))
  model_df$cumltvWt <- cumsum(model_df$AICwt)
  if(!is.null(digits)){
    model_df$AIC <- round(model_df$AIC, digits = digits)
    model_df$delta <- round(model_df$delta, digits = digits)
    model_df$AICwt <- round(model_df$AICwt, digits = digits)
    model_df$cumltvWt <- round(model_df$cumltvWt, digits = digits)
  }
  row.names(model_df) <- NULL
  return(model_df)

}
