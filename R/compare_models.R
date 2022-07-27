#' @title Compare auto_occ models via AIC
#'
#' @description Given a (preferabbly named) list of models, calculate AIC,
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
#'
#' @export

compare_models <- function(
    model_list,
    add_formula = FALSE,
    digits = NULL){

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
