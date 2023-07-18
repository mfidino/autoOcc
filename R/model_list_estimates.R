#' @title Collect all parameter estimates from a model list
#'
#' @description Takes a list of auto_occ_fit objects and retrieves all parameter estimates from them.
#'
#'
#' @param model_list A list of auto_occ_fit objects. If list is named, those names will be retained.
#' See details for additional information.
#'
#'
#' @details
#'
#' This function will return a named list object with two elements. The first
#' element, estimates, is a data.frame with parameter estimates from each model.
#' The number of rows will be equal to the total number of unique parameters across
#' models in the model set. The second element, ses, is a data.frame of standard error
#' estimates from each model. Like the first element, the number of rows are euqal to the
#' total number of unique parameters, with \code{NA} values present if that specific
#' parameter is not in a model.
#'
#' For the \code{model_list}, the data set must be the same across each model or else
#' this function will fail.
#'
#' @examples
#'  data(opossum_det_hist)
#'  data(opossum_covariates)
#'
#'  opossum_y <- autoOcc::format_y(
#'  x = opossum_det_hist,
#'  site_column = "Site",
#'  time_column = "Season",
#'  history_columns = "^Week" # start with Week
#'  )
#'
#'  # going down to just 30 sites for this example
#'  opossum_y <- opossum_y[1:30,,]
#'
#'  # scale the covariates for analysis
#'  oc_scaled <- as.data.frame(
#'    lapply(
#'      opossum_covariates[1:30,],
#'      function(x){
#'        if(is.numeric(x)){
#'          scale(x)
#'        }else{
#'          x
#'        }
#'      }
#'    )
#'  )
#'  m1 <- auto_occ(
#'    ~1~1,
#'    y = opossum_y,
#'    det_covs = oc_scaled,
#'    occ_covs = oc_scaled
#'  )
#'
#'  m2 <- auto_occ(
#'    ~Income~Income,
#'    y = opossum_y,
#'    det_covs = oc_scaled,
#'    occ_covs = oc_scaled
#'  )
#'  m3 <- auto_occ(
#'    ~1~Income,
#'    y = opossum_y,
#'    det_covs = oc_scaled,
#'    occ_covs = oc_scaled
#'  )
#'
#'  m4 <- auto_occ(
#'    ~Income~1,
#'    y = opossum_y,
#'    det_covs = oc_scaled,
#'    occ_covs = oc_scaled
#'  )
#'
#'  all_ests <- model_list_estimates(
#'    list(m1, m2, m3, m4)
#'  )
#'
#'  # you can also name the list
#'  all_ests <- model_list_estimates(
#'    list(m1 = m1, m2 = m2, m3 = m3, m4 = m4)
#'  )
#'
#'
#' @export

model_list_estimates <- function(model_list){
  # check if all are auto occ objects
  if(!is.list(model_list)){
    stop("model_list must be a list object.")
  }
  if(
    !all(
      sapply(
        model_list,
        class
      ) == "auto_occ_fit"
    )
  ){
    stop("All objects in model_list must be of class auto_occ_fit")
  }
  if(
    !all(
      sapply(
        model_list,
        function(x)x@opt$convergence == 0
      )
    )
  ){
    to_go <-  sapply(
      model_list,
      function(x)x@opt$convergence
    )
    to_go <- which(to_go != 0)
    if(!is.null(names(model_list))){
      conv_names <- names(model_list)
    } else {
      conv_names <-paste0("model ", seq_len(length(model_list)))
    }
    warning(
      paste0(
        "Removing models that did not converge.\n",
        "Removed: ", toString(conv_names[to_go])
      )
    )
    model_list <- model_list[-to_go]
  }
  if(length(model_list) == 1){
    stop("There is only one model, use summary(<insert your model object here>) to retrieve parameter estimates instead.")
  }
  for(i in 2:length(model_list)){
    data_check <- identical(
      model_list[[i-1]]@y,
      model_list[[i]]@y
    )
    if(!data_check){
      stop("Data sets differ across model objects.")
    }
  }
  # get names of all parameters
  parms <- lapply(
    model_list,
    function(x){
      x@estimates$parameter
    }
  )
  # get all unique parameters
  all_parms <- unique(
    unlist(
      parms
    )
  )

  model_df <- matrix(
    ncol = length(all_parms),
    nrow = length(model_list)
  )
  colnames(model_df) <- all_parms
  row.names(model_df) <- names(model_list)
  if(is.null(row.names(model_df))){
    row.names(model_df) <- seq_len(nrow(model_df))
  }
  ests <- list(
    estimates = as.data.frame(model_df),
    ses = as.data.frame(model_df)
  )
  for(i in 1:length(model_list)){
    ests$estimates[i,parms[[i]]] <- model_list[[i]]@estimates$Est
    ests$ses[i,parms[[i]]] <- model_list[[i]]@estimates$SE
  }
  return(ests)
}
