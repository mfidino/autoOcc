#'
#' Generate design matrix for fitting autologistic occuapncy models (internal).
#'
#' @param x the occupancy or detection list / data.frame.
#'
#' @param formula the specific formula for either occupancy or detection.
#'
#' @param type Either \code{"psi"} for occupancy or \code{"rho"} for detection.
#'
#' @param y A three-dimensional array of species detections. The first dimension
#' is sites, the second dimension denotes primary sampling periods, and the third
#' dimension denotes the secondary sampling periods within each primary sampling
#' period. If the species was detected on a given survey, that element would receive
#' a 1, otherwise it is 0. If sampling did not occur for a given sampling period, those
#' elements should be NA.
#'
#' @param to_drop which row indexes need to be dropped from the analysis,
#' based on sites that do not have any data.
#'
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#'
#' @noRd


get_dm <- function(x, my_formula, type = c("psi","rho"), y, to_drop = NULL){
  # if just a site by nparam matrix for psi
  nsite <- dim(y)[1]
  nseason <- dim(y)[2]
  nrep <- dim(y)[3]

  # if psi does not temporally vary

  if(
    is.data.frame(x) &
    type == "psi"
  ){
    if(length(to_drop)>0){
      x <- x[-to_drop, ,drop = FALSE]
    }
    x <- factor_df_cols(
      x,
      temp_var = FALSE,
      type = type
    )

    to_return <- vector(
      "list",
      length = nseason
    )
    tmp <- try(
        model.frame(
        formula = my_formula,
        data = x,
        na.action = NULL
      ),
      silent = TRUE
    )
    if(inherits(tmp,"try-error")){
      tmp <- gsub(
        "not found\n",
        "not found. This can happen if you did not input temporally varying covariates as a named list or if you left this covariate out of the covariates supplied to the 'occ_covs' argument of auto_occ().\n",
        tmp
      )
      stop(tmp)
    }
    tmp <- model.matrix(
      my_formula,
      tmp
    )
    for(i in 1:nseason){
      to_return[[i]] <- tmp
    }
    return(to_return)
  }
  # if temporally varying psi
  if(
    is.list(x) &
    !is.data.frame(x) &
    type == "psi"
  ){
    # check if there is any temporal variation, and if so
    #  it has the correct number of columns.
    ncols <- sapply(
      x,
      function(k){
        cc <- ncol(k)
        ifelse(is.null(cc), -1, cc)
      }
    )

    for(i in 1:length(ncols)){
      ncols[[i]] <- ncols[[i]] == -1 |
        ncols[[i]] == nseason
    }
    if(sum(ncols) != length(x)){
      baddies <- which(ncols!=1)
      error_report <- paste0(
        "Occupancy covariates must either be a vector,\n",
        "or a matrix/data.frame with a number of columns equal to the\n",
        "number of primary sampling periods.",
        "These covariates are not properly\n",
        "set up: ", paste0(names(x)[baddies] , collapse = ", ")
      )
      stop(error_report)
    }

    if(length(to_drop)>0){
      x <- lapply(
        x,
        function(k){
          if(length(ncol(k)) == 0){
            k[-to_drop]
          }else{
            k[-to_drop, , drop = FALSE]
          }
        }
      )
    }
    for(i in 1:length(x)){
      x[[i]] <- factor_df_cols(
        df = as.data.frame(x[[i]]),
        name = names(x)[i],
        temp_var = TRUE,
        type = type
      )
    }

    to_return <- vector("list", length = nseason)
    for(i in 1:nseason){
      to_return[[i]] <- lapply(
        x,
        function(k){
          if(is.null(ncol(k))| ncol(k) == 1){
            return(k)
          } else {
            return(k[,i])
          }
        }
      )
      to_return[[i]] <- do.call(
        "cbind.data.frame",
        to_return[[i]]
      )

      colnames(to_return[[i]]) <- names(x)
      to_return[[i]] <- model.frame(
        formula = my_formula,
        data = to_return[[i]],
        na.action = NULL
      )
      to_return[[i]] <- model.matrix(
        my_formula,
        to_return[[i]]
      )
    }
    return(to_return)
  }
  if(
    type == "rho"
  ){
    if(length(to_drop)>0){
      x <- lapply(
        x,
        function(k){
          if(is.null(ncol(k))){
            k[-to_drop]
          }else{
            k[-to_drop,]
          }
        }
      )
    }
    #quick checks, give -1 if no temporal variation
    ncols <- sapply(
        x,
        function(k){
          cc <- ncol(k)
          ifelse(is.null(cc), -1, cc)
        }
      )

    for(i in 1:length(ncols)){
      ncols[[i]] <- ncols[[i]] == -1 |
        ncols[[i]] == nseason |
        ncols[[i]] == (nseason * nrep)
    }
    if(sum(ncols) != length(x)){
      baddies <- which(ncols!=1)
      error_report <- paste0(
        "Detection covariates must either be a vector,\n",
        "a matrix/data.frame with a number of columns equal to the\n",
        "number of primary sampling periods, or a matrix/data.frame with\n",
        "a number of columns equal to the number of primary\n",
        "sampling periods times the number of secondary observations\n",
        "within time periods. These covariates are not properly\n",
        "set up: ", paste0(names(x)[baddies] , collapse = ", ")
      )
      stop(error_report)
    }
    # create factors as needed
    for(i in 1:length(x)){
      x[[i]] <- factor_df_cols(
        df = as.data.frame(x[[i]]),
        name = names(x)[i],
        temp_var = TRUE,
        type = type
      )
    }

    to_return <- vector("list", length = nseason)
    names(to_return) <- paste0(
      "time_", 1:nseason
    )
    for(i in 1:nseason){
      for(j in 1:nrep){
        if(j == 1){
          to_return[[i]] <- vector("list", length = nrep)
          names(to_return[[i]]) <- paste0(
            "observation_", 1:nrep
          )
        }
        tmp <-  lapply(
          x,
          function(k){
            # if no temporal variation
            if(is.null(ncol(k)) | ncol(k) == 1){
              return(k)
            }
            # if variation across each sampling period
            if(ncol(k) == nseason){
              return(k[,i])
            } else{
              nx <- ((1:ncol(k) - 0.001) %/% nrep) + 1
              nx <- which(nx == i)[j]
              return(k[,nx])
            }
          }
        )
        tmp <- do.call(
          "cbind.data.frame",
          tmp
        )
        colnames(tmp) <- names(x)

        tmp <- try(
          model.frame(
            formula = my_formula,
            data = tmp,
            na.action = NULL
          ),
          silent = TRUE
        )
        if(inherits(tmp, "try-error")){
          tmp <- gsub(
            "not found\n",
            "not found. This can happen if you did not input temporally varying covariates as a named list or if you left this covariate out of the covariates supplied to the 'rho_covs' argument of auto_occ().\n",
            tmp
          )
          stop(tmp)
        }
        tmp <- model.matrix(
          my_formula,
          tmp
        )
        to_return[[i]][[j]] <- tmp
      }
    }
    return(to_return)
  }
}



#'
#' Make data.frame columns into factors (internal). Mostly taken from unmarked.
#'
#' @param df The data.frame
#'
#' @param name the column name being converted
#'
#' @param temp_var whether or not the data.frame contains temporally varying covariates. Defaults to \code{FALSE}
#'
#' @param type Either 'psi' for occupancy or 'rho' for detection.
#'
#' @noRd

factor_df_cols <- function(df, name=NULL, temp_var = FALSE, type = type){
  stopifnot(
    inherits(
      df, "data.frame"
    )
  )
  char_cols <- sapply(
    df,
    is.character
  )
  dm <- ifelse(
    type == "psi",
    "occ_covs",
    "det_covs"
  )
  if(any(char_cols)){
    if(is.null(name)){
      name <- names(char_cols)[char_cols]
    }

    warning(
      paste0(
        "\n",dm,": '",  name,
        "' column is a character object. Converted it to a factor."
      ),
      call.=FALSE,
      immediate. = TRUE
    )
  }
  if(temp_var){
    # If temporally varying we need unique categories across ALL columns
    df_unq <- unique(
      unlist(
        df
      )
    )
    to_change <- df[,char_cols, drop=FALSE]
    df[,char_cols] <- lapply(
      to_change,
      function(x) factor(x, levels = df_unq)
    )
    return(df)
  }
  to_change <- df[,char_cols, drop=FALSE]
  df[,char_cols] <- lapply(
    to_change,
    function(k){
      my_lvl <- unique(k)
      to_return <- factor(k, levels = my_lvl)
    }
  )
  return(df)
}
