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
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#'
#' @noRd


get_dm <- function(x, my_formula, type = c("psi","rho"), y){
  # if just a site by nparam matrix for psi
  nsite <- dim(y)[1]
  nseason <- dim(y)[2]
  nrep <- dim(y)[3]
  # if psi does not temporally vary
  if(
    is.data.frame(x) &
    type == "psi"
  ){
    to_return <- vector(
      "list",
      length = nseason
    )
    tmp <- model.frame(
      formula = my_formula,
      data = x,
      na.action = NULL
    )
    tmp <- model.matrix(
      my_formula,
      tmp
    )
    for(i in 1:nseason){
      to_return[[i]] <- tmp
    }
    return(to_return)
    temp_var_psi <- FALSE
  }
  # if temporally varying psi
  if(
    is.list(x) &
    !is.data.frame(x) &
    type == "psi"
  ){
    to_return <- vector("list", length = nseason)
    for(i in 1:nseason){
      to_return[[i]] <- lapply(
        x,
        function(k){
          if(ncol(k) == 1){
            return(k[,1])
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
      return(to_return)
    }
    temp_var_psi <- TRUE
  }
  if(
    type == "rho"
  ){
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
        "a matrix with a number of columns equal to the\n",
        "number of primary sampling periods, or a matrix with\n",
        "a number of columns equal to the number of primary\n",
        "sampling periods times the number of secondary observations\n",
        "within time periods. These covariates are not properly\n",
        "set up: ", paste0(names(x)[baddies] , collapse = ", ")
      )
      stop(error_report)
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
            if(is.null(ncol(k))){
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

        tmp <- model.frame(
          formula = my_formula,
          data = tmp,
          na.action = NULL
        )
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
