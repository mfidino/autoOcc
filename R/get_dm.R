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
      to_return[[i]] <- data.frame(
        dplyr::bind_cols(
          to_return[[i]]
        )
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
    to_return <- vector("list", length = nseason)
    for(i in 1:nseason){
      for(j in 1:nrep){
        if(j == 1){
          to_return[[i]] <- vector("list", length = nrep)
        }
        tmp <-  lapply(
          x,
          function(k){
            if(ncol(k) == nseason){
              return(k[,i])
            } else{
              nx <- ((1:ncol(k) - 0.001) %/% nrep) + 1
              nx <- which(nx == i)[j]
              return(k[,nx])
            }
          }
        )
        tmp <- data.frame(
          dplyr::bind_cols(
            tmp
          )
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


