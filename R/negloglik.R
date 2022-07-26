#' Likelihood function for auto_occ (internal).
#'
#' @param parms the parameters to estimate
#'
#' @param y A three-dimensional array of species detections. The first dimension
#' is sites, the second dimension denotes primary sampling periods, and the third
#' dimension denotes the secondary sampling periods within each primary sampling
#' period. If the species was detected on a given survey, that element would receive
#' a 1, otherwise it is 0. If sampling did not occur for a given sampling period, those
#' elements should be NA.
#'
#' @param cov_list A list of covariates supplied by \code{auto_occ}.
#' @importFrom stats plogis
#'
#' @noRd
negloglik <- function(parms,y,cov_list){

  # get names of psi parameters
  psi_parm_names <- colnames(cov_list$psi[[1]])

  # get names of rho parameters
  rho_parm_names <- colnames(cov_list$rho[[1]][[1]])

  psi_parms <- parms[1:length(psi_parm_names)]
  names(psi_parms) <- psi_parm_names
  theta <- parms[length(psi_parm_names)+1]
  rho_parms <- parms[(length(psi_parm_names)+2):length(parms)]
  names(rho_parms) <- rho_parm_names

  # some indices
  nsite <- dim(y)[1]
  nseason <- dim(y)[2]
  nrep <- dim(y)[3]

  # get probabilities
  psi <- psi_theta <- matrix(
    NA,
    ncol = nseason,
    nrow = nsite
  )
  for(i in 1:nseason){
      psi[,i] <- cov_list$psi[[i]] %*% psi_parms
  }
  psi_theta <- psi + theta
  psi_prob <- plogis(psi)
  psi_theta_prob <- plogis(psi_theta)

  rho <- array(
    NA,
    dim = c(nsite, nseason, nrep)
  )

  for(i in 1:nseason){
    for(j in 1:nrep){
      rho[,i,j] <- cov_list$rho[[i]][[j]] %*% rho_parms
    }
  }
  rho_prob <- plogis(rho)
  # get initial occupancy
  #init_psi <- psi[,1] / (psi[,1] + (1 - psi_theta[,1]))
  #psi_mat <- cbind(
  #  init_psi,
  #  1 - init_psi
  #)
  psi_mat <- cbind(
    psi_prob[,1],
    1 - psi_prob[,1]
  )
  tpm <- array(
    c(
      psi_theta_prob,psi_prob, 1 - psi_theta_prob, 1 - psi_prob
    ),
    dim = c(nsite, nseason, 2,2)
  )

  rho_mat <- function(my_y,my_rho){
    j <- sum(!is.na(my_y))
    ndet <- sum(my_y, na.rm = TRUE)
    my_y <- my_y[!is.na(my_y)]
    to_return <-
      diag(
        c(
          prod(
            ifelse(my_y>0, my_rho, 1-my_rho)
          ),
          ifelse(ndet>0, 0, 1)
        )
      )
    return(to_return)
  }

  likelihood <- rep(NA, nsite)
  for(i in 1:nsite){
    for(t in 1:nseason){
      if(t == 1){
        tmp_mat <- psi_mat[i,] %*% rho_mat(y[i,t,], my_rho = rho_prob[i,t,])
        next
      }
      if(t != nseason){
        tmp_mat <- tmp_mat %*% tpm[i,t,,] %*% rho_mat(y[i,t,], my_rho = rho_prob[i,t,])
        next
      }
      if(t == nseason){
        tmp_mat <- tmp_mat %*% tpm[i,t,,] %*% matrix(
          diag(
            rho_mat(y[i,t,], my_rho = rho_prob[i,t,])
          ), ncol = 1
        )
      }
    }
    likelihood[i] <- tmp_mat
  }
  nll <- -sum(log(likelihood))
  return(nll)
}
