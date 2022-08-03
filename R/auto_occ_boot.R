setMethod("parboot", "auto_occ_fit",
          function(object, statistic=auto_occ_sse, nsim=10, report, seed = NULL, parallel = TRUE, ncores, ...)
          {
            dots <- list(...)
            statistic <- match.fun(statistic)
            call <- match.call(call = sys.call(-1))
            formula <- object@formula
            umf <- object@occcovs
            umf <- as.data.frame(umf)
            y <- object@y
            est <- object@estimates
            est <- as.numeric(est$Est[grep("psi",est$parameter)])

            starts <- est
            t0 <- statistic(object)
            lt0 <- length(t0)
            t.star <- matrix(NA, nsim, lt0)
            if(!missing(report))
              cat("t0 =", t0, "\n")
            simdata <- umf
            if (!is.null(seed)) set.seed(seed)
            simList <- simulate(object, nsim = nsim, na.rm = FALSE)
            availcores <- detectCores()
            if(missing(ncores)) ncores <- availcores - 1
            if(ncores > availcores) ncores <- availcores

            no_par <- ncores < 2 || nsim < 100 || !parallel

            if (no_par) {
              if (!missing(report)) {
                for(i in 1:nsim) {
                  simdata <- replaceY(simdata, simList[[i]])
                  fit <- update(object, data=simdata, starts=starts, se=FALSE)
                  t.star[i,] <- statistic(fit, ...)
                  if(!missing(report)) {
                    if (nsim > report && i %in% seq(report, nsim, by=report))
                      cat("iter", i, ": ", t.star[i, ], "\n")
                  }
                }
              } else {
                t.star <- pbsapply(1:nsim, function(i) {
                  simdata <- replaceY(simdata, simList[[i]])
                  fit <- update(object, data=simdata, starts=starts, se=FALSE)
                  t.star.tmp <- statistic(fit, ...)
                })
                if (lt0 > 1)
                  t.star <- t(t.star)
                else
                  t.star <- matrix(t.star, ncol = lt0)
              }
            } else {
              message("Running parametric bootstrap in parallel on ", ncores, " cores.")
              if (!missing(report)) message("Bootstrapped statistics not reported during parallel processing.")
              cl <- makeCluster(ncores)
              if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
              on.exit(stopCluster(cl))
              varList <- c("simList", "y", "object", "simdata", "starts", "statistic", "dots")
              # If call formula is an object, include it too
              fm.nms <- all.names(object@call)
              if (!any(grepl("~", fm.nms))) varList <- c(varList, fm.nms[2])
              ## Hack to get piFun for unmarkedFitGMM and unmarkedFitMPois
              if(.hasSlot(umf, "piFun")) varList <- c(varList, umf@piFun)
              clusterExport(cl, varList, envir = environment())
              clusterEvalQ(cl, library(unmarked))
              clusterEvalQ(cl, list2env(dots))
              t.star.parallel <- pblapply(1:nsim, function(i) {
                simdata <- replaceY(simdata, simList[[i]])
                fit <- update(object, data = simdata, starts = starts, se = FALSE)
                t.star <- statistic(fit, ...)
              }, cl = cl)
              t.star <- matrix(unlist(t.star.parallel), nrow = length(t.star.parallel), byrow = TRUE)
            }
            if (!is.null(names(t0)))
              colnames(t.star) <- names(t0)
            else colnames(t.star) <- paste("t*", 1:lt0, sep="")
            out <- new("parboot", call = call, t0 = t0, t.star = t.star)
            return(out)
          })


auto_occ_sse <- function(fit,...){
  sse <- sum(residuals(fit)^2, na.rm=TRUE)
  return(c(SSE=sse))
}

auto_occ_fitted <- function(object, na.rm = FALSE){
    char <- lapply(object@formula, function(x) {
      paste(deparse(x), collapse = "")
    })
    rho_formula <- as.formula(char[[2]])
    psi_formula <- as.formula(paste("~", char[[3]]))

    data <- object@occcovs
    if(all(sapply(data,function(x) is.null(dim(x))))){
      data <- as.data.frame(data)
    }

    # get the design matrix for psi
    mf <- get_dm(x = data, my_formula  = psi_formula, type = "psi",y = object@y)

    # any offset if it is present
    x.offset <- lapply(
      mf,
      function(x){
        tmp <- grep("(offset)", colnames(x))
        if(length(tmp) == 0){
          tmp <- rep(0, nrow(x))
        }else{
          tmp <- x[,tmp]
        }
        return(tmp)
      }
    )
    est <- object@estimates
    est <- est$Est[grep("psi",est$parameter)]
    state <- vector("list", length= length(mf))
    for(i in 1:length(mf)){
      e1 <- plogis(cbind(mf[[i]],0) %*% est + x.offset[[i]])
      e2 <- plogis(cbind(mf[[i]],1) %*% est + x.offset[[i]])
      state[[i]] <- e1 / (e1 + (1 - e2))
    }
    # get probability of detection as well.
    data <- object@detcovs
    mf <- get_dm(x = data, my_formula  = rho_formula, type = "rho",y = object@y)

    # any offset if it is present
    x.offset <- lapply(
      mf,
      function(x){
        lapply(
          x,
          function(y){
            tmp <- grep("(offset)", colnames(y))
            if(length(tmp) == 0){
              tmp <- rep(0, nrow(y))
            }else{
              tmp <- y[,tmp]
            }
            return(tmp)
          }
        )
      }
    )
    est <- object@estimates
    est <- est$Est[grep("rho",est$parameter)]
    p <- vector("list", length= length(mf))
    for(i in 1:length(mf)){
      p[[i]] <- vector("list", length(mf[[i]]))
      for(j in 1:length(mf[[i]])){
        p[[i]][[j]] <- plogis(mf[[i]][[j]] %*% est + x.offset[[i]][[j]])
      }
    }
    # true for models with E[Y] = p * E[X]
    fitted <- object@y
    for(i in 1:length(mf)){
      for(j in 1:length(mf[[i]])){
        fitted[,i,j] <- state[[i]] * p[[i]][[j]]
      }
    }
    return(fitted)
}

longshot <- auto_occ_fitted(m2)

auto_occ_residuals <- function(object){
  y <- object@y
  e <- auto_occ_fitted(object)
  r <- y - e
  return(r)
}

auto_occ_sse <- function(object){
  sse <- sum(auto_occ_residuals(object)^2, na.rm = TRUE)
  return(sse)
}

ack <- auto_occ_sse(m2)

yo <- auto_occ_residuals(m2)




setMethod("nonparboot", "unmarkedFit",



        auto_occ_boot <-   function(object, nsim = 10, seed = NULL) {
          y <- object@y

          if(!is.null(seed)){
            set.seed(seed)
          }
          nsite <- dim(y)[1]

          parm_matrix <- matrix(
            NA,
            ncol = nrow(object@estimates),
            nrow = nsim
          )
          colnames(parm_matrix) <- object@estimates$parameter
          cat(
            paste0("\nBootstrapping ", nsim, " times...\n")
          )
          pb <- txtProgressBar(max = nsim)
          for(i in 1:nsim){
            new_y <- y
            my_resample <- sample(
                1:nsite,
                nsite,
                replace = TRUE
            )
            new_y <- y[my_resample,,]
            new_occ_covs <- object@occcovs[my_resample,]
            new_det_covs <- lapply(
              object@detcovs,
              function(x, mr = my_resample)
                x[mr,]
            )


            new_fit <- auto_occ(
              object@formula, y = new_y,det_covs = new_det_covs,
              occ_covs = new_occ_covs
            )
            if(new_fit@opt$convergence == 0){
              parm_matrix[i,] <- new_fit@estimates$Est
            }
            setTxtProgressBar(pb, i)
          }
          parm_matrix <- parm_matrix[complete.cases(parm_matrix),]
          return(parm_matrix)
          })


hm2 <- auto_occ_boot(m2, nsim = 500)

all_hm <- rbind(hm, hm2)

ack <- rbind(pmcmc, hm2[,1:5])
# do some predictions
