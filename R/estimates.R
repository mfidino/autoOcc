setMethod("vcov", "auto_occ_fit",
          function (object, type , ...)
          {
            if (is.null(object@opt$hessian)) {
              stop("Hessian was not computed for this model.")
            }
            v <- solve(object@opt$hessian)
            pnames <- object@estimates$parameter
            pnames <- substr(pnames,1,3)
            tmp <- rle(pnames)
            tmp <- sapply(tmp$lengths, function(x) 1:x)
            tmp <- unlist(tmp)
            pnames <- paste0(
              substr(
                pnames,1,3
              ),"_",tmp
            )
            rownames(v) <- colnames(v) <- pnames
            if (missing(type)) {
              return (v)
            } else {
              inds <- grep(type,colnames(v))
              return (v[inds, inds, drop = FALSE])
            }
          })

setMethod("linearComb",
          signature(obj = "unmarkedEstimate", coefficients = "matrixOrVector"),
          function(obj, coefficients, offset = NULL, re.form = NULL)
          {
            if(!is(coefficients, "matrix"))
              coefficients <- t(as.matrix(coefficients))
            est <- obj@estimates
            covMat <- obj@covMat
            if(!is.null(re.form) & .hasSlot(obj, "fixed")){
              est <- est[obj@fixed]
              covMat <- covMat[obj@fixed, obj@fixed, drop=FALSE]
            }
            stopifnot(ncol(coefficients) == length(est))
            if (is.null(offset))
              offset <- rep(0, nrow(coefficients))
            e <- as.vector(coefficients %*% est) + offset
            v <- coefficients %*% covMat %*% t(coefficients)
            if (!is.null(obj@covMatBS)) {
              v.bs <- coefficients %*% obj@covMatBS %*% t(coefficients)
            } else {
              v.bs <- NULL
            }
            umelc <- new("unmarkedLinComb", parentEstimate = obj,
                         estimate = e, covMat = v, covMatBS = v.bs,
                         coefficients = coefficients)
            umelc
          })
