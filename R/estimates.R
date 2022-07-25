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

setMethod("SE", "auto_occ_fit", function(obj,...)
{
  v <- vcov(obj,...)
  sqrt(diag(v))
})

setMethod("confint","auto_occ_fit",
          function(object,parm, level = 0.95,type){

  if(missing(type)){
    stop("Must specify type as either 'psi' or 'rho'")
  }
  if(missing(parm)){
    parm <- grep(type, object@estimates$parameter)
  }else{
    to_grep <- paste0("psi - ", parm, collapse = "|")
    parm <- which(object@estimates$parameter %in% to_grep)
    if(length(parm)<1){
      stop(
        paste0(
          "One or more of parm supplied are not in model.",
        )
      )
    }
  }

  ests <- object@estimates[parm,]
  lwr <- qnorm(
    (1 - level)/2,
    ests$Est,
    ests$SE
  )
  upr <- qnorm(
    (1 - (1 - level)/2),
    ests$Est,
    ests$SE
  )
  return(
    data.frame(
      parameter = ests$parameter,
      Est = ests$Est,
      lower = lwr,
      upper = upr
    )
  )

})

