newdata <- data.frame(
  x = seq(-3,3,length.out = 50),
  j = factor("B", levels = c("A","B")))

# things to do, get the dm stuff to work.

setMethod(
  "predict",
  "auto_occ_fit",
  function(model,type, newdata = NULL,backTransform = TRUE,level = 0.95){
    if(missing(type)){
      stop("Supply either 'psi' or 'rho' as type for predictions.")
    }
    if(!type %in% c("psi","rho")){
      stop("type must be either' psi' for occupancy or 'rho' for detection." )
    }
    char <- lapply(model@formula, function(x) {
      paste(deparse(x), collapse = "")
    })
    if(type == "psi"){
      my_formula <- as.formula(paste("~", char[[3]]))
    }else{
      my_formula <- as.formula(char[[2]])
    }
    if(
      is.null(newdata)
    ){
     if(type == "psi"){
       newdata <- model@occcovs
     }else{
       newdata <- model@detcovs

     }
    }

    # get the old dm
    data <- ifelse(
      type == "psi",
      model@occcovs,
      model@deccovs)
    mf <- model.frame(my_formula, model@, na.action=stats::na.pass)
    X.terms <- stats::terms(data[[1]])
    fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
    xlevs <- lapply(fac_cols, levels)
    xlevs <- xlevs[names(xlevs) %in% names(mf)]
    nmf <- model.frame(my_formula, newdata, na.action=stats::na.pass, xlev = c("A","B"))
    #X <- model.matrix(X.terms, newdata, xlev=xlevs)
    X <- model.matrix(form_nobars, nmf)
    offset <- model.offset(nmf)

    # get variance covariance matrix
    covMat <- vcov(
      model,
      type = type
    )

    est <- model@estimates
    est <- est$Est[grep(type,est$parameter)]
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

  }
)
