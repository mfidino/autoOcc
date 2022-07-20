setClass("auto_occ_fit",
         representation(fitType = "character",
                        call = "call",
                        formula = "formula",
                        y = "array",
                        estimates = "data.frame",
                        AIC = "numeric",
                        opt = "list",
                        negLogLike = "numeric",
                        nllFun = "function",
                        detcovs = "list",
                        occcovs = "list"))

# constructor for auto_occ_fit objects
auto_occ_fit <- function(fitType, call, formula, y,
                        estimates, AIC, opt, negLogLike, nllFun,
                        detcovs, occcovs)
{
  aofit <- new("auto_occ_fit", fitType = fitType, call = call,
               formula = formula, y = y,
               estimates = estimates, AIC = AIC, opt = opt,
               negLogLike = negLogLike,
               nllFun = nllFun,
               detcovs = detcovs, occcovs = occcovs)
  return(aofit)
}



setMethod("summary", "auto_occ_fit", function(object)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  summaryOut <- summary(object@estimates)
  cat("AIC:", object@AIC,"\n")
  cat("\noptim convergence code:", object@opt$convergence)
  cat("\noptim iterations:", object@opt$counts[1], "\n")
  if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
  invisible(summaryOut)
  cat("\nParameter estimates:\n\n")
  print(object@estimates)
})


setMethod("show", "auto_occ_fit", function(object)
{
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  show(object@estimates)
  cat("\nAIC:", object@AIC,"\n")
  if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
})
