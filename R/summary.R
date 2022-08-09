#' @title An S4 class to represent the summary of a fitted autologistic occupancy model
#
#' @slot call The formula call.
#' @slot optim_convergence_code Convergence code from \code{\link[stats]{optim}}
#' @slot optim_iterations Number of iterations from \code{\link[stats]{optim}}
#' @slot psi The occupancy parameters in a data.frame
#' @slot rho The detection parameters in a data.frame
#' @slot AIC the AIC of the fitted model.
#'
#' @name auto_occ_summary-class
#' @rdname auto_occ_summary-class
#' @export
#'

#' @importFrom methods new
setClass("auto_occ_summary",
         representation(call = "call",
                        optim_convergence_code = "integer",
                        optim_iterations = "integer",
                        psi = "data.frame",
                        rho = "data.frame",
                        AIC = "numeric"))

  # constructor for summary.auto_occ_fit
  auto_occ_summary <- function(call, optim_convergence_code,
                                   optim_iterations,
                                   psi,rho,AIC){
    my_summary <- new(
      "auto_occ_summary",
      call = call,
      optim_convergence_code = optim_convergence_code,
      optim_iterations = optim_iterations,
      psi = psi,
      rho = rho,
      AIC = AIC
    )
    return(my_summary)
  }

#' @importFrom methods show
#' @noRd
  setMethod("show", "auto_occ_summary", function(object)
  {
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    cat("\noptim convergence code:", object@optim_convergence_code)
    cat("\noptim iterations:", object@optim_iterations, "\n")
    if(!identical(object@optim_convergence_code, 0L))
      warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
    cat("\nOccupancy estimates:\n\n")
    print(object@psi, digits = 3)
    cat("\nNote: psi - theta is the autologistic term\n")
    cat("\nDetection estimates:\n\n")
    print(object@rho, digits = 3)
    cat("\n")
    cat("AIC:", object@AIC,"\n")

  })

