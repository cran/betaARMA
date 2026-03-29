#' Print Method for a barma Summary
#'
#' @description
#' S3 method for printing the detailed summary of a `"barma"` model object.
#'
#' @details
#' This function formats and displays the summary list created by
#' `summary.barma()`, including the call, coefficient table, and
#' information criteria.
#'
#' @param x A fitted model summary object of class `"summary.barma"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the original object `x`.
#'
#' @export
#' @method print summary.barma
print.summary.barma <- function(x, ...) {
  
  # ---------------------------------------------------------------------- #
  # 1. Print Model and Call 
  # ---------------------------------------------------------------------- #
  cat("Beta Autoregressive Moving Average Model\n")
  
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", 
        sep = "")
  }
  
  # ---------------------------------------------------------------------- #
  # 2. Print Link Function 
  # ---------------------------------------------------------------------- #
  if (!is.null(x$link)) {
    cat("\nLink function:", x$link, "\n")
  }
  
  # ---------------------------------------------------------------------- #
  # 3. Print Coefficients Table 
  # ---------------------------------------------------------------------- #
  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    print(x$coefficients)
  } else {
    cat("\nNo coefficients table found.\n")
  }
  
  # ---------------------------------------------------------------------- #
  # 4. Print Convergence Status 
  # ---------------------------------------------------------------------- #
  if (!is.null(x$conv) && x$conv != 0) {
    cat("\nWarning: Optimization algorithm did not converge.\n")
  }
  
  # ---------------------------------------------------------------------- #
  # 5. Print Log-Likelihood and Information Criteria 
  # ---------------------------------------------------------------------- #
  if (!is.null(x$loglik)) {
    cat("\n---")
    cat("\nLog-likelihood:", round(x$loglik, 4))
    cat("\nAIC:", round(x$aic, 4), 
        "| BIC:", round(x$bic, 4), 
        "| HQ:", round(x$hq, 4), "\n")
  }
  
  invisible(x)
}