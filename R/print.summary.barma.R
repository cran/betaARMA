# source_code.txt

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
      sep = ""
    )
  }

  # ---------------------------------------------------------------------- #
  # 2. Print Link Function
  # ---------------------------------------------------------------------- #
  if (!is.null(x$link)) {
    cat("\nLink function:", x$link, "\n")
  }

  # ---------------------------------------------------------------------- #
  # 2.5 Print Penalization Information
  # ---------------------------------------------------------------------- #
  if (isTRUE(x$penalty)) {
    lambda_str <- if (!is.null(x$lambda)) format(x$lambda, digits = 4) else "Unknown"
    cat(sprintf("Penalization : Ridge (lambda = %s)\n", lambda_str))
  }

  # ---------------------------------------------------------------------- #
  # 3. Print Coefficients Table
  # ---------------------------------------------------------------------- #
  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    stats::printCoefmat(x$coefficients,
      digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"),
      na.print = "NA"
    )
  } else {
    cat("\nNo coefficients table found.\n")
  }

  # ---------------------------------------------------------------------- #
  # 4. Print Optimization and Convergence Status
  # ---------------------------------------------------------------------- #
  if (!is.null(x$conv)) {
    # 1. Extract raw method
    raw_method <- if (!is.null(x$opt_method)) x$opt_method else "Unknown"

    # 2. Format the string with package information
    method_str <- switch(raw_method,
      "BFGS"     = "BFGS (stats::optim)",
      "L-BFGS-B" = "L-BFGS-B (stats::optim)",
      "lbfgs"    = "L-BFGS (lbfgs::lbfgs)",
      raw_method
    ) # Fallback to raw name if no match

    # 3. Append bounds information if applicable
    if (raw_method == "L-BFGS-B" && isTRUE(x$opt_bounds)) {
      method_str <- paste0(method_str, " with box constraints")
    }

    # 4. Print the final diagnostic message
    if (x$conv == 0) {
      if (is.null(x$opt$counts)) {
        cat(sprintf(
          "\nOptimization: %s converged.\n",
          method_str
        ))
      } else {
        cat(sprintf(
          "\nOptimization: %s converged after %s objective function evaluations.\n",
          method_str, x$opt$counts[1]
        ))
      }
    } else {
      cat(sprintf(
        "\nWarning: Optimization algorithm [%s] did not converge (Code: %s).\n",
        method_str, x$conv
      ))
    }
  }

  # ---------------------------------------------------------------------- #
  # 5. Print Log-Likelihood and Information Criteria
  # ---------------------------------------------------------------------- #
  if (!is.null(x$loglik)) {
    cat("---\n")
    cat("Log-likelihood:", formatC(x$loglik, digits = 4), "\n")
    cat(
      "AIC:", formatC(x$aic, digits = 4),
      "  BIC:", formatC(x$bic, digits = 4),
      "  HQ:", formatC(x$hq, digits = 4), "\n"
    )
  }

  invisible(x)
}
