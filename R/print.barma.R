#' Print Method for a barma Model
#'
#' @description
#' S3 method for printing a summary of a fitted `"barma"` model object.
#'
#' @details
#' This function provides a concise summary, showing the function call 
#' that generated the model, the link function used, and the final 
#' estimated coefficients.
#'
#' @param x A fitted model object of class `"barma"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the original object `x`.
#'
#' @export
#' @method print barma
print.barma <- function(x, ...) {
  
  # Print the title
  cat("Beta Autoregressive Moving Average Model\n")
  
  # Print the call, handling line breaks
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", 
        sep = "")
  }
  
  # Print the link function
  if (!is.null(x$link)) {
    cat("\nLink function:", x$link, "\n")
  }
  
  # Print the coefficients
  if (!is.null(x$coef)) {
    cat("\nCoefficients:\n")
    print(x$coef)
  } else {
    cat("\nNo coefficients found in object.\n")
  }
  
  # Add a warning if convergence was not reached
  if (!is.null(x$conv) && x$conv != 0) {
    cat("\nWarning: Optimization algorithm did not converge.\n")
  }
  
  # Return the object invisibly, which is standard for print methods
  invisible(x)
}