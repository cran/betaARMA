#' Extract Fitted Values from a barma Model
#'
#' @description
#' S3 method for extracting the fitted mean values (mu-hat) from a fitted 
#' model object of class `"barma"`.
#'
#' @details
#' The fitted values are returned as a time series (`ts`) object, 
#' matching the time properties of the original input `y`. The first
#' `max_lag` observations are `NA`, as they cannot be fitted.
#'
#' @param object A fitted model object of class `"barma"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A `ts` object of the fitted mean values.
#'
#' @importFrom stats fitted
#' @export
#' @method fitted barma
fitted.barma <- function(object, ...) {
  return(object$fitted)
}