#' Extract Coefficients from a barma Model
#'
#' @description
#' S3 method for extracting the vector of estimated coefficients from a fitted 
#' model object of class `"barma"`.
#'
#' @param object A fitted model object of class `"barma"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A named numeric vector of all estimated coefficients 
#' (e.g., alpha, varphi, theta, phi).
#'
#' @importFrom stats coef
#' @export
#' @method coef barma
coef.barma <- function(object, ...) {
  return(object$coef)
}