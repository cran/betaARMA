#' @title Create Link Function Structure for BARMA Models
#'
#' @description
#' A helper function that constructs a list containing the link function,
#' its inverse, and the derivative of the mean function. It extends the
#' standard \code{make.link} by adding support for the \code{"loglog"} link.
#'
#' @details
#' This function is used by \code{\link{barma}}, \code{\link{loglik_barma}},
#' \code{\link{score_vector_barma}}, and \code{\link{fim_barma}} to handle
#' the \code{link} argument in a standardized way. It is also exported for
#' users who implement custom workflows.
#'
#' For the \code{"logit"}, \code{"probit"}, and \code{"cloglog"} links,
#' the function acts as a wrapper around the base R \code{make.link}.
#'
#' For the \code{"loglog"} link, which is not available in \code{make.link},
#' the necessary components are defined explicitly:
#' \itemize{
#'   \item Link function: \eqn{g(\mu) = -\log(-\log(\mu))}
#'   \item Inverse link: \eqn{g^{-1}(\eta) = \exp(-\exp(-\eta))}
#'   \item Derivative \eqn{\frac{d\mu}{d\eta}}:
#'     \eqn{\exp(-\exp(-\eta) - \eta)}
#' }
#'
#' If an unsupported link is provided, the function stops with an error
#' message listing the available options.
#'
#' @param link A character string specifying the link function.
#'   Defaults to \code{"logit"}. Accepted values are \code{"logit"},
#'   \code{"probit"}, \code{"cloglog"}, and \code{"loglog"}.
#'
#' @importFrom stats make.link
#'
#' @return A list with three components:
#' \item{linkfun}{The link function \eqn{g(\mu)}, transforming
#'   \eqn{\mu \in (0,1)} to \eqn{\eta \in (-\infty, \infty)}.}
#' \item{linkinv}{The inverse link function \eqn{g^{-1}(\eta)},
#'   transforming \eqn{\eta} back to \eqn{\mu}.}
#' \item{mu.eta}{The derivative \eqn{d\mu/d\eta} of the inverse
#'   link function.}
#'
#' @author
#' Original R code by Fabio M. Bayer
#' (Federal University of Santa Maria, \email{bayer@ufsm.br}).
#' Modified and improved by Everton da Costa
#' (Federal University of Pernambuco, \email{everto.cost@gmail.com}).
#'
#' @seealso \code{\link{barma}}, \code{\link[stats]{make.link}}
#'
#' @examples
#' # --- Create a logit link structure ---
#' logit_link <- make_link_structure(link = "logit")
#'
#' # Apply the link function
#' mu <- 0.5
#' eta <- logit_link$linkfun(mu)
#' print(eta) # Should be 0
#'
#' # Apply the inverse link function
#' mu_restored <- logit_link$linkinv(eta)
#' print(mu_restored) # Should be 0.5
#'
#' # --- Create a loglog link structure ---
#' loglog_link <- make_link_structure(link = "loglog")
#'
#' # Apply the loglog link function
#' mu_loglog <- 0.8
#' eta_loglog <- loglog_link$linkfun(mu_loglog)
#' print(eta_loglog)
#'
#' # Apply the inverse loglog link function
#' mu_restored_loglog <- loglog_link$linkinv(eta_loglog)
#' print(mu_restored_loglog) # Should be ~0.8
#'
#' @export
make_link_structure <- function(link = "logit") {
  
  # ======================================================================== #
  # 1. Validate Link Argument
  # ======================================================================== #
  
  if (!is.character(link) || length(link) != 1) {
    stop("Argument 'link' must be a single character string.")
  }
  
  # ======================================================================== #
  # 2. Select Link Structure
  # ======================================================================== #
  
  if (link == "logit") {
    stats <- make.link("logit")
  } else if (link == "probit") {
    stats <- make.link("probit")
  } else if (link == "cloglog") {
    stats <- make.link("cloglog")
  } else if (link == "loglog") {
    stats <- list(
      linkfun = function(mu)  -log(-log(mu)),
      linkinv = function(eta) exp(-exp(-eta)),
      mu.eta  = function(eta) exp(-exp(-eta) - eta)
    )
  } else {
    stop(
      "'", link, "' link not available. ",
      "Available links are: 'logit', 'probit', 'cloglog', 'loglog'."
    )
  }
  
  # ======================================================================== #
  # 3. Return Structure
  # ======================================================================== #
  
  return(list(
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta  = stats$mu.eta
  ))
  
}