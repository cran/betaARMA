#' Simulate a Beta Autoregressive Moving Average (BARMA) Time Series
#'
#' @description
#' Generates a random time series from a Beta Autoregressive Moving Average
#' (BARMA) model. The function can simulate BARMA(p,q), BAR(p), or BMA(q)
#' processes.
#'
#' @details
#' The model type is determined by the `varphi` and `theta` parameters:
#' \itemize{
#'   \item **BARMA(p,q):** Both `varphi` and `theta` are provided.
#'   \item **BAR(p):** Only `varphi` is provided (`theta` is `NA`).
#'   \item **BMA(q):** Only `theta` is provided (`varphi` is `NA`).
#' }
#'
#' @param n The desired length of the final time series (after burn-in).
#' @param varphi A numeric vector of autoregressive (AR) parameters (\eqn{\varphi}).
#'   Default is `NA` for models without an AR component.
#' @param theta A numeric vector of moving average (MA) parameters (\eqn{\theta}).
#'   Default is `NA` for models without an MA component.
#' @param alpha The intercept term (\eqn{\alpha}) in the linear predictor.
#'   Default is `0.0`.
#' @param phi The precision parameter (\eqn{\phi > 0}) of the Beta distribution.
#'   Higher values result in less variance. Default is `20`.
#' @param freq The frequency of the resulting `ts` object (e.g., 12 for monthly,
#'   4 for quarterly). Default is `12`.
#' @param link The link function to connect the mean \eqn{\mu_t} to the linear
#'   predictor. Must be one of `"logit"` (default), `"probit"`, `"cloglog"`, or
#'   `"loglog"`.
#'
#' @importFrom stats rbeta
#'
#' @return A `ts` object of length `n` containing the simulated Beta-distributed
#'   time series.
#'
#' @author
#' Original R code by: Fabio M. Bayer (Federal University of Santa Maria, <bayer@ufsm.br>)
#' Modified and improved by: Everton da Costa (Federal University of Pernambuco, <everto.cost@gmail.com>)
#'
#' @seealso \code{\link{barma}}, \code{\link{make_link_structure}}
#'
#' @examples
#' # Set seed for reproducibility
#' # --- Example 1: Simulate a BAR(1) process ---
#' # y_t depends on y_{t-1}
#' set.seed(2025)
#' bar1_series <- simu_barma(n = 250, alpha = 0.0, varphi = 0.5, phi = 20,
#' link = "logit")
#' plot(bar1_series, main = "Simulated BAR(1) Process", ylab = "Value")
#'
#' # --- Example 2: Simulate a BMA(1) process ---
#' # y_t depends on the previous error term
#' set.seed(2025)
#' bma1_series <- simu_barma(n = 250, alpha = 0.0, theta = -0.2, phi = 20)
#' plot(bma1_series, main = "Simulated BMA(1) Process", ylab = "Value")
#'
#' # --- Example 3: Simulate a BARMA(2,1) process with a cloglog link ---
#' set.seed(2025)
#' barma21_series <- simu_barma(
#'   n = 200,
#'   alpha = 0.0,
#'   varphi = c(0.4, 0.2), # AR(2) components
#'   theta = -0.3,         # MA(1) component
#'   phi = 20,
#'   link = "cloglog"
#' )
#' plot(barma21_series, main = "Simulated BARMA(2,1) Process", ylab = "Value")
#'
#' @export
simu_barma <- function(n,
                       alpha = 0.0,
                       varphi = NA, theta = NA,
                       phi = 20,
                       freq = 12, link = "logit") {


  ar <- NA
  ma <- NA

  if (any(is.na(varphi) == F)) ar <- 1:length(varphi)
  if (any(is.na(theta) == F)) ma <- 1:length(theta)

  # Link functions
  # The link function structure provides the necessary components
  # to perform the \beta ARMA model simulation with the specified link function.
  # Available link functions: "logit", "loglog", "cloglog"
  link_structure <- make_link_structure(link)

  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta


  # ===========================================================================
  # BARMA model
  # ===========================================================================
  if (any(is.na(varphi) == F) && any(is.na(theta) == F)) {
    # This section simulates a \beta ARMA(p ,q) process
    # where both AR and MA orders are specified

    # Compute the AR and MA orders
    p <- max(ar)
    q <- max(ma)
    m <- max(p, q)

    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)

    # E(error) = 0
    error <- rep(0, n + m)
    eta <- y <- rep(NA, n + m)

    for (i in (m + 1):(n + m)) {

      eta[i] <- alpha + varphi %*% ynew[i - ar] + theta %*% error[i - ma]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])
      error[i] <- ynew[i] - eta[i]

    }

    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }


  # ===========================================================================
  # BAR model
  # ===========================================================================
  if (any(is.na(varphi) == F) && any(is.na(theta) == T)) {
    # This section simulates a \beta AR(p) process
    # where only the AR order is specified

    # Compute the AR order
    p <- max(ar)
    m <- p

    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)

    eta <- y <- rep(NA, n + m)

    for (i in (m + 1):(n + m)) {

      eta[i] <- alpha + varphi %*% ynew[i - ar]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])

    }

    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }


  # ===========================================================================
  # BMA model
  # ===========================================================================
  if (any(is.na(varphi) == T) && any(is.na(theta) == F)) {
    # This section simulates a \beta MA(q) process
    # where only the MA order is specified

    # Compute the MA order
    q <- max(ma)
    m <- q

    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)

    # E(error)=0
    eta <- y <- error <- rep(0, n + m)

    for (i in (m + 1):(n + m)) {

      eta[i] <- alpha + theta %*% error[i - ma]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])
      error[i] <- ynew[i] - eta[i]

    }

    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }


}



