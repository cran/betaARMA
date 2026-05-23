#' Fit Beta Autoregressive Moving Average (BARMA) Models
#'
#' @title Fit Beta Autoregressive Moving Average (BARMA) Models via
#' Conditional Maximum Likelihood
#'
#' @description
#' Fits a Beta Autoregressive Moving Average (BARMA) model to time series data
#' valued in (0, 1) using Conditional Maximum Likelihood Estimation (CMLE).
#' The function performs complete model estimation including parameter
#' estimation, hypothesis testing infrastructure, and model diagnostics.
#'
#' @details
#' This function fits the BARMA(p,q) model as proposed by Rocha & Cribari-Neto
#' (2009, with erratum 2017). It serves as the main wrapper for the optimization
#' process, calling specialized helper functions for likelihood computation,
#' gradient calculation, and Fisher Information Matrix estimation.
#'
#' **Model Specification**: The BARMA model is defined as:
#'
#' \deqn{g(\mu_t) = \alpha + X_t\beta + \sum_{i=1}^{p} \varphi_i
#' (g(y_{t-i}) - X_{t-i}\beta) + \sum_{j=1}^{q} \theta_j \epsilon_{t-j}}
#'
#' where \eqn{y_t | F_{t-1} \sim Beta(\mu_t\phi, (1-\mu_t)\phi)},
#' \eqn{g(\cdot)} is the link function, and \eqn{F_{t-1}} is the information
#' set at time t-1.
#'
#' **Model Types** (specified via `ar` and `ma` arguments):
#' \itemize{
#'   \item **BARMA(p,q):** Both `ar` and `ma` are specified.
#'   \item **BAR(p):** Only `ar` is specified; omit `ma` or pass `integer(0)`.
#'   \item **BMA(q):** Only `ma` is specified; omit `ar` or pass `integer(0)`.
#' }
#'
#' **External Regressors**: Covariates can be included via `xreg`. The model
#' becomes a regression BARMA, where the mean depends on current covariates
#' and lagged responses/errors.
#'
#' **Optimization**: The function uses quasi-Newton algorithms
#' (BFGS or L-BFGS-B) via \code{\link[stats]{optim}}, or the
#' \code{\link[lbfgs]{lbfgs}} function, utilizing analytic gradients for
#' efficiency. Initial values are obtained from the internal `start_values()`
#' function.
#'
#' **Implementation Notes**:
#' - The conditional log-likelihood is computed by conditioning on the first
#'   \eqn{m = \max(p,q)} observations, which are required to initialize the
#'   recursive structure of the model.
#' - The 2017 Erratum corrections are implemented for correct handling of
#'   moving average components in the score vector and Fisher Information Matrix.
#' - All computations are vectorized where possible for efficiency.
#'
#' @author
#' Everton da Costa (Federal University of Pernambuco,
#' \email{everto.cost@gmail.com});
#' Francisco Cribari-Neto (Federal University of Pernambuco,
#' \email{francisco.cribari@ufpe.br})
#'
#' @note
#' The original version of this function was developed by Fabio M. Bayer
#' (Federal University of Santa Maria, \email{bayer@ufsm.br}). It has been
#' substantially modified and improved by Everton da Costa, with suggestions
#' and contributions from Francisco Cribari-Neto.
#'
#' @references
#' Rocha, A.V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. \emph{TEST}, 18(3), 529-545.
#' \doi{10.1007/s11749-008-0112-z}
#'
#' Rocha, A.V., & Cribari-Neto, F. (2017). Erratum to: Beta autoregressive
#' moving average models. \emph{TEST}, 26, 451-459.
#' \doi{10.1007/s11749-017-0528-4}
#'
#' Cribari-Neto, F., Costa, E., & Fonseca, R.V. (2025). Numerical stability
#' enhancements in beta autoregressive moving average model estimation.
#' \emph{Brazilian Journal of Probability and Statistics}, 39(4), 410-437.
#' \doi{10.1214/25-BJPS645}
#'
#' @param y
#' A time series object (`ts`) with values strictly in the open interval (0, 1).
#' Must have at least \eqn{\max(p,q) + 1} observations.
#'
#' @param ar
#' A numeric vector specifying autoregressive (AR) lags
#' (e.g., `c(1, 2)` for AR(2)). Defaults to `integer(0)`, which omits the
#' AR component entirely. Absence should be expressed by omitting this
#' argument or passing `integer(0)`.
#'
#' @param ma
#' A numeric vector specifying moving average (MA) lags
#' (e.g., `1` for MA(1)). Defaults to `integer(0)`, which omits the
#' MA component entirely. Absence should be expressed by omitting this
#' argument or passing `integer(0)`.
#'
#' @param link
#' The link function connecting the mean \eqn{\mu_t} to the linear predictor
#' \eqn{\eta_t}. One of:
#' \itemize{
#'   \item `"logit"` (default): \eqn{g(x) = \log(x/(1-x))}
#'   \item `"probit"`: \eqn{g(x) = \Phi^{-1}(x)}
#'   \item `"cloglog"`: Complementary log-log link
#'   \item `"loglog"`: Log-log link
#' }
#'
#' @param xreg
#' A matrix or data frame of external regressors (covariates), optional.
#' Must have the same number of rows as `y`. If provided, its columns
#' are included in the linear predictor with associated coefficients in `beta`.
#'
#' @param optimization
#' A named list controlling the optimization procedure. Recognized fields:
#' \describe{
#'   \item{\code{method}}{Character string specifying the algorithm. One of
#'     \code{"BFGS"} (default), \code{"L-BFGS-B"}, or \code{"lbfgs"}.}
#'   \item{\code{lower}}{Numeric vector or scalar of lower bounds, used only
#'     by \code{"L-BFGS-B"}. Defaults to \code{-Inf}.}
#'   \item{\code{upper}}{Numeric vector or scalar of upper bounds, used only
#'     by \code{"L-BFGS-B"}. Defaults to \code{Inf}.}
#' }
#'
#' @param penalty
#' Logical. If \code{TRUE}, the ridge penalization scheme of Cribari-Neto,
#' Costa and Fonseca (2025) is applied during conditional maximum likelihood
#' estimation. The penalty term \eqn{(n - a)\lambda_n \|\boldsymbol{\nu}\|^2}
#' is subtracted from the conditional log-likelihood, where
#' \eqn{\lambda_n = 1/(n - a)^{0.9}} and \eqn{\boldsymbol{\nu}} collects
#' \eqn{\alpha}, \eqn{\boldsymbol{\varphi}}, and \eqn{\boldsymbol{\theta}}.
#' The precision parameter \eqn{\phi} and regression coefficients
#' \eqn{\boldsymbol{\beta}} are excluded from penalization.
#' Ridge penalization enhances the curvature of the conditional log-likelihood
#' surface, reducing convergence failures and the occurrence of implausible
#' estimates. Defaults to \code{FALSE}.
#'
#' @return
#' An object of class `"barma"` containing:
#' \item{coef}{Named vector of all estimated parameters, ordered as:
#'   alpha, AR parameters, MA parameters, beta parameters, phi.}
#' \item{vcov}{The variance-covariance matrix of the estimators, computed as
#'   the inverse of the observed Fisher Information Matrix.}
#' \item{model}{A summary table with coefficients, standard errors, z-statistics,
#'   and p-values for hypothesis tests \eqn{H_0: \theta_i = 0}.}
#' \item{fitted}{Fitted conditional mean values as a `ts` object
#'   (NA-padded for the first \eqn{m = \max(p,q)} observations).}
#' \item{muhat}{Alias for `fitted` (fitted mean values).}
#' \item{etahat}{Estimated linear predictor values (full vector, NA-padded).}
#' \item{errorhat}{Estimated errors on predictor scale (full vector, 0-padded).}
#' \item{loglik}{The conditional log-likelihood at the CMLE.}
#' \item{fisher_info_mat}{The expected Fisher Information Matrix
#' (Rocha & Cribari-Neto, 2009).}
#' \item{conv}{Convergence code from `optim` (0 = success).}
#' \item{alpha, beta, varphi, theta, phi}{Individual parameter estimates.}
#' \item{start_values}{Initial parameter values used in optimization.}
#' \item{call}{The original function call.}
#' \item{opt}{Cleaned output list from the optimization procedure.}
#' \item{opt_raw}{Raw output object returned directly by `optim()` or `lbfgs()`.}
#'
#' @seealso
#' \code{\link{simu_barma}} for simulation,
#' \code{\link{loglik_barma}} for likelihood computation,
#' \code{\link{score_vector_barma}} for gradient computation,
#' \code{\link{fim_barma}} for Fisher Information Matrix
#'
#' @examples
#' \donttest{
#' # Example 1: Fit a BAR(1) model (no MA component)
#' set.seed(2025)
#' y_sim_bar1 <- simu_barma(
#'   n      = 250,
#'   alpha  = 0.0,
#'   varphi = 0.6,
#'   phi    = 25.0,
#'   link   = "logit",
#'   freq   = 12
#' )
#'
#' fit_bar1 <- barma(y_sim_bar1, ar = 1, link = "logit")
#' summary(fit_bar1)
#' coef(fit_bar1)
#'
#' # Example 2: Fit a BARMA(1, 1) model
#' set.seed(2025)
#' y_sim_barma11 <- simu_barma(
#'   n      = 250,
#'   alpha  = 0.0,
#'   varphi = 0.6,
#'   theta  = 0.3,
#'   phi    = 25.0,
#'   link   = "logit",
#'   freq   = 12
#' )
#'
#' fit_barma11 <- barma(y_sim_barma11, ar = 1, ma = 1, link = "logit")
#' summary(fit_barma11)
#'
#' # Example 3: Fit a BMA(1) model (no AR component)
#' set.seed(2025)
#' y_sim_bma1 <- simu_barma(
#'   n      = 250,
#'   alpha  = 0.0,
#'   theta  = 0.3,
#'   phi    = 20.0,
#'   link   = "logit",
#'   freq   = 12
#' )
#'
#' fit_bma1 <- barma(y_sim_bma1, ma = 1, link = "logit")
#' summary(fit_bma1)
#'
#' # Example 4: BARMA(1, 1) model with harmonic seasonal regressors
#' hs <- sin(2 * pi * seq_along(y_sim_barma11) / 12)
#' hc <- cos(2 * pi * seq_along(y_sim_barma11) / 12)
#' X <- cbind(hs = hs, hc = hc)
#'
#' fit_barma11_xreg <- barma(
#'   y_sim_barma11,
#'   ar   = 1,
#'   ma   = 1,
#'   link = "logit",
#'   xreg = X
#' )
#' summary(fit_barma11_xreg)
#'
#' # Example 5: Fit a BARMA(1, 1) model using L-BFGS-B
#' set.seed(2025)
#' y_sim_barma11_LBFGSB <- simu_barma(
#'   n      = 250,
#'   alpha  = 0.0,
#'   varphi = 0.6,
#'   theta  = 0.3,
#'   phi    = 25.0,
#'   link   = "logit",
#'   freq   = 12
#' )
#'
#' fit_barma11_LBFGSB <- barma(
#'   y_sim_barma11_LBFGSB,
#'   ar = 1,
#'   ma = 1,
#'   link = "logit",
#'   optimization = list(
#'     method = "L-BFGS-B",
#'     upper = c(Inf, 1, 1, Inf)
#'   )
#' )
#'
#' summary(fit_barma11_LBFGSB)
#' }
#'
#' @importFrom stats is.ts optim dbeta frequency pnorm start ts
#'
#' @export
barma <- function(
    y,
    ar = integer(0),
    ma = integer(0),
    link = "logit",
    xreg = NULL,
    optimization = list(method = "BFGS", lower = -Inf, upper = Inf),
    penalty = FALSE) {
  # 0. CAPTURE CALL ==========================================================
  z <- list(call = match.call())

  # 1. VALIDATE INPUT PARAMETERS =============================================

  # Check time series object -------------------------------------------------
  if (!is.ts(y)) {
    stop(
      "'y' must be a time series object (ts). ",
      "Use ts(data, frequency = ...) to convert."
    )
  }

  # Check for missing values -------------------------------------------------
  if (any(is.na(y))) {
    stop("'y' contains missing values (NA).")
  }

  # Check unit interval bounds -----------------------------------------------
  if (min(y) <= 0 || max(y) >= 1) {
    stop(
      "'y' values must be strictly in (0, 1). ",
      "Current range: [", round(min(y), 4), ", ", round(max(y), 4), "]."
    )
  }

  # Check link function ------------------------------------------------------
  valid_links <- c("logit", "probit", "cloglog", "loglog")
  if (!is.character(link) || !(link %in% valid_links)) {
    stop(
      "'link' must be one of: ",
      paste(valid_links, collapse = ", "), ". ",
      "Current: '", link, "'."
    )
  }

  # Check optimization argument ----------------------------------------------
  if (!is.list(optimization)) {
    stop("'optimization' must be a named list, e.g. list(method = 'BFGS').")
  }

  valid_methods <- c("BFGS", "L-BFGS-B", "lbfgs")

  if (is.null(optimization$method)) {
    stop(
      "'optimization$method' is missing. Must be one of: ",
      paste(valid_methods, collapse = ", "), "."
    )
  }

  if (!optimization$method %in% valid_methods) {
    stop(
      "Unknown method '", optimization$method, "'. Must be one of: ",
      paste(valid_methods, collapse = ", "), "."
    )
  }

  # Check xreg ----------------------------------------------------------------
  if (!is.null(xreg)) {
    if (!is.matrix(xreg)) {
      xreg <- as.matrix(xreg)
    }

    if (nrow(xreg) != length(y)) {
      stop(
        "'xreg' must have ", length(y), " rows to match length of 'y'. ",
        "Current: ", nrow(xreg), "."
      )
    }

    if (any(is.na(xreg))) {
      stop("'xreg' contains missing values (NA).")
    }
  }

  # 2. DETERMINE MODEL STRUCTURE =============================================

  # Backward compatibility ---------------------------------------------------
  # in v1.0.1, absence was expressed as ar = NA or ma = NA.
  # Convert silently to integer(0) and warn the user.

  if (length(ar) == 1 && is.na(ar)) {
    warning(
      "Passing 'ar = NA' is deprecated. ",
      "Use 'ar = integer(0)' or omit the argument to indicate no AR component."
    )
    ar <- integer(0)
  }
  if (length(ma) == 1 && is.na(ma)) {
    warning(
      "Passing 'ma = NA' is deprecated. ",
      "Use 'ma = integer(0)' or omit the argument to indicate no MA component."
    )
    ma <- integer(0)
  }

  # Coerce to integer — enforces correct type if user passes numeric ---------
  ar_lags <- as.integer(ar)
  ma_lags <- as.integer(ma)

  has_xreg <- !is.null(xreg)

  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  n_beta_params <- if (has_xreg) ncol(xreg) else 0L

  has_ar <- n_ar_params > 0
  has_ma <- n_ma_params > 0

  # Maximum lag required to initialize the recursive structure
  ar_order <- if (has_ar) max(ar_lags) else 0L
  ma_order <- if (has_ma) max(ma_lags) else 0L
  max_lag <- max(ar_order, ma_order)

  n_obs <- length(y)

  if (n_obs <= max_lag) {
    stop(
      "Insufficient observations (", n_obs, ") ",
      "for specified lag structure (max_lag = ", max_lag, "). ",
      "Need at least ", max_lag + 1, " observations."
    )
  }

  # Parameter names for output -----------------------------------------------
  names_varphi <- paste0("varphi", ar_lags)
  names_theta <- paste0("theta", ma_lags)

  beta_names <- if (has_xreg) {
    if (!is.null(colnames(xreg))) {
      colnames(xreg)
    } else {
      paste0("beta", seq_len(n_beta_params))
    }
  } else {
    character(0)
  }

  # 3. SETUP LINK FUNCTION AND TRANSFORMED DATA ==============================

  # Get link function structure
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  # Transform response variable to predictor scale
  y_transformed <- linkfun(y)

  # 4. SETUP PARAMETER INDICES ===============================================
  par_idx <- make_par_indices(n_ar_params, n_ma_params, n_beta_params)

  idx_alpha <- par_idx$idx_alpha
  idx_varphi <- par_idx$idx_varphi
  idx_theta <- par_idx$idx_theta
  idx_beta <- par_idx$idx_beta
  idx_phi <- par_idx$idx_phi
  n_params <- par_idx$n_params

  # 5. GET INITIAL PARAMETER VALUES ==========================================

  # Call start_values function to get initial parameter estimates
  init_pars <- start_values(
    y,
    ar = ar_lags,
    ma = ma_lags,
    xreg = xreg,
    link = link
  )

  if (is.null(init_pars) || length(init_pars) != n_params) {
    warning(
      "start_values() returned unexpected length or NULL. ",
      "Optimization may fail or take longer."
    )
  }

  # 6. OPTIMIZATION (CMLE) ===================================================

  z$opt_method <- optimization$method
  z$opt_bounds <- !(all(optimization$lower == -Inf) &&
                      all(optimization$upper == Inf))

  # Call optim with BFGS algorithm using analytic gradient
  opt_list <- run_optimizer(
    y = y,
    init_pars = init_pars,
    ar_lags = ar_lags,
    ma_lags = ma_lags,
    xreg = xreg,
    link = link,
    par_idx = par_idx,
    optimization = optimization,
    penalty = penalty
  )

  # Check convergence
  if (opt_list$opt$convergence != 0) {
    warning(
      "Optimization did not converge (method = '", optimization$method, "'). ",
      "Convergence code: ", opt_list$opt$convergence, ". ",
      "Results may be unreliable. ",
      "Consider checking input data."
    )
  }

  z$opt <- opt_list$opt
  z$opt_raw <- opt_list$opt_raw
  z$conv <- opt_list$opt$convergence
  z$loglik <- opt_list$opt$loglik

  coef_raw <- opt_list$opt$par

  # Extract estimates
  z$alpha <- coef_raw[idx_alpha]
  z$varphi <- coef_raw[idx_varphi]
  z$theta <- coef_raw[idx_theta]
  z$beta <- coef_raw[idx_beta]
  z$phi <- coef_raw[idx_phi]

  z$coef <- coef_raw

  # 7. COMPUTE FISHER INFORMATION MATRIX & FITTED VALUES =====================

  # Call fim_barma to compute Fisher Information Matrix and fitted values
  fim_results <- fim_barma(
    y = y,
    ar = ar_lags,
    ma = ma_lags,
    alpha = z$alpha,
    varphi = z$varphi,
    theta = z$theta,
    phi = z$phi,
    link = link,
    xreg = xreg,
    beta = z$beta,
    penalty = penalty
  )

  # Extract FIM and model fit statistics
  z$fisher_info_mat <- fim_results$fisher_info_mat
  z$fitted <- fim_results$fitted_ts # ts object with NAs
  z$muhat <- fim_results$fitted_ts # Alias for fitted
  z$etahat <- fim_results$etahat_full # Full vector, NA-padded
  z$errorhat <- fim_results$errorhat_full # Full vector, 0-padded

  # 8. STORE ADDITIONAL INFORMATION FOR METHODS ==============================

  z$y <- y
  z$link <- link
  z$start_values <- init_pars
  z$n_params <- n_params
  z$n_obs <- n_obs
  z$max_lag <- max_lag
  z$ar_lags <- ar_lags
  z$ma_lags <- ma_lags
  z$xreg <- xreg # Store for predict() and residuals() methods
  z$penalty <- penalty
  z$lambda <- if (penalty) .barma_ridge_lambda(n_obs, max_lag)

  # 9. ASSIGN CLASS AND RETURN ===============================================

  class(z) <- "barma"

  return(z)
}