#' Fit Beta Autoregressive Moving Average (BARMA) Models
#'
#' @title
#' Fit Beta Autoregressive Moving Average (BARMA) Models via Maximum Likelihood
#'
#' @description
#' Fits a Beta Autoregressive Moving Average (BARMA) model to time series data
#' valued in (0, 1) using Maximum Likelihood Estimation (MLE). The function
#' performs complete model estimation including parameter estimation, hypothesis
#' testing infrastructure, and model diagnostics.
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
#' **Optimization**: The function uses the BFGS quasi-Newton algorithm via
#' \code{\link{optim}} with analytic gradients for efficiency. Initial values
#' are obtained from the `start_values()` function.
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
#' \item{loglik}{The conditional log-likelihood at the MLE.}
#' \item{fisher_info_mat}{The observed Fisher Information Matrix.}
#' \item{conv}{Convergence code from `optim` (0 = success).}
#' \item{alpha, beta, varphi, theta, phi}{Individual parameter estimates.}
#' \item{start_values}{Initial parameter values used in optimization.}
#' \item{call}{The original function call.}
#' \item{opt}{Raw output object from `optim()` call.}
#'
#' @seealso
#' \code{\link{simu_barma}} for simulation,
#' \code{\link{loglik_barma}} for likelihood computation,
#' \code{\link{score_vector_barma}} for gradient computation,
#' \code{\link{fim_barma}} for Fisher Information Matrix
#'
#' @examples
#' \donttest{
#'   # Example 1: Fit a BAR(1) model (no MA component)
#'   set.seed(2025)
#'   y_sim_bar1 <- simu_barma(
#'     n      = 250,
#'     alpha  = 0.0,
#'     varphi = 0.6,
#'     phi    = 25.0,
#'     link   = "logit",
#'     freq   = 12
#'   )
#'
#'   fit_bar1 <- barma(y_sim_bar1, ar = 1, link = "logit")
#'   summary(fit_bar1)
#'   coef(fit_bar1)
#'
#'   # Example 2: Fit a BARMA(1, 1) model
#'   set.seed(2025)
#'   y_sim_barma11 <- simu_barma(
#'     n      = 250,
#'     alpha  = 0.0,
#'     varphi = 0.6,
#'     theta  = 0.3,
#'     phi    = 25.0,
#'     link   = "logit",
#'     freq   = 12
#'   )
#'
#'   fit_barma11 <- barma(y_sim_barma11, ar = 1, ma = 1, link = "logit")
#'   summary(fit_barma11)
#'
#'   # Example 3: Fit a BMA(1) model (no AR component)
#'   set.seed(2025)
#'   y_sim_bma1 <- simu_barma(
#'     n      = 250,
#'     alpha  = 0.0,
#'     theta  = 0.3,
#'     phi    = 20.0,
#'     link   = "logit",
#'     freq   = 12
#'   )
#'
#'   fit_bma1 <- barma(y_sim_bma1, ma = 1, link = "logit")
#'   summary(fit_bma1)
#'
#'   # Example 4: BARMA(1, 1) model with harmonic seasonal regressors
#'   hs <- sin(2 * pi * seq_along(y_sim_barma11) / 12)
#'   hc <- cos(2 * pi * seq_along(y_sim_barma11) / 12)
#'   X  <- cbind(hs = hs, hc = hc)
#'
#'   fit_barma11_xreg <- barma(
#'     y_sim_barma11,
#'     ar   = 1,
#'     ma   = 1,
#'     link = "logit",
#'     xreg = X
#'   )
#'   summary(fit_barma11_xreg)
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
    xreg = NULL
) {
  
  # --------------------------------------------------------------------------
  # 1. VALIDATE INPUT PARAMETERS
  # --------------------------------------------------------------------------
  
  # Check time series object
  if (!is.ts(y)) {
    stop(
      "y must be a time series object (ts class). ",
      "Use ts(data, frequency = ...) to convert."
    )
  }
  
  # Check bounds (0, 1)
  if (any(is.na(y))) {
    stop("y contains missing values (NA).")
  }
  
  if (min(y) <= 0 || max(y) >= 1) {
    stop(
      "Response values must be strictly in (0, 1). ",
      "Current range: [", round(min(y), 4), ", ", round(max(y), 4), "]."
    )
  }
  
  # Check link function
  if (!is.character(link) || !(link %in% c("logit", 
                                           "probit", 
                                           "cloglog", 
                                           "loglog"))) {
    stop(
      "link must be one of: 'logit', 'probit', 'cloglog', or 'loglog'. ",
      "Current: '", link, "'"
    )
  }
  
  # Store original call
  z <- list(call = match.call())
  
  # --------------------------------------------------------------------------
  # 2. DETERMINE MODEL STRUCTURE
  # --------------------------------------------------------------------------
  
  # Backward compatibility: in v1.0.1, absence was expressed as ar = NA or
  # ma = NA. Convert silently to integer(0) and warn the user.
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
  
  # Resolve lag vectors to integer(0) if absent.
  ar_lags <- if (length(ar) > 0) ar else integer(0)
  ma_lags <- if (length(ma) > 0) ma else integer(0)
  
  has_xreg <- !is.null(xreg)
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  # Named flags for readability throughout the function
  has_ar <- n_ar_params > 0
  has_ma <- n_ma_params > 0
  
  # Setup and Validate Regressors
  if (has_xreg) {
    if (!is.matrix(xreg)) {
      xreg <- as.matrix(xreg)
    }
    
    if (nrow(xreg) != length(y)) {
      stop(
        "Number of rows in xreg (", nrow(xreg), ") ",
        "must match length of y (", length(y), ")."
      )
    }
    
    if (any(is.na(xreg))) {
      stop("xreg contains missing values (NA).")
    }
    
    n_beta_params <- ncol(xreg)
    beta_names <- colnames(xreg)
    if (is.null(beta_names)) {
      beta_names <- paste0("beta", 1:n_beta_params)
    }
  } else {
    n_beta_params <- 0
    beta_names <- character(0)
  }
  
  # Create parameter names for output
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta  <- if (has_ma) paste0("theta",  ma_lags) else character(0)
  
  # --------------------------------------------------------------------------
  # 3. SETUP TIME SERIES PROPERTIES
  # --------------------------------------------------------------------------
  
  # Determine maximum lag required to initialize the recursive structure
  ar_order <- if (has_ar) max(ar_lags) else 0L
  ma_order <- if (has_ma) max(ma_lags) else 0L
  max_lag  <- max(ar_order, ma_order)
  
  n_obs <- length(y)
  
  # Check for sufficient observations
  if (n_obs <= max_lag) {
    stop(
      "Insufficient observations (", n_obs, ") ",
      "for specified lag structure (max_lag = ", max_lag, "). ",
      "Need at least ", max_lag + 1, " observations."
    )
  }
  
  # --------------------------------------------------------------------------
  # 4. SETUP LINK FUNCTION AND TRANSFORMED DATA
  # --------------------------------------------------------------------------
  
  # Get link function structure
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  
  # Transform response variable to predictor scale
  y_transformed <- linkfun(y)
  
  # --------------------------------------------------------------------------
  # 5. SETUP PARAMETER INDICES (FOR OPTIMIZATION)
  # --------------------------------------------------------------------------
  
  # Parameter order (matching start_values order):
  # 1. alpha (intercept)
  # 2. AR parameters (if present)
  # 3. MA parameters (if present)
  # 4. beta (regression coefficients, if present)
  # 5. phi (precision parameter)
  
  idx_alpha <- 1
  
  idx_varphi <- if (has_ar) {
    2:(1 + n_ar_params)
  } else {
    integer(0)
  }
  
  last_idx <- if (has_ar) max(idx_varphi) else idx_alpha
  
  idx_theta <- if (has_ma) {
    (last_idx + 1):(last_idx + n_ma_params)
  } else {
    integer(0)
  }
  
  if (has_ma) last_idx <- max(idx_theta)
  
  # Regression coefficients (if present)
  idx_beta <- if (has_xreg) {
    (last_idx + 1):(last_idx + n_beta_params)
  } else {
    integer(0)
  }
  
  if (has_xreg) last_idx <- max(idx_beta)
  
  # Precision parameter (always last)
  idx_phi <- last_idx + 1
  last_idx <- idx_phi
  
  n_params <- last_idx
  
  # --------------------------------------------------------------------------
  # 6. GET INITIAL PARAMETER VALUES
  # --------------------------------------------------------------------------
  
  # Call start_values function to get initial parameter estimates
  init_pars <- start_values(
    y,
    link = link,
    ar = ar_lags,
    ma = ma_lags,
    xreg = xreg
  )
  
  if (is.null(init_pars) || length(init_pars) != n_params) {
    warning(
      "start_values() returned unexpected length or NULL. ",
      "Optimization may fail or take longer."
    )
  }
  
  # --------------------------------------------------------------------------
  # 7. OPTIMIZATION (MLE)
  # --------------------------------------------------------------------------
  
  # Call optim with BFGS algorithm using analytic gradient
  opt <- optim(
    par = init_pars,
    fn = function(x) {
      # Negative log-likelihood (for minimization)
      (-1) * loglik_barma(
        y = y,
        ar = ar_lags,
        ma = ma_lags,
        alpha = x[idx_alpha],
        varphi = x[idx_varphi],
        theta = x[idx_theta],
        phi = x[idx_phi],
        link = link,
        xreg = xreg,
        beta = x[idx_beta]
      )
    },
    gr = function(x) {
      # Negative score vector (for minimization)
      (-1) * score_vector_barma(
        y = y,
        ar = ar_lags,
        ma = ma_lags,
        alpha = x[idx_alpha],
        varphi = x[idx_varphi],
        theta = x[idx_theta],
        phi = x[idx_phi],
        link = link,
        xreg = xreg,
        beta = x[idx_beta]
      )
    },
    method = "BFGS"
  )
  
  # Check convergence
  if (opt$convergence != 0) {
    warning(
      "BFGS optimization did not converge. ",
      "Convergence code: ", opt$convergence, ". ",
      "Results may be unreliable. ",
      "Consider checking input data or initial values."
    )
  }
  
  z$conv <- opt$convergence
  z$opt <- opt
  z$loglik <- -1 * opt$value
  
  # --------------------------------------------------------------------------
  # 8. EXTRACT AND ORGANIZE FINAL PARAMETERS
  # --------------------------------------------------------------------------
  
  # Extract raw parameters from optimizer
  coef_raw <- opt$par
  
  # Extract individual parameter vectors
  z$alpha  <- coef_raw[idx_alpha]
  z$varphi <- if (has_ar) coef_raw[idx_varphi] else numeric(0)
  z$theta  <- if (has_ma) coef_raw[idx_theta] else numeric(0)
  z$phi    <- coef_raw[idx_phi]
  z$beta   <- if (has_xreg) coef_raw[idx_beta] else numeric(0)
  
  # Create named coefficient vector for user output
  # Order: alpha, AR, MA, beta, phi
  coef_final <- c(z$alpha, z$varphi, z$theta, z$beta, z$phi)
  coef_names <- c("alpha", names_varphi, names_theta, beta_names, "phi")
  names(coef_final) <- coef_names
  
  z$coef <- coef_final
  
  # --------------------------------------------------------------------------
  # 9. COMPUTE FISHER INFORMATION MATRIX & FITTED VALUES
  # --------------------------------------------------------------------------
  
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
    beta = z$beta
  )
  
  # Extract FIM and model fit statistics
  z$fisher_info_mat <- fim_results$fisher_info_mat
  z$fitted   <- fim_results$fitted_ts        # ts object with NAs
  z$muhat    <- fim_results$fitted_ts        # Alias for fitted
  z$etahat   <- fim_results$etahat_full      # Full vector, NA-padded
  z$errorhat <- fim_results$errorhat_full    # Full vector, 0-padded
  
  # --------------------------------------------------------------------------
  # 10. STORE ADDITIONAL INFORMATION FOR METHODS
  # --------------------------------------------------------------------------
  
  z$y <- y
  z$link <- link
  z$start_values <- init_pars
  z$n_params <- n_params
  z$n_obs <- n_obs
  z$max_lag <- max_lag
  z$ar_lags <- ar_lags
  z$ma_lags <- ma_lags
  z$xreg <- xreg  # Store for predict() and residuals() methods
  
  # --------------------------------------------------------------------------
  # 11. ASSIGN CLASS AND RETURN
  # --------------------------------------------------------------------------
  
  class(z) <- "barma"
  
  return(z)
}