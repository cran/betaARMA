#' Log-Likelihood for BARMA Models
#'
#' @title
#' Compute Conditional Log-Likelihood for Beta Autoregressive
#' Moving Average Models
#'
#' @description
#' Computes the conditional log-likelihood of a Beta Autoregressive
#' Moving Average (BARMA) model. This function is designed for users who:
#' \itemize{
#'   \item Implement custom optimization algorithms
#'   \item Verify theoretical properties
#'   \item Conduct simulation studies
#'   \item Debug model fitting
#'   \item Integrate BARMA models into their own workflows
#' }
#'
#' @details
#' The log-likelihood is computed as:
#' \deqn{\ell = \sum_{t=m+1}^{n} \log f(y_t | \mu_t, \phi)}
#' 
#' where \eqn{f} is the Beta density with shape parameters
#' \eqn{shape1 = \mu_t * \phi} and \eqn{shape2 = (1-\mu_t)*\phi}.
#'
#' The linear predictor is constructed as:
#' \deqn{\eta_t = \alpha + X_t \beta + \sum_{i=1}^{p} \varphi_i 
#' (y_{t-i} - X_{t-i}\beta) + \sum_{j=1}^{q} \theta_j \epsilon_{t-j}}
#'
#' **Important**: This function implements the corrections from the 2017
#' Erratum (Rocha & Cribari-Neto, 2017) for moving average components.
#' See References section for details.
#'
#' **Parameter Order**: Parameters should be supplied in the order:
#' alpha, varphi (AR), theta (MA), phi, beta (regressors).
#' This matches the parameter order used by \code{\link{barma}}.
#'
#' @param y
#' A numeric vector representing the time series data, with values
#' strictly in (0, 1).
#'
#' @param ar
#' A numeric vector specifying the autoregressive (AR) lags
#' (e.g., \code{c(1, 2)}). Can be \code{NA} or \code{NULL} if no AR component.
#'
#' @param ma
#' A numeric vector specifying the moving average (MA) lags
#' (e.g., \code{1}). Can be \code{NA} or \code{NULL} if no MA component.
#'
#' @param alpha
#' The intercept term (numeric scalar).
#'
#' @param varphi
#' A numeric vector of autoregressive (AR) parameters.
#' Use \code{numeric(0)} or empty vector if no AR component.
#' Length must match \code{ar} specification.
#'
#' @param theta
#' A numeric vector of moving average (MA) parameters.
#' Use \code{numeric(0)} or empty vector if no MA component.
#' Length must match \code{ma} specification.
#'
#' @param phi
#' The precision parameter of the BARMA model (must be positive and finite).
#' Larger values indicate less variance for a given mean.
#'
#' @param link
#' A character string specifying the link function:
#' \code{"logit"} (default), \code{"probit"}, or \code{"cloglog"}.
#'
#' @param xreg
#' A matrix or data frame of static regressors (optional).
#' Must have the same number of rows as length of \code{y}.
#'
#' @param beta
#' A numeric vector of regression coefficients for \code{xreg} (optional).
#' Length must match number of columns in \code{xreg}.
#'
#' @return
#' A numeric scalar representing the conditional log-likelihood value.
#' Returns \code{-Inf} if:
#' \itemize{
#'   \item \code{phi} is non-positive or non-finite
#'   \item Insufficient observations for specified lag structure
#'   \item Fitted values are outside (0, 1)
#'   \item Any numerical issues occur
#' }
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
#' @seealso
#' \code{\link{barma}} for model fitting,
#' \code{\link{score_vector_barma}} for gradient computation,
#' \code{\link{fim_barma}} for Fisher Information Matrix
#'
#' @examples
#' \donttest{
#'   # Example 1: Log-likelihood for a BAR(1) model
#'   set.seed(2025)
#'   y_sim_bar <- simu_barma(
#'     n = 250,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     phi = 25.0,
#'     link = "logit",
#'     freq = 12
#'   )
#'
#'   loglik_barma(
#'     y = y_sim_bar,
#'     ar = 1,
#'     ma = NA,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = numeric(0),
#'     phi = 25.0,
#'     link = "logit"
#'   )
#'
#'   # Example 2: Log-likelihood for a BARMA(1,1) model
#'   set.seed(2025)
#'   y_sim_barma <- simu_barma(
#'     n = 250,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = 0.3,
#'     phi = 25.0,
#'     link = "logit",
#'     freq = 12
#'   )
#'
#'   loglik_barma(
#'     y = y_sim_barma,
#'     ar = 1,
#'     ma = 1,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = 0.3,
#'     phi = 25.0,
#'     link = "logit"
#'   )
#'
#' }
#'
#' @export
loglik_barma <- function(
    y,
    ar,
    ma,
    alpha,
    varphi,
    theta,
    phi,
    link,
    xreg = NULL,
    beta = NULL
) {
  
  # --------------------------------------------------------------------------
  # 1. VALIDATE INPUT PARAMETERS 
  # --------------------------------------------------------------------------
  
  # Validate response variable
  if (!is.numeric(y)) {
    stop("y must be a numeric vector.")
  }
  if (length(y) < 2) {
    stop("y must have at least 2 observations.")
  }
  if (anyNA(y)) {
    stop("y contains missing values (NA).")
  }
  if (min(y) <= 0 || max(y) >= 1) {
    stop(
      "Response values must be strictly in (0, 1). ",
      "Current range: [", round(min(y), 4), ", ", round(max(y), 4), "]."
    )
  }
  
  # Validate precision parameter
  if (!is.numeric(phi) || length(phi) != 1) {
    stop("phi must be a single numeric value.")
  }
  if (phi <= 0 || !is.finite(phi)) {
    return(-Inf)
  }
  
  # --------------------------------------------------------------------------
  # 2. DETERMINE MODEL STRUCTURE 
  # --------------------------------------------------------------------------
  
  # Check for presence of AR, MA, and external regressors
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  has_xreg <- !is.null(xreg)
  
  # Handle empty lag cases
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  n_ar_params <- length(varphi)
  n_ma_params <- length(theta)
  
  # Validate parameter-lag consistency ---
  if (has_ar && n_ar_params != length(ar_lags)) {
    stop(
      "Length of 'varphi' (", n_ar_params, ") ",
      "must match length of 'ar' (", length(ar_lags), ")."
    )
  }
  if (has_ma && n_ma_params != length(ma_lags)) {
    stop(
      "Length of 'theta' (", n_ma_params, ") ",
      "must match length of 'ma' (", length(ma_lags), ")."
    )
  }
  
  # Setup Regressors ---
  if (has_xreg) {
    if (!is.matrix(xreg)) xreg <- as.matrix(xreg)
    
    if (is.null(beta)) {
      stop("If 'xreg' is provided, 'beta' must also be provided.")
    }
    if (!is.numeric(beta)) {
      stop("beta must be numeric.")
    }
    if (ncol(xreg) != length(beta)) {
      stop(
        "Length of 'beta' (", length(beta), ") ",
        "must match columns of 'xreg' (", ncol(xreg), ")."
      )
    }
    if (nrow(xreg) != length(y)) {
      stop(
        "Number of rows in 'xreg' (", nrow(xreg), ") ",
        "must match length of 'y' (", length(y), ")."
      )
    }
    if (anyNA(xreg)) {
      stop("xreg contains missing values (NA).")
    }
    
    # Pre-compute X * beta for efficiency (vectorized operation)
    xb <- as.vector(xreg %*% beta)
  } else {
    # Zero contribution from regressors if none provided
    xb <- numeric(length(y))
  }
  
  # --------------------------------------------------------------------------
  # 3. SETUP LINK FUNCTIONS AND TIME SERIES PROPERTIES
  # --------------------------------------------------------------------------
  
  # Get link function structures
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  
  # Transform response using link function
  y_transformed <- linkfun(y)
  
  # Determine maximum lag for burn-in period
  ar_order <- if (has_ar) max(ar_lags) else 0
  ma_order <- if (has_ma) max(ma_lags) else 0
  max_lag  <- max(ar_order, ma_order)
  n_obs <- length(y)
  
  # Check for sufficient observations
  if (n_obs <= max_lag) {
    warning(
      "Insufficient observations: have ", n_obs, 
      " observations but need > ", max_lag, " (max lag order)."
    )
    return(-Inf)
  }
  
  # --------------------------------------------------------------------------
  # 4. CALCULATE ERROR AND PREDICTOR ITERATIVELY
  # --------------------------------------------------------------------------
  
  # Initialize containers for linear predictor and errors
  # Observations 1:max_lag will have 0/NA values (burn-in period)
  error <- rep(0, n_obs)
  eta   <- rep(NA_real_, n_obs)
  
  # Recursively compute predictor and errors
  for (t in (max_lag + 1):n_obs) {
    
    # Part A: Base Linear Predictor ---
    # Initialize with intercept and regressor contribution
    eta[t] <- alpha + xb[t]
    
    # Part B: Autoregressive Component ---
    # AR(p) term: acts on lagged transformed response, adjusted for regression
    if (has_ar) {
      # Extract lagged transformed response, adjusted for regression effect
      prev_terms <- y_transformed[t - ar_lags] - xb[t - ar_lags]
      # Add AR contribution: sum(varphi_i * (y_{t-i} - X_{t-i}*beta))
      eta[t] <- eta[t] + drop(crossprod(varphi, prev_terms))
    }
    
    # Part C: Moving Average Component ---
    # MA(q) term: acts on lagged prediction errors
    if (has_ma) {
      eta[t] <- eta[t] + drop(crossprod(theta, error[t - ma_lags]))
    }
    
    # Part D: Compute Prediction Error ---
    # Error is deviation of transformed response from predictor
    error[t] <- y_transformed[t] - eta[t]
  }
  
  # --------------------------------------------------------------------------
  # 5. CALCULATE THE FINAL LOG-LIKELIHOOD
  # --------------------------------------------------------------------------
  
  # Extract observations used in likelihood (excluding burn-in period)
  idx_effective <- (max_lag + 1):n_obs
  eta_eff <- eta[idx_effective]
  y_eff   <- y[idx_effective]
  
  # Transform linear predictor to mean scale
  mu_eff <- linkinv(eta = eta_eff)
  
  # Check for valid mu values (must be in (0,1))
  if (any(mu_eff <= 0 | mu_eff >= 1 | !is.finite(mu_eff))) {
    return(-Inf)
  }
  
  # Compute log-likelihood using Beta density
  # dbeta parameterization: shape1 = mu*phi, shape2 = (1-mu)*phi
  ll_terms <- dbeta(
    y_eff,
    shape1 = mu_eff * phi,
    shape2 = (1 - mu_eff) * phi,
    log = TRUE
  )
  
  # Check for numerical issues (non-finite values)
  if (any(!is.finite(ll_terms))) {
    return(-Inf)
  }
  
  # Return sum of log-likelihood terms
  final_loglik <- sum(ll_terms)
  return(final_loglik)
}
