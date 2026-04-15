#' @title Score Vector for the BARMA Model
#' @description Computes the score vector (gradient of the log-likelihood)
#'   for the Beta Autoregressive Moving Average (BARMA) model at a given
#'   parameter vector. This function is designed for users who:
#' \itemize{
#'   \item Implement custom optimization algorithms
#'   \item Verify theoretical properties
#'   \item Conduct simulation studies
#'   \item Debug model fitting
#'   \item Integrate BARMA models into their own workflows
#' }
#'
#' @param y A time series object (\code{ts}) with values strictly in (0, 1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#'   Defaults to \code{integer(0)}, which omits the AR component entirely.
#'   Absence should be expressed by omitting this argument or passing
#'   \code{integer(0)}.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#'   Defaults to \code{integer(0)}, which omits the MA component entirely.
#'   Absence should be expressed by omitting this argument or passing
#'   \code{integer(0)}.
#' @param alpha The intercept parameter.
#' @param varphi A numeric vector of AR parameters.
#'   Absence should be expressed by omitting this argument or passing
#'   \code{numeric(0)}.
#' @param theta A numeric vector of MA parameters.
#'   Absence should be expressed by omitting this argument or passing
#'   \code{numeric(0)}.
#' @param phi The precision parameter (must be positive).
#' @param link A character string specifying the link function.
#'   One of \code{"logit"}, \code{"probit"}, \code{"cloglog"},
#'   or \code{"loglog"}.
#' @param xreg An optional matrix of external regressors.
#' @param beta An optional numeric vector of regression coefficients
#'   corresponding to \code{xreg}.
#'
#' @return A numeric vector of the same length as the parameter vector,
#'   giving the partial derivatives of the log-likelihood with respect
#'   to each parameter. The order of the components is:
#'   \code{(alpha, varphi, theta, beta, phi)}.
#'
#' @seealso \code{\link{barma}}, \code{\link{loglik_barma}},
#'   \code{\link{fim_barma}}
#'
#' @examples
#' \donttest{
#'   # Example 1: Score vector for a BAR(1) model (no MA component)
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
#'   score_vector_barma(
#'     y      = y_sim_bar1,
#'     ar     = 1,
#'     alpha  = 0.0,
#'     varphi = 0.6,
#'     theta  = numeric(0),
#'     phi    = 25.0,
#'     link   = "logit"
#'   )
#'
#'   # Example 2: Score vector for a BARMA(1, 1) model
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
#'   score_vector_barma(
#'     y      = y_sim_barma11,
#'     ar     = 1,
#'     ma     = 1,
#'     alpha  = 0.0,
#'     varphi = 0.6,
#'     theta  = 0.3,
#'     phi    = 25.0,
#'     link   = "logit"
#'   )
#'
#'   # Example 3: Score vector for a BMA(1) model (no AR component)
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
#'   score_vector_barma(
#'     y      = y_sim_bma1,
#'     ma     = 1,
#'     alpha  = 0.0,
#'     varphi = numeric(0),
#'     theta  = 0.6,
#'     phi    = 25.0,
#'     link   = "logit"
#'   )
#' }
#'
#' @export
score_vector_barma <- function(y, ar = integer(0), ma = integer(0),
                               alpha, varphi, theta, phi, link,
                               xreg = NULL, beta = NULL) {
  
  # ------------------------------------------------------------------------
  # 1. Validate Precision Parameter
  # ------------------------------------------------------------------------
  if (phi <= 0 || !is.finite(phi)) {
    warning("phi must be positive and finite; returning zero gradient")
    n_params <- 1 + length(varphi) + length(theta) + 1 + length(beta)
    return(rep(0, n_params))
  }
  
  # ------------------------------------------------------------------------
  # 2. Determine Model Structure
  # ------------------------------------------------------------------------
  
  # Resolve lag vectors to integer(0) if absent.
  ar_lags <- if (length(ar) > 0) ar else integer(0)
  ma_lags <- if (length(ma) > 0) ma else integer(0)
  
  # Force varphi / theta to numeric(0) when the corresponding component
  # is absent, so that drop(crossprod(numeric(0), numeric(0))) == 0.
  varphi <- if (length(ar_lags) > 0) varphi else numeric(0)
  theta  <- if (length(ma_lags) > 0) theta  else numeric(0)
  
  # Setup Regressors
  has_xreg <- !is.null(xreg)
  
  if (has_xreg) {
    if (!is.matrix(xreg)) xreg <- as.matrix(xreg)
    if (is.null(beta)) stop("If 'xreg' is provided, 'beta' is required.")
    if (ncol(xreg) != length(beta)) {
      stop("Length of 'beta' must match columns of 'xreg'.")
    }
    if (nrow(xreg) != length(y)) {
      stop("Number of rows in 'xreg' must match length of 'y'.")
    }
    # Pre-compute X * beta for efficiency
    xb <- as.vector(xreg %*% beta)
    n_beta_params <- length(beta)
  } else {
    xb <- numeric(length(y))
    n_beta_params <- 0
  }
  
  n_ar_params <- length(varphi)
  n_ma_params <- length(theta)
  
  # Validate parameter-lag consistency
  if (n_ar_params != length(ar_lags)) {
    stop("Mismatch between 'ar' lags and 'varphi' parameters.")
  }
  if (n_ma_params != length(ma_lags)) {
    stop("Mismatch between 'ma' lags and 'theta' parameters.")
  }
  
  # ------------------------------------------------------------------------
  # 3. Setup Link Functions and Time Series Properties
  # ------------------------------------------------------------------------
  
  # Get link function structures
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu_eta_fun <- link_structure$mu.eta
  
  # Transform response using link function
  y_transformed <- linkfun(y)
  
  # Determine maximum lag
  ar_order <- if (length(ar_lags) > 0) max(ar_lags) else 0
  ma_order <- if (length(ma_lags) > 0) max(ma_lags) else 0
  max_lag  <- max(ar_order, ma_order)
  n_obs <- length(y)
  
  if (n_obs <= max_lag) {
    warning("Insufficient observations for the specified lag structure")
    n_params <- 1 + n_ar_params + n_ma_params + 1 + n_beta_params
    return(rep(0, n_params))
  }
  
  # ------------------------------------------------------------------------
  # 4. Recursive Calculation: Predictor, Error, and Derivatives
  # ------------------------------------------------------------------------
  error <- rep(0, n_obs)
  eta   <- rep(NA_real_, n_obs)
  
  # Initialize derivative matrices
  d_eta_d_alpha  <- rep(0, n_obs)
  d_eta_d_varphi <- matrix(0, nrow = n_obs, ncol = n_ar_params)
  d_eta_d_theta  <- matrix(0, nrow = n_obs, ncol = n_ma_params)
  d_eta_d_beta   <- matrix(0, nrow = n_obs, ncol = n_beta_params)
  
  for (t in (max_lag + 1):n_obs) {
    
    # Part A: Compute Linear Predictor (eta)
    eta[t] <- alpha + xb[t] +
      drop(crossprod(varphi, y_transformed[t - ar_lags] - xb[t - ar_lags])) +
      drop(crossprod(theta,  error[t - ma_lags]))
    
    error[t] <- y_transformed[t] - eta[t]
    
    # Part B: Compute Recursive Derivatives
    
    # 1. Derivative w.r.t. alpha
    d_eta_d_alpha[t] <- 1 -
      drop(crossprod(theta, d_eta_d_alpha[t - ma_lags]))
    
    # 2. Derivative w.r.t. AR (varphi)
    if (n_ar_params > 0) {
      d_eta_d_varphi[t, ] <- y_transformed[t - ar_lags] - xb[t - ar_lags] -
        drop(crossprod(theta, d_eta_d_varphi[t - ma_lags, , drop = FALSE]))
    }
    
    # 3. Derivative w.r.t. MA (theta)
    if (n_ma_params > 0) {
      d_eta_d_theta[t, ] <- error[t - ma_lags] -
        drop(crossprod(theta, d_eta_d_theta[t - ma_lags, , drop = FALSE]))
    }
    
    # 4. Derivative w.r.t. beta (regressors)
    if (n_beta_params > 0) {
      d_eta_d_beta[t, ] <- xreg[t, ] -
        drop(crossprod(varphi, xreg[t - ar_lags, , drop = FALSE])) -
        drop(crossprod(theta,  d_eta_d_beta[t - ma_lags, , drop = FALSE]))
    }
  }
  
  # ------------------------------------------------------------------------
  # 5. Calculate Score Vector Components
  # ------------------------------------------------------------------------
  eta_effective <- eta[(max_lag + 1):n_obs]
  y_effective   <- y[(max_lag + 1):n_obs]
  mu_effective  <- linkinv(eta = eta_effective)
  
  # Numerical stability check
  if (any(mu_effective <= 0 | mu_effective >= 1 | !is.finite(mu_effective))) {
    warning("mu_effective out of bounds; returning zero gradient")
    n_params <- 1 + n_ar_params + n_ma_params + 1 + n_beta_params
    return(rep(0, n_params))
  }
  
  mu_eta_val <- mu_eta_fun(eta = eta_effective)
  ystar  <- log(y_effective / (1 - y_effective))
  mustar <- digamma(mu_effective * phi) - digamma((1 - mu_effective) * phi)
  
  # Chain rule common term: dL/deta = dL/dmu * dmu/deta
  common_term <- mu_eta_val * (ystar - mustar)
  
  # Compute final scores. Using crossprod for sum(vector * vector)
  
  idx <- (max_lag + 1):n_obs
  
  score_alpha <- as.numeric(
    phi * crossprod(d_eta_d_alpha[idx], common_term)
  )
  
  score_varphi <- as.numeric(
    phi * crossprod(d_eta_d_varphi[idx, , drop = FALSE], common_term)
  )
  
  score_theta <- as.numeric(
    phi * crossprod(d_eta_d_theta[idx, , drop = FALSE], common_term)
  )
  
  score_beta <- as.numeric(
    phi * crossprod(d_eta_d_beta[idx, , drop = FALSE], common_term)
  )
  
  score_phi <- sum(
    mu_effective * (ystar - mustar) +
      log(1 - y_effective) -
      digamma((1 - mu_effective) * phi) +
      digamma(phi)
  )
  
  # ------------------------------------------------------------------------
  # 6. Assemble and Return
  # ------------------------------------------------------------------------
  final_score <- c(score_alpha,
                   score_varphi, 
                   score_theta, 
                   score_beta,
                   score_phi)
  
  if (any(!is.finite(final_score))) {
    warning("Non-finite values in score vector; returning zeros")
    return(rep(0, length(final_score)))
  }
  
  return(as.numeric(final_score))
}