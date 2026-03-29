#' @title Score Vector for the BARMA Model
#' @description Computes the score vector (gradient of the log-likelihood)
#'   for the Beta Autoregressive Moving Average (BARMA) model at a given
#'   parameter vector. Used internally by \code{\link{barma}} during
#'   optimization via the BFGS algorithm.
#'
#' @param y A time series object (\code{ts}) with values strictly in (0, 1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#'   Use \code{NA} or \code{NULL} if no AR component.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#'   Use \code{NA} or \code{NULL} if no MA component.
#' @param alpha The intercept parameter.
#' @param varphi A numeric vector of AR parameters.
#'   Use \code{numeric(0)} if no AR component.
#' @param theta A numeric vector of MA parameters.
#'   Use \code{numeric(0)} if no MA component.
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
#'   to each parameter.
#'
#' @seealso \code{\link{barma}}, \code{\link{loglik_barma}},
#'   \code{\link{fim_barma}}
#'
#' @keywords internal
score_vector_barma <- function(y, ar, ma, alpha, varphi, theta, phi, link,
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
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
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
  
  # Get lag specifications
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  # Validate parameter-lag consistency
  if (has_ar && n_ar_params != length(ar_lags)) {
    stop("Mismatch between 'ar' lags and 'varphi' parameters.")
  }
  if (has_ma && n_ma_params != length(ma_lags)) {
    stop("Mismatch between 'ma' lags and 'theta' parameters.")
  }
  
  # ------------------------------------------------------------------------
  # 3. Setup Link Functions and Time Series Properties
  # ------------------------------------------------------------------------
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu_eta_fun <- link_structure$mu.eta
  ynew <- linkfun(y)
  
  # Determine maximum lag
  ar_order <- if (has_ar) max(ar_lags) else 0
  ma_order <- if (has_ma) max(ma_lags) else 0
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
  
  d_eta_d_varphi <- if (has_ar) {
    matrix(0, nrow = n_obs, ncol = n_ar_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  d_eta_d_theta <- if (has_ma) {
    matrix(0, nrow = n_obs, ncol = n_ma_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  d_eta_d_beta <- if (has_xreg) {
    matrix(0, nrow = n_obs, ncol = n_beta_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  for (t in (max_lag + 1):n_obs) {
    
    # Part A: Compute Linear Predictor (Eta) ---
    # Model: alpha + X*beta + AR(y - X*beta) + MA(error)
    eta[t] <- alpha + xb[t]
    
    if (has_ar) {
      # AR term acts on the "regression residual" (ynew - xb)
      prev_terms <- ynew[t - ar_lags] - xb[t - ar_lags]
      eta[t] <- eta[t] + as.numeric(crossprod(varphi, prev_terms))
    }
    
    if (has_ma) {
      eta[t] <- eta[t] + as.numeric(crossprod(theta, error[t - ma_lags]))
    }
    
    error[t] <- ynew[t] - eta[t]
    
    # Part B: Compute Recursive Derivatives ---
    
    # 1. Derivative w.r.t. Alpha
    d_eta_d_alpha[t] <- 1
    if (has_ma) {
      d_eta_d_alpha[t] <- 1 - as.numeric(
        crossprod(theta, d_eta_d_alpha[t - ma_lags])
      )
    }
    
    # 2. Derivative w.r.t. AR (Varphi)
    if (has_ar) {
      # Base: ynew - X*beta at lags
      d_eta_d_varphi[t, ] <- ynew[t - ar_lags] - xb[t - ar_lags]
      if (has_ma) {
        ma_effect <- crossprod(theta, d_eta_d_varphi[t - ma_lags, , drop = FALSE])
        d_eta_d_varphi[t, ] <- d_eta_d_varphi[t, ] - as.numeric(ma_effect)
      }
    }
    
    # 3. Derivative w.r.t. MA (Theta)
    if (has_ma) {
      d_eta_d_theta[t, ] <- error[t - ma_lags]
      ma_effect <- crossprod(theta, d_eta_d_theta[t - ma_lags, , drop = FALSE])
      d_eta_d_theta[t, ] <- d_eta_d_theta[t, ] - as.numeric(ma_effect)
    }
    
    # 4. Derivative w.r.t. Beta (Regressors)
    if (has_xreg) {
      # Base: x_t - sum(varphi * x_{t-k})
      base_grad <- xreg[t, ]
      if (has_ar) {
        # Calculate sum(varphi_k * x_{t-k})
        # crossprod: (1 x p) * (p x k) -> (1 x k)
        ar_adjustment <- crossprod(
          varphi,
          xreg[t - ar_lags, , drop = FALSE]
        )
        base_grad <- base_grad - as.numeric(ar_adjustment)
      }
      
      d_eta_d_beta[t, ] <- base_grad
      
      if (has_ma) {
        ma_effect <- crossprod(theta, d_eta_d_beta[t - ma_lags, , drop = FALSE])
        d_eta_d_beta[t, ] <- d_eta_d_beta[t, ] - as.numeric(ma_effect)
      }
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
  ystar <- linkfun(y_effective)
  mustar <- digamma(mu_effective * phi) - digamma((1 - mu_effective) * phi)
  
  # Chain rule common term: dL/deta = dL/dmu * dmu/deta
  common_term <- mu_eta_val * (ystar - mustar)
  
  # Compute final scores (gradient = sum(dL/deta * deta/dparam))
  # Using crossprod for sum(vector * vector)
  
  idx <- (max_lag + 1):n_obs
  
  score_alpha <- as.numeric(
    phi * crossprod(d_eta_d_alpha[idx], common_term)
  )
  
  score_varphi <- if (has_ar) {
    as.numeric(phi * crossprod(d_eta_d_varphi[idx, , drop = FALSE], common_term))
  } else {
    numeric(0)
  }
  
  score_theta <- if (has_ma) {
    as.numeric(phi * crossprod(d_eta_d_theta[idx, , drop = FALSE], common_term))
  } else {
    numeric(0)
  }
  
  score_beta <- if (has_xreg) {
    as.numeric(phi * crossprod(d_eta_d_beta[idx, , drop = FALSE], common_term))
  } else {
    numeric(0)
  }
  
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
                   score_phi, 
                   score_beta)
  
  if (any(!is.finite(final_score))) {
    warning("Non-finite values in score vector; returning zeros")
    return(rep(0, length(final_score)))
  }
  
  return(as.numeric(final_score))
}