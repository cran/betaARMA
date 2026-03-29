#' Forecast a barma Model
#'
#' @description
#' S3 method for producing forecasts from a fitted `barma` model object.
#'
#' @details
#' This function computes dynamic, multi-step-ahead point forecasts.
#' It implements a "Regression with ARMA errors" logic, where the AR components
#' are applied to the deviations from the regression line:
#' \eqn{AR(g(y_{t-k}) - x_{t-k}^\top \beta)}.
#'
#' @param object A fitted model object of class `barma`. Must contain
#'   `object$xreg` if regressors were used.
#' @param h The number of steps to forecast ahead (forecast horizon).
#'   Default is 6.
#' @param xreg A matrix of future regressor values for the forecast horizon.
#'   Should have h rows and the same number of columns as xreg used in model
#'   fitting.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A `ts` object containing the point forecasts for h steps ahead.
#'
#' @importFrom forecast forecast
#' @export
#' @method forecast barma
forecast.barma <- function(object, h = 6, xreg = NULL, ...) {
  
  # --------------------------------------------------------------------------
  # 1. Extract Model Components
  # --------------------------------------------------------------------------
  
  # --- Get parameters ---
  alpha <- object$alpha
  varphi <- object$varphi
  theta <- object$theta
  beta <- object$beta # Regression coefficients
  
  # --- Get link function ---
  link_structure <- make_link_structure(object$link)
  link_inv <- link_structure$linkinv
  link_fun <- link_structure$linkfun
  
  # --- Get lags ---
  ar_lags <- object$ar_lags
  ma_lags <- object$ma_lags
  has_ar <- length(ar_lags) > 0
  has_ma <- length(ma_lags) > 0
  has_regressors <- !is.null(beta) && length(beta) > 0
  
  # --- Validate regressors ---
  if (has_regressors) {
    if (is.null(xreg)) {
      stop(
        "Model was fit with regressors, but `xreg` is NULL. ",
        "Please provide `xreg` for forecasting."
      )
    }
    # Ensure object has the training xreg to calculate historical deviations
    if (is.null(object$xreg)) {
      stop(
        "The fitted model object does not contain `xreg`. ",
        "Historical regressors are required for AR terms in this specification."
      )
    }
    
    # Check dimensions
    if (nrow(xreg) != h) {
      stop("`xreg` must have ", h, " rows, but has ", nrow(xreg), " rows.")
    }
    if (ncol(xreg) != length(beta)) {
      stop(
        "`xreg` must have ", length(beta), " columns, but has ",
        ncol(xreg), " columns."
      )
    }
  } else {
    if (!is.null(xreg)) {
      warning(
        "`xreg` provided but model has no regressors. ",
        "`xreg` will be ignored."
      )
      xreg <- NULL
    }
  }
  
  # --------------------------------------------------------------------------
  # 2. Prepare Data Matrices
  # --------------------------------------------------------------------------
  
  # --- Get historical data ---
  y <- object$y
  n_obs <- object$n_obs
  
  # Get g(y)
  y_new <- link_fun(y)
  
  # Get historical errors (predictor scale)
  error_hat <- object$errorhat
  
  # --- Setup Padded Vectors ---
  # y_new_padded holds g(y_t) for t <= n_obs and eta_t for t > n_obs
  y_new_padded <- c(y_new, rep(NA_real_, h))
  error_padded <- c(error_hat, rep(NA_real_, h))
  
  # --- Setup Regressor Matrix (History + Future) ---
  if (has_regressors) {
    xreg_train <- as.matrix(object$xreg)
    xreg_test <- as.matrix(xreg)
    
    # Combine training and testing regressors into one large matrix.
    # This allows easy indexing for [t - ar_lags] regardless of whether the
    # lag falls in the past or future.
    xreg_combined <- rbind(xreg_train, xreg_test)
  }
  
  # Vector to store final forecasts
  forecast_values <- rep(NA_real_, h)
  
  # --------------------------------------------------------------------------
  # 3. Iterate and Compute Forecasts
  # --------------------------------------------------------------------------
  
  for (i in seq_len(h)) {
    t <- n_obs + i # The time step we are forecasting
    
    # --- 1. Base Level (Intercept) ---
    eta_forecast <- alpha
    
    # --- 2. Add Current Regression Component ---
    if (has_regressors) {
      # Add X_t * beta
      eta_forecast <- eta_forecast +
        as.numeric(crossprod(beta, xreg_combined[t, ]))
    }
    
    # --- 3. Add AR Component (on Deviations) ---
    if (has_ar) {
      # Calculate indices for the lags
      lag_indices <- t - ar_lags
      
      # Retrieve lagged values of Y (on link scale)
      y_lags <- y_new_padded[lag_indices]
      
      if (has_regressors) {
        # Retrieve lagged values of X * beta
        # We use xreg_combined to handle lags safely
        x_lags_beta <- xreg_combined[lag_indices, , drop = FALSE] %*% beta
        
        # AR term applies to the residual: (Y_{t-k} - X_{t-k} * beta)
        ar_term <- as.numeric(crossprod(varphi, y_lags - x_lags_beta))
      } else {
        # Standard AR if no regressors
        ar_term <- as.numeric(crossprod(varphi, y_lags))
      }
      
      eta_forecast <- eta_forecast + ar_term
    }
    
    # --- 4. Add MA Component ---
    if (has_ma) {
      # MA term applies to errors
      eta_forecast <- eta_forecast +
        as.numeric(crossprod(theta, error_padded[t - ma_lags]))
    }
    
    # --- Store Forecast ---
    mu_forecast <- link_inv(eta_forecast)
    forecast_values[i] <- mu_forecast
    
    # --- Update padded vectors for next iteration ---
    # The expected value of g(y_t) is the forecasted eta_t
    y_new_padded[t] <- eta_forecast
    
    # The expected value of future errors is 0
    error_padded[t] <- 0
  }
  
  # --------------------------------------------------------------------------
  # 4. Format and Return as ts Object
  # --------------------------------------------------------------------------
  
  y_ts <- stats::ts(y)
  ts_start <- stats::tsp(y_ts)[2] + (1 / stats::frequency(y_ts))
  ts_freq <- stats::frequency(y_ts)
  
  forecast_ts <- stats::ts(
    forecast_values,
    start = ts_start,
    frequency = ts_freq
  )
  
  return(forecast_ts)
}