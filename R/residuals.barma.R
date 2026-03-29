#' Calculate Residuals for a barma Model Object
#'
#' @description
#' Computes various types of residuals for a fitted Beta Autoregressive Moving 
#' Average (BARMA) model object of class `"barma"`.
#'
#' @details
#' This function is an S3 method for the generic `residuals` function, 
#' tailored for objects returned by the `barma()` function. 
#' 
#' By default (`type = "standardized"`), it calculates standardized residuals 
#' on the predictor scale, defined as:
#' \deqn{r_t = \frac{g(y_t) - \hat{\eta}_t}{\sqrt{Var(g(y_t) - \hat{\eta}_t)}}}
#' where \eqn{Var(g(y_t) - \hat{\eta}_t) \approx (g'(\hat{\mu}_t))^2 \cdot}
#' \deqn{\frac{\hat{\mu}_t(1-\hat{\mu}_t)}{1+\hat{\phi}}}
#'. These residuals are useful for diagnostic checking, 
#' as they are approximately standard normal if the model is correctly 
#' specified.
#'
#' @param object A fitted model object of class `"barma"`, typically the 
#'   result of a call to `barma()`.
#' @param type The type of residuals to compute. Currently, only 
#'   `"standardized"` (default) is implemented, returning standardized 
#'   residuals on the predictor scale.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector or `ts` object containing the requested residuals. 
#'   For standardized residuals, the first `max_lag` values will be `NA`.
#'
#' @export
#' @method residuals barma
residuals.barma <- function(object, type = "standardized", ...) {
  
  if (type != "standardized") {
    warning("Only type = 'standardized' is currently implemented. ",
            "Returning standardized residuals.")
    type <- "standardized" 
  }
  
  if (type == "standardized") {
    
    # --- Extract necessary components directly from the object ---
    y <- object$y 
    if (is.null(y)) {
      stop("Original time series 'y' not found in the 'barma' object.")
    }
    
    etahat <- object$etahat   # Full NA-padded vector
    fitted_mus <- object$fitted # Full NA-padded ts object
    phi <- object$phi
    link <- object$link
    max_lag <- object$max_lag # Use stored max_lag
    
    if (is.null(max_lag)) {
      stop("'max_lag' not found in the 'barma' object.")
    }
    
    # Get link function derivatives
    link_structure <- make_link_structure(link)
    linkfun <- link_structure$linkfun
    mu.eta_fun <- link_structure$mu.eta
    
    # --- Calculate residuals only for effective observations ---
    n_obs <- length(y)
    if (n_obs <= max_lag) {
      # This warning should now only appear if truly applicable
      warning("Insufficient observations relative to max_lag.")
      return(rep(NA_real_, n_obs))
    }
    
    idx_effective <- (max_lag + 1):n_obs
    
    y_effective <- y[idx_effective]
    etahat_effective <- etahat[idx_effective]
    # Ensure fitted_mus is treated as a simple vector if needed
    muhat_effective <- as.numeric(fitted_mus[idx_effective]) 
    
    # g(y_t)
    ynew_effective <- linkfun(y_effective)
    
    # g(y_t) - eta_hat_t
    raw_resids_predictor <- ynew_effective - etahat_effective
    
    # Calculate standard deviation component
    Vmu_hat <- muhat_effective * (1 - muhat_effective)
    g_prime_sq_inv <- mu.eta_fun(eta = etahat_effective)^2
    g_prime_sq <- 1 / g_prime_sq_inv
    var_resids_predictor <- g_prime_sq * Vmu_hat / (1 + phi)
    sd_resids_predictor <- sqrt(var_resids_predictor)
    
    # Standardized residuals
    standardized_resids <- raw_resids_predictor / sd_resids_predictor
    
    # --- Pad with NAs and return as ts object ---
    final_resids <- c(rep(NA_real_, max_lag), standardized_resids)
    
    if (stats::is.ts(y)) {
      final_resids <- stats::ts(final_resids, 
                                start = stats::start(y), 
                                frequency = stats::frequency(y))
    }
    
    return(final_resids)
    
  } else {
    stop("Residual type '", type, "' not implemented.")
  }
}