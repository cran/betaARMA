#' Summarize a barma Model Fit
#'
#' (Full documentation block...)
#'
#' @export
#' @method summary barma
#' @param object A fitted model object of class `"barma"`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A list object of class `"summary.barma"`...
summary.barma <- function(object, ...) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Calculate vcov and model table ---
  # ------------------------------------------------------------------------- #
  
  # Inverse of Fisher information matrix
  vcov <- try(solve(object$fisher_info_mat, tol = 1e-20), silent = TRUE)
  
  # Check if the inverse matrix is possible
  if (inherits(vcov, "try-error")) {
    
    warning("FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
    
    # Create a coefficient table with NAs for std. errors
    n_params <- length(object$coef)
    model_table <- cbind(
      Estimate = round(object$coef, 4),
      `Std. Error` = rep(NA_real_, n_params),
      `z value` = rep(NA_real_, n_params),
      `Pr(>|z|)` = rep(NA_real_, n_params)
    )
    
  } else {
    
    # Calculate Std. Errors and p-values
    stderror <- sqrt(diag(vcov))
    z_zstat <- abs(object$coef / stderror)
    z_pvalues <- 2 * (1 - pnorm(z_zstat))
    
    model_table <- cbind(
      "Estimate" = round(object$coef, 4),
      "Std. Error" = round(stderror, 4),
      "z value" = round(z_zstat, 4),
      "Pr(>|z|)" = round(z_pvalues, 4)
    )
  }
  
  # ------------------------------------------------------------------------- #
  # --- 2. Calculate Information Criteria ---
  # ------------------------------------------------------------------------- #
  aux_info1 <- -2 * object$loglik
  aux_info2 <- object$n_params 
  
  # Use effective number of observations
  # log_n <- log(object$n_obs - object$max_lag) 
  
  # Use the number of observations
  log_n <- log(object$n_obs)
  
  aic <- aux_info1 + 2 * aux_info2
  bic <- aux_info1 + log_n * aux_info2
  hq <- aux_info1 + log(log_n) * 2 * aux_info2
  
  # ------------------------------------------------------------------------- #
  # --- 3. Assemble the summary object ---
  # ------------------------------------------------------------------------- #
  summary_list <- list(
    call = object$call,
    link = object$link,
    coefficients = model_table,
    vcov = vcov,
    loglik = object$loglik,
    aic = aic,
    bic = bic,
    hq = hq,
    conv = object$conv
  )
  
  # Assign the new class
  class(summary_list) <- "summary.barma"
  
  return(summary_list)
}