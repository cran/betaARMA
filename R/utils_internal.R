# Internal utility functions — not exported
# These helpers are shared across loglik_barma(), score_vector_barma(),
# and fim_barma() to avoid duplication of critical constants.

# .barma_ridge_lambda: data-driven L2 penalty scalar per Cribari-Neto,
# Costa and Fonseca (2025, BJPS, Section 2): lambda_n = 1/(n-a)^ridge_rate
.barma_ridge_lambda <- function(n_obs, max_lag, ridge_rate = 0.9) {
  1 / (n_obs - max_lag)^ridge_rate
}