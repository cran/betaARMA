## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 3.5,
  fig.align = "center",
  warning   = FALSE,
  message   = FALSE
)

## ----libraries----------------------------------------------------------------
library(betaARMA)
library(forecast)
library(ggplot2)
library(gridExtra)
library(zoo)
library(scales)
library(moments)

## ----config_plot, echo=FALSE--------------------------------------------------
# Shared ggplot configuration
sample_size     <- length(brasilia_ts)
end_time_series <- time(brasilia_ts)[sample_size] + 1

ggplot_size_font <- 9

ggplot_theme <- theme_minimal() + 
  theme(
  legend.position = "bottom",
  title        = element_text(size = ggplot_size_font),
  axis.text    = element_text(size = ggplot_size_font),
  axis.title   = element_text(size = ggplot_size_font),
  legend.text  = element_text(size = ggplot_size_font),
  legend.title = element_text(size = ggplot_size_font)
)

ggplot_scale_y <- scale_y_continuous(
  breaks = seq(0.0, 1.00, 0.10)
)

ggplot_scale_x <- scale_x_continuous(
  breaks = seq(1999, end_time_series, 2),
  limits = c(1999, end_time_series)
)

## ----subsample----------------------------------------------------------------
y_subsample_ts <- window(
  brasilia_ts,
  start = c(2006, 02),
  end   = c(2019, 07)
)

y_train <- window(
  y_subsample_ts,
  end = c(2018, 07)
)

y_test <- window(
  y_subsample_ts,
  start = c(2018, 08)
)

## ----ggplot_subsample, echo=FALSE---------------------------------------------
len_y_ts        <- length(y_train)
start_subsample <- time(y_train)[1]
end_subsample   <- time(y_train)[len_y_ts]

ggplot_subsample <- ggplot(brasilia_df, aes(time, y)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = 1)) +
  geom_vline(xintercept = start_subsample,
             linetype = "dashed", linewidth = 0.5, color = "red") +
  geom_vline(xintercept = end_subsample,
             linetype = "dashed", linewidth = 0.5, color = "red") +
  labs(x = "Year", y = "Relative humidity (proportion)") +
  ggplot_scale_x +
  ggplot_scale_y +
  ggplot_theme

## ----print_ggplot_subsample, echo=FALSE, fig.cap="Relative humidity in Brasília. Dashed red lines delimit the training period used in the instability analysis."----
print(ggplot_subsample)

## ----descriptive_subsample, echo=FALSE----------------------------------------
descriptive_sub_df <- data.frame(
  Min             = min(y_subsample_ts),
  Max             = max(y_subsample_ts),
  Median          = median(y_subsample_ts),
  Mean            = mean(y_subsample_ts),
  SD              = sd(y_subsample_ts),
  Skewness        = moments::skewness(y_subsample_ts),
  Excess_kurtosis = moments::kurtosis(y_subsample_ts)
)
descriptive_sub_df <- round(descriptive_sub_df, 2)

## ----print_descriptive_subsample, echo=FALSE----------------------------------
knitr::kable(
  descriptive_sub_df,
  caption = "Descriptive statistics of the monthly relative humidity in
  Brasília — subsample (February 2006 to July 2019)."
)

## ----regressors---------------------------------------------------------------
freq              <- frequency(y_train)
n_test            <- length(y_test)
trend_index_train <- seq_along(y_train)
trend_index_test  <- (max(trend_index_train) + 1):
                     (max(trend_index_train) + n_test)

hs_train <- sin(2 * pi * trend_index_train / freq)
hc_train <- cos(2 * pi * trend_index_train / freq)
hs_test  <- sin(2 * pi * trend_index_test  / freq)
hc_test  <- cos(2 * pi * trend_index_test  / freq)

X_train <- cbind(hs = hs_train, hc = hc_train)
X_test  <- cbind(hs = hs_test,  hc = hc_test)

## ----fit_cmle-----------------------------------------------------------------
ar_lags <- c(10, 18)
ma_lags <- c(1, 11, 13)

fit_cmle <- barma(
  y       = y_train,
  ar      = ar_lags,
  ma      = ma_lags,
  penalty = FALSE,
  xreg    = X_train
)

## ----summary_cmle-------------------------------------------------------------
summary(fit_cmle)

## ----fit_pcmle----------------------------------------------------------------
fit_pcmle <- barma(
  y       = y_train,
  ar      = ar_lags,
  ma      = ma_lags,
  penalty = TRUE,
  xreg    = X_train
)

## ----summary_pcmle------------------------------------------------------------
summary(fit_pcmle)

## ----diagnostics_pcmle, fig.height=7.0, fig.cap="Diagnostic plots for the PCMLE fit: AR=(10,18), MA=(1,11,13)."----
plot(fit_pcmle, which = "all")

## ----fit_cmle_reduced---------------------------------------------------------
ar_lags_reduced <- c(10, 18)
ma_lags_reduced <- c(1, 13)

fit_cmle_reduced <- barma(
  y       = y_train,
  ar      = ar_lags_reduced,
  ma      = ma_lags_reduced,
  penalty = FALSE,
  xreg    = X_train
)

## ----summary_cmle_reduced-----------------------------------------------------
summary(fit_cmle_reduced)

## ----diagnostics_all, fig.height=7.0, fig.cap="Diagnostic plots for the reduced CMLE model AR=(10,18), MA=(1,13)."----
plot(fit_cmle_reduced, which = "all")

## ----diagnostics_qq, fig.height=4.0, fig.width=4.5, fig.cap="Normal Q-Q plot of quantile residuals for the reduced CMLE model."----
plot(fit_cmle_reduced, which = "qq")

## ----forecast_setup-----------------------------------------------------------
# Define the translation dictionary
month_map <- c(
  "Jan" = "Jan", "Fev" = "Feb", "Mar" = "Mar", "Abr" = "Apr",
  "Mai" = "May", "Jun" = "Jun", "Jul" = "Jul", "Ago" = "Aug",
  "Set" = "Sep", "Out" = "Oct", "Nov" = "Nov", "Dez" = "Dec"
)

dates     <- as.Date(zoo::as.yearmon(time(y_test)))
pt_months <- tools::toTitleCase(format(dates, "%b"))
months    <- unname(month_map[pt_months])

forecasts_df <- data.frame(
  observed = as.numeric(y_test),
  forecast = forecast(fit_cmle_reduced, xreg = X_test, h = n_test)$mean,
  month    = months
)
forecasts_df$time <- seq_len(n_test)

## ----ggplot_forecast, echo=FALSE----------------------------------------------
ggplot_forecast <- ggplot(forecasts_df, aes(x = time)) +
  geom_line(aes(y = forecast, colour = "Forecast"), linewidth = 0.8) +
  geom_point(aes(y = observed, shape = "Observed"),
             color = "#00BFC4", size = 3) +
  scale_colour_manual(name   = " ",
                      values = c("Forecast" = "#C77CFF")) +
  scale_shape_manual(name   = " ",
                     values = c("Observed" = 16)) +
  scale_x_continuous(breaks = seq_len(n_test),
                     limits = c(1, n_test),
                     labels = forecasts_df$month) +
  scale_y_continuous(breaks = seq(0.3, 1.0, 0.1),
                     limits = c(0.3, 1.0),
                     labels = number_format(accuracy = 0.1)) +
  labs(x = "Time (months)", y = "Relative humidity (proportion)") +
  ggplot_theme

## ----print_forecast, echo=FALSE, fig.cap="Out-of-sample forecasts from the reduced CMLE model against the observed test data (teal dots)."----
print(ggplot_forecast)

## ----forecast_accuracy--------------------------------------------------------
mae  <- function(obs, pred) mean(abs(obs - pred))
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

obs <- forecasts_df$observed

accuracy_table <- data.frame(
  Model = "betaARMA CMLE (reduced)",
  MAE   = round(mae(obs,  forecasts_df$forecast) * 100, 2),
  RMSE  = round(rmse(obs, forecasts_df$forecast) * 100, 2)
)

knitr::kable(
  accuracy_table,
  caption = paste0(
    "MAE and RMSE ($\\times 100$) over the ",
    n_test,
    "-month forecast horizon."
  )
)

## ----session_info, echo=FALSE-------------------------------------------------
cat("=================================================================\n")
cat("Session Information\n")
cat("=================================================================\n")
cat("This report was generated at:",
    format(Sys.time(), "%B %d, %Y at %I:%M %p"), "\n\n")
print(sessionInfo())

