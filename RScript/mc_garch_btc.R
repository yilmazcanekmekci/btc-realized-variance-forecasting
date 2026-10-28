# =============================================================================
# MC-GARCH -- Bitcoin hourly realized variance
#
# Realized-variance adaptation of the Engle-Sokalska multiplicative-component
# framework, applied to log RV:
#   RV_t = h_t * s_hod(t) * z_t ,   log z_t ~ AR(1)-GARCH(1,1), zero-mean
#   h_t   : mean hourly RV of the previous complete UTC calendar day (known at t)
#   s_hod : hour-of-day factor = exp( mean over estimation hours of log(RV/h) ).
#           The geometric mean centres log z at zero within every UTC hour of the
#           estimation sample, which is consistent with include.mean = FALSE.


library(data.table)
library(rugarch)

# Configuration 
CSV_PATH <- "C:/Users/YILMAZ/Desktop/BTCUSD_5min_2015_2025.csv"
OUTPUT_DIR <- "C:/Users/YILMAZ/Desktop"
stopifnot(dir.exists(OUTPUT_DIR), file.exists(CSV_PATH))

SAMPLE_START <- as.POSIXct("2016-12-01 15:00:00", tz = "UTC")   # first full hour
DATA_END <- as.POSIXct("2025-12-31 23:55:00", tz = "UTC")
TRAIN_END <- as.POSIXct("2021-12-31 23:59:59", tz = "UTC")
VAL_END <- as.POSIXct("2023-12-31 23:59:59", tz = "UTC")
TEST_START <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
WFV_N_SPLITS <- 5L
ARCH_LAGS <- 12L
LB_LAGS <- 24L

cat("MC-GARCH v6 -- Bitcoin hourly realized variance\n")

# [1] Load data, continuous returns, hourly RV 
cat("\n[1] Loading data and computing hourly RV\n")

dt <- fread(CSV_PATH, select = c("timestamp", "close"),
            colClasses = list(character = "timestamp"))
dt[, timestamp := as.POSIXct(timestamp, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")]
stopifnot(!anyNA(dt$timestamp))
setkey(dt, timestamp)
stopifnot(!anyDuplicated(dt$timestamp), all(diff(as.numeric(dt$timestamp)) == 300))
stopifnot(all(is.finite(dt$close)), all(dt$close > 0))
cat(sprintf(" Full file: %d bars, complete 5-min grid\n", nrow(dt)))

# Bar timestamp = bar OPEN (verified against raw trades). The return in row t is
# log C_t - log C_{t-1} = price(t+5min) - price(t), i.e. it spans [t, t+5min).
# Each return is assigned to the clock hour CONTAINING its interval: the RV
# labelled H:00 holds the twelve returns spanning H:00 - (H+1):00 UTC.
# Returns are computed over the WHOLE series before any window cut.
dt[, r := c(NA_real_, diff(log(close)))]
dt[, hour := as.POSIXct(floor(as.numeric(timestamp) / 3600) * 3600,
                        origin = "1970-01-01", tz = "UTC")]
dt <- dt[timestamp >= SAMPLE_START & timestamp <= DATA_END]

rv_hourly <- dt[, .(RV = sum(r^2), n_obs = sum(!is.na(r))), by = hour]
setkey(rv_hourly, hour)
stopifnot(all(rv_hourly$n_obs == 12L),
          all(diff(as.numeric(rv_hourly$hour)) == 3600),
          isTRUE(all.equal(sum(dt$r^2), sum(rv_hourly$RV))))
rv_hourly[, log_RV  := fifelse(RV > 0, log(RV), NA_real_)]
rv_hourly[, prev_RV := shift(RV, 1L)]
# Common evaluation rule for ALL models
rv_hourly[, eval_ok := RV > 0 & !is.na(prev_RV) & prev_RV > 0]
cat(sprintf("  Bars in window: %d | hourly obs: %d | zero-RV hours: %d | hours after a zero-RV hour: %d\n",
            nrow(dt), nrow(rv_hourly), sum(rv_hourly$RV == 0),
            sum(rv_hourly$RV > 0 & !is.na(rv_hourly$prev_RV) & rv_hourly$prev_RV == 0)))

# [2] Splits 
cat("\n[2] Splits\n")
rv_train <- rv_hourly[hour <= TRAIN_END]
rv_val <- rv_hourly[hour > TRAIN_END & hour <= VAL_END]
rv_test <- rv_hourly[hour > VAL_END]
rv_all <- rv_hourly[hour <= VAL_END]
stopifnot(nrow(rv_val) == 17520L, nrow(rv_test) == 17544L)
for (nm in c("Train", "Validation", "Test")) {
  s <- switch(nm, Train = rv_train, Validation = rv_val, Test = rv_test)
  cat(sprintf("  %-10s: %6d hours (%s -> %s) | evaluable: %d\n", nm, nrow(s),
              format(min(s$hour), "%Y-%m-%d %H:%M", tz = "UTC"),
              format(max(s$hour), "%Y-%m-%d %H:%M", tz = "UTC"), sum(s$eval_ok)))
}
QLIKE_FLOOR <- as.numeric(quantile(rv_train$RV[rv_train$RV > 0], 0.001))

# [3] Functions 
hod_of <- function(h) as.integer(format(h, "%H", tz = "UTC"))

# h_prev = mean RV of the previous COMPLETE calendar day (NA otherwise)
add_daily_level <- function(rv) {
  d <- copy(rv)
  d[, date := as.Date(hour, tz = "UTC")]
  daily <- d[, .(daily_RV = if (.N == 24L) mean(RV) else NA_real_), by = date]
  setkey(daily, date)
  stopifnot(all(diff(as.numeric(daily$date)) == 1))
  daily[, h_prev := shift(daily_RV, 1L)]
  d <- merge(d, daily[, .(date, h_prev)], by = "date", all.x = TRUE)
  setkey(d, hour)
  d[, hod := hod_of(hour)]
  d
}

# A fit is accepted only if it converged, parameters are finite and valid, and
# the GARCH(1,1) is covariance-stationary (alpha + beta < 1)
valid_fit <- function(fit) {
  if (is.null(fit) || fit@fit$convergence != 0L) return(FALSE)
  p <- coef(fit)
  all(is.finite(p)) && p["omega"] > 0 && p["alpha1"] >= 0 &&
    p["beta1"] >= 0 && p["alpha1"] + p["beta1"] < 1
}
fit_garch <- function(spec, x) {
  for (sv in c("hybrid", "nlminb", "solnp")) {
    out <- tryCatch(ugarchfit(spec, data = x, solver = sv, out.sample = 0),
                    error = function(e) NULL)
    if (valid_fit(out)) return(list(fit = out, solver = sv))
  }
  stop("All GARCH solvers failed or returned non-converged / invalid fits")
}

# Estimate s_hod (geometric) + AR(1)-GARCH(1,1) + homoskedastic AR(1) on a sample
estimate_mcgarch <- function(rv_data, final = FALSE) {
  d <- add_daily_level(rv_data)
  d <- d[!is.na(h_prev) & h_prev > 0 & RV > 0]
  d[, log_y := log(RV / h_prev)]
  s_i <- d[, .(s_i = exp(mean(log_y))), by = hod]
  setkey(s_i, hod)
  d <- merge(d, s_i, by = "hod", all.x = TRUE)
  setkey(d, hour)
  d[, log_z := log_y - log(s_i)]
  
  # Identification check: mean(log z) must be ~0
  m_lz <- mean(d$log_z)
  stopifnot(abs(m_lz) < 1e-8)
  cat(sprintf(" log z: mean=%.2e  median=%+.4f  sd=%.4f  (n=%d)\n",
              m_lz, median(d$log_z), sd(d$log_z), nrow(d)))
  
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0), include.mean = FALSE),
    distribution.model = "norm")
  ff <- fit_garch(spec, d$log_z)
  p <- coef(ff$fit)
  par <- list(phi = unname(p["ar1"]), omega = unname(p["omega"]),
              alpha = unname(p["alpha1"]), beta = unname(p["beta1"]))
  par$persistence <- par$alpha + par$beta
  stopifnot(par$persistence < 1)
  # Long-run unconditional variance of the AR residual process
  par$resid_lr_var <- par$omega / (1 - par$persistence)
  
  ar1 <- arima(d$log_z, order = c(1, 0, 0), include.mean = FALSE)
  par$phi_ar1 <- unname(coef(ar1)["ar1"])
  
  list(s_i = s_i, params = par, fit = ff$fit, solver = ff$solver,
       convergence = ff$fit@fit$convergence, n_obs_fit = nrow(d), log_z = d$log_z)
}

# One-step-ahead forecasts for hours in [from, to]. rv_hist = FULL history up to
# `to`. Hour t uses only: h_t (previous day), s_hod(t), log z_{t-1}.
mcgarch_forecast <- function(rv_hist, s_i, params, from, to) {
  d <- add_daily_level(rv_hist)
  d <- merge(d, s_i, by = "hod", all.x = TRUE)
  setkey(d, hour)
  d[, log_z := fifelse(RV > 0 & h_prev > 0, log(RV / (h_prev * s_i)), NA_real_)]
  n <- nrow(d)
  lz0 <- fifelse(is.na(d$log_z), 0, d$log_z)    # 0 = model's unconditional mean
  
  base <- log(d$h_prev) + log(d$s_i)
  pred  <- rep(NA_real_, n); pred_ar1 <- rep(NA_real_, n); sigma2 <- rep(NA_real_, n)
  pred[2:n] <- base[2:n] + params$phi     * lz0[1:(n - 1)]
  pred_ar1[2:n] <- base[2:n] + params$phi_ar1 * lz0[1:(n - 1)]
  
  s2 <- params$resid_lr_var
  for (j in 2:n) {
    eps_prev <- lz0[j - 1] - params$phi * (if (j - 1 == 1) 0 else lz0[j - 2])
    s2 <- params$omega + params$alpha * eps_prev^2 + params$beta * s2
    sigma2[j] <- s2
  }
  out <- data.table(datetime_utc = d$hour, log_RV_hat = pred,
                    log_RV_hat_ar1 = pred_ar1, sigma2_t = sigma2)
  out[datetime_utc >= from & datetime_utc <= to]
}

actual_eval <- rv_hourly[eval_ok == TRUE, .(datetime_utc = hour, log_RV)]

eval_forecasts <- function(preds, col = "log_RV_hat") {
  ev <- merge(preds, actual_eval, by = "datetime_utc")
  yhat <- ev[[col]]
  stopifnot(all(is.finite(yhat)))
  e  <- ev$log_RV - yhat
  ratio <- exp(ev$log_RV) / pmax(exp(yhat), QLIKE_FLOOR)
  ql <- mean(ratio - log(ratio) - 1)
  list(rmse = sqrt(mean(e^2)), mae = mean(abs(e)), bias = mean(e), qlike = ql, n = nrow(ev))
}

excl_report <- function(preds, label) {
  w <- merge(preds[, .(datetime_utc)],
             rv_hourly[, .(datetime_utc = hour, RV, prev_RV, eval_ok)],
             by = "datetime_utc", all.x = TRUE)
  n_zero  <- sum(w$RV == 0)
  n_after <- sum(w$RV > 0 & !is.na(w$prev_RV) & w$prev_RV == 0)
  cat(sprintf("    %s: forecasts=%d  evaluated=%d  excluded: zero-RV=%d, hour-after-zero=%d\n",
              label, nrow(w), sum(w$eval_ok), n_zero, n_after))
  c(n_forecast = nrow(w), n_eval = sum(w$eval_ok), n_zero = n_zero, n_after_zero = n_after)
}

# [4] Walk-forward validation (calendar-defined folds)
cat("\n[4] Walk-forward validation (", WFV_N_SPLITS, " folds, expanding window)\n", sep = "")
chunk <- nrow(rv_val) %/% WFV_N_SPLITS
stopifnot(chunk * WFV_N_SPLITS == nrow(rv_val), chunk == 3504L)
val_hours <- rv_val$hour

wfv_rows <- list(); wfv_fc <- list()
for (k in seq_len(WFV_N_SPLITS)) {
  f_start <- val_hours[(k - 1L) * chunk + 1L]
  f_end <- val_hours[k * chunk]
  tr_f  <- rv_hourly[hour < f_start]
  cat(sprintf("  Fold %d: %s -> %s (train hours: %d)\n", k,
              format(f_start, "%Y-%m-%d %H:%M", tz = "UTC"),
              format(f_end, "%Y-%m-%d %H:%M", tz = "UTC"), nrow(tr_f)))
  fit_k <- estimate_mcgarch(tr_f)
  preds_k <- mcgarch_forecast(rv_hourly[hour <= f_end], fit_k$s_i, fit_k$params, f_start, f_end)
  stopifnot(nrow(preds_k) == chunk, !anyNA(preds_k$log_RV_hat))
  ex <- excl_report(preds_k, sprintf("fold %d", k))
  g <- eval_forecasts(preds_k, "log_RV_hat"); a <- eval_forecasts(preds_k, "log_RV_hat_ar1")
  cat(sprintf("    GARCH-phi: RMSE=%.4f MAE=%.4f bias=%+.4f QLIKE=%.4f | AR1-ML: RMSE=%.4f QLIKE=%.4f | persistence=%.4f solver=%s conv=%d\n",
              g$rmse, g$mae, g$bias, g$qlike, a$rmse, a$qlike,
              fit_k$params$persistence, fit_k$solver, fit_k$convergence))
  wfv_fc[[k]] <- cbind(preds_k[, .(datetime_utc, log_RV_hat, log_RV_hat_ar1)], fold = k)
  wfv_rows[[k]] <- data.table(
    fold = k,
    from_utc = format(f_start, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    to_utc  = format(f_end,   "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    train_hrs = nrow(tr_f), n_forecast = ex["n_forecast"], n_eval = ex["n_eval"],
    n_excl_zero = ex["n_zero"], n_excl_after_zero = ex["n_after_zero"],
    rmse = g$rmse, mae = g$mae, bias = g$bias, qlike = g$qlike,
    rmse_ar1 = a$rmse, qlike_ar1 = a$qlike,
    phi = fit_k$params$phi, phi_ar1 = fit_k$params$phi_ar1,
    persistence = fit_k$params$persistence,
    convergence = fit_k$convergence, solver = fit_k$solver)
}
wfv_df <- rbindlist(wfv_rows); wfv_fc_df <- rbindlist(wfv_fc)

# [5] Final estimation on TRAIN + VAL 
cat("\n[5] Final estimation on TRAIN + VAL\n")
fit_final <- estimate_mcgarch(rv_all, final = TRUE)
P <- fit_final$params
cat(sprintf(" phi=%.6f  phi_ar1=%.6f  omega=%.8f  alpha=%.6f  beta=%.6f  persistence=%.6f\n",
            P$phi, P$phi_ar1, P$omega, P$alpha, P$beta, P$persistence))
cat(sprintf(" phi difference (GARCH - AR1): %+.6f\n", P$phi - P$phi_ar1))
cat(sprintf(" residual long-run variance: %.4f | n_obs_fit: %d | solver: %s | convergence: %d\n",
            P$resid_lr_var, fit_final$n_obs_fit, fit_final$solver, fit_final$convergence))
cat(" Hour-of-day factors s_hod (geometric):\n"); print(fit_final$s_i[order(hod)])

# Residual diagnostics (in-sample, final model)
zres <- as.numeric(residuals(fit_final$fit, standardize = TRUE))
arch_lm <- function(x, lags) {
  x2 <- x^2; X <- embed(x2, lags + 1)
  r2 <- summary(lm(X[, 1] ~ X[, -1]))$r.squared
  stat <- nrow(X) * r2
  c(stat = stat, p = pchisq(stat, lags, lower.tail = FALSE))
}
lb1 <- Box.test(zres,   lag = LB_LAGS, type = "Ljung-Box", fitdf = 1)
lb2 <- Box.test(zres^2, lag = LB_LAGS, type = "Ljung-Box", fitdf = 2)
al <- arch_lm(zres, ARCH_LAGS)
ac1 <- acf(zres, lag.max = 1, plot = FALSE)$acf[2]
ac2 <- acf(zres^2, lag.max = 1, plot = FALSE)$acf[2]
cat("\n  Residual diagnostics (standardized, in-sample; large n => tests reject for tiny dependence, read effect sizes):\n")
cat(sprintf("    mean=%.4f  var=%.4f | lag-1 autocorr: resid=%+.4f, resid^2=%+.4f\n",
            mean(zres), var(zres), ac1, ac2))
cat(sprintf(" Ljung-Box(%d) resid: Q=%.1f p=%.3g | resid^2: Q=%.1f p=%.3g | ARCH-LM(%d): stat=%.1f p=%.3g\n",
            LB_LAGS, lb1$statistic, lb1$p.value, lb2$statistic, lb2$p.value,
            ARCH_LAGS, al["stat"], al["p"]))

# [6] Test forecasts 
cat("\n[6] Test forecasts (2024-2025)\n")
test_end <- max(rv_test$hour)
tp <- mcgarch_forecast(rv_hourly, fit_final$s_i, P, TEST_START, test_end)
stopifnot(nrow(tp) == 17544L, !anyNA(tp$log_RV_hat), !anyNA(tp$log_RV_hat_ar1), !anyNA(tp$sigma2_t))
ex_t <- excl_report(tp, "test")
g <- eval_forecasts(tp, "log_RV_hat"); a <- eval_forecasts(tp, "log_RV_hat_ar1")
cat(sprintf(" GARCH-phi : RMSE=%.4f MAE=%.4f bias=%+.4f QLIKE=%.4f\n", g$rmse, g$mae, g$bias, g$qlike))
cat(sprintf(" AR1-ML    : RMSE=%.4f MAE=%.4f bias=%+.4f QLIKE=%.4f\n", a$rmse, a$mae, a$bias, a$qlike))

# No-look-ahead check 1: truncating the future must not change past forecasts
cut <- as.POSIXct("2024-06-30 23:00:00", tz = "UTC")
tp_cut <- mcgarch_forecast(rv_hourly[hour <= cut], fit_final$s_i, P, TEST_START, cut)
tp_ref <- tp[datetime_utc <= cut]
stopifnot(nrow(tp_cut) == nrow(tp_ref),
          isTRUE(all.equal(tp_cut$log_RV_hat, tp_ref$log_RV_hat)),
          isTRUE(all.equal(tp_cut$sigma2_t, tp_ref$sigma2_t)))
# No-look-ahead check 2: perturbing RV at t0 must not change forecasts up to t0,
# but must change the forecast for t0+1
t0 <- as.POSIXct("2024-03-15 12:00:00", tz = "UTC")
rv2 <- copy(rv_hourly); rv2[hour == t0, RV := RV * 5]
tp2 <- mcgarch_forecast(rv2, fit_final$s_i, P, TEST_START, test_end)
stopifnot(isTRUE(all.equal(tp2[datetime_utc <= t0]$log_RV_hat, tp[datetime_utc <= t0]$log_RV_hat)),
          abs(tp2[datetime_utc == t0 + 3600]$log_RV_hat - tp[datetime_utc == t0 + 3600]$log_RV_hat) > 1e-6)
cat(" No-look-ahead checks passed (truncation + perturbation)\n")

# [7] Write outputs 
cat("\n[7] Writing CSV outputs to:", OUTPUT_DIR, "\n")
fwrite(tp, file.path(OUTPUT_DIR, "mcgarch_forecasts.csv"), dateTimeAs = "write.csv")
fwrite(wfv_fc_df[, .(datetime_utc, fold, log_RV_hat, log_RV_hat_ar1)],
       file.path(OUTPUT_DIR, "mcgarch_wfv_forecasts.csv"), dateTimeAs = "write.csv")
fwrite(wfv_df, file.path(OUTPUT_DIR, "mcgarch_wfv.csv"))
fwrite(rv_hourly[, .(datetime_utc = hour, RV, log_RV, n_obs, prev_RV, eval_ok)],
       file.path(OUTPUT_DIR, "rv_hourly_check.csv"), dateTimeAs = "write.csv")
fwrite(fit_final$s_i[order(hod)], file.path(OUTPUT_DIR, "mcgarch_diurnal.csv"))
fwrite(data.table(
  phi = P$phi, phi_ar1 = P$phi_ar1, omega = P$omega, alpha = P$alpha, beta = P$beta,
  persistence = P$persistence, resid_long_run_variance = P$resid_lr_var,
  mean_log_z_fit = mean(fit_final$log_z), n_obs_fit = fit_final$n_obs_fit,
  convergence = fit_final$convergence, solver = fit_final$solver,
  resid_mean = mean(zres), resid_var = var(zres),
  resid_acf1 = ac1, resid2_acf1 = ac2,
  lb_resid_Q = unname(lb1$statistic), lb_resid_p = lb1$p.value,
  lb_resid2_Q = unname(lb2$statistic), lb_resid2_p = lb2$p.value,
  archlm_stat = unname(al["stat"]), archlm_p = unname(al["p"]),
  sample_start = format(SAMPLE_START, tz = "UTC"), data_end = format(DATA_END, tz = "UTC"),
  fit_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  file.path(OUTPUT_DIR, "mcgarch_metadata.csv"))
cat("DONE (6 files: forecasts, wfv_forecasts, wfv, rv_hourly_check, diurnal, metadata).\n")