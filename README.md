# Bitcoin Realized Variance Forecasting

Seven models are compared out of sample for one-hour-ahead forecasts of Bitcoin's hourly realized variance: Random Walk, Historical Mean, HAR-Classic, HAR-Extended, HAR-Calendar, MC-GARCH and an LSTM network.

## Repository structure

```
r_script/
  final_mc_garch_btc.R                               MC-GARCH estimation and forecasts (R)
notebooks/
  retrivering_btc_data.ipynb                         Downloads the 5-minute BTC/USD data
  final_main_canonical_training_RETRAIN_TRUE.ipynb   Full analysis, trains the LSTM (GPU)
  final_main_reporting_RETRAIN_FALSE.ipynb           Same analysis, loads the saved LSTM forecasts
data/
  Raw/                                               5-minute BTC/USD data (1 CSV file)
  mcgarch/                                           Output files of the R script (6 CSV files)
  lstm/                                              Saved LSTM forecasts (2 CSV files)
results/                                             Tables and forecast series behind the thesis results
```

## Data and setup

| | |
|---|---|
| Source | BTC/USD, Bitstamp, 5-minute OHLCV (via `ccxt`) |
| Study window | 2016-12-01 15:00 to 2025-12-31 23:55 UTC (79,617 hourly observations) |
| Target | log of hourly realized variance (sum of twelve squared 5-minute returns) |
| Horizon | one hour ahead |
| Split | Train 2016-2021, Validation 2022-2023 (5 walk-forward folds), Test 2024-2025 |
| Evaluation hours | 17,514 (validation) and 17,539 (test), identical for every model |
| Metrics | RMSE, MAE, Bias, R², QLIKE |
| Model comparison | Diebold-Mariano tests (Newey-West variance), squared log error and QLIKE, Holm correction over 21 pairs |

## Models

| Model | Description |
|---|---|
| Random Walk | previous hour's log RV |
| Historical Mean | mean log RV of the estimation sample |
| HAR-Classic | Corsi (2009) log-linear HAR: hourly, daily and weekly lags of log RV (OLS, HAC errors) |
| HAR-Extended | HAR-Classic plus daily jump and downside-semivariance terms |
| HAR-Calendar | HAR-Extended plus sine/cosine terms for hour of day and day of week (the same nine inputs as the LSTM, used linearly) |
| MC-GARCH | Multiplicative-component model adapted to log RV (R): daily level, hour-of-day factor, AR(1)-GARCH(1,1) |
| LSTM | Two-layer LSTM on 168-hour windows of the nine inputs, three seeds averaged |

## Requirements

- **Python** (notebooks, developed and run on Kaggle, GPU T4 x2): see `requirements.txt` (TensorFlow 2.21.0).
- **R** (MC-GARCH): packages `data.table` and `rugarch`; the script prints the versions it ran with. `rstudioapi` is optional and only needed to run the script from inside RStudio.

## How to run

1. **Download the data**: run `notebooks/retrivering_btc_data.ipynb`, or use the file in `data/Raw/`.
2. **MC-GARCH**: run `r_script/final_mc_garch_btc.R` from the repository root (`Rscript r_script/final_mc_garch_btc.R`) or open it in RStudio. Paths are resolved relative to the repository: it reads `data/Raw/BTCUSD_5min_2015_2025.csv` and writes the six CSV files to `data/mcgarch/`.
3. **Kaggle dataset**: upload the 5-minute CSV and the six MC-GARCH CSV files together to one Kaggle dataset. The notebooks find the files automatically.
4. **Train the LSTM**: run `final_main_canonical_training_RETRAIN_TRUE.ipynb` on Kaggle with a GPU. It trains the network, saves the two LSTM forecast files, computes SHAP importances and writes all result tables.
5. **Reproduce the results**: add the two files from `data/lstm/` to the Kaggle dataset and run `final_main_reporting_RETRAIN_FALSE.ipynb`. It recomputes all metrics, tables and tests from the saved LSTM forecasts in a few minutes. SHAP is only computed in the training run.

The two notebooks are two run modes of the same analysis (the `RETRAIN_LSTM` setting), not two different experiments. LSTM training on a GPU is not guaranteed to give identical numbers on other hardware or library versions, which is why the final LSTM forecasts are saved in `data/lstm/`. Both notebooks produced identical forecast metrics.

## Results files

`results/` holds the output tables of the final run:

| File | Content |
|---|---|
| `table_metrics_validation.csv`, `table_metrics_test.csv` | RMSE, MAE, Bias, R², QLIKE per model |
| `table_dm_validation.csv`, `table_dm_test.csv` | pairwise Diebold-Mariano tests (lags 12, 24, 168; Holm-adjusted p-values) |
| `table_dm_mcgarch_vs_ar1.csv` | diagnostic test MC-GARCH vs the homoskedastic MC-AR(1) variant |
| `table_mz_hac.csv` | Mincer-Zarnowitz calibration regressions with HAC Wald test |
| `table_wfv_folds.csv` | metrics per walk-forward fold |
| `table_descriptive_stats.csv`, `table_vif.csv` | data description and multicollinearity check |
| `all_forecasts_validation.csv`, `all_forecasts_test.csv` | forecast series of every model |
| `run_manifest.json` | settings of the run |

## Main results

Forecast errors on log RV (lower is better). Validation = 2022-2023 walk-forward folds, Test = 2024-2025.

| Model | Val RMSE | Val QLIKE | Test RMSE | Test QLIKE |
|---|---|---|---|---|
| Random Walk | 1.0565 | 1.0393 | 0.9122 | 0.6573 |
| Historical Mean | 1.7185 | 1.0556 | 1.4054 | 0.8038 |
| HAR-Classic | 0.9327 | 0.7524 | 0.8300 | 0.5292 |
| HAR-Extended | 0.9339 | 0.7556 | 0.8305 | 0.5173 |
| HAR-Calendar | 0.9180 | 0.7283 | 0.8099 | 0.4933 |
| MC-GARCH | 0.9542 | 0.8446 | 0.8541 | 0.5597 |
| LSTM | 0.9426 | 0.7107 | **0.7841** | **0.4654** |

- **Test period:** the LSTM is the most accurate model under both losses, including against HAR-Calendar (Diebold-Mariano, Holm-adjusted).
- **Validation period:** against HAR-Classic and HAR-Extended the LSTM is better under QLIKE and not significantly different under squared log error. HAR-Calendar is more accurate than the LSTM under squared log error and not significantly different under QLIKE.
- **Calendar information:** adding it to the HAR model improves out-of-sample accuracy in both periods, so part of the LSTM's advantage over HAR-Extended comes from calendar information rather than from its nonlinear structure.
- **Calibration:** the Mincer-Zarnowitz hypothesis (intercept 0, slope 1) is rejected for every model with a non-degenerate forecast; with about 17,500 observations the deviations are small for the HAR models and the LSTM on the test period (slopes about 1.02-1.05).
- The results apply to one asset, one exchange and a one-hour horizon. The validation walk-forward forecasts and the test forecasts come from models trained on different amounts of data.
