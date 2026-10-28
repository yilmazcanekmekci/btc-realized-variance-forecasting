# Bitcoin Realized Variance Forecasting

Six models are compared out of sample: Random Walk, Historical Mean, HAR-Classic, HAR-Extended, MC-GARCH and an LSTM network.

## Repository structure

```
RScript/
  mc_garch_btc.R                               MC-GARCH estimation and forecasts (R)
notebooks/
  retrivering_btc_data.ipynb                   Downloads the 5-minute BTC/USD data
  main_canonical_training_RETRAIN_TRUE.ipynb   Full analysis, trains the LSTM (GPU)
  main_reporting_RETRAIN_FALSE.ipynb           Same analysis, loads the saved LSTM forecasts
data/
  Raw/                                         5-minute BTC/USD data
  mcgarch/                                     Output files of the R script (6 CSV files)
  lstm/                                        Saved LSTM forecasts (2 CSV files)
```

## Data and setup

| | |
|---|---|
| Source | BTC/USD, Bitstamp, 5-minute OHLCV (via `ccxt`) |
| Study window | 2016-12-01 15:00 to 2025-12-31 23:55 UTC (79,617 hourly observations) |
| Target | log of hourly realized variance (sum of twelve squared 5-minute returns) |
| Horizon | one hour ahead |
| Split | Train 2016-2021, Validation 2022-2023 (5 walk-forward folds), Test 2024-2025 |
| Metrics | RMSE, MAE, Bias, R², QLIKE |
| Model comparison | Diebold-Mariano tests (Newey-West variance), squared log error and QLIKE, Holm correction |

## How to run

1. **Download the data**: run `notebooks/retrivering_btc_data.ipynb`, or use the file in `data/Raw/`.
2. **MC-GARCH**: run `RScript/mc_garch_btc.R` in RStudio. It writes the six CSV files in `data/mcgarch/`.
3. **Kaggle dataset**: upload the 5-minute CSV and the six MC-GARCH CSV files together to one Kaggle dataset. The notebooks find the files automatically.
4. **Train the LSTM**: run `main_canonical_training_RETRAIN_TRUE.ipynb` on Kaggle with a GPU. It saves the two LSTM forecast files.
5. **Reproduce the results**: add the two files from `data/lstm/` to the Kaggle dataset and run `main_reporting_RETRAIN_FALSE.ipynb`. It recomputes all tables and tests without retraining.

The two notebooks are two run modes of the same analysis (the `RETRAIN_LSTM` setting), not two different experiments. LSTM training on a GPU is not guaranteed to give identical numbers on other hardware, which is why the final LSTM forecasts are saved in `data/lstm/`.

## Main results

Forecast errors on log RV (lower is better). Validation = 2022-2023 walk-forward folds, Test = 2024-2025.

| Model | Val RMSE | Val QLIKE | Test RMSE | Test QLIKE |
|---|---|---|---|---|
| Random Walk | 1.0565 | 1.0393 | 0.9122 | 0.6573 |
| Historical Mean | 1.7185 | 1.0556 | 1.4054 | 0.8038 |
| HAR-Classic | 0.9327 | 0.7524 | 0.8300 | 0.5292 |
| HAR-Extended | 0.9339 | 0.7556 | 0.8305 | 0.5173 |
| MC-GARCH | 0.9542 | 0.8446 | 0.8541 | 0.5597 |
| LSTM | 0.9426 | 0.7107 | **0.7841** | **0.4654** |

- On the test period the LSTM is the most accurate model under both losses (Diebold-Mariano, Holm-adjusted).
- On the validation period the LSTM is better under QLIKE, and not significantly different from the HAR models under squared log error.
- The results apply to one asset, one exchange and a one-hour horizon.
