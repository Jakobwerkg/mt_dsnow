# Calibration

ΔSNOW (`nixmass::swe.delta.snow`) is driven by **snow depth only**; the calibration tunes its parameters so
that the modelled SWE matches a reference SWE series. The pipeline is

```
prepare observations  →  optimise parameters  →  collect results  →  validate on Mag25  →  overview notebooks
```

Two model variants are calibrated in separate folders:

| Variant | Parameters | Folder | Optimisers |
|---|---|---|---|
| ΔSNOW, static ρ_max (original model) | 7 (`rho.max`, `rho.null`, `c.ov`, `k.ov`, `k`, `tau`, `eta.null`) | `calibration_win21/`, `calibration_snowpack/` | Nelder–Mead `_nm.R`, Differential Evolution `_de.R` |
| ΔSNOW2.0, dynamic ρ_max(age) (Winkler et al. 2021) | 10 (`sigma`, `mu`, `rho_h`, `rho_l` replace `rho.max`) | `calibration_rho_dyn/` | Nelder–Mead `_nm_dsnow2.0.R` |

All scripts share the objective `score = Σ w_i · metric_i` over SWE and bulk density ρ = SWE / HS
(NRMSE, NBIAS, 1 − KGE for each), weights given on the command line. A weight combination is a
**phase** (`1A` … `6D`, table in `calibration_rho_dyn/collect_rho_dyn_results.R`).

## Data requirements

The model needs a **gap-free daily HS series** per hydrological year (1 Aug – 31 Jul): `nixmass` stops on
any `NA` snow depth, so in the optimisers a winter block with a single missing HS day is silently dropped
from the fit (`tryCatch` → `NA`). The reference SWE may be sparse — the score is evaluated only on days that
have a SWE value (ρ additionally needs HS > 0).

| Dataset | Prepared by → file read by the calibration | HS (drives the model) | SWE (scored against) | Used for |
|---|---|---|---|---|
| **Win21** — 17 stations of Winkler et al. (2021) | `calibration_data/raw_data/win21/source_data/H_SWE_obs.Rda` (list of `zoo` series, `Hobs`, `SWEobs`) — copy it by hand to `optimisation_output/win21/data/H_SWE_obs.Rda` | daily, **cm** (scripts scale ×1/100), complete | manual measurements, ~every **14 days** | calibration (`calibration_win21/`, `calibration_rho_dyn/` with `DSNOW_SUBSET=win21`) |
| **SNOWPACK** — `sp_all` (18 stations), `sp_rg` (rain-gauge stations), `sp_b2000` (< 2000 m) | hourly SNOWPACK `.smet` under `calibration_data/snowpack_data/{data_18_all, data_rain_gauge, station_below_2000}/` → `calibration_snowpack/prepare_snowpack_data.R` → `optimisation_output/<subset>/data/d_obs_SNOWPACK.rda` | daily (first record of each day of `HS_mod`), m | **daily** SNOWPACK `SWE` (model-to-model calibration) | calibration (`calibration_snowpack/`, `calibration_rho_dyn/`) |
| **Mag25** — 41 SLF stations, Sep 2016 – Aug 2022 | SLF text files → `calibration_data/raw_data/mag25/raw2csv_mag25.py` → `calibration_data/raw_data/mag25/slf_dataset/Mag25_all.nc` (`HS` m, `SWE` mm, `HNW` mm, `HN` m, `altitude`) | daily, m, complete (validation fills gaps with 0 and forces 0 at both ends of a season) | manual **pits ~every 15 days** at 23 stations; plus **daily HNW** (Nov–Apr, ~80 % of days) | independent validation only: `calibration_rho_dyn/run_full_validation_rho_dyn.R`, `hnw_validation/`, `peak_SWE/` |

Rules applied by the optimisers to every station: winters with fewer than 200 days, or incomplete winters
with snow (> 5 cm) at either end, are skipped; the remaining winters alternate between the **validation set**
(1st, 3rd, …) and the **fit set** (2nd, 4th, …). Win21 only: `kuehtai` and `Weissfluhjoch` are dropped,
`Sta.Maria` is validation-only.

## Running

```bash
# 1. observations (once per SNOWPACK subset; Win21: copy H_SWE_obs.Rda, see table)
DSNOW_SUBSET=sp_rg Rscript calibration/calibration_snowpack/prepare_snowpack_data.R

# 2. one weight combination = one phase           SWE_NRMSE RHO_NRMSE SWE_NBIAS RHO_NBIAS [SWE_KGE] [RHO_KGE]
./calibration/run_calibration.sh                  1.0       0.0       0.0       0.0

# 3. the phase sweep (several hours; the phases to run are listed in the script)
./calibration/run_all_phases.sh
```

`run_calibration.sh` runs the optimiser scripts that are enabled (uncommented) inside it — currently only
the ΔSNOW2.0 Nelder–Mead script; `run_all_phases.sh` likewise lists the active phases. The SNOWPACK subset
comes from `DSNOW_SUBSET` (`sp_all`, `sp_rg`, `sp_b2000`; the ΔSNOW2.0 script also accepts `win21`).
Defaults: `sp_all` for `prepare_snowpack_data.R` and the ΔSNOW2.0 script, `sp_rg` for the static
SNOWPACK optimisers. Results are written as `opt_results__<weights>.rds` (NM) /
`opt_results_DE__<weights>.rds` (DE), tagged by the weight combination.

Static variant: results land in `optimisation_output/<subset>/data/R_opt_logs[_DE]/`, are collected by
`optimisation_output/helpers/collect_opt_results.R` (→ `<subset>/res/opt_results_summary.csv`) and
`collect_summaries.py` (→ `optimisation_output/combinded_res/all_summaries.csv`), and validated by
`hnw_validation/full_validation/run_full_validation.R`.

## `calibration_rho_dyn/` — ΔSNOW2.0 (dynamic ρ_max)

Run in this order; every step reads the previous one's output under `dyn_rho_max_res/` (git-ignored):

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `dsnow_parameter_optimization_nm_dsnow2.0.R <weights>` | `optimisation_output/<subset>/data/` (table above) | `dyn_rho_max_res/<subset>/nm_res/opt_results__<weights>.rds` |
| 2 | `collect_rho_dyn_results.R` | all `opt_results*.rds` under `dyn_rho_max_res/` | `dyn_rho_max_res/all_summaries_dyn_rho_max.csv` (+ `.rds`); one row per run, phase label from the weights, plus a control row with the nixmass defaults |
| 3 | `run_full_validation_rho_dyn.R` (needs the `future` package) | step 2 + `Mag25_all.nc` | `dyn_rho_max_res/validation/all_summaries_validated_dyn_rho_max.csv` (SWE & HNW metrics appended) and one NetCDF per run in `validation/nc/` (daily `SWE_obs/mod`, `HNW_obs/mod` per station) |
| 4 | `rho_dyn_results_overview.ipynb` | step 3 | ranking by SWE RMSE, \|HNW rel. bias\| and their combination → `dyn_rho_max_res/rho_dyn_results_combined_rank_phases.csv`; figures in `dyn_rho_max_res/plots/` |
| 5 | `../../peak_SWE/peak_swe_comparison.ipynb` | step 3 NetCDFs + step 4 CSV + `Mag25_all.nc` | peak-SWE / bulk-density evaluation of every run |

`hnw_validation_helper.py` holds the metric and plotting helpers used by the notebook.

---

*README created by Claude Opus 5 (Anthropic) on 2026-09-18 from the scripts and data in this folder,
under the direction of the author.*
