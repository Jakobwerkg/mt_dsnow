# DeltaSnow Parameter Calibration

Calibrates the 7 free parameters of the **deltasnow** model (`nixmass::swe.delta.snow`) against observed HS and SWE data from two independent datasets. The optimizer minimises a weighted combination of normalised SWE error, density error, and SWE bias.

---

## Pipeline Overview

```
Raw data                Data prep              Calibration              Results
─────────────────────   ────────────────────   ──────────────────────   ─────────────────────
SNOWPACK .smet files ──▶ prepare_SNOWPACK_      run_all_phases.sh
                         data.R                  └─ run_calibration.sh
                          │                           ├─ SNOWPACK NM ──▶ R_opt_logs/
                          ▼                           ├─ SNOWPACK DE ──▶ R_opt_logs_DE/
                     d_obs_SNOWPACK.rda               ├─ Win21 NM    ──▶ R_opt_logs/
                                                      └─ Win21 DE    ──▶ R_opt_logs_DE/
SLF Mag25 txt files ───▶ raw2csv_Mag25.py                                      │
Win21 H_SWE_obs.Rda ───▶ raw2csv_Win21.R       ◀────────────────────────────────┘
                          │                    collect_opt_results.R
                          ▼                      └─▶ opt_results_summary_grid_search.csv
                  HS_SWE_by_station/*.csv
                          │
                          ▼
                  csv2rda_calib.R
                          │
                          ▼
                   d_obs_WIN_MAG.rda
                   (used by Win21 scripts)
```

---

## File Descriptions

### Shell scripts (entry points)

| File | Purpose |
|------|---------|
| `run_all_phases.sh` | Top-level script. Runs a sequence of weight combinations (Phases 1–5). Each phase calls `run_calibration.sh`. Comment/uncomment phases to control what runs. |
| `run_calibration.sh` | Runs all 4 optimisation scripts sequentially for one weight triple `(w_swe, w_rho, w_bias)`. Usage: `./run_calibration.sh 0.4 0.6 0.0` |

### Data preparation

| File | Input | Output |
|------|-------|--------|
| `calibration_SNOWPACK/prepare_SNOWPACK_data.R` | `.smet` files from SNOWPACK simulation output (`par_sens/SNOWPACK_data/data_rain_gauge/raw_alpsolut/`) | `calibration_SNOWPACK/data/d_obs_SNOWPACK.rda` — named list of zoo objects, one per station, with columns `Hobs` (m) and `SWEobs` (mm) |
| `calibration_data/raw_data/Mag25/raw2csv_Mag25.py` | SLF raw text files: `OBS-HN.txt`, `OBS-HNW.txt`, `OBS-HS-STAKE.txt`, `OBS-SWE-PROFILE.txt`, station list | Per-station CSVs in `calibration_data/output/HS_SWE_by_station/` with columns `date, hs [m], swe_obs [mm]`. Also saves `Mag25_all.nc`. |
| `calibration_data/raw_data/dsnow/raw2csv_Win21.R` | `calibration_Win21/data/H_SWE_obs.Rda` (Win et al. 2021 dataset) | Per-station CSVs in `calibration_data/output/HS_SWE_by_station/` (same format as above) |
| `calibration_data/combining_data/csv2rda_calib.R` | Per-station CSVs from `HS_SWE_by_station/` | `calibration_data/output/calibration_rda_files/d_obs_WIN_MAG.rda` — named list of zoo objects used by the Win21 calibration scripts |

### Calibration scripts (4 optimisers)

All four scripts share the same logic — they differ only in dataset and optimisation algorithm.

| File | Dataset | Algorithm | Output dir |
|------|---------|-----------|-----------|
| `calibration_SNOWPACK/dsnow_parameter_optimization.R` | `d_obs_SNOWPACK.rda` | Nelder-Mead (`optimx`) | `calibration_SNOWPACK/data/R_opt_logs/` |
| `calibration_SNOWPACK/dsnow_parameter_optimization_DE.R` | `d_obs_SNOWPACK.rda` | Differential Evolution (`DEoptim`) | `calibration_SNOWPACK/data/R_opt_logs_DE/` |
| `calibration_Win21/WIN21_dsnow_paprameter_optimization.R` | `H_SWE_obs.Rda` | Nelder-Mead (`optimx`) | `calibration_Win21/data/R_opt_logs/` |
| `calibration_Win21/WIN21_dsnow_parameter_opitmization_DE.R` | `H_SWE_obs.Rda` | Differential Evolution (`DEoptim`) | `calibration_Win21/data/R_opt_logs_DE/` |

Each script saves one `.rds` file per run, named:
```
opt_results__SWE_NRMSE_<w>__RHO_NRMSE_<w>__SWE_NBIAS_<w>.rds
```

### Result collection

| File | Purpose |
|------|---------|
| `collect_opt_results.R` | Scans all `opt_results*.rds` files, extracts best parameters and scores, writes `opt_results_summary_grid_search.csv` and `.rds`. Run from the repo root: `Rscript calibration/collect_opt_results.R` |

---

## Objective Function

All scripts minimise the same weighted score:

```
score = w1 · NRMSE_SWE + w2 · NRMSE_rho + w3 · NBIAS_SWE

NRMSE_SWE = RMSE(SWE_mod, SWE_obs) / mean(SWE_obs)
NRMSE_rho = RMSE(rho_mod, rho_obs) / mean(rho_obs)     # bulk density = SWE / HS
NBIAS_SWE = |mean(SWE_mod - SWE_obs)| / mean(SWE_obs)
```

Weights are passed via command line and must sum to 1. The filename encodes them (e.g. `SWE_NRMSE_0p4__RHO_NRMSE_0p6__SWE_NBIAS_0p0`).

---

## Parameters Being Optimised

| Parameter | Physical meaning | Bounds |
|-----------|-----------------|--------|
| `rho.max` | Maximum bulk snow density [kg/m³] | 300 – 600 |
| `rho.null` | Fresh snow density [kg/m³] | 60 – 150 |
| `c.ov` | Overburden compaction coefficient | 1e-6 – 0.01 |
| `k.ov` | Overburden compaction exponent | 0.01 – 1.0 |
| `k` | Settling rate | 0.001 – 0.1 |
| `tau` | Time constant | 0.001 – 0.1 |
| `eta.null` | Snow viscosity [Pa·s] | 1e5 – 1e8 |

---

## Train / Validation Split

Per station, winters are split **alternating by season**:
- Even winters (counter 0, 2, 4, …) → **fit set**
- Odd winters (counter 1, 3, 5, …) → **validation set**

Winters shorter than 200 days or with snow remaining at the season boundary (incomplete record) are skipped entirely. The hydrological year runs **1 Aug → 31 Jul**.

---

## Running the Pipeline

### Full Phase 5 run (currently active)

```bash
cd ~/code/mt_dsnow/calibration
./run_all_phases.sh
```

This runs the 4 weight combinations defined in Phase 5 of `run_all_phases.sh`:

| Phase | w_SWE | w_RHO | w_BIAS |
|-------|-------|-------|--------|
| 5A | 0.8 | 0.1 | 0.1 |
| 5B | 0.1 | 0.8 | 0.1 |
| 5C | 0.4 | 0.4 | 0.2 |
| 5D | 0.25 | 0.25 | 0.5 |

Each combination runs all 4 optimisers → **4 × 4 = 16 `.rds` files** produced.

### Single weight combination

```bash
./run_calibration.sh 0.5 0.5 0.0
```

### Collect and compare all results

```bash
Rscript calibration/collect_opt_results.R
# → opt_results_summary_grid_search.csv
```

---

## Dependencies (R)

```r
install.packages(c("optimx", "DEoptim", "zoo", "foreach",
                   "doParallel", "lubridate", "nixmass", "tidyverse"))
```

The `nixmass` package provides `swe.delta.snow()`, the deltasnow forward model used inside the objective function.






Stora: 





from pathlib import Path
import os
import sys
import io
import warnings
from contextlib import redirect_stdout
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from joblib import Parallel, delayed
from tqdm.auto import tqdm


def _find_repo_root():
    candidates = [
        Path.cwd(),
        Path.cwd().parent,
        Path('/Users/jakobwerkgarner/code/mt_dsnow'),
    ]
    for c in candidates:
        if (c / 'snow_to_swe_master' / 'main.py').exists():
            return c.resolve()
    raise FileNotFoundError('Could not locate repository root containing snow_to_swe_master/main.py')


repo_root = _find_repo_root()
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'snow_to_swe_master') not in sys.path:
    sys.path.insert(0, str(repo_root / 'snow_to_swe_master'))

import HNW_validation.HNW_validation_helper as val_helper
from main import SnowToSwe

SHOW_SNOWTOSWE_BANNER = False


def _plot_validation_on_ax(ax, x, y, stats, title, lim, xlabel, ylabel):
    vmax = max(1, len(x) / 10)

    h = ax.hist2d(
        x, y,
        bins=50,
        range=[lim, lim],
        norm=LogNorm(vmin=1, vmax=vmax),
        cmap='viridis',
    )

    cb = plt.colorbar(h[3], ax=ax, label='Number of observations')
    ticks = [t for t in [1, 10, 100, 1000, 10000] if t <= vmax]
    if vmax not in ticks:
        ticks.append(vmax)
    cb.set_ticks(ticks)
    cb.set_ticklabels([f'{int(t)}' for t in ticks[:-1]] + [f'{int(vmax)}'])

    ax.plot(lim, lim, '--', color='gray', linewidth=1.2)
    ticks_xy = np.linspace(lim[0], lim[1], 5)
    ax.set_xticks(ticks_xy)
    ax.set_yticks(ticks_xy)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)

    textstr = (
        f"R2: {stats['R2']:.2f}\n"
        f"Bias: {stats['Bias']:.2f}\n"
        f"RMSE: {stats['RMSE']:.2f}\n"
        f"Rel_BIAS: {stats['Rel_BIAS']:.1%}\n"
        f"N: {stats['N']}"
    )
    ax.text(
        0.03, 0.97, textstr,
        transform=ax.transAxes,
        fontsize=10,
        va='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.85),
    )

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(False)


def _process_station_winter(station_name, winter_year, time_idx, hs_values, snow_to_swe):
    try:
        hs = pd.Series(hs_values).fillna(0).clip(lower=0).astype(float)
        if len(hs) == 0:
            return None

        # SnowToSwe expects snow-free boundaries
        if hs.iloc[0] != 0:
            hs.iloc[0] = 0.0
        if hs.iloc[-1] != 0:
            hs.iloc[-1] = 0.0

        swe_list = snow_to_swe.convert_list(hs.tolist(), timestep=24, verbose=False)
        if swe_list is None:
            return None

        swe_arr = np.asarray(swe_list, dtype=float)
        if swe_arr.shape[0] != len(hs):
            return None

        return (station_name, winter_year, time_idx, swe_arr)
    except Exception as e:
        return ('__ERROR__', station_name, f'winter {winter_year}/{winter_year + 1}: {e}')


def _run_and_plot_one_calibration(row, base_data, tasks, out_dir, run_idx, total_runs, SnowToSwe):
    if SHOW_SNOWTOSWE_BANNER:
        snow_to_swe = SnowToSwe(
            rho_max=float(row['rho.max_phys']),
            rho_null=float(row['rho.null_phys']),
            c_ov=float(row['c.ov_phys']),
            k_ov=float(row['k.ov_phys']),
            k=float(row['k_phys']),
            tau=float(row['tau_phys']),
            eta_null=float(row['eta.null_phys']),
        )
    else:
        with redirect_stdout(io.StringIO()):
            snow_to_swe = SnowToSwe(
                rho_max=float(row['rho.max_phys']),
                rho_null=float(row['rho.null_phys']),
                c_ov=float(row['c.ov_phys']),
                k_ov=float(row['k.ov_phys']),
                k=float(row['k_phys']),
                tau=float(row['tau_phys']),
                eta_null=float(row['eta.null_phys']),
            )

    mag = base_data.copy()
    mag['SWE_mod'] = xr.full_like(mag['HS'], np.nan)

    results = Parallel(n_jobs=-1, return_as='generator')(
        delayed(_process_station_winter)(stn, y, tidx, hsv, snow_to_swe)
        for (stn, y, tidx, hsv) in tasks
    )

    errors = 0
    ok_tasks = 0
    for result in results:
        if result is None:
            continue
        if isinstance(result, tuple) and result and result[0] == '__ERROR__':
            errors += 1
            continue

        station_name, _, time_idx, swe_arr = result
        mag['SWE_mod'].loc[dict(station=station_name, time=time_idx)] = swe_arr
        ok_tasks += 1

    hnw_mod = mag['SWE_mod'].diff(dim='time').clip(min=0)
    hnw_mod = hnw_mod.reindex(time=mag['time'])
    mag['HNW_mod'] = hnw_mod

    all_df = (
        mag[['HNW', 'HNW_mod', 'SWE', 'SWE_mod', 'HS']]
        .to_dataframe()
        .reset_index()
        .rename(columns={'HNW': 'HNW_obs', 'SWE': 'SWE_obs'})
    )

    # Match the original validation filtering strategy
    swe_val = all_df[all_df['SWE_obs'].notna()].copy()
    swe_val.index = pd.to_datetime(swe_val['time'])
    swe_val = val_helper._filter_season(swe_val, full_season=True)
    swe_val = swe_val.dropna(subset=['SWE_obs', 'SWE_mod'])
    swe_val = swe_val[swe_val['SWE_obs'] >= 0]
    swe_val = swe_val[np.isfinite(swe_val['SWE_obs']) & np.isfinite(swe_val['SWE_mod'])]

    hnw_val = all_df[all_df['HNW_obs'].notna()].copy()
    hnw_val.index = pd.to_datetime(hnw_val['time'])
    hnw_val = val_helper._filter_season(hnw_val, full_season=False)
    hnw_val = hnw_val.dropna(subset=['HNW_obs', 'HNW_mod'])
    hnw_val = hnw_val[hnw_val['HNW_obs'] >= 0]
    hnw_val = hnw_val[np.isfinite(hnw_val['HNW_obs']) & np.isfinite(hnw_val['HNW_mod'])]

    hs = val_helper._calculate_metrics(hnw_val['HNW_obs'].values, hnw_val['HNW_mod'].values)
    hs['N'] = len(hnw_val)
    ss = val_helper._calculate_metrics(swe_val['SWE_obs'].values, swe_val['SWE_mod'].values)
    ss['N'] = len(swe_val)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    _plot_validation_on_ax(
        axes[0],
        hnw_val['HNW_obs'].values,
        hnw_val['HNW_mod'].values,
        hs,
        title='HNW validation',
        lim=[0, 100],
        xlabel='Observed HNW (mm)',
        ylabel='Modeled HNW (mm)',
    )
    _plot_validation_on_ax(
        axes[1],
        swe_val['SWE_obs'].values,
        swe_val['SWE_mod'].values,
        ss,
        title='SWE validation',
        lim=[0, 1000],
        xlabel='Observed SWE (mm)',
        ylabel='Modeled SWE (mm)',
    )

    meta_title = (
        f"Run {run_idx}/{total_runs} | Dataset: {row['dataset']} | Algorithm: {row['algorithm']} | "
        f"Weights (SWE/RHO/BIAS)=({row['w_SWE_NRMSE']:.1f}/{row['w_RHO_NRMSE']:.1f}/{row['w_SWE_NBIAS']:.1f}) | "
        f"Score={row['best_value']:.4f}"
    )
    fig.suptitle(meta_title, fontsize=12, fontweight='bold', y=1.03)

    param_text = (
        f"rho_max={row['rho.max_phys']:.4g}, rho_null={row['rho.null_phys']:.4g}, c_ov={row['c.ov_phys']:.4g}, "
        f"k_ov={row['k.ov_phys']:.4g}, k={row['k_phys']:.4g}, tau={row['tau_phys']:.4g}, eta_null={row['eta.null_phys']:.4g}"
    )
    fig.text(0.5, -0.02, param_text, ha='center', va='top', fontsize=9)
    plt.tight_layout()

    file_stub = (
        f"run_{run_idx:02d}__{row['dataset']}__{row['algorithm']}__"
        f"SWE_{row['w_SWE_NRMSE']:.1f}_RHO_{row['w_RHO_NRMSE']:.1f}_BIAS_{row['w_SWE_NBIAS']:.1f}"
    ).replace(' ', '_').replace('/', '-')
    out_file = out_dir / f"{file_stub}.png"
    fig.savefig(out_file, dpi=250, bbox_inches='tight')
    plt.show()

    return {
        'run_idx': run_idx,
        'dataset': row['dataset'],
        'algorithm': row['algorithm'],
        'w_SWE': row['w_SWE_NRMSE'],
        'w_RHO': row['w_RHO_NRMSE'],
        'w_BIAS': row['w_SWE_NBIAS'],
        'score': row['best_value'],
        'HNW_RMSE': hs['RMSE'],
        'HNW_R2': hs['R2'],
        'SWE_RMSE': ss['RMSE'],
        'SWE_R2': ss['R2'],
        'ok_tasks': ok_tasks,
        'errors': errors,
        'figure': str(out_file),
    }


warnings.filterwarnings('ignore', category=FutureWarning, module=r'main')

mag25_candidates = [
    repo_root / 'calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc',
    Path('/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc'),
]
mag25_path = next((p for p in mag25_candidates if p.exists()), None)
if mag25_path is None:
    raise FileNotFoundError('Could not find Mag25_all.nc')

base_data = xr.open_dataset(mag25_path).drop_sel(station='Weisfluh_Joch')
times = pd.to_datetime(base_data['time'].values)
hyd_year = np.where(times.month > 8, times.year, times.year - 1)
winter_years = np.unique(hyd_year)
station_list = base_data['station'].values
hs_by_station = {s: base_data['HS'].sel(station=s).values for s in station_list}

tasks = []
for stn in station_list:
    hs_full = hs_by_station[stn]
    for y in winter_years:
        mask = hyd_year == y
        time_idx = times[mask].values
        hs_vals = hs_full[mask]
        tasks.append((stn, int(y), time_idx, hs_vals))

print(f'Prepared {len(tasks)} station×winter tasks for validation.')

out_dir = repo_root / 'calibration/plots/hnw_swe_validation_per_run'
out_dir.mkdir(parents=True, exist_ok=True)

runs = df_filtered.copy().reset_index(drop=True)
print(f'Generating side-by-side HNW/SWE plots for {len(runs)} calibration runs...')

summaries = []
for i, (_, row) in enumerate(tqdm(runs.iterrows(), total=len(runs), desc='Calibration runs'), start=1):
    try:
        summary = _run_and_plot_one_calibration(
            row=row,
            base_data=base_data,
            tasks=tasks,
            out_dir=out_dir,
            run_idx=i,
            total_runs=len(runs),
            SnowToSwe=SnowToSwe,
        )
        summaries.append(summary)
    except Exception as e:
        print(f'[WARNING] Run {i} failed ({row["dataset"]} / {row["algorithm"]}): {e}')

validation_summary = pd.DataFrame(summaries)
display(validation_summary)
print(f'Finished. Saved figures to: {out_dir}')