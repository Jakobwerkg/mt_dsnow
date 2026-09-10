# mt_dsnow — Master's Thesis Data Analysis

Data analysis code for my Master's thesis in the **Master of Atmospheric and Cryospheric
Sciences (ACINN)** at the **University of Innsbruck**.

The thesis compares and recalibrates two empirical snow models that convert snow-depth (HS)
measurements into snow water equivalent (SWE):

- **ΔSnow** (`nixmass::swe.delta.snow`, Winkler et al. 2021) — R
- **HS2SWE** (Magnusson et al.) — Python

Work covers three datasets: **Win21** (Winkler et al. 2021 manual observations),
**Mag25** (Magnusson et al. 2025, SLF automatic stations) and **SNOWPACK**
(18 Alpsolut stations, SNOWPACK-modelled SWE). The SNOWPACK subsets are labelled
`SP_all` (all stations), `SP_RG` (rain-gauge subset) and `SP_b2000` (below 2000 m);
their directories use the lowercase equivalents `sp_all`, `sp_rg`, `sp_b2000`.

---

## Getting started

```bash
pip install -r requirements.txt      # Python
Rscript r_requirements.R             # R
```

Every script and notebook locates the repository through the **`.projectroot`** marker
file in this directory, so the repository can be cloned anywhere and each notebook runs
from its own folder. Do not delete or move that file.

```python
# Python / notebooks
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / ".projectroot").exists())
```

```r
# R — see the find_project_root() helper at the top of each script
ROOT <- find_project_root()
```

### External dependency

The Morris sensitivity notebooks in `par_sens/` import a Python port of the ΔSnow model
from a `snow_to_swe_master/` directory in this folder. That code is **not** part of this
repository (it is GPL-3 licensed third-party code, kept out of git deliberately); obtain
it separately and place it at `snow_to_swe_master/` before running those notebooks.

### Data

Raw and derived data (`.nc`, `.csv`, `.rda`, `.rds`, `.npy`, `.smet`) and all generated
figures are git-ignored — see `.gitignore`. The repository holds **code only**; every
dataset and figure is reproduced by running the scripts and notebooks in order.

---

## Folder structure

| Folder | What happens there |
|---|---|
| `calibration/calibration_data/` | Raw-data ingestion, quality checks and conversion of all three datasets into common `.nc` / `.rda` / `.csv` formats. Per-dataset subfolders under `raw_data/` (`win21`, `mag25`) and `snowpack_data/`. |
| `calibration/calibration_win21/`, `calibration/calibration_snowpack/` | ΔSnow parameter optimisation in R — `_nm.R` = Nelder–Mead, `_de.R` = Differential Evolution. |
| `calibration/optimisation_output/` | One directory per subset (`win21`, `sp_all`, `sp_rg`, `sp_b2000`, `dyn_rho_max`) holding that run's input data, optimiser logs and result summaries. `helpers/` collects and tabulates them. |
| `hnw_validation/` | Independent validation of the calibrated parameter sets against observed new-snow water equivalent (HNW) and SWE. `full_validation/` runs all parameter sets; results are plotted in `plot_validation_results.ipynb`. |
| `par_sens/` | Morris parameter-sensitivity analysis, run separately per dataset and compared jointly in `morris_comparison.ipynb`. |
| `plot_style.py` | Project-wide plot style — colours, linestyles, subset labels, subplot lettering. Imported by every notebook. |

## Running the calibration

```bash
./calibration/run_all_phases.sh          # full phase sweep (several hours)
./calibration/run_calibration.sh 0.5 0.5 0.0 0.0   # a single weight combination
```

The SNOWPACK scripts pick their subset from the `DSNOW_SUBSET` environment variable
(`sp_all`, `sp_rg`, `sp_b2000`, `dyn_rho_max`; default `sp_rg`):

```bash
DSNOW_SUBSET=sp_b2000 Rscript calibration/calibration_snowpack/prepare_snowpack_data.R
DSNOW_SUBSET=sp_b2000 Rscript calibration/calibration_snowpack/dsnow_parameter_optimization_de.R 0.5 0.5 0 0
```

Each subset reads and writes `calibration/optimisation_output/<subset>/data/`.

## Conventions

- Analysis lives in **Jupyter notebooks** (Python); model runs and calibration are **R scripts**.
- Files and folders are lowercase `snake_case`; `.R` is the only capitalised extension.
- Notebooks are committed **with executed outputs**; re-render from the notebook's own
  directory.
- All figures use `plot_style.py` so colours and labels stay consistent across the thesis.
- `Archive/` folders hold superseded work and are **git-ignored** (local only), as are
  data files, figures and other large outputs.
