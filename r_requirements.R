# R dependencies for the calibration, validation and model-run scripts.
# Run once:  Rscript r_requirements.R

packages <- c(
  "nixmass",      # DeltaSnow model (Winkler et al. 2021)
  "optimx",       # Nelder-Mead parameter optimisation
  "DEoptim",      # Differential Evolution parameter optimisation
  "ncdf4",        # NetCDF I/O
  "zoo",          # irregular time series
  "lubridate",    # date handling
  "tidyverse",    # data wrangling
  "foreach",      # parallel loops
  "doParallel"    # parallel backend for foreach
)

missing <- packages[!packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}
