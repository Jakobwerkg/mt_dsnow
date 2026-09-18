# run_dsnow_defaults.R
# -----------------------------------------------------------------------------
# Reference runs for peak_swe_comparison.ipynb: nixmass::swe.delta.snow with
# its built-in defaults on the Mag25 HS data (no calibration):
#
#   reference_runs/dsnow_default.nc     original ΔSNOW  (dyn_rho_max = FALSE,
#                                       constant rho.max = 401.26 kg m-3)
#   reference_runs/dsnow2.0_default.nc  ΔSNOW2.0        (dyn_rho_max = TRUE,
#                                       rho_max(age) S-curve, Winkler et al. 2021)
#
# Model driving (one run per station x hydrological year Sep-Aug, hs forced to
# 0 at the segment ends, Weisfluh_Joch skipped) and the NetCDF layout
# (station_name, SWE_mod, HNW_mod, SWE_obs, HNW_obs) are the same as in
# calibration/calibration_rho_dyn/run_full_validation_rho_dyn.R, so the two
# files are read by the notebook exactly like the calibrated runs.
#
# Usage:
#   Rscript peak_SWE/run_dsnow_defaults.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})

# ─────────────────────────────────────────────────────────────────────────────
# Project root (.projectroot marker), works for Rscript and interactive use
# ─────────────────────────────────────────────────────────────────────────────
find_project_root <- function(start = NULL) {
  if (is.null(start)) {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    start <- if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
  p <- normalizePath(start, mustWork = TRUE)
  while (!file.exists(file.path(p, ".projectroot"))) {
    if (dirname(p) == p) stop("Could not locate the project root (.projectroot marker not found).")
    p <- dirname(p)
  }
  p
}
ROOT <- find_project_root()

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
MAG25_NC         <- file.path(ROOT, "calibration/calibration_data/raw_data/mag25/slf_dataset/Mag25_all.nc")
OUT_DIR          <- file.path(ROOT, "peak_SWE/reference_runs")
EXCLUDE_STATIONS <- c("Weisfluh_Joch")

# label -> (output file, dyn_rho_max). model_opts = list() -> nixmass defaults.
RUNS <- list(
  "dsnow_default"    = list(file = "dsnow_default.nc",    dyn_rho_max = FALSE, label = "dSNOW default"),
  "dsnow2.0_default" = list(file = "dsnow2.0_default.nc", dyn_rho_max = TRUE,  label = "dSNOW2.0 default")
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers (as in run_full_validation_rho_dyn.R)
# ─────────────────────────────────────────────────────────────────────────────
hyd_year_of <- function(d) ifelse(month(d) > 8L, year(d), year(d) - 1L)

# ΔSnow on one station x one hydrological year; hs forced to 0 at both ends so
# every winter starts and ends snow-free.
run_dsnow <- function(dates, hs_m, dyn_rho_max) {
  hs <- pmax(as.numeric(hs_m), 0)
  hs[is.na(hs)] <- 0
  if (length(hs) == 0) return(NULL)
  hs[1]          <- 0
  hs[length(hs)] <- 0

  df <- data.frame(date = as.character(dates), hs = hs, stringsAsFactors = FALSE)
  out <- tryCatch(
    nixmass::swe.delta.snow(df,
                            model_opts  = list(),        # package defaults
                            dyn_rho_max = dyn_rho_max,
                            layers      = FALSE,
                            strict_mode = FALSE,
                            verbose     = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  if (is.list(out)) as.numeric(out$SWE) else as.numeric(out)
}

calc_metrics <- function(obs, mod) {
  ok  <- is.finite(obs) & is.finite(mod) & obs >= 0
  obs <- obs[ok]; mod <- mod[ok]
  res <- mod - obs
  c(RMSE = sqrt(mean(res^2)), Bias = mean(res), Rel_BIAS = sum(res) / sum(obs),
    R2 = 1 - sum(res^2) / sum((obs - mean(obs))^2), N = length(obs))
}

# ─────────────────────────────────────────────────────────────────────────────
# Load Mag25
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(MAG25_NC))
nc_in <- nc_open(MAG25_NC)
time_dim   <- nc_in$dim[["time"]]
time_raw   <- time_dim$vals
time_units <- time_dim$units
time_cal   <- if (!is.null(time_dim$calendar) && nchar(time_dim$calendar) > 0) time_dim$calendar else "standard"
dates_all  <- as.Date(substr(sub("days since ", "", time_units), 1, 10)) + as.integer(time_raw)

station_names <- nc_in$dim[["station"]]$vals
HS_all        <- ncvar_get(nc_in, "HS")    # [station x time], m
SWE_obs_all   <- ncvar_get(nc_in, "SWE")   # [station x time], mm
HNW_obs_all   <- if ("HNW" %in% names(nc_in$var)) ncvar_get(nc_in, "HNW") else
                 matrix(NA_real_, length(station_names), length(dates_all))
nc_close(nc_in)

Ns <- length(station_names); Nt <- length(dates_all)
hyd_years_all <- hyd_year_of(dates_all)
winter_years  <- sort(unique(hyd_years_all))
keep_station  <- !(station_names %in% EXCLUDE_STATIONS)

message(sprintf("Mag25: %d stations x %d days (%s - %s)", Ns, Nt, dates_all[1], dates_all[Nt]))

# ─────────────────────────────────────────────────────────────────────────────
# NetCDF writer (same layout as the calibrated validation runs)
# ─────────────────────────────────────────────────────────────────────────────
write_validation_nc <- function(nc_file, SWE_mod, HNW_mod, dyn_rho_max, label) {
  dir.create(dirname(nc_file), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(nc_file)) file.remove(nc_file)

  dim_station <- ncdim_def("station", units = "", vals = seq_len(Ns), create_dimvar = FALSE)
  dim_time    <- ncdim_def("time", units = time_units, vals = as.integer(time_raw), create_dimvar = TRUE)
  dim_nchar   <- ncdim_def("nchar", units = "", vals = seq_len(max(nchar(station_names))), create_dimvar = FALSE)

  v_stn     <- ncvar_def("station_name", "",   list(dim_nchar, dim_station), prec = "char")
  v_swe_mod <- ncvar_def("SWE_mod",      "mm", list(dim_station, dim_time), missval = NA_real_, prec = "double")
  v_hnw_mod <- ncvar_def("HNW_mod",      "mm", list(dim_station, dim_time), missval = NA_real_, prec = "double")
  v_swe_obs <- ncvar_def("SWE_obs",      "mm", list(dim_station, dim_time), missval = NA_real_, prec = "double")
  v_hnw_obs <- ncvar_def("HNW_obs",      "mm", list(dim_station, dim_time), missval = NA_real_, prec = "double")

  nc_out <- nc_create(nc_file, vars = list(v_stn, v_swe_mod, v_hnw_mod, v_swe_obs, v_hnw_obs))
  on.exit(nc_close(nc_out), add = TRUE)

  ncvar_put(nc_out, v_stn,     station_names)
  ncvar_put(nc_out, v_swe_mod, SWE_mod)
  ncvar_put(nc_out, v_hnw_mod, HNW_mod)
  ncvar_put(nc_out, v_swe_obs, SWE_obs_all)
  ncvar_put(nc_out, v_hnw_obs, HNW_obs_all)
  ncatt_put(nc_out, "time", "calendar", time_cal)

  defaults <- if (dyn_rho_max) {
    list(sigma = 0.03, mu = 80, rho_h = 600, rho_l = 380, rho.null = 80.73706,
         c.ov = 0.0005170964, k.ov = 0.3782312, k = 0.029297, tau = 0.02356521,
         eta.null = 8543502, timestep = 24)
  } else {
    list(rho.max = 401.2588, rho.null = 81.19417, c.ov = 0.0005104722, k.ov = 0.37856737,
         k = 0.02993175, tau = 0.02362476, eta.null = 8523356, timestep = 24)
  }
  ncatt_put(nc_out, 0, "dsnow_parameters",
            paste(sprintf("%s=%g", names(defaults), unlist(defaults)), collapse = "; "))
  ncatt_put(nc_out, 0, "dyn_rho_max", as.integer(dyn_rho_max))
  ncatt_put(nc_out, 0, "source",      "nixmass::swe.delta.snow")
  ncatt_put(nc_out, 0, "subset",      "reference")
  ncatt_put(nc_out, 0, "phase",       "default")
  ncatt_put(nc_out, 0, "algorithm",   "none")
  ncatt_put(nc_out, 0, "run_label",   label)
  ncatt_put(nc_out, 0, "source_path", sprintf("nixmass %s defaults (dyn_rho_max = %s) on %s",
                                              as.character(packageVersion("nixmass")),
                                              dyn_rho_max, basename(MAG25_NC)))
  invisible(nc_file)
}

# ─────────────────────────────────────────────────────────────────────────────
# Run both defaults
# ─────────────────────────────────────────────────────────────────────────────
for (key in names(RUNS)) {
  r  <- RUNS[[key]]
  t0 <- Sys.time()
  message(sprintf("\n%s (dyn_rho_max = %s) ...", r$label, r$dyn_rho_max))

  SWE_mod <- matrix(NA_real_, Ns, Nt)
  for (si in seq_len(Ns)) {
    if (!keep_station[si]) next
    for (wy in winter_years) {
      idx <- which(hyd_years_all == wy)
      if (length(idx) < 10L) next
      seg <- run_dsnow(dates_all[idx], HS_all[si, idx], r$dyn_rho_max)
      if (!is.null(seg) && length(seg) == length(idx)) SWE_mod[si, idx] <- seg
    }
  }

  # HNW_mod = daily increase of SWE_mod, melt clipped to 0
  HNW_mod <- matrix(NA_real_, Ns, Nt)
  for (si in seq_len(Ns)) {
    d <- diff(SWE_mod[si, ]); d[d < 0] <- 0
    HNW_mod[si, 2:Nt] <- d
  }

  nc_file <- file.path(OUT_DIR, r$file)
  write_validation_nc(nc_file, SWE_mod, HNW_mod, r$dyn_rho_max, r$label)

  # Sanity check: full-year SWE metrics as in the R validation
  m <- calc_metrics(as.vector(SWE_obs_all[keep_station, ]), as.vector(SWE_mod[keep_station, ]))
  message(sprintf("  SWE: RMSE %.2f | Bias %+.2f | Rel_BIAS %+.3f | R2 %.3f | N %d   (%.1f min)",
                  m[["RMSE"]], m[["Bias"]], m[["Rel_BIAS"]], m[["R2"]], as.integer(m[["N"]]),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  message("  wrote ", sub(paste0("^", ROOT, "/?"), "", nc_file))
}
