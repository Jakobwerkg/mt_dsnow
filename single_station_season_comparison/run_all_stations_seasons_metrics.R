# run_all_stations_seasons_metrics.R
#
# Batch version of single_station_season_comparison.R:
# runs ΔSnow (default + best_new) for EVERY Mag25 station × season
# (hydrological years, Sep–Aug), loads the precomputed HS2SWE output, and
# writes per-run metrics (RMSE, rel.bias, R², n — separately for SWE and HNW)
# plus the underlying daily timeseries to CSV.
#
# Outputs (in this directory):
#   metrics_all_stations_seasons.csv     station × season × model metrics
#   timeseries_all_stations_seasons.csv  daily obs + modelled SWE per run
#
# Both CSVs are consumed by best_runs_ranking.ipynb.

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})

# ─────────────────────────────────────────────────────────────────────────────
# Project root
# Located through the .projectroot marker file, so the repository can live
# anywhere.  Works both for `Rscript path/to/script.R` and interactive use.
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
    if (dirname(p) == p) {
      stop("Could not locate the project root (.projectroot marker not found).")
    }
    p <- dirname(p)
  }
  p
}

ROOT <- find_project_root()


if (!l10n_info()$`UTF-8`)
  suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
exclude_stations <- c("Weisfluh_Joch")   # project-wide exclusion
dyn_rho_max      <- FALSE

model_opts_default <- list(
  rho.null = 81,    rho.max = 401,
  eta.null = 8.5e6, k = 0.0300,
  tau = 0.024,      c.ov = 5.1e-4,
  k.ov = 0.38
)

# "best_new" calibration — c.ov = 6.0e-4 (see single_station_season_comparison.R)
model_opts_best_new <- list(
  rho.null = 101.17,   rho.max = 380.01,
  eta.null = 8.233e6,  k = 0.0269,
  tau = 0.0227,        c.ov = 6.0e-4,
  k.ov = 0.4117
)

Mag25_nc_file  <- file.path(ROOT, "calibration/calibration_data/raw_data/mag25/slf_dataset/Mag25_all.nc")
HS2SWE_nc_file <- file.path(ROOT, "model_diff/layerwise_data/combined_layerwise_default_Mag25.nc")
out_dir        <- file.path(ROOT, "single_station_season_comparison")

metrics_csv    <- file.path(out_dir, "metrics_all_stations_seasons.csv")
timeseries_csv <- file.path(out_dir, "timeseries_all_stations_seasons.csv")

# ─────────────────────────────────────────────────────────────────────────────
# Helpers (identical conventions to single_station_season_comparison.R)
# ─────────────────────────────────────────────────────────────────────────────
hyd_year_of <- function(d) ifelse(month(d) > 8L, year(d), year(d) - 1L)

nc_dates <- function(nc) {
  time_dim  <- nc$dim[["time"]]
  origin_dt <- as.Date(substr(sub("days since ", "", time_dim$units), 1, 10))
  origin_dt + as.integer(time_dim$vals)
}

nc_station_names <- function(nc) {
  vals <- nc$dim[["station"]]$vals
  if (is.character(vals)) return(vals)
  as.character(ncvar_get(nc, "station"))
}

run_dsnow <- function(dates, hs_m, model_opts, dyn_rho_max) {
  hs <- pmax(as.numeric(hs_m), 0)
  hs[is.na(hs)] <- 0
  if (length(hs) == 0)       return(NULL)
  hs[1]          <- 0
  hs[length(hs)] <- 0

  df  <- data.frame(date = as.character(dates), hs = hs, stringsAsFactors = FALSE)
  out <- tryCatch(
    nixmass::swe.delta.snow(df,
                            model_opts  = model_opts,
                            dyn_rho_max = dyn_rho_max,
                            layers      = FALSE,
                            strict_mode = FALSE,
                            verbose     = FALSE),
    error = function(e) { warning("swe.delta.snow error: ", e$message); NULL }
  )
  if (is.null(out)) return(NULL)
  if (is.list(out)) as.numeric(out$SWE) else as.numeric(out)
}

hnw_from_swe <- function(swe) {
  hnw <- c(NA_real_, diff(swe))
  hnw[hnw < 0] <- 0
  hnw
}

metrics <- function(mod, obs) {
  ok <- is.finite(mod) & is.finite(obs)
  n  <- sum(ok)
  if (n < 2L) return(list(rmse = NA_real_, rel_bias = NA_real_, r2 = NA_real_, n = n))
  m <- mod[ok]; o <- obs[ok]
  list(
    rmse     = sqrt(mean((m - o)^2)),
    rel_bias = 100 * sum(m - o) / sum(o),
    r2       = suppressWarnings(cor(m, o))^2,
    n        = n
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Load Mag25 observations
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(Mag25_nc_file), file.exists(HS2SWE_nc_file))

nc_mag       <- nc_open(Mag25_nc_file)
dates_mag    <- nc_dates(nc_mag)
stations_mag <- nc_station_names(nc_mag)
HS_all       <- ncvar_get(nc_mag, "HS")    # [station × time], m
SWE_all      <- ncvar_get(nc_mag, "SWE")   # [station × time], mm
HNW_all      <- ncvar_get(nc_mag, "HNW")   # [station × time], mm
nc_close(nc_mag)

nc_h2s       <- nc_open(HS2SWE_nc_file)
dates_h2s    <- nc_dates(nc_h2s)
stations_h2s <- nc_station_names(nc_h2s)
SWE_h2s_all  <- ncvar_get(nc_h2s, "hs2swe_swe_total")   # [station × time], mm
nc_close(nc_h2s)

t_match <- match(dates_mag, dates_h2s)   # align HS2SWE onto the Mag25 time axis

stations <- setdiff(stations_mag, exclude_stations)

hyd_all <- hyd_year_of(dates_mag)
seasons <- sort(unique(hyd_all))
seasons <- seasons[vapply(seasons, function(y) sum(hyd_all == y) >= 300L, logical(1))]
message(sprintf("%d stations × %d seasons (%s)", length(stations), length(seasons),
                paste(sprintf("%d_%d", seasons, seasons + 1), collapse = ", ")))

# ─────────────────────────────────────────────────────────────────────────────
# Batch loop
# ─────────────────────────────────────────────────────────────────────────────
metric_rows <- list()
ts_rows     <- list()

for (stn in stations) {
  si_mag <- match(stn, stations_mag)
  si_h2s <- match(stn, stations_h2s)

  for (hy in seasons) {
    season_lbl <- sprintf("%d_%d", hy, hy + 1)
    idx        <- which(hyd_all == hy)

    dates   <- dates_mag[idx]
    hs_obs  <- HS_all[si_mag, idx]
    swe_obs <- SWE_all[si_mag, idx]
    hnw_obs <- HNW_all[si_mag, idx]
    swe_h2s <- if (is.na(si_h2s)) rep(NA_real_, length(idx))
               else SWE_h2s_all[si_h2s, t_match[idx]]

    swe_ds_def <- run_dsnow(dates, hs_obs, model_opts_default,  dyn_rho_max)
    swe_ds_new <- run_dsnow(dates, hs_obs, model_opts_best_new, dyn_rho_max)
    if (is.null(swe_ds_def)) swe_ds_def <- rep(NA_real_, length(idx))
    if (is.null(swe_ds_new)) swe_ds_new <- rep(NA_real_, length(idx))

    runs <- list(
      dSnow_default  = swe_ds_def,
      dSnow_best_new = swe_ds_new,
      HS2SWE         = swe_h2s
    )

    for (mod_name in names(runs)) {
      swe_mod <- runs[[mod_name]]
      m_swe   <- metrics(swe_mod,               swe_obs)
      m_hnw   <- metrics(hnw_from_swe(swe_mod), hnw_obs)
      metric_rows[[length(metric_rows) + 1L]] <- data.frame(
        station      = stn,
        season       = season_lbl,
        model        = mod_name,
        swe_rmse     = m_swe$rmse,
        swe_rel_bias = m_swe$rel_bias,
        swe_r2       = m_swe$r2,
        swe_n        = m_swe$n,
        hnw_rmse     = m_hnw$rmse,
        hnw_rel_bias = m_hnw$rel_bias,
        hnw_r2       = m_hnw$r2,
        hnw_n        = m_hnw$n,
        stringsAsFactors = FALSE
      )
    }

    ts_rows[[length(ts_rows) + 1L]] <- data.frame(
      date                = dates,
      station             = stn,
      season              = season_lbl,
      hs_obs_m            = hs_obs,
      swe_obs_mm          = swe_obs,
      hnw_obs_mm          = hnw_obs,
      swe_dsnow_default   = swe_ds_def,
      swe_dsnow_best_new  = swe_ds_new,
      swe_hs2swe          = swe_h2s,
      stringsAsFactors    = FALSE
    )
  }
  message("Done: ", stn)
}

metrics_df <- do.call(rbind, metric_rows)
ts_df      <- do.call(rbind, ts_rows)

write.csv(metrics_df, metrics_csv,    row.names = FALSE)
write.csv(ts_df,      timeseries_csv, row.names = FALSE)

message(sprintf("Wrote %s  (%d rows)", metrics_csv,    nrow(metrics_df)))
message(sprintf("Wrote %s  (%d rows)", timeseries_csv, nrow(ts_df)))
