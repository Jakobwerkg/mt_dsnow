# run_full_validation_rho_dyn.R
# -----------------------------------------------------------------------------
# Independent HNW & SWE validation of every calibrated parameter set of the
# dynamic-rho.max deltaSNOW variant (dyn_rho_max = TRUE, 10 parameters).
#
# Reads all_summaries_dyn_rho_max.csv (one row per calibration run, written by
# collect_rho_dyn_results.R), runs nixmass::swe.delta.snow on the full Mag25
# multi-station dataset for each row, derives HNW_mod = clip(diff(SWE_mod), 0),
# and computes HNW and SWE validation metrics.  The input table is written back
# out with the metric columns appended, and for every parameter set a
# validation NetCDF (daily HNW_mod, SWE_mod, HNW_obs, SWE_obs; same layout as
# hnw_validation/run_single_dsnow_validation.R, readable by
# hnw_validation/plot_validation_results.ipynb) is written to NC_DIR.
#
# Metric / filtering logic mirrors hnw_validation/full_validation/
# run_full_validation.R (and hnw_validation_helper.compute_metrics_hnw_swe):
#   * obs = x, mod = y, residual = mod - obs
#   * RMSE      = sqrt(mean(res^2))
#   * Bias      = mean(res)                     (mean model - obs)
#   * Rel_BIAS  = sum(res) / sum(obs)           (PBIAS)
#   * R2        = 1 - SS_res / SS_tot           (Nash-Sutcliffe)
#   * N         = number of valid obs-mod pairs
#   * Weisfluh_Joch excluded; obs >= 0 required; pairs must be finite.
#   * HNW restricted to the snow season (Nov 1 - Apr 30); SWE uses the full year.
#
# Output columns appended (per variable SWE / HNW):
#   <VAR>_RMSE | <VAR>_Bias | <VAR>_Rel_BIAS | <VAR>_R2 | <VAR>_N | nc_file
#
# Usage:
#   Rscript calibration/calibration_rho_dyn/collect_rho_dyn_results.R      # first
#   Rscript calibration/calibration_rho_dyn/run_full_validation_rho_dyn.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})
# The parallel loop below drives the cluster through future:: (not attached).
if (!requireNamespace("future", quietly = TRUE))
  stop("Package 'future' is needed for the parallel loop: install.packages('future')")

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


# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
RES_DIR     <- file.path(ROOT, "calibration/calibration_rho_dyn/dyn_rho_max_res")
SUMMARY_CSV <- file.path(RES_DIR, "/all_summaries_dyn_rho_max.csv")
MAG25_NC    <- file.path(ROOT, "calibration/calibration_data/raw_data/mag25/slf_dataset/Mag25_all.nc")
EXCLUDE_STATIONS <- c("Weisfluh_Joch")

# The 10 parameters of the dynamic-rho.max model expected in the summary CSV
# (underscore notation).  sigma/mu/rho_h/rho_l replace rho_max.
PAR_COLS <- c("sigma", "mu", "rho_h", "rho_l",
              "rho_null", "eta_null", "k", "tau", "c_ov", "k_ov")

# This script validates the dynamic variant only.
DYN_RHO_MAX <- TRUE

OUT_DIR <- file.path(RES_DIR, "validation")
OUT_CSV <- file.path(OUT_DIR, "all_summaries_validated_dyn_rho_max.csv")

# One validation NetCDF per parameter set (~3 MB each):
#   <NC_DIR>/<subset>_phase=<phase>_algorithm=<algorithm>.nc
# Set SAVE_NC <- FALSE to only compute the metrics table.
SAVE_NC <- TRUE
NC_DIR  <- file.path(OUT_DIR, "nc")

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
hyd_year_of <- function(d) ifelse(month(d) > 8L, year(d), year(d) - 1L)

# Translate underscore parameter names to nixmass dot-notation.
to_model_opts <- function(row) {
  list(
    sigma    = as.numeric(row[["sigma"]]),
    mu       = as.numeric(row[["mu"]]),
    rho_h    = as.numeric(row[["rho_h"]]),
    rho_l    = as.numeric(row[["rho_l"]]),
    rho.null = as.numeric(row[["rho_null"]]),
    c.ov     = as.numeric(row[["c_ov"]]),
    k.ov     = as.numeric(row[["k_ov"]]),
    k        = as.numeric(row[["k"]]),
    tau      = as.numeric(row[["tau"]]),
    eta.null = as.numeric(row[["eta_null"]])
  )
}

# Run ΔSnow on one station x one hydrological year (Sep-Aug). hs forced to 0
# at the segment ends so each winter starts and ends snow-free.
run_dsnow <- function(dates, hs_m, model_opts, dyn_rho_max) {
  hs <- pmax(as.numeric(hs_m), 0)
  hs[is.na(hs)] <- 0
  if (length(hs) == 0) return(NULL)
  hs[1]          <- 0
  hs[length(hs)] <- 0

  df <- data.frame(date = as.character(dates), hs = hs, stringsAsFactors = FALSE)
  out <- tryCatch(
    nixmass::swe.delta.snow(df,
                            model_opts  = model_opts,
                            dyn_rho_max = dyn_rho_max,
                            layers      = FALSE,
                            strict_mode = FALSE,
                            verbose     = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  if (is.list(out)) as.numeric(out$SWE) else as.numeric(out)
}

# RMSE / Bias / Rel_BIAS / R2 / N for one obs (x) vs mod (y) vector pair.
calc_metrics <- function(obs, mod) {
  ok  <- is.finite(obs) & is.finite(mod)
  obs <- obs[ok]; mod <- mod[ok]
  n   <- length(obs)
  if (n < 1) return(c(RMSE = NA_real_, Bias = NA_real_, Rel_BIAS = NA_real_,
                      R2 = NA_real_, N = 0))
  res    <- mod - obs
  rmse   <- sqrt(mean(res^2))
  bias   <- mean(res)
  relb   <- if (sum(obs) != 0) sum(res) / sum(obs) else NA_real_
  ss_res <- sum((obs - mod)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2     <- if (ss_tot != 0) 1 - ss_res / ss_tot else NA_real_
  c(RMSE = rmse, Bias = bias, Rel_BIAS = relb, R2 = r2, N = n)
}

# Console progress bar (base R only): percentage, done/total, elapsed, ETA.
# Redrawn in place when stderr is a terminal (or R is interactive); otherwise
# one line per update so redirected logs stay readable.
make_progress_bar <- function(total, width = 30L) {
  t_start  <- Sys.time()
  in_place <- interactive() || isatty(stderr())
  fmt_dur  <- function(secs) if (secs < 60) sprintf("%.0f s", secs) else
                                             sprintf("%.1f min", secs / 60)
  function(done) {
    if (total < 1L) return(invisible())
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    eta     <- if (done > 0L) elapsed / done * (total - done) else NA_real_
    n_fill  <- round(width * done / total)
    cat(sprintf("%s[%s%s] %3.0f%% (%d/%d) | elapsed %s | finished in %s%s",
                if (in_place) "\r" else "",
                strrep("=", n_fill), strrep(" ", width - n_fill),
                100 * done / total, done, total,
                fmt_dur(elapsed), if (is.na(eta)) "--" else fmt_dur(eta),
                if (in_place && done < total) "" else "\n"),
        file = stderr())
    flush(stderr())
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Load Mag25 (stations as a dimension, time embedded as numeric dim values)
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(MAG25_NC), file.exists(SUMMARY_CSV))

nc_in <- nc_open(MAG25_NC)
time_dim   <- nc_in$dim[["time"]]
time_raw   <- time_dim$vals
time_units <- time_dim$units
time_cal   <- if (!is.null(time_dim$calendar) && nchar(time_dim$calendar) > 0)
                time_dim$calendar else "standard"
origin_dt  <- as.Date(substr(sub("days since ", "", time_units), 1, 10))
dates_all  <- origin_dt + as.integer(time_raw)

station_names <- nc_in$dim[["station"]]$vals
avail_vars    <- names(nc_in$var)

HS_all      <- ncvar_get(nc_in, "HS")    # [station x time], m
SWE_obs_all <- ncvar_get(nc_in, "SWE")   # [station x time], mm
HNW_obs_all <- if ("HNW" %in% avail_vars) ncvar_get(nc_in, "HNW") else
               matrix(NA_real_, length(station_names), length(dates_all))
nc_close(nc_in)

Ns <- length(station_names)
Nt <- length(dates_all)

hyd_years_all <- hyd_year_of(dates_all)
winter_years  <- sort(unique(hyd_years_all))

# Pre-computed masks shared by every parameter set.
keep_station <- !(station_names %in% EXCLUDE_STATIONS)        # drop Weisfluh_Joch
season_mask  <- month(dates_all) >= 11 | month(dates_all) <= 4 # HNW: Nov-Apr

message(sprintf("Mag25: %d stations x %d days (%s - %s)",
                Ns, Nt, dates_all[1], dates_all[Nt]))

# ─────────────────────────────────────────────────────────────────────────────
# Write ONE validation NetCDF (same layout as run_single_dsnow_validation.R)
# ─────────────────────────────────────────────────────────────────────────────
write_validation_nc <- function(nc_file, SWE_mod, HNW_mod, model_opts, dyn_rho_max,
                                global_atts = list()) {
  dir.create(dirname(nc_file), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(nc_file)) file.remove(nc_file)

  nchar_max   <- max(nchar(station_names))
  dim_station <- ncdim_def("station", units = "", vals = seq_len(Ns),                  create_dimvar = FALSE)
  dim_time    <- ncdim_def("time",    units = time_units, vals = as.integer(time_raw), create_dimvar = TRUE)
  dim_nchar   <- ncdim_def("nchar",   units = "", vals = seq_len(nchar_max),           create_dimvar = FALSE)

  # "station_name" avoids a name clash with the "station" dimension (NC convention
  # reserves same-name 1-D variables as coordinate variables).
  v_stn     <- ncvar_def("station_name", "",   list(dim_nchar, dim_station), prec = "char")
  v_swe_mod <- ncvar_def("SWE_mod",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
  v_hnw_mod <- ncvar_def("HNW_mod",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
  v_swe_obs <- ncvar_def("SWE_obs",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
  v_hnw_obs <- ncvar_def("HNW_obs",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")

  nc_out <- nc_create(nc_file, vars = list(v_stn, v_swe_mod, v_hnw_mod, v_swe_obs, v_hnw_obs))
  on.exit(nc_close(nc_out), add = TRUE)   # inside a function on.exit is safe

  ncvar_put(nc_out, v_stn,     station_names)
  ncvar_put(nc_out, v_swe_mod, SWE_mod)
  ncvar_put(nc_out, v_hnw_mod, HNW_mod)
  ncvar_put(nc_out, v_swe_obs, SWE_obs_all)
  ncvar_put(nc_out, v_hnw_obs, HNW_obs_all)

  ncatt_put(nc_out, "time", "calendar", time_cal)

  defaults <- if (dyn_rho_max) {
    list(sigma = 0.03, mu = 80, rho_h = 600, rho_l = 380,
         rho.null = 80.73706, c.ov = 0.0005170964, k.ov = 0.3782312,
         k = 0.029297, tau = 0.02356521, eta.null = 8543502, timestep = 24)
  } else {
    list(rho.max = 401.2588, rho.null = 81.19417,
         c.ov = 0.0005104722, k.ov = 0.37856737, k = 0.02993175,
         tau = 0.02362476, eta.null = 8523356, timestep = 24)
  }
  effective_opts <- utils::modifyList(defaults, model_opts)
  param_str <- paste(mapply(function(k, v) sprintf("%s=%g", k, v),
                            names(effective_opts), unlist(effective_opts)),
                     collapse = "; ")
  ncatt_put(nc_out, 0, "dsnow_parameters", param_str)
  ncatt_put(nc_out, 0, "dyn_rho_max",      as.integer(dyn_rho_max))
  ncatt_put(nc_out, 0, "source",           "nixmass::swe.delta.snow")
  for (nm in names(global_atts)) {
    v <- global_atts[[nm]]
    if (length(v) == 1 && !is.na(v)) ncatt_put(nc_out, 0, nm, as.character(v))
  }
  invisible(nc_file)
}

# ─────────────────────────────────────────────────────────────────────────────
# Validate ONE parameter set -> named metric vector (+ optional NetCDF)
# ─────────────────────────────────────────────────────────────────────────────
validate_one <- function(model_opts, dyn_rho_max, nc_file = NA_character_,
                         global_atts = list()) {
  SWE_mod <- matrix(NA_real_, Ns, Nt)

  for (si in seq_len(Ns)) {
    if (!keep_station[si]) next
    hs_stn <- as.numeric(HS_all[si, ])
    for (wy in winter_years) {
      idx <- which(hyd_years_all == wy)
      if (length(idx) < 10L) next
      seg <- run_dsnow(dates_all[idx], hs_stn[idx], model_opts, dyn_rho_max)
      if (!is.null(seg) && length(seg) == length(idx)) SWE_mod[si, idx] <- seg
    }
  }

  # HNW_mod = diff(SWE_mod) along time, melt (negative) clipped to 0.
  HNW_mod <- matrix(NA_real_, Ns, Nt)
  for (si in seq_len(Ns)) {
    d <- diff(SWE_mod[si, ]); d[d < 0] <- 0
    HNW_mod[si, 2:Nt] <- d
  }

  if (!is.na(nc_file)) {
    write_validation_nc(nc_file, SWE_mod, HNW_mod, model_opts, dyn_rho_max, global_atts)
  }

  rk <- which(keep_station)

  # SWE: full year, obs >= 0
  swe_obs <- as.vector(SWE_obs_all[rk, ]); swe_mod <- as.vector(SWE_mod[rk, ])
  s_sel   <- is.finite(swe_obs) & is.finite(swe_mod) & swe_obs >= 0
  swe_m   <- calc_metrics(swe_obs[s_sel], swe_mod[s_sel])

  # HNW: snow season (Nov-Apr), obs >= 0
  cs      <- which(season_mask)
  hnw_obs <- as.vector(HNW_obs_all[rk, cs]); hnw_mod <- as.vector(HNW_mod[rk, cs])
  h_sel   <- is.finite(hnw_obs) & is.finite(hnw_mod) & hnw_obs >= 0
  hnw_m   <- calc_metrics(hnw_obs[h_sel], hnw_mod[h_sel])

  c(SWE_RMSE = swe_m[["RMSE"]], SWE_Bias = swe_m[["Bias"]],
    SWE_Rel_BIAS = swe_m[["Rel_BIAS"]], SWE_R2 = swe_m[["R2"]], SWE_N = swe_m[["N"]],
    HNW_RMSE = hnw_m[["RMSE"]], HNW_Bias = hnw_m[["Bias"]],
    HNW_Rel_BIAS = hnw_m[["Rel_BIAS"]], HNW_R2 = hnw_m[["R2"]], HNW_N = hnw_m[["N"]])
}

# ─────────────────────────────────────────────────────────────────────────────
# Parameter sets + one NetCDF name per row
# ─────────────────────────────────────────────────────────────────────────────
opt <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE, check.names = FALSE)
message(sprintf("Loaded %d parameter sets from %s", nrow(opt), basename(SUMMARY_CSV)))

missing_cols <- setdiff(PAR_COLS, names(opt))
if (length(missing_cols) > 0)
  stop("Summary CSV lacks parameter columns: ", paste(missing_cols, collapse = ", "))

has_par <- stats::complete.cases(opt[, PAR_COLS, drop = FALSE])
message(sprintf("Rows with full parameters: %d / %d", sum(has_par), nrow(opt)))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# <subset>_phase=<phase>_algorithm=<algorithm>.nc — the naming already used in
# hnw_validation/validation_data/.  Rows without a phase label get their
# weight combination instead.
weight_cols <- grep("^w_", names(opt), value = TRUE)
phase_tag <- ifelse(
  !is.na(opt$phase %||% NA) & nzchar(opt$phase %||% ""),
  opt$phase,
  apply(opt[, weight_cols, drop = FALSE], 1, function(w)
        paste(sub("^w_", "", weight_cols), gsub("\\.", "p", format(w, trim = TRUE)),
              sep = "-", collapse = "_"))
)
nc_stem <- sprintf("%s_phase=%s_algorithm=%s",
                   opt$subset %||% "unknown", phase_tag, opt$algorithm %||% "unknown")
nc_stem <- gsub("[^A-Za-z0-9._=+-]", "_", nc_stem)
nc_stem <- make.unique(nc_stem, sep = "_dup")
nc_files <- if (SAVE_NC) file.path(NC_DIR, paste0(nc_stem, ".nc")) else
            rep(NA_character_, nrow(opt))
nc_files[!has_par] <- NA_character_

# ─────────────────────────────────────────────────────────────────────────────
# Loop over all parameter sets (parallel over rows, with a live progress bar)
# ─────────────────────────────────────────────────────────────────────────────
METRIC_NAMES <- c("SWE_RMSE", "SWE_Bias", "SWE_Rel_BIAS", "SWE_R2", "SWE_N",
                  "HNW_RMSE", "HNW_Bias", "HNW_Rel_BIAS", "HNW_R2", "HNW_N")

if (SAVE_NC) dir.create(NC_DIR, showWarnings = FALSE, recursive = TRUE)

# Rows without a full parameter set keep an all-NA metric row and are skipped.
todo         <- which(has_par)
metrics_list <- rep(list(setNames(rep(NA_real_, length(METRIC_NAMES)), METRIC_NAMES)),
                    nrow(opt))

nc_workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
message(sprintf("Running validation of %d parameter sets on %d cores%s...",
                length(todo), nc_workers,
                if (SAVE_NC) sprintf(", NetCDFs -> %s", NC_DIR) else ""))
t0 <- Sys.time()

# Rows are handed to the PSOCK cluster one at a time through future::future()
# rather than one blocking foreach %dopar% call, so the master can collect
# finished rows as they come in and advance the progress bar after each one.
# future ships validate_one() and the ~2 MB of Mag25 arrays it needs with every
# task; negligible next to the minutes each task takes.
# NOTE: no on.exit(stopCluster) here. At the top level of a source()d script
# on.exit fires right after its own line (killing the cluster before use ->
# "invalid connection"); the cluster is stopped in the `finally` below instead.
cl       <- NULL
progress <- make_progress_bar(length(todo))
n_done   <- 0L
tryCatch({
  cl <- parallel::makeCluster(nc_workers)
  future::plan(future::cluster, workers = cl)
  progress(n_done)

  queue   <- todo
  running <- list()   # row index (as name) -> future still being evaluated
  while (length(queue) > 0L || length(running) > 0L) {
    # Hand out rows while a worker is idle.
    while (length(queue) > 0L && future::nbrOfFreeWorkers() > 0L) {
      i <- queue[1L]; queue <- queue[-1L]
      model_opts_i <- to_model_opts(opt[i, ])
      nc_file_i    <- nc_files[i]
      atts_i       <- list(subset      = opt$subset[i]      %||% NA,
                           phase       = opt$phase[i]       %||% NA,
                           algorithm   = opt$algorithm[i]   %||% NA,
                           source_path = opt$source_path[i] %||% NA)
      running[[as.character(i)]] <- future::future(
        validate_one(model_opts_i, DYN_RHO_MAX, nc_file = nc_file_i, global_atts = atts_i),
        packages = c("nixmass", "ncdf4"),
        seed = NULL,                 # deterministic model, skip the RNG check
        conditions = character(0)    # like foreach: drop worker warnings, keep errors
      )
    }
    # Collect whatever has finished and advance the bar.
    for (key in names(running)) {
      if (future::resolved(running[[key]])) {
        metrics_list[[as.integer(key)]] <- future::value(running[[key]])
        running[[key]] <- NULL
        n_done <- n_done + 1L
        progress(n_done)
      }
    }
    if (length(running) > 0L) Sys.sleep(0.5)   # poll interval
  }
}, finally = {
  if (n_done < length(todo)) cat("\n", file = stderr())   # bar left mid-line
  future::plan(future::sequential)
  if (!is.null(cl)) parallel::stopCluster(cl)
})

metrics_mat <- do.call(rbind, metrics_list)
metrics_df  <- as.data.frame(metrics_mat, row.names = FALSE)
names(metrics_df) <- METRIC_NAMES
metrics_df$SWE_N <- as.integer(metrics_df$SWE_N)
metrics_df$HNW_N <- as.integer(metrics_df$HNW_N)
# stored relative to the project root, so the table stays portable
metrics_df$nc_file <- ifelse(is.na(nc_files), NA_character_,
                             sub(paste0("^", ROOT, "/?"), "", nc_files))

message(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ─────────────────────────────────────────────────────────────────────────────
# Merge + save
# ─────────────────────────────────────────────────────────────────────────────
opt_out <- cbind(opt[, setdiff(names(opt), c(METRIC_NAMES, "nc_file")), drop = FALSE],
                 metrics_df)

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(opt_out, OUT_CSV, row.names = FALSE)

message("\n=== Validation summary (head) ===")
show_cols <- intersect(c("subset", "phase", "algorithm",
                         "SWE_RMSE", "SWE_Rel_BIAS", "SWE_R2",
                         "HNW_RMSE", "HNW_Rel_BIAS", "HNW_R2"), names(opt_out))
print(utils::head(opt_out[, show_cols], 10), row.names = FALSE)

message("\nWrote: ", OUT_CSV)
if (SAVE_NC) message(sprintf("Wrote %d NetCDF files to: %s", sum(!is.na(nc_files)), NC_DIR))
