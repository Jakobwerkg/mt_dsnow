# run_full_validation.R
# -----------------------------------------------------------------------------
# Independent HNW & SWE validation of every calibrated parameter set.
#
# Reads all_summaries.csv (one row per calibration run, with the 7 ΔSnow
# parameters), runs nixmass::swe.delta.snow on the full Mag25 multi-station
# dataset for each row, derives HNW_mod = clip(diff(SWE_mod), 0), and computes
# HNW and SWE validation metrics. The input table is written back out with the
# metric columns appended.
#
# Metric / filtering logic mirrors HNW_validation_helper.compute_metrics_hnw_swe
# (used by Archive/HNW_validation_stats.ipynb):
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
#   <VAR>_RMSE | <VAR>_Bias | <VAR>_Rel_BIAS | <VAR>_R2 | <VAR>_N
#
# Usage:
#   Rscript HNW_validation/full_R_Validation/run_full_validation.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
  library(foreach)
  library(doParallel)
})

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR    <- "/Users/jakobwerkgarner/code/mt_dsnow"
SUMMARY_CSV <- file.path(BASE_DIR, "calibration/AA_opt_out/workflow_helper/all_summaries.csv")
MAG25_NC    <- file.path(BASE_DIR, "calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc")
EXCLUDE_STATIONS <- c("Weisfluh_Joch")

# The 7 ΔSnow parameters expected in the summary CSV (underscore notation).
PAR_COLS <- c("rho_max", "rho_null", "eta_null", "k", "tau", "c_ov", "k_ov")

# Dynamic rho_max parameterisation?  FALSE = static (Winkler 2021 style, default).
# When TRUE, "_dyn_rho_max" is appended to the output filename automatically.
# NOTE: dynamic mode also expects sigma/mu/rho_h/rho_l, which are NOT stored in
# the summary CSV — only enable this if those parameters are supplied elsewhere.
DYN_RHO_MAX <- FALSE

OUT_CSV <- file.path(
  BASE_DIR, "HNW_validation/full_R_Validation",
  paste0("all_summaries_validated_R",
         if (DYN_RHO_MAX) "_dyn_rho_max" else "", ".csv")
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
hyd_year_of <- function(d) ifelse(month(d) > 8L, year(d), year(d) - 1L)

# Translate underscore parameter names to nixmass dot-notation.
to_model_opts <- function(row) {
  list(
    rho.max  = as.numeric(row[["rho_max"]]),
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

# ─────────────────────────────────────────────────────────────────────────────
# Load Mag25 (stations as a dimension, time embedded as numeric dim values)
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(MAG25_NC), file.exists(SUMMARY_CSV))

nc_in <- nc_open(MAG25_NC)
time_dim   <- nc_in$dim[["time"]]
time_units <- time_dim$units
origin_dt  <- as.Date(substr(sub("days since ", "", time_units), 1, 10))
dates_all  <- origin_dt + as.integer(time_dim$vals)

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
# Validate ONE parameter set -> named metric vector
# ─────────────────────────────────────────────────────────────────────────────
validate_one <- function(model_opts, dyn_rho_max) {
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
# Loop over all parameter sets (parallel over rows)
# ─────────────────────────────────────────────────────────────────────────────
opt <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE, check.names = FALSE)
message(sprintf("Loaded %d parameter sets from %s", nrow(opt), basename(SUMMARY_CSV)))

has_par <- !is.na(opt$rho_max) & !is.na(opt$rho_null) & !is.na(opt$eta_null) &
           !is.na(opt$k) & !is.na(opt$tau) & !is.na(opt$c_ov) & !is.na(opt$k_ov)
message(sprintf("Rows with full parameters: %d / %d", sum(has_par), nrow(opt)))

METRIC_NAMES <- c("SWE_RMSE", "SWE_Bias", "SWE_Rel_BIAS", "SWE_R2", "SWE_N",
                  "HNW_RMSE", "HNW_Bias", "HNW_Rel_BIAS", "HNW_R2", "HNW_N")

nc_workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(nc_workers)
registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(cl, c(
  "HS_all", "SWE_obs_all", "HNW_obs_all", "dates_all", "station_names",
  "hyd_years_all", "winter_years", "Ns", "Nt", "keep_station", "season_mask",
  "run_dsnow", "calc_metrics", "validate_one", "to_model_opts", "DYN_RHO_MAX"
), envir = environment())
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(nixmass))))

message(sprintf("Running validation on %d cores...", nc_workers))
t0 <- Sys.time()

metrics_mat <- foreach(i = seq_len(nrow(opt)), .combine = rbind,
                       .packages = "nixmass") %dopar% {
  if (!has_par[i]) return(setNames(rep(NA_real_, length(METRIC_NAMES)), METRIC_NAMES))

  validate_one(to_model_opts(opt[i, ]), DYN_RHO_MAX)
}

metrics_df <- as.data.frame(metrics_mat, row.names = FALSE)
names(metrics_df) <- METRIC_NAMES
metrics_df$SWE_N <- as.integer(metrics_df$SWE_N)
metrics_df$HNW_N <- as.integer(metrics_df$HNW_N)

message(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ─────────────────────────────────────────────────────────────────────────────
# Merge + save
# ─────────────────────────────────────────────────────────────────────────────
opt_out <- cbind(opt[, setdiff(names(opt), METRIC_NAMES), drop = FALSE], metrics_df)

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(opt_out, OUT_CSV, row.names = FALSE)

message("\n=== Validation summary (head) ===")
show_cols <- intersect(c("subset", "dataset", "phase", "algorithm",
                         "SWE_RMSE", "SWE_Rel_BIAS", "SWE_R2",
                         "HNW_RMSE", "HNW_Rel_BIAS", "HNW_R2"), names(opt_out))
print(utils::head(opt_out[, show_cols], 10), row.names = FALSE)

message("\nWrote: ", OUT_CSV)
