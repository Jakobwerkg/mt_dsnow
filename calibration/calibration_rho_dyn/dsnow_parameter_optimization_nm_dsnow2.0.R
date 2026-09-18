##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION — NELDER-MEAD (dynamic rho.max)
##############################################################################
#
# Objective (minimized):
#   score = w_swe_nrmse * NRMSE_SWE
#         + w_rho_nrmse * NRMSE_rho
#         + w_swe_bias  * NBIAS_SWE
#         + w_rho_bias  * NBIAS_rho
#         + w_kge_swe   * (1 - KGE_SWE)
#         + w_kge_rho   * (1 - KGE_rho)
#
# Metrics:
#   NRMSE = RMSE / mean(obs)
#   NBIAS = |mean(mod - obs)| / mean(obs)
#   KGE   = 1 - sqrt((r-1)^2 + (beta-1)^2 + (alpha-1)^2)   [Gupta et al. 2009]
#
# Model variant: dyn_rho_max = TRUE.  The constant rho.max of the original
# deltaSNOW model is replaced by an age-dependent maximum density
#   rho_max(age) = atan(sigma * (age - mu)) * (rho_h - rho_l) / pi + (rho_h + rho_l) / 2
# so ten parameters are calibrated instead of seven.
#
# The optimiser searches the scaled parameter space (par / par_scale) so that
# all parameters are of order 1; minimize_score() multiplies back by
# par_scale to recover physical units.
#
# Command-line usage (all weights optional, positional):
#   Rscript script.R <w_swe_nrmse> <w_rho_nrmse> <w_swe_bias> <w_rho_bias>
#                    <w_kge_swe> <w_kge_rho>
#
# This script shares its data handling, objective function, parameter bounds
# and output format with the Differential Evolution script next to it and with
# the corresponding pair for the other dataset; only the blocks marked
# "DATASET-SPECIFIC" and "OPTIMIZER-SPECIFIC" differ.
#
##############################################################################

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)


# -----------------------------------------------------------------------------
# Project root
# Located through the .projectroot marker file, so the repository can live
# anywhere.  Works both for `Rscript path/to/script.R` and interactive use.
# -----------------------------------------------------------------------------
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

# =============================================================================
# DATASET-SPECIFIC SETTINGS
# =============================================================================

# -----------------------------------------------------------------------------
# Which observations to calibrate on.  Change the default here, or override
# without editing this file:  DSNOW_SUBSET=sp_b2000 Rscript <script>
#
#   sp_all      all 18 Alpsolut SNOWPACK stations
#   sp_rg       SNOWPACK rain-gauge subset
#   sp_b2000    SNOWPACK stations below 2000 m a.s.l.
#   win21       Winkler et al. (2021) observation stations
#
# The prepared observations are read from the subset's existing working
# directory optimisation_output/<subset>/data (SNOWPACK subsets are written
# there by calibration_snowpack/prepare_snowpack_data.R).  Results of this
# dynamic-rho.max variant go to calibration_rho_dyn/dyn_rho_max_res/<subset>/
# (de_res/ and nm_res/), separate from the static-rho.max calibrations.
# -----------------------------------------------------------------------------

SUBSET <- Sys.getenv("DSNOW_SUBSET", "sp_all")

SUBSETS <- list(
  sp_all   = list(obs_file = "d_obs_SNOWPACK.rda", hs_scale = 1,
                  drop = character(0),                  fit_only_exclude = character(0)),
  sp_rg    = list(obs_file = "d_obs_SNOWPACK.rda", hs_scale = 1,
                  drop = character(0),                  fit_only_exclude = character(0)),
  sp_b2000 = list(obs_file = "d_obs_SNOWPACK.rda", hs_scale = 1,
                  drop = character(0),                  fit_only_exclude = character(0)),
  win21    = list(obs_file = "H_SWE_obs.Rda",      hs_scale = 1 / 100,   # Hobs in cm -> m
                  drop = c("kuehtai", "Weissfluhjoch"), fit_only_exclude = "Sta.Maria")
)
if (!SUBSET %in% names(SUBSETS)) {
  stop("Unknown DSNOW_SUBSET '", SUBSET, "'. Valid: ",
       paste(names(SUBSETS), collapse = ", "))
}
DS <- SUBSETS[[SUBSET]]

OBS_DIR  <- file.path(ROOT, "calibration", "optimisation_output", SUBSET, "data")
RES_DIR <- file.path(ROOT, "calibration", "calibration_rho_dyn", "dyn_rho_max_res", SUBSET)

OBS_FILE <- DS$obs_file       # prepared observations inside OBS_DIR

HS_SCALE <- DS$hs_scale       # factor bringing Hobs to metres

DYN_RHO_MAX <- TRUE           # this variant: maximum bulk density varies with layer age

# Stations held out of the calibration (only Win21 holds any out):
#   DROP_STATIONS    removed from the data entirely (never fitted, never validated)
#   FIT_ONLY_EXCLUDE kept as independent validation stations, never fitted
DROP_STATIONS    <- DS$drop
FIT_ONLY_EXCLUDE <- DS$fit_only_exclude

cat(sprintf("Subset: %s\n  observations: %s\n  results     : %s\n", SUBSET, OBS_DIR, RES_DIR))

# =============================================================================
# CONFIGURATION
# =============================================================================

# --- Default weights (overridden by command-line args if provided) ---
WEIGHT_SWE_NRMSE <- 1.0
WEIGHT_RHO_NRMSE <- 0.0
WEIGHT_SWE_NBIAS <- 0.0
WEIGHT_RHO_NBIAS <- 0.0
WEIGHT_SWE_KGE   <- 0.0
WEIGHT_RHO_KGE   <- 0.0

# --- Command-line override ---
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) WEIGHT_SWE_NRMSE <- as.numeric(.args[1])
if (length(.args) >= 2) WEIGHT_RHO_NRMSE <- as.numeric(.args[2])
if (length(.args) >= 3) WEIGHT_SWE_NBIAS <- as.numeric(.args[3])
if (length(.args) >= 4) WEIGHT_RHO_NBIAS <- as.numeric(.args[4])
if (length(.args) >= 5) WEIGHT_SWE_KGE   <- as.numeric(.args[5])
if (length(.args) >= 6) WEIGHT_RHO_KGE   <- as.numeric(.args[6])

cat(sprintf(
  "Weights: SWE_NRMSE=%.3f  RHO_NRMSE=%.3f  SWE_NBIAS=%.3f  RHO_NBIAS=%.3f  SWE_KGE=%.3f  RHO_KGE=%.3f\n",
  WEIGHT_SWE_NRMSE, WEIGHT_RHO_NRMSE, WEIGHT_SWE_NBIAS,
  WEIGHT_RHO_NBIAS, WEIGHT_SWE_KGE,   WEIGHT_RHO_KGE
))

# --- Season settings ---
SEASON_START <- "-08-01"   # hydrological year start (Aug 1)
SEASON_END   <- "-07-31"   # hydrological year end   (Jul 31)
SEASON_MONTH <- 8          # month used to assign season label

# --- Model settings ---
EPS <- 1e-6                # minimum snow depth [m] for density calculation

# --- deltaSNOW parameters (dynamic rho.max variant), in the order the model expects them ---
# The S-curve parameters sigma, mu, rho_h, rho_l replace the constant rho.max.
PAR_NAMES <- c("sigma", "mu", "rho_h", "rho_l",
               "rho.null", "c.ov", "k.ov", "k", "tau", "eta.null")

# =============================================================================
# OPTIMIZER-SPECIFIC SETTINGS
# =============================================================================

# Nelder-Mead (via optimx, followed by BFGS) runs unconstrained, like the
# reference implementation; the parameter bounds below are therefore only a
# sanity check on the start values, not a constraint on the search.

# --- Where results are written (relative to RES_DIR) ---
LOG_SUBDIR    <- "nm_res"
RESULT_PREFIX <- "opt_results__"

# =============================================================================
# DATA
# =============================================================================

d_obs <- get(load(
  file.path(OBS_DIR, OBS_FILE)
))

# Stations removed from the data entirely (see DATASET-SPECIFIC SETTINGS)
d_obs <- d_obs[setdiff(names(d_obs), DROP_STATIONS)]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Assign a hydrological season label (year in which Aug 1 falls)
set_season <- function(date, start_month = SEASON_MONTH) {
  date <- as.Date(date)
  ifelse(month(date) < start_month, year(date) - 1L, year(date))
}

# Root Mean Square Error
rmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((mod[ok] - obs[ok])^2))
}

# Normalized RMSE: RMSE / mean(obs)
nrmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((mod[ok] - obs[ok])^2)) / mean(obs[ok])
}

# Normalized Absolute Bias: |mean(mod - obs)| / mean(obs)
nbias <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  abs(mean(mod[ok] - obs[ok])) / mean(obs[ok])
}

# Kling-Gupta Efficiency (Gupta et al. 2009); perfect = 1
# Enters the score as (1 - KGE) so that 0 = perfect, larger = worse
kge <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (sum(ok) < 3) return(NA_real_)
  r     <- cor(obs[ok], mod[ok])
  beta  <- mean(mod[ok]) / mean(obs[ok])
  alpha <- sd(mod[ok])   / sd(obs[ok])
  1 - sqrt((r - 1)^2 + (beta - 1)^2 + (alpha - 1)^2)
}

# =============================================================================
# COMBINED SCORE
# =============================================================================

combined_score <- function(df) {
  swe_obs <- df$swe_obs
  swe_mod <- df$swe_mod
  hs      <- df$hs

  # SWE metrics
  nrmse_swe <- nrmse(swe_obs, swe_mod)
  nbias_swe <- nbias(swe_obs, swe_mod)
  kge_swe   <- kge(swe_obs, swe_mod)

  # Bulk density metrics — only where snow is present
  snow      <- is.finite(hs) & hs > EPS
  rho_obs   <- ifelse(snow, swe_obs / hs, NA_real_)
  rho_mod   <- ifelse(snow, swe_mod / hs, NA_real_)
  nrmse_rho <- nrmse(rho_obs, rho_mod)
  nbias_rho <- nbias(rho_obs, rho_mod)
  kge_rho   <- kge(rho_obs, rho_mod)

  # Weighted score (KGE enters as 1-KGE: 0 = perfect)
  score <- WEIGHT_SWE_NRMSE * nrmse_swe      +
           WEIGHT_RHO_NRMSE * nrmse_rho      +
           WEIGHT_SWE_NBIAS * nbias_swe      +
           WEIGHT_RHO_NBIAS * nbias_rho      +
           WEIGHT_SWE_KGE   * (1 - kge_swe) +
           WEIGHT_RHO_KGE   * (1 - kge_rho)

  attr(score, "metrics") <- list(
    rmse_swe  = rmse(swe_obs, swe_mod),
    rmse_rho  = rmse(rho_obs, rho_mod),
    nrmse_swe = nrmse_swe,
    nrmse_rho = nrmse_rho,
    bias_swe  = mean(swe_mod - swe_obs, na.rm = TRUE),
    nbias_swe = nbias_swe,
    nbias_rho = nbias_rho,
    kge_swe   = kge_swe,
    kge_rho   = kge_rho,
    mean_swe  = mean(swe_obs[is.finite(swe_obs) & is.finite(swe_mod)], na.rm = TRUE),
    mean_rho  = mean(rho_obs[is.finite(rho_obs) & is.finite(rho_mod)], na.rm = TRUE)
  )

  return(score)
}

# =============================================================================
# TRAIN / VALIDATION SPLIT  (even counter -> fit, odd counter -> validation)
# =============================================================================

d_obs_fit <- list()
d_obs_val <- list()

for (station in names(d_obs)) {
  cat("Processing station:", station, "\n")

  d     <- d_obs[[station]]
  years <- unique(year(index(d)))

  fit_list <- list()
  val_list <- list()
  counter  <- 1L

  for (y in years[seq_len(length(years) - 1L)]) {

    winter <- subset(
      d,
      index(d) >= as.Date(paste0(y,     SEASON_START)) &
      index(d) <  as.Date(paste0(y + 1, SEASON_END)) + 1
    )

    if (nrow(winter) < 200) {
      cat("  skipping", y, "— fewer than 200 days\n"); next
    }
    if (nrow(winter) < 365 &&
        (as.numeric(winter$Hobs[1])            > 0.05 ||
         as.numeric(winter$Hobs[nrow(winter)]) > 0.05)) {
      cat("  skipping", y, "— incomplete winter with snow at edge\n"); next
    }

    if (counter %% 2 == 0) fit_list[[length(fit_list) + 1]] <- winter
    else                    val_list[[length(val_list) + 1]] <- winter
    counter <- counter + 1L
  }

  d_obs_fit[[station]] <- if (length(fit_list) > 0) do.call(rbind, fit_list) else NULL
  d_obs_val[[station]] <- if (length(val_list) > 0) do.call(rbind, val_list) else NULL
}

# =============================================================================
# CONVERT ZOO TO TIBBLE
# =============================================================================

zoo_to_tibble <- function(obs_list, exclude = character(0)) {
  out <- list()
  for (station in names(obs_list)) {
    if (station %in% exclude) {
      cat("Skipping", station, "— held out of this set\n"); next
    }

    x <- obs_list[[station]]
    if (is.null(x) || length(x) == 0) next

    dat <- tryCatch(as_tibble(coredata(x)), error = function(e) NULL)
    if (is.null(dat) || ncol(dat) == 0) next
    if (ncol(dat) == 2 && all(is.na(colnames(dat)))) colnames(dat) <- c("Hobs", "SWEobs")
    if (!all(c("Hobs", "SWEobs") %in% colnames(dat))) {
      cat("Skipping", station, "— missing Hobs or SWEobs\n"); next
    }

    out[[station]] <- tibble(
      date    = as.Date(index(x)),
      name    = station,
      hs      = dat$Hobs * HS_SCALE,   # -> m
      swe_obs = dat$SWEobs,
      block   = set_season(index(x))
    )
  }
  bind_rows(out)
}

# FIT_ONLY_EXCLUDE stations stay in the validation set but are never fitted.
d_obs_fit_tibble <- zoo_to_tibble(d_obs_fit, exclude = FIT_ONLY_EXCLUDE)
d_obs_val_tibble <- zoo_to_tibble(d_obs_val)

cat("\nFit set       :", length(unique(d_obs_fit_tibble$name)), "stations,",
    nrow(d_obs_fit_tibble), "days\n")
cat("Validation set:", length(unique(d_obs_val_tibble$name)), "stations,",
    nrow(d_obs_val_tibble), "days\n")

# =============================================================================
# OBJECTIVE FUNCTION
# =============================================================================

minimize_score <- function(par, data, scale, verbose = FALSE) {

  result <- tryCatch({

    par_real   <- setNames(par * scale, PAR_NAMES)
    model_opts <- as.list(par_real)

    if (verbose) {
      cat("scaled   =", paste(paste0(PAR_NAMES, "=", round(par,      6)), collapse = ", "), "\n")
      cat("unscaled =", paste(paste0(PAR_NAMES, "=", round(par_real, 6)), collapse = ", "), "\n")
    }

    # The S-curve must rise with layer age, i.e. rho_l < rho_h
    if (par_real[["rho_l"]] >= par_real[["rho_h"]]) {
      cat("rho_l >= rho_h — returning penalty.\n"); return(1e12)
    }

    # Run model per station and season block in parallel
    station_results <- foreach(
      station   = unique(data$name),
      .packages = c("dplyr", "tidyr", "nixmass"),
      .combine  = bind_rows
    ) %dopar% {

      data_station <- filter(data, name == station)

      block_results <- lapply(unique(data_station$block), function(b) {

        df_block   <- filter(data_station, block == b)
        full_dates <- tibble(date = seq(min(df_block$date), max(df_block$date), by = "1 day"))

        joined <- df_block %>%
          select(date, hs, swe_obs) %>%
          right_join(full_dates, by = "date") %>%
          arrange(date)

        swe_mod <- tryCatch(
          joined %>%
            select(date, hs) %>%
            mutate(date = as.character(date)) %>%
            nixmass::swe.delta.snow(
              model_opts  = model_opts,
              dyn_rho_max = DYN_RHO_MAX
            ),
          error = function(e) rep(NA_real_, nrow(joined))
        )

        drop_na(mutate(joined, swe_mod = swe_mod))
      })

      bind_rows(block_results)
    }

    if (nrow(station_results) == 0) {
      cat("No valid output — returning penalty.\n"); return(1e12)
    }

    score_obj <- combined_score(station_results)
    metrics   <- attr(score_obj, "metrics")
    score     <- as.numeric(score_obj)

    if (!is.finite(score)) {
      cat("Non-finite score — returning penalty.\n"); return(1e12)
    }

    if (verbose && !is.null(metrics)) {
      cat(sprintf(
        "| bias=%.4f | RMSE_SWE=%.4f | NRMSE_SWE=%.4f | NBIAS_SWE=%.4f | NBIAS_rho=%.4f | KGE_SWE=%.4f | KGE_rho=%.4f | NRMSE_rho=%.4f | score=%.4f\n",
        metrics$bias_swe,  metrics$rmse_swe,  metrics$nrmse_swe,
        metrics$nbias_swe, metrics$nbias_rho,
        metrics$kge_swe,   metrics$kge_rho,
        metrics$nrmse_rho, score
      ))
    } else {
      cat("Score =", round(score, 4), "\n")
    }

    return(score)

  }, error = function(e) { cat("Error:", e$message, "\n"); 1e12 })

  if (length(result) != 1 || !is.numeric(result)) return(1e12)
  result
}

# =============================================================================
# PARAMETER BOUNDS AND START VALUES
# =============================================================================

# Physical (unscaled) parameter bounds.  The parameters shared with the static
# model keep the Schellander reference ranges of the other calibration scripts;
# the S-curve parameters get ranges around the nixmass defaults.
par_lower <- c(sigma = 0.005, mu = 10,  rho_h = 400, rho_l = 200,
               rho.null = 50,  c.ov = 1e-6, k.ov = 0.01, k = 0.01, tau = 0.01, eta.null = 1e6)

par_upper <- c(sigma = 0.3,   mu = 300, rho_h = 700, rho_l = 500,
               rho.null = 200, c.ov = 1e-3, k.ov = 1,    k = 0.2,  tau = 0.2,  eta.null = 2e7)
# k.ov: nixmass stops with "'k.ov' must be < 1", so an upper bound above 1
# (as in the static scripts) only wastes evaluations on the penalty value.

# Physical start values: nixmass defaults for dyn_rho_max = TRUE
par_delta <- c(
  sigma    = 0.03,
  mu       = 80,
  rho_h    = 600,
  rho_l    = 380,
  rho.null = 80.73706,
  c.ov     = 0.0005170964,
  k.ov     = 0.3782312,
  k        = 0.029297,
  tau      = 0.02356521,
  eta.null = 8543502
)

# Per-parameter scaling so the optimizer works on values of order 1.
# minimize_score multiplies back by `scale` to recover physical units.
par_scale <- c(sigma = 0.1, mu = 100, rho_h = 1000, rho_l = 1000,
               rho.null = 1000, c.ov = 0.001, k.ov = 1, k = 0.1, tau = 0.1, eta.null = 1e7)

stopifnot(
  identical(names(par_lower), PAR_NAMES),
  identical(names(par_upper), PAR_NAMES),
  identical(names(par_delta), PAR_NAMES),
  identical(names(par_scale), PAR_NAMES)
)

par_start <- par_delta / par_scale

cat("\nParameter bounds (physical):\n")
print(rbind(lower = par_lower, upper = par_upper))
cat("In-bounds check:\n")
print(ifelse(par_lower <= par_delta & par_delta <= par_upper, "in bounds", "ERROR"))
cat("\nStarting parameters (scaled):\n");   print(par_start)
cat("\nStarting parameters (unscaled):\n"); print(par_start * par_scale)

# =============================================================================
# PARALLEL SETUP
# =============================================================================

nc <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

parallel::clusterExport(cl, c(
  "WEIGHT_SWE_NRMSE", "WEIGHT_RHO_NRMSE",
  "WEIGHT_SWE_NBIAS", "WEIGHT_RHO_NBIAS",
  "WEIGHT_SWE_KGE",   "WEIGHT_RHO_KGE",
  "EPS", "DYN_RHO_MAX", "PAR_NAMES",
  "rmse", "nrmse", "nbias", "kge", "combined_score"
))

# =============================================================================
# TEST AT INITIAL GUESS
# =============================================================================

cat("\nTesting objective function at start values...\n")
cat("Initial score =",
    minimize_score(par_start, d_obs_fit_tibble, scale = par_scale, verbose = TRUE), "\n")

# =============================================================================
# OPTIMIZATION
# =============================================================================

cat("\nStarting optimization (optimx default: Nelder-Mead -> BFGS, follow.on)...\n")
opt <- optimx(
  par     = par_start,
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  scale   = par_scale,
  verbose = TRUE,
  control = list(trace = 1, follow.on = TRUE)
)

opt_df      <- as.data.frame(opt)
best_idx    <- which.min(opt_df$value)
best_scaled <- setNames(as.numeric(opt_df[best_idx, names(par_start), drop = TRUE]), PAR_NAMES)
best_value  <- opt_df$value[best_idx]

# =============================================================================
# RESULTS
# =============================================================================

best_par <- best_scaled * par_scale   # back to physical units

cat("\nOptimization complete!\n")
cat("Best score:", best_value, "\n")
cat("Best parameters (scaled):\n");   print(best_scaled)
cat("Best parameters (unscaled):\n"); print(best_par)

# Print as R call (pydeltasnow has no dynamic-rho.max variant)
fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 10)
cat(sprintf(
  "\nswe_mod <- nixmass::swe.delta.snow(\n  data,\n  model_opts = list(\n%s\n  ),\n  dyn_rho_max = %s\n)\n",
  paste0("    ", PAR_NAMES, " = ", vapply(best_par, fmt, ""), collapse = ",\n"),
  DYN_RHO_MAX
))

# =============================================================================# SAVE RESULTS

fmt_tag <- function(x) gsub("-", "m", gsub("\\.", "p",
  format(x, scientific = FALSE, trim = TRUE, digits = 6)))

weight_vals <- c(SWE_NRMSE = WEIGHT_SWE_NRMSE, RHO_NRMSE = WEIGHT_RHO_NRMSE,
                 SWE_NBIAS = WEIGHT_SWE_NBIAS,  RHO_NBIAS = WEIGHT_RHO_NBIAS,
                 SWE_KGE   = WEIGHT_SWE_KGE,    RHO_KGE   = WEIGHT_RHO_KGE)

weight_tag <- paste0(names(weight_vals), "_", fmt_tag(weight_vals), collapse = "__")

save_dir <- file.path(RES_DIR, LOG_SUBDIR)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0(RESULT_PREFIX, weight_tag, ".rds"))

saveRDS(list(
  opt        = opt,
  best_par   = best_par,
  best_value = best_value,
  weights    = weight_vals,
  subset      = SUBSET,
  dyn_rho_max = DYN_RHO_MAX,
  fit_data   = d_obs_fit_tibble,
  val_data   = d_obs_val_tibble
), file = save_file)

cat("\nOptimization finished. Results saved to:\n", save_file, "\n", sep = "")

# =============================================================================
# CLEANUP
# =============================================================================

parallel::stopCluster(cl)
