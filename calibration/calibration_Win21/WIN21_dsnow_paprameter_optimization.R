##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION — NELDER-MEAD
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
# Command-line usage (all weights optional, positional):
#   Rscript script.R <w_swe_nrmse> <w_rho_nrmse> <w_swe_bias> <w_rho_bias> <w_kge_swe> <w_kge_rho>
#
##############################################################################

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)

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
SEASON_START    <- "-08-01"   # hydrological year start (Aug 1)
SEASON_END      <- "-07-31"   # hydrological year end   (Jul 31)
SEASON_MONTH    <- 8          # month used to assign season label

# --- Model settings ---
EPS <- 1e-6                   # minimum snow depth [m] for density calculation

# =============================================================================
# DATA
# =============================================================================

d_obs <- get(load(
  "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_Win21/data/H_SWE_obs.Rda"
))

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

zoo_to_tibble <- function(obs_list) {
  out <- list()
  for (station in names(obs_list)) {
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
      hs      = dat$Hobs,
      swe_obs = dat$SWEobs,
      block   = set_season(index(x))
    )
  }
  bind_rows(out)
}

d_obs_fit_tibble <- zoo_to_tibble(d_obs_fit)
d_obs_val_tibble <- zoo_to_tibble(d_obs_val)

# =============================================================================
# OBJECTIVE FUNCTION
# =============================================================================

minimize_score <- function(par, data, verbose = FALSE) {

  if (any(par < par_lower | par > par_upper)) return(1e12)

  result <- tryCatch({

    if (verbose) {
      nms <- if (!is.null(names(par))) names(par) else paste0("p", seq_along(par))
      cat("params =", paste(paste0(nms, "=", round(par, 6)), collapse = ", "), "\n")
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
              model_opts = list(
                rho.max  = par[1], rho.null = par[2],
                c.ov     = par[3], k.ov     = par[4],
                k        = par[5], tau      = par[6],
                eta.null = par[7]
              ),
              dyn_rho_max = FALSE
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

par_lower <- c(rho.max  = 200,  rho.null = 20,  c.ov = 1e-6,
               k.ov     = 0,    k        = 0.01,   tau = 0,   eta.null = 1e6)

par_upper <- c(rho.max  = 600,  rho.null = 200, c.ov = 0.01,
               k.ov     = 1,    k        = 0.5, tau = 1,   eta.null = 2e7)

par_start <- c(
  rho.max  = 401.2588,
  rho.null = 81.19417,
  c.ov     = 0.0005104722,
  k.ov     = 0.37856737,
  k        = 0.02993175,
  tau      = 0.02362476,
  eta.null = 8523356
)

cat("\nParameter bounds:\n")
print(rbind(lower = par_lower, upper = par_upper))
cat("\nStarting parameters:\n"); print(par_start)

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
  "EPS", "rmse", "nrmse", "nbias", "kge", "combined_score"
))

# =============================================================================
# TEST AT INITIAL GUESS
# =============================================================================

cat("\nTesting objective function at start values...\n")
cat("Initial score =", minimize_score(par_start, d_obs_fit_tibble, verbose = TRUE), "\n")

# =============================================================================
# OPTIMIZATION
# =============================================================================

cat("\nStarting optimization (Nelder-Mead)...\n")
opt <- optimx(
  par     = par_start,
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  verbose = TRUE,
  method  = "Nelder-Mead",
  control = list(trace = 1, maxit = 200)
)

opt_df     <- as.data.frame(opt)
best_idx   <- which.min(opt_df$value)
best_par   <- as.numeric(opt_df[best_idx, names(par_start), drop = TRUE])
names(best_par) <- names(par_start)
best_value <- opt_df$value[best_idx]

cat("\nBest score:", best_value, "\n")
cat("Best parameters:\n"); print(best_par)

fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 10)
cat(sprintf(
  "\nswe_results = swe_deltasnow(hs,\n  rho_max=%s, rho_null=%s, c_ov=%s,\n  k_ov=%s, k=%s, tau=%s, eta_null=%s\n)\n",
  fmt(best_par["rho.max"]),  fmt(best_par["rho.null"]), fmt(best_par["c.ov"]),
  fmt(best_par["k.ov"]),     fmt(best_par["k"]),        fmt(best_par["tau"]),
  fmt(best_par["eta.null"])
))

# =============================================================================
# SAVE RESULTS
# =============================================================================

fmt_tag <- function(x) gsub("-", "m", gsub("\\.", "p",
  format(x, scientific = FALSE, trim = TRUE, digits = 6)))

weight_vals <- c(SWE_NRMSE = WEIGHT_SWE_NRMSE, RHO_NRMSE = WEIGHT_RHO_NRMSE,
                 SWE_NBIAS = WEIGHT_SWE_NBIAS,  RHO_NBIAS = WEIGHT_RHO_NBIAS,
                 SWE_KGE   = WEIGHT_SWE_KGE,    RHO_KGE   = WEIGHT_RHO_KGE)

weight_tag <- paste0(names(weight_vals), "_", fmt_tag(weight_vals), collapse = "__")

save_dir <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_WIN21/data/R_opt_logs"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0("opt_results__", weight_tag, ".rds"))

saveRDS(list(
  opt        = opt,
  best_par   = best_par,
  best_value = best_value,
  weights    = weight_vals,
  fit_data   = d_obs_fit_tibble,
  val_data   = d_obs_val_tibble
), file = save_file)

cat("\nOptimization finished. Results saved to:\n", save_file, "\n", sep = "")

# =============================================================================
# CLEANUP
# =============================================================================

parallel::stopCluster(cl)