# single_station_season_comparison.R
#
# Compares, for ONE station and ONE season (hydrological year, Sep–Aug):
#   * ΔSnow (nixmass::swe.delta.snow) with default parameters
#   * ΔSnow with the "best_new" calibrated parameters
#   * HS2SWE (precomputed Python output from
#     model_diff/layerwise_data/combined_layerwise_default_Mag25.nc)
# against the Magnusson (Mag25) SWE / HNW observations.
#
# Produces a single-panel SWE plot with a centred metrics table below,
# styled to match plot_style.py.

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


# Rscript may start in the C locale, which mangles Δ / ² / — in plot text
if (!l10n_info()$`UTF-8`)
  suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))


# ─────────────────────────────────────────────────────────────────────────────
# CONFIG  (user-editable)
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: "Zermatt" is NOT in the Mag25 dataset — the script prints all valid
# station names if the one below is not found. Mag25 stations include e.g.
# Adelboden, Montana, Saas_Fee, Davos_Flueelastr, Ulrichen, Zuoz, ...
# (Maloja has good snow-course coverage in 2019/20: 11 SWE obs.)

station  <- "Saas_Fee"    # set target station
season   <- "2020_2021"   # set target season (hydrological year Sep–Aug)
save_fig <- TRUE

# ΔSnow parameter sets -------------------------------------------------------
# Units follow nixmass::swe.delta.snow: tau [m], eta.null [Pa s], c.ov [-].
model_opts_default <- list(
  rho.null = 81,      rho.max = 401,
  eta.null = 8.5e6,   k = 0.0300,
  tau      = 0.024,   c.ov = 5.1e-4,
  k.ov     = 0.38
)

# "best_new" calibration (Rain_Gauge phase=2C, Nelder-Mead) — identical to
# the set used in hnw_validation/run_single_dsnow_validation.R.
# c.ov: the mantissa 6.0 converts to 6.0e-4 (NOT 6.0e-6) — verified against
# run_single_dsnow_validation.R (c.ov = 0.0006) and the nixmass default scale
# (~5.1e-4).
model_opts_best_new <- list(
  rho.null = 101.17,    rho.max = 380.01,
  eta.null = 8.233e6,   k = 0.0269,
  tau      = 0.0227,    c.ov = 6.0e-4,
  k.ov     = 0.4117
)

# Static rho_max parameterisation (both param sets carry a fixed rho.max)
dyn_rho_max <- FALSE

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
Mag25_nc_file  <- file.path(ROOT, "calibration/calibration_data/raw_data/mag25/slf_dataset/Mag25_all.nc")
HS2SWE_nc_file <- file.path(ROOT, "model_diff/layerwise_data/combined_layerwise_default_Mag25.nc")
out_dir        <- file.path(ROOT, "single_station_season_comparison")

# ─────────────────────────────────────────────────────────────────────────────
# Style — mirrors plot_style.py (Okabe-Ito palette, FIG.SINGLE, 350 dpi)
# ─────────────────────────────────────────────────────────────────────────────
C_OBS    <- "#222222"   # near-black   — observed reference
C_DSNOW  <- "#E69F00"   # amber        — ΔSnow model
C_HS2SWE <- "#009E73"   # blue-green   — HS2SWE model
C_HS     <- "#999999"   # mid-grey     — snow height (hs), auxiliary
C_GRID   <- "#EBEBEB"   # subtle horizontal grid
C_EDGE   <- "#BBBBBB"   # light axes borders
C_TICK   <- "#777777"   # muted tick colour
C_LABEL  <- "#444444"   # axis-label colour
C_TITLE  <- "#333333"   # title colour

FIG_SIZE <- c(10, 7.2)  # slightly taller to accommodate centred table
FIG_DPI  <- 350
LWD_MAIN <- 2.1         # ≈ matplotlib lines.linewidth 1.6 pt
LWD_HS   <- 1.4         # thinner — HS is auxiliary

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
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

# ΔSnow wrapper — identical conventions to run_single_dsnow_validation.R
run_dsnow <- function(dates, hs_m, model_opts, dyn_rho_max) {
  hs <- pmax(as.numeric(hs_m), 0)
  hs[is.na(hs)] <- 0
  if (length(hs) == 0) return(NULL)
  hs[1]          <- 0
  hs[length(hs)] <- 0

  df  <- data.frame(date = as.character(dates), hs = hs,
                    stringsAsFactors = FALSE)
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

# HNW_mod = diff(SWE_mod), melt (negative) clipped to 0 — project convention
hnw_from_swe <- function(swe) {
  hnw <- c(NA_real_, diff(swe))
  hnw[hnw < 0] <- 0
  hnw
}

metrics <- function(mod, obs) {
  ok <- is.finite(mod) & is.finite(obs)
  n  <- sum(ok)
  if (n < 2L)
    return(list(rmse = NA_real_, rel_bias = NA_real_, r2 = NA_real_, n = n))
  m <- mod[ok]; o <- obs[ok]
  list(
    rmse     = sqrt(mean((m - o)^2)),
    rel_bias = 100 * sum(m - o) / sum(o),
    r2       = suppressWarnings(cor(m, o))^2,
    n        = n
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Parse season → hydrological-year window
# ─────────────────────────────────────────────────────────────────────────────
yrs <- as.integer(strsplit(season, "_")[[1]])
stopifnot(length(yrs) == 2L, yrs[2] == yrs[1] + 1L)
hyd_year <- yrs[1]

# ─────────────────────────────────────────────────────────────────────────────
# Load Magnusson observations (SWE, HNW, HS forcing)
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(Mag25_nc_file))
nc_mag <- nc_open(Mag25_nc_file)

dates_mag    <- nc_dates(nc_mag)
stations_mag <- nc_station_names(nc_mag)

if (!station %in% stations_mag) {
  nc_close(nc_mag)
  stop("Station '", station, "' not found in Mag25. Available stations:\n  ",
       paste(sort(stations_mag), collapse = ", "))
}
si_mag <- match(station, stations_mag)

# ncdf4 returns [station x time] for NC variables declared (time, station)
HS_obs  <- ncvar_get(nc_mag, "HS")[si_mag, ]    # m
SWE_obs <- ncvar_get(nc_mag, "SWE")[si_mag, ]   # mm (biweekly snow course)
HNW_obs <- ncvar_get(nc_mag, "HNW")[si_mag, ]   # mm
nc_close(nc_mag)

idx <- which(hyd_year_of(dates_mag) == hyd_year)
if (length(idx) < 10L)
  stop("Season ", season, " not (sufficiently) covered by Mag25 (",
       format(min(dates_mag)), " - ", format(max(dates_mag)), ").")

dates   <- dates_mag[idx]
hs_obs  <- HS_obs[idx]
swe_obs <- SWE_obs[idx]
hnw_obs <- HNW_obs[idx]

message(sprintf(
  "Station %s, season %s: %d days (%s - %s), %d SWE obs, %d HNW obs",
  station, season, length(dates), min(dates), max(dates),
  sum(is.finite(swe_obs)), sum(is.finite(hnw_obs))))

# ─────────────────────────────────────────────────────────────────────────────
# Load precomputed HS2SWE output (Python model — do not recompute)
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(HS2SWE_nc_file))
nc_h2s <- nc_open(HS2SWE_nc_file)

dates_h2s    <- nc_dates(nc_h2s)
stations_h2s <- nc_station_names(nc_h2s)
stopifnot(station %in% stations_h2s)
si_h2s <- match(station, stations_h2s)

SWE_h2s_all <- ncvar_get(nc_h2s, "hs2swe_swe_total")[si_h2s, ]  # mm
nc_close(nc_h2s)

# Align by date (both files share the Mag25 time axis, but don't assume it)
swe_h2s <- SWE_h2s_all[match(dates, dates_h2s)]
if (all(!is.finite(swe_h2s)))
  stop("HS2SWE output has no data for ", station, " in season ", season, ".")

# ─────────────────────────────────────────────────────────────────────────────
# Run ΔSnow for both parameter sets
# ─────────────────────────────────────────────────────────────────────────────
message("Running \u0394Snow (default parameters) ...")
swe_ds_def <- run_dsnow(dates, hs_obs, model_opts_default, dyn_rho_max)

message("Running \u0394Snow (best_new parameters) ...")
swe_ds_new <- run_dsnow(dates, hs_obs, model_opts_best_new, dyn_rho_max)

if (is.null(swe_ds_def) || is.null(swe_ds_new))
  stop("ΔSnow run failed — see warnings above.")

# Modelled HNW: positive daily SWE increments
hnw_ds_def <- hnw_from_swe(swe_ds_def)
hnw_ds_new <- hnw_from_swe(swe_ds_new)
hnw_h2s    <- hnw_from_swe(swe_h2s)

# ─────────────────────────────────────────────────────────────────────────────
# Metrics (per model, vs. obs, separately for SWE and HNW)
# ─────────────────────────────────────────────────────────────────────────────
runs <- list(
  "\u0394Snow default"   = list(swe = swe_ds_def, hnw = hnw_ds_def,
                                col = C_DSNOW,    lty = 2),
  "\u0394Snow best_new"  = list(swe = swe_ds_new, hnw = hnw_ds_new,
                                col = C_DSNOW,    lty = 1),
  "HS2SWE"               = list(swe = swe_h2s,    hnw = hnw_h2s,
                                col = C_HS2SWE,   lty = 1)
)

for (nm in names(runs)) {
  runs[[nm]]$m_swe <- metrics(runs[[nm]]$swe, swe_obs)
  runs[[nm]]$m_hnw <- metrics(runs[[nm]]$hnw, hnw_obs)
}

# ─────────────────────────────────────────────────────────────────────────────
# Build centred metrics table (console + figure panel)
# ─────────────────────────────────────────────────────────────────────────────
mono_fam <- if (capabilities("aqua")) "Menlo" else "mono"

# Column widths
W_MOD  <- 20   # model name
W_NUM  <-  9   # numeric columns
W_N    <-  6   # count column

mk_rule <- function(cross = "+", h = "-") {
  paste0(
    cross, strrep(h, W_MOD  + 2),
    cross, strrep(h, W_NUM  + 2),
    cross, strrep(h, W_NUM  + 2),
    cross, strrep(h, W_NUM  + 2),  # R2 same width
    cross, strrep(h, W_N    + 2),
    cross, strrep(h, W_NUM  + 2),  # double rule between SWE / HNW blocks
    cross, strrep(h, W_NUM  + 2),
    cross, strrep(h, W_NUM  + 2),
    cross, strrep(h, W_N    + 2),
    cross
  )
}

rule_outer <- mk_rule("+", "=")   # heavy rule (top / bottom / mid-header)
rule_inner <- mk_rule("+", "-")   # light rule (between data rows)

# Group header
grp_swe  <- "SWE (vs. biweekly obs.)"
grp_hnw  <- "HNW (vs. daily obs.)"
grp_span <- (W_NUM + 2) * 3 + (W_N + 2) + 2   # 4 columns incl. separators

fmt_grp <- function(label, span) {
  pad <- span - nchar(label)
  paste0(strrep(" ", floor(pad / 2)), label, strrep(" ", ceiling(pad / 2)))
}

row_group <- paste0(
  "| ", strrep(" ", W_MOD), " | ",
  fmt_grp(grp_swe, grp_span), " | ",
  fmt_grp(grp_hnw, grp_span), " |"
)

# Column header
row_colhdr <- sprintf(
  "| %-*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s |",
  W_MOD, "Model",
  W_NUM, "RMSE[mm]",  W_NUM, "rel.bias%",
  W_NUM, "R\u00b2",  W_N,   "n",
  W_NUM, "RMSE[mm]",  W_NUM, "rel.bias%",
  W_NUM, "R\u00b2",  W_N,   "n"
)

fmt_row <- function(nm, r) {
  fmt_val <- function(x, fmt) if (is.na(x)) sprintf(paste0("%", W_NUM, "s"), "NA") else sprintf(fmt, x)
  sprintf(
    "| %-*s | %*s | %*s | %*s | %*d | %*s | %*s | %*s | %*d |",
    W_MOD, nm,
    W_NUM, fmt_val(r$m_swe$rmse,     sprintf("%%%d.1f", W_NUM)),
    W_NUM, fmt_val(r$m_swe$rel_bias, sprintf("%%%d.1f", W_NUM)),
    W_NUM, fmt_val(r$m_swe$r2,       sprintf("%%%d.2f", W_NUM)),
    W_N,   r$m_swe$n,
    W_NUM, fmt_val(r$m_hnw$rmse,     sprintf("%%%d.2f", W_NUM)),
    W_NUM, fmt_val(r$m_hnw$rel_bias, sprintf("%%%d.1f", W_NUM)),
    W_NUM, fmt_val(r$m_hnw$r2,       sprintf("%%%d.2f", W_NUM)),
    W_N,   r$m_hnw$n
  )
}

data_rows <- vapply(names(runs), function(nm) fmt_row(nm, runs[[nm]]),
                    character(1))

table_lines <- c(
  rule_outer,
  row_group,
  rule_outer,
  row_colhdr,
  rule_inner,
  data_rows[1],
  rule_inner,
  data_rows[2],
  rule_inner,
  data_rows[3],
  rule_outer
)

cat("\n", paste(table_lines, collapse = "\n"), "\n\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────────────
fig_file <- file.path(out_dir, sprintf(
  "Model_intercomp_%s_%s_RMSE(SWE)_rel_bias(HNW).png", station, season))

if (save_fig) {
  png(fig_file,
      width  = FIG_SIZE[1], height = FIG_SIZE[2],
      units  = "in", res = FIG_DPI, bg = "white",
      type   = if (capabilities("aqua")) "quartz" else "cairo")
}

op <- par(no.readonly = TRUE)
layout(matrix(1:2, ncol = 1), heights = c(4.8, 2.4))   # plot / table panel

# ── Upper panel: SWE time series ─────────────────────────────────────────────
par(
  mar      = c(3.2, 3.6, 2.4, 3.6),
  mgp      = c(2.1, 0.6, 0),
  family   = "sans",
  col.axis = C_TICK,
  col.lab  = C_LABEL,
  fg       = C_EDGE,
  tcl      = -0.3,
  las      = 1
)

ylim <- range(0, swe_obs, swe_ds_def, swe_ds_new, swe_h2s, na.rm = TRUE) *
        c(1, 1.05)

plot(dates, swe_ds_new, type = "n", ylim = ylim,
     xlab = "", ylab = "SWE [mm]", xaxt = "n",
     cex.lab = 0.95, cex.axis = 0.85)

abline(h = pretty(ylim), col = C_GRID, lwd = 0.7)

month_ticks <- seq(as.Date(sprintf("%d-09-01", yrs[1])),
                   as.Date(sprintf("%d-09-01", yrs[2])), by = "1 month")
axis.Date(1, at = month_ticks, format = "%b", cex.axis = 0.85)

# Snow height — secondary right axis
hs_max   <- max(hs_obs, 0.1, na.rm = TRUE) * 1.05
hs_scale <- ylim[2] / hs_max
lines(dates, hs_obs * hs_scale, col = C_HS, lty = 2, lwd = LWD_HS)
hs_ticks <- pretty(c(0, hs_max))
hs_ticks <- hs_ticks[hs_ticks <= hs_max]
axis(4, at = hs_ticks * hs_scale, labels = hs_ticks, cex.axis = 0.85)
mtext("HS [m]", side = 4, line = 2.1, cex = 0.95, col = C_LABEL, las = 0)

box(col = C_EDGE)

for (nm in names(runs)) {
  r <- runs[[nm]]
  lines(dates, r$swe, col = r$col, lty = r$lty, lwd = LWD_MAIN)
}

ok_obs <- is.finite(swe_obs)
points(dates[ok_obs], swe_obs[ok_obs], col = C_OBS, pch = 16, cex = 0.7)

title(
  main     = sprintf("%s \u2014 season %s: \u0394Snow vs. HS2SWE vs. Mag25 obs",
                     gsub("_", " ", station), gsub("_", "/", season)),
  cex.main = 1.05, font.main = 1, col.main = C_TITLE
)

legend("topleft",
       legend  = c("Mag25 SWE obs", names(runs), "HS obs (right axis)"),
       col     = c(C_OBS,
                   vapply(runs, function(r) r$col, character(1)),
                   C_HS),
       lty     = c(NA,
                   vapply(runs, function(r) r$lty, numeric(1)),
                   2),
       lwd     = c(NA, rep(LWD_MAIN, length(runs)), LWD_HS),
       pch     = c(16, rep(NA, length(runs)),        NA),
       cex     = 0.8, seg.len = 2.4,
       bty     = "o", box.col = "#DDDDDD",
       bg      = adjustcolor("white", alpha.f = 0.85))

# ── Lower panel: centred metrics table ───────────────────────────────────────
par(mar    = c(0.2, 0, 0.4, 0),
    family = mono_fam)
plot.new()

n_lines <- length(table_lines)
line_h  <- 1 / (n_lines + 1)

for (i in seq_len(n_lines)) {
  text(x      = 0.5,
       y      = 1 - i * line_h,
       labels = table_lines[i],
       adj    = c(0.5, 0.5),      # centred anchor
       cex    = 0.70,
       col    = C_TITLE,
       xpd    = NA,
       family = mono_fam)
}

par(op)
layout(1)

if (save_fig) {
  dev.off()
  message("Saved: ", fig_file)
}