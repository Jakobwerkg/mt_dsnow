# run_dsnow_validation.R
# Runs nixmass::swe.delta.snow on the Mag25 dataset and writes a validation
# NetCDF containing daily HNW_mod, SWE_mod and the original HNW_obs / SWE_obs.
# The output NC is consumed by HNW_validation.ipynb.

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
Mag25_nc_file    <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc"
out_nc_file      <- "/Users/jakobwerkgarner/code/mt_dsnow/HNW_validation/validation_data/Mag25_2B_WIN21_DE.nc"

exclude_stations <- c("Weisfluh_Joch")

# Use dynamic rho_max parameterisation?  FALSE = static (Winkler 2021 style)
dyn_rho_max <- FALSE

# Model parameters — uncomment and edit to override nixmass defaults
model_opts <- list(
        rho.max  = 345.6,
        rho.null = 73.8,
        c.ov     = 0.008001,
        k.ov     = 1.0,
        k        = 0.1463,
        tau      = 0.03,
        eta.null = 10510000
)

# nixmass uses dot-notation; translate any underscore keys
model_opts <- local({
  nms <- names(model_opts)
  nms <- gsub("rho_max",  "rho.max",  nms, fixed = TRUE)
  nms <- gsub("rho_null", "rho.null", nms, fixed = TRUE)
  nms <- gsub("eta_null", "eta.null", nms, fixed = TRUE)
  nms <- gsub("c_ov",     "c.ov",     nms, fixed = TRUE)
  nms <- gsub("k_ov",     "k.ov",     nms, fixed = TRUE)
  names(model_opts) <- nms
  model_opts
})


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
hyd_year_of <- function(d) {
  ifelse(month(d) > 8L, year(d), year(d) - 1L)
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


# ─────────────────────────────────────────────────────────────────────────────
# Load Mag25
# Mag25_all.nc stores station names as dimension values (not a variable),
# and time as a dimension with embedded numeric values (not a coordinate variable).
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(file.exists(Mag25_nc_file))
nc_in <- nc_open(Mag25_nc_file)
on.exit(nc_close(nc_in), add = TRUE)

time_dim   <- nc_in$dim[["time"]]
time_raw   <- time_dim$vals
time_units <- time_dim$units
time_cal   <- if (!is.null(time_dim$calendar) && nchar(time_dim$calendar) > 0)
                time_dim$calendar else "standard"
origin_dt  <- as.Date(substr(sub("days since ", "", time_units), 1, 10))
dates_all  <- origin_dt + as.integer(time_raw)

station_names <- nc_in$dim[["station"]]$vals

avail_vars <- names(nc_in$var)
message("Variables in NC: ", paste(avail_vars, collapse = ", "))

HS_all      <- ncvar_get(nc_in, "HS")   # [station × time], m
SWE_obs_all <- ncvar_get(nc_in, "SWE")  # [station × time], mm
HNW_obs_all <- if ("HNW" %in% avail_vars) {
  ncvar_get(nc_in, "HNW")
} else {
  message("HNW not in NC — HNW_obs will be all NA")
  matrix(NA_real_, nrow = length(station_names), ncol = length(dates_all))
}

Nt <- length(dates_all)
Ns <- length(station_names)
message(sprintf("Loaded: %d stations × %d days  (%s – %s)",
                Ns, Nt, dates_all[1], dates_all[Nt]))


# ─────────────────────────────────────────────────────────────────────────────
# Run ΔSnow per station × hydrological year
# ─────────────────────────────────────────────────────────────────────────────
hyd_years_all <- hyd_year_of(dates_all)
winter_years  <- sort(unique(hyd_years_all))
SWE_mod_all   <- matrix(NA_real_, nrow = Ns, ncol = Nt)

for (si in seq_len(Ns)) {
  stn <- station_names[si]
  if (stn %in% exclude_stations) { message("Skipping: ", stn); next }

  hs_stn <- as.numeric(HS_all[si, ])
  for (wy in winter_years) {
    idx <- which(hyd_years_all == wy)
    if (length(idx) < 10L) next
    swe_seg <- run_dsnow(dates_all[idx], hs_stn[idx], model_opts, dyn_rho_max)
    if (!is.null(swe_seg) && length(swe_seg) == length(idx))
      SWE_mod_all[si, idx] <- swe_seg
  }
  message("Done: ", stn)
}


# ─────────────────────────────────────────────────────────────────────────────
# Derive HNW_mod = diff(SWE_mod), clip melt (negative) to 0
# ─────────────────────────────────────────────────────────────────────────────
HNW_mod_all <- matrix(NA_real_, nrow = Ns, ncol = Nt)
for (si in seq_len(Ns)) {
  d <- diff(SWE_mod_all[si, ])
  d[d < 0] <- 0
  HNW_mod_all[si, 2:Nt] <- d
}


# ─────────────────────────────────────────────────────────────────────────────
# Write output NetCDF
# ─────────────────────────────────────────────────────────────────────────────
dir.create(dirname(out_nc_file), showWarnings = FALSE, recursive = TRUE)
# Close any leaked nc_out handle from a previous failed run before removing the file.
if (exists("nc_out")) { tryCatch(nc_close(nc_out), error = function(e) NULL); rm(nc_out) }
if (file.exists(out_nc_file)) file.remove(out_nc_file)

nchar_max   <- max(nchar(station_names))
dim_station <- ncdim_def("station", units = "", vals = seq_len(Ns),                    create_dimvar = FALSE)
dim_time    <- ncdim_def("time",    units = time_units, vals = as.integer(time_raw),   create_dimvar = TRUE)
dim_nchar   <- ncdim_def("nchar",   units = "", vals = seq_len(nchar_max),             create_dimvar = FALSE)

# "station_name" avoids a name clash with the "station" dimension (NC convention
# reserves same-name 1-D variables as coordinate variables).
v_stn     <- ncvar_def("station_name", "",   list(dim_nchar, dim_station), prec = "char")
v_swe_mod <- ncvar_def("SWE_mod",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
v_hnw_mod <- ncvar_def("HNW_mod",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
v_swe_obs <- ncvar_def("SWE_obs",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")
v_hnw_obs <- ncvar_def("HNW_obs",      "mm", list(dim_station, dim_time),  missval = NA_real_, prec = "double")

nc_out <- nc_create(out_nc_file, vars = list(v_stn, v_swe_mod, v_hnw_mod, v_swe_obs, v_hnw_obs))
on.exit(nc_close(nc_out), add = TRUE)   # guarantees close even if a put fails

ncvar_put(nc_out, v_stn,     station_names)
ncvar_put(nc_out, v_swe_mod, SWE_mod_all)
ncvar_put(nc_out, v_hnw_mod, HNW_mod_all)
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

message("Wrote: ", out_nc_file)