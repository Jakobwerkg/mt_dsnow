# rda2nc.R
# Convert Win21 calibration data (H_SWE_obs.Rda) to NetCDF4.
# Output: Win21_all.nc  (read by morris_sensitivity_Win21.ipynb)

suppressPackageStartupMessages({
  library(zoo)
  library(ncdf4)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
rda_path <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/dsnow/Win21_calib/H_SWE_obs.Rda"
nc_path  <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/dsnow/Win21_calib/Win21_all.nc"

# ── Station metadata ──────────────────────────────────────────────────────────
station_meta <- data.frame(
  Station = c("Holzgau", "Ladis", "Obernberg", "Koessen", "Felbertauern",
              "Innervillgraten", "Muerren", "Truebsee", "Ulrichen", "Zermatt",
              "Davos Flueelastr.", "Klosters KW", "San Bernardino", "Sta.Maria", "Zuoz"),
  Lon = c(10.3333, 10.6492, 11.4292, 12.4028, 12.5056,
          12.3750, 7.890193, 8.395291, 8.308283, 7.751165,
          9.848163, 9.895973, 9.184634, 10.419344, 9.962676),
  Lat = c(47.2500, 47.0969, 47.0194, 47.6717, 47.1181,
          46.8083, 46.558180, 46.791210, 46.504610, 46.023400,
          46.812550, 46.860580, 46.463260, 46.599810, 46.604330),
  Alt = c(1100, 1350, 1360, 590, 1650,
          1400, 1650, 1780, 1350, 1600,
          1560, 1200, 1640, 1415, 1710),
  stringsAsFactors = FALSE
)

# ── Load Rda into a clean environment (avoid polluting / clobbering globals) ──
stopifnot(file.exists(rda_path))
local_env <- new.env(parent = emptyenv())
load(rda_path, envir = local_env)
obj_names <- ls(local_env)
stopifnot(length(obj_names) >= 1L)
data_list <- get(obj_names[1L], envir = local_env)
cat("Loaded object  :", obj_names[1L], "  class:", class(data_list), "\n")
cat("Stations in Rda:", paste(names(data_list), collapse = ", "), "\n")

# ── Filter to stations with metadata ─────────────────────────────────────────
valid_stations <- names(data_list)[names(data_list) %in% station_meta$Station]
cat("Stations matched:", length(valid_stations), "\n")
dropped <- names(data_list)[!names(data_list) %in% station_meta$Station]
if (length(dropped) > 0L)
  cat("Dropped         :", paste(dropped, collapse = ", "), "\n")
stopifnot(length(valid_stations) > 0L)
data_list  <- data_list[valid_stations]
n_stations <- length(valid_stations)

# ── Build unified daily time axis ─────────────────────────────────────────────
# Normalise any Date / POSIXct / numeric index to Date before combining;
# keeps do.call(c, ...) from silently stripping the Date class.
extract_dates <- function(z) {
  idx <- tryCatch(index(z), error = function(e) NULL)
  if (is.null(idx) || length(idx) == 0L) return(as.Date(character(0)))
  as.Date(idx)   # handles Date, POSIXct, character "YYYY-MM-DD", numeric
}
all_dates <- sort(unique(do.call(c, lapply(data_list, extract_dates))))
n_time    <- length(all_dates)
cat(sprintf("Time steps      : %d  (%s – %s)\n",
            n_time, format(all_dates[1L]), format(all_dates[n_time])))
stopifnot(n_time > 0L)

# ── Fill data matrices [time × station] ───────────────────────────────────────
mat_HS  <- matrix(NA_real_, nrow = n_time, ncol = n_stations)
mat_SWE <- matrix(NA_real_, nrow = n_time, ncol = n_stations)

for (i in seq_along(valid_stations)) {
  st   <- valid_stations[i]
  z    <- data_list[[st]]
  rows <- match(as.Date(index(z)), all_dates)   # integer positions in all_dates
  ok   <- !is.na(rows)                          # guard against unmatched dates
  cd   <- coredata(z)
  if ("Hobs"   %in% colnames(cd)) mat_HS [rows[ok], i] <- as.numeric(cd[ok, "Hobs"])
  if ("SWEobs" %in% colnames(cd)) mat_SWE[rows[ok], i] <- as.numeric(cd[ok, "SWEobs"])
  cat(sprintf("  %-22s  %d HS obs  %d SWE obs\n",
              st, sum(!is.na(mat_HS[, i])), sum(!is.na(mat_SWE[, i]))))
}

# ── Write NetCDF4 ─────────────────────────────────────────────────────────────
if (file.exists(nc_path)) file.remove(nc_path)

time_vals <- as.integer(all_dates)   # days since 1970-01-01
nchar_max <- max(nchar(valid_stations))

dim_time    <- ncdim_def("time",    "days since 1970-01-01", time_vals,
                         unlim = TRUE, calendar = "standard")
dim_station <- ncdim_def("station", "", seq_len(n_stations), create_dimvar = FALSE)
dim_nchar   <- ncdim_def("nchar",   "", seq_len(nchar_max),  create_dimvar = FALSE)

# station_name: char variable (nchar × station) — decoded by the Python notebook
v_stn_name <- ncvar_def("station_name", "",  list(dim_nchar, dim_station), prec = "char")
v_lon      <- ncvar_def("lon", "degrees_east",  list(dim_station), NA_real_, prec = "double")
v_lat      <- ncvar_def("lat", "degrees_north", list(dim_station), NA_real_, prec = "double")
v_alt      <- ncvar_def("alt", "m",             list(dim_station), NA_real_, prec = "double")
# HS in metres (model and notebook threshold hs.max() < 0.05 assume metres)
v_HS  <- ncvar_def("HS",  "m",  list(dim_time, dim_station), NA_real_,
                   longname = "Observed snow depth",           prec = "float")
v_SWE <- ncvar_def("SWE", "mm", list(dim_time, dim_station), NA_real_,
                   longname = "Observed snow water equivalent", prec = "float")

nc <- nc_create(nc_path,
                list(v_stn_name, v_lon, v_lat, v_alt, v_HS, v_SWE),
                force_v4 = TRUE)
on.exit(nc_close(nc), add = TRUE)   # guarantees close even if a put fails

ncvar_put(nc, v_stn_name, valid_stations)
meta_ordered <- station_meta[match(valid_stations, station_meta$Station), ]
ncvar_put(nc, v_lon, meta_ordered$Lon)
ncvar_put(nc, v_lat, meta_ordered$Lat)
ncvar_put(nc, v_alt, meta_ordered$Alt)
ncvar_put(nc, v_HS,  mat_HS / 100)   # Hobs is in cm → convert to m for SnowToSwe
ncvar_put(nc, v_SWE, mat_SWE)

# Tag HS and SWE with auxiliary coordinates so xarray promotes them automatically
ncatt_put(nc, "HS",  "coordinates", "station_name lon lat alt")
ncatt_put(nc, "SWE", "coordinates", "station_name lon lat alt")

ncatt_put(nc, 0,      "title",       "Win21 calibration dataset – HS and SWE observations")
ncatt_put(nc, 0,      "source",      basename(rda_path))
ncatt_put(nc, 0,      "created",     format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
ncatt_put(nc, 0,      "Conventions", "CF-1.8")
ncatt_put(nc, "time", "calendar",    "standard")

cat("Written:", nc_path, "\n")
