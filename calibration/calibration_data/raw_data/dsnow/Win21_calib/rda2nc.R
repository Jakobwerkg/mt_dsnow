# Clear workspace
rm(list = ls())
suppressPackageStartupMessages({
  library(zoo)
  library(ncdf4)
})

# -------- Paths --------
rda_path <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/dsnow/Win21_calib/H_SWE_obs.Rda"
nc_path  <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/dsnow/Win21_calib/Win21_calibration.nc"

# -------- Station metadata --------
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

# -------- Load Rda --------
load(rda_path)
loaded_obj <- ls()[!ls() %in% c("rda_path", "nc_path", "station_meta")]
data_list  <- get(loaded_obj[1])

# -------- Filter to stations with metadata --------
valid_stations <- names(data_list)[names(data_list) %in% station_meta$Station]
cat("Stations matched:", length(valid_stations), "\n")
cat("Dropped:", paste(names(data_list)[!names(data_list) %in% station_meta$Station], collapse = ", "), "\n")
data_list <- data_list[valid_stations]
n_stations <- length(valid_stations)

# -------- Unified time axis --------
all_dates <- sort(unique(do.call(c, lapply(data_list, index))))
n_time    <- length(all_dates)
cat("Time steps:", n_time, "\n")

# -------- Fill matrices --------
mat_HS  <- matrix(NA_real_, nrow = n_time, ncol = n_stations)
mat_SWE <- matrix(NA_real_, nrow = n_time, ncol = n_stations)

for (i in seq_along(valid_stations)) {
  st  <- valid_stations[i]
  z   <- data_list[[st]]
  idx <- match(index(z), all_dates)
  if ("Hobs"   %in% colnames(z)) mat_HS [idx, i] <- as.numeric(coredata(z)[, "Hobs"])
  if ("SWEobs" %in% colnames(z)) mat_SWE[idx, i] <- as.numeric(coredata(z)[, "SWEobs"])
}

# -------- NC dimensions --------
dim_time    <- ncdim_def("time",    "days since 1970-01-01", as.numeric(all_dates),
                         unlim = TRUE, calendar = "standard")
dim_station <- ncdim_def("station", "", seq_len(n_stations), create_dimvar = FALSE)

# -------- NC variables --------
var_lon <- ncvar_def("lon", "degrees_east",  list(dim_station), NA_real_, longname = "longitude",                    prec = "double")
var_lat <- ncvar_def("lat", "degrees_north", list(dim_station), NA_real_, longname = "latitude",                     prec = "double")
var_alt <- ncvar_def("alt", "m",             list(dim_station), NA_real_, longname = "altitude above sea level",     prec = "double")
var_HS  <- ncvar_def("HS",  "cm",            list(dim_time, dim_station), NA_real_, longname = "Observed snow depth",           prec = "float")
var_SWE <- ncvar_def("SWE", "mm",            list(dim_time, dim_station), NA_real_, longname = "Observed snow water equivalent", prec = "float")

# -------- Write NC --------
nc <- nc_create(nc_path, list(var_lon, var_lat, var_alt, var_HS, var_SWE))

ncatt_put(nc, 0, "title",        "Win21 calibration dataset – HS and SWE observations")
ncatt_put(nc, 0, "source",       basename(rda_path))
ncatt_put(nc, 0, "created",      format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
ncatt_put(nc, 0, "Conventions",  "CF-1.8")
ncatt_put(nc, 0, "station_names", paste(valid_stations, collapse = ","))

meta_ordered <- station_meta[match(valid_stations, station_meta$Station), ]
ncvar_put(nc, var_lon, meta_ordered$Lon)
ncvar_put(nc, var_lat, meta_ordered$Lat)
ncvar_put(nc, var_alt, meta_ordered$Alt)
ncvar_put(nc, var_HS,  mat_HS)
ncvar_put(nc, var_SWE, mat_SWE)

nc_close(nc)
cat("Written:", nc_path, "\n")