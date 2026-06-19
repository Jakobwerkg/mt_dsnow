# run_dsnow_layerwise.R
# Runs nixmass::swe.delta.snow with layers=TRUE on all Mag25 stations
# (per hydrological year Sep–Aug) and writes a NetCDF with layerwise
# variables (h, swe, rho, age) for all stations.
#
# Output: dsnow_layerwise_Mag25.nc
# Dimensions: layer × time × station
# Run from the mt_dsnow project root.

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
Mag25_nc <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc"
out_nc   <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/layerwise_data/dyn_rho_max/dsnow_layerwise_Mag25_default.nc"

# ── Model settings (Winkler 2021 static defaults) ──────────────────────────────
dyn_rho_max <- TRUE
model_opts  <- list()

# ── Process encoding ───────────────────────────────────────────────────────────
# nixmass::swe.delta.snow (layers=TRUE) returns out$processes: a character vector
# with one entry per timestep describing what the model did.
# We encode it as integers for compact NetCDF storage.
#   0  = no snow (NA in proc)
#  -1  = runoff (all snow melted, HS dropped to 0)
#   1  = produced first layer (start of snow season)
#   2  = created new layer (snowfall exceeded tau threshold)
#   3  = drenching / melt (HS decreased by more than tau)
#  99  = unrecognised label
encode_proc <- function(p) {
  code <- integer(length(p))
  code[!is.na(p) & p == "runoff"]              <- -1L
  code[!is.na(p) & startsWith(p, "produce")]  <-  1L
  code[!is.na(p) & startsWith(p, "create")]   <-  2L
  code[!is.na(p) & startsWith(p, "drenching")]<-  3L
  code[!is.na(p) & code == 0L & !is.na(p) &
         !p %in% c("runoff") &
         !startsWith(p, "produce") &
         !startsWith(p, "create") &
         !startsWith(p, "drenching")]          <- 99L
  code
}

# ── Load Mag25 ─────────────────────────────────────────────────────────────────
nc_in         <- nc_open(Mag25_nc)
time_units    <- ncatt_get(nc_in, "time", "units")$value
time_cal      <- ncatt_get(nc_in, "time", "calendar")$value
time_raw      <- ncvar_get(nc_in, "time")
origin_dt     <- as.Date(substr(sub("days since ", "", time_units), 1, 10))
dates_all     <- origin_dt + as.integer(time_raw)
station_names <- ncvar_get(nc_in, "station")
HS_all        <- ncvar_get(nc_in, "HS")   # [station × time], metres
nc_close(nc_in)

Nt           <- length(dates_all)
Ns           <- length(station_names)
hyd_year     <- ifelse(month(dates_all) > 8L, year(dates_all), year(dates_all) - 1L)
winter_years <- sort(unique(hyd_year))

message(sprintf("Loaded: %d stations × %d time steps (%s – %s)",
                Ns, Nt, dates_all[1], dates_all[Nt]))

# ── First pass: run model, collect results in lists ────────────────────────────
# Each entry: list(idx, h_mat, swe_mat, rho_mat, age_mat)
# where h_mat is [nlayers × ntime_seg]
seg_data  <- vector("list", Ns)
swe_total <- matrix(NA_real_, nrow = Ns, ncol = Nt)
max_layers <- 1L

EXCLUDE_STATIONS <- "Weisfluh_Joch"

for (si in seq_len(Ns)) {
  if (station_names[si] %in% EXCLUDE_STATIONS) {
    message(sprintf("Skipping: %s (excluded from evaluation)", station_names[si]))
    next
  }
  hs_stn    <- as.numeric(HS_all[si, ])
  seg_data[[si]] <- list()

  for (wy in winter_years) {
    idx <- which(hyd_year == wy)
    if (length(idx) < 10L) next

    hs           <- pmax(as.numeric(hs_stn[idx]), 0)
    hs[is.na(hs)] <- 0
    if (hs[1]          != 0) hs[1]          <- 0
    if (hs[length(hs)] != 0) hs[length(hs)] <- 0

    df  <- data.frame(date = as.character(dates_all[idx]), hs = hs,
                      stringsAsFactors = FALSE)

    out <- tryCatch(
      nixmass::swe.delta.snow(df, model_opts = model_opts,
                              dyn_rho_max = dyn_rho_max, layers = TRUE,
                              strict_mode = TRUE, verbose = FALSE),
      error = function(e) {
        message(sprintf("  Error [%s wy=%d]: %s", station_names[si], wy, e$message))
        NULL
      }
    )
    if (is.null(out)) next

    swe_total[si, idx] <- as.numeric(out$SWE)

    # out$h is [nlayers × ntime_seg] — confirmed by existing run_lyer_delta_snow.R
    h_mat   <- as.matrix(out$h)
    swe_mat <- as.matrix(out$swe)
    age_mat <- as.matrix(out$age)
    nl      <- nrow(h_mat)

    # Compute density: rho [kg/m³] = swe [mm] / h [m]
    rho_mat <- matrix(NA_real_, nrow = nl, ncol = ncol(h_mat))
    pos     <- !is.na(h_mat) & h_mat > 0
    rho_mat[pos] <- swe_mat[pos] / h_mat[pos]

    # Diagnostics: encode process label and count active layers per timestep
    proc_vec  <- encode_proc(out$processes)
    nlyr_vec  <- as.integer(colSums(!is.na(h_mat) & h_mat > 0))

    if (nl > max_layers) max_layers <- nl

    seg_data[[si]] <- c(seg_data[[si]],
                        list(list(idx = idx, h = h_mat, swe = swe_mat,
                                  rho = rho_mat, age = age_mat,
                                  proc = proc_vec, nlayers = nlyr_vec)))
  }
  message(sprintf("Done: %s", station_names[si]))
}

message(sprintf("Max layers across all station-winters: %d", max_layers))

# ── Allocate 3D arrays: [layer × time × station] ──────────────────────────────
# ncdf4 column-major: first listed dim varies fastest in memory.
# Variable defined with list(dim_layer, dim_time, dim_station) → R array [layer, time, station].
h_arr   <- array(NA_real_, dim = c(max_layers, Nt, Ns))
swe_arr <- array(NA_real_, dim = c(max_layers, Nt, Ns))
rho_arr <- array(NA_real_, dim = c(max_layers, Nt, Ns))
age_arr <- array(NA_real_, dim = c(max_layers, Nt, Ns))
# 2D diagnostic arrays: [station × time]
proc_arr <- matrix(0L,       nrow = Ns, ncol = Nt)
nlyr_arr <- matrix(0L,       nrow = Ns, ncol = Nt)

for (si in seq_len(Ns)) {
  for (seg in seg_data[[si]]) {
    idx <- seg$idx
    nl  <- nrow(seg$h)
    # seg$h is [nl × ntime_seg]; assign to [1:nl, idx, si]
    h_arr  [1:nl, idx, si] <- seg$h
    swe_arr[1:nl, idx, si] <- seg$swe
    rho_arr[1:nl, idx, si] <- seg$rho
    age_arr[1:nl, idx, si] <- seg$age
    proc_arr[si, idx]      <- seg$proc
    nlyr_arr[si, idx]      <- seg$nlayers
  }
}

# ── Write output NetCDF ────────────────────────────────────────────────────────
dir.create(dirname(out_nc), recursive = TRUE, showWarnings = FALSE)
if (file.exists(out_nc)) file.remove(out_nc)

nchar_max  <- max(nchar(station_names))
dim_layer  <- ncdim_def("layer",   "layer_number", vals = seq_len(max_layers), create_dimvar = TRUE)
dim_time   <- ncdim_def("time",    time_units,     vals = as.integer(time_raw), create_dimvar = TRUE)
dim_station <- ncdim_def("station", "",            vals = seq_len(Ns), create_dimvar = FALSE)
dim_nchar  <- ncdim_def("nchar",   "",             vals = seq_len(nchar_max), create_dimvar = FALSE)

mv <- NA_real_

v_stn     <- ncvar_def("station",          "",      list(dim_nchar, dim_station), prec = "char")
v_swe_tot <- ncvar_def("dsnow_swe_total",  "mm",    list(dim_station, dim_time),               missval = mv, prec = "double")
v_h       <- ncvar_def("dsnow_h",          "m",     list(dim_layer, dim_time, dim_station),    missval = mv, prec = "float")
v_swe     <- ncvar_def("dsnow_swe",        "mm",    list(dim_layer, dim_time, dim_station),    missval = mv, prec = "float")
v_rho     <- ncvar_def("dsnow_rho",        "kg/m3", list(dim_layer, dim_time, dim_station),    missval = mv, prec = "float")
v_age     <- ncvar_def("dsnow_age",        "days",  list(dim_layer, dim_time, dim_station),    missval = mv, prec = "float")
v_proc    <- ncvar_def("dsnow_process",    "1",     list(dim_station, dim_time),               missval = NA_integer_, prec = "integer")
v_nlyr    <- ncvar_def("dsnow_n_layers",   "1",     list(dim_station, dim_time),               missval = NA_integer_, prec = "integer")

nc_out <- nc_create(out_nc, vars = list(v_stn, v_swe_tot, v_h, v_swe, v_rho, v_age, v_proc, v_nlyr))

ncvar_put(nc_out, v_stn,     station_names)
ncvar_put(nc_out, v_swe_tot, swe_total)
ncvar_put(nc_out, v_h,       h_arr)
ncvar_put(nc_out, v_swe,     swe_arr)
ncvar_put(nc_out, v_rho,     rho_arr)
ncvar_put(nc_out, v_age,     age_arr)
ncvar_put(nc_out, v_proc,    proc_arr)
ncvar_put(nc_out, v_nlyr,    nlyr_arr)

ncatt_put(nc_out, "time",            "calendar",      time_cal)
ncatt_put(nc_out, "layer",           "long_name",     "Snow layer index (1=bottom/oldest, N=top/newest)")
ncatt_put(nc_out, "dsnow_process",   "long_name",     "DeltaSnow process code per timestep")
ncatt_put(nc_out, "dsnow_process",   "flag_values",   "-1 0 1 2 3 99")
ncatt_put(nc_out, "dsnow_process",   "flag_meanings",
          "runoff no_snow first_layer new_layer drenching unknown")
ncatt_put(nc_out, "dsnow_process",   "note",
          paste("-1=runoff(season end), 0=no snow, 1=first layer produced (season start),",
                "2=new layer created (fresh snow > tau), 3=drenching/melt (HS decrease > tau)"))
ncatt_put(nc_out, "dsnow_n_layers",  "long_name",     "Number of active snow layers (h > 0)")
ncatt_put(nc_out, 0, "source",            "nixmass::swe.delta.snow, layers=TRUE")
ncatt_put(nc_out, 0, "forcing",           Mag25_nc)
ncatt_put(nc_out, 0, "dyn_rho_max",       as.integer(dyn_rho_max))
ncatt_put(nc_out, 0, "n_stations",        Ns)
ncatt_put(nc_out, 0, "n_timesteps",       Nt)
ncatt_put(nc_out, 0, "max_layers",        max_layers)
ncatt_put(nc_out, 0, "creation_date",     as.character(Sys.time()))
ncatt_put(nc_out, 0, "description",
          "Layer-wise DeltaSnow output for all Mag25 stations, per hydrological year (Sep-Aug)")

nc_close(nc_out)
message(sprintf("Wrote: %s  [%d layers × %d time × %d stations]",
                out_nc, max_layers, Nt, Ns))
