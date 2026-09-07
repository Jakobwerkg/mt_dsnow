# Clear workspace
rm(list = ls())

suppressPackageStartupMessages({
  library(zoo)
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


# -------- Paths --------
rda_path <- file.path(ROOT, "calibration/calibration_data/raw_data/win21/source_data/H_SWE_obs.Rda")
out_dir  <- file.path(ROOT, "calibration/calibration_data/output/HS_SWE_by_station")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Load RDA into obs_list --------
tmp_env <- new.env()
loaded_names <- load(rda_path, envir = tmp_env)
if (length(loaded_names) < 1) stop("Rda file did not contain any objects.")
obs_list <- tmp_env[[ loaded_names[1] ]]
rm(tmp_env)

stopifnot(is.list(obs_list))

# -------- Helper: extract date, hs, swe_obs from a zoo series --------
extract_station <- function(z) {
  if (!inherits(z, "zoo")) return(NULL)
  
  core <- coredata(z)
  cn   <- colnames(core)
  
  # Simpler: require exactly these columns, no fallback logic
  if (!all(c("Hobs", "SWEobs") %in% cn)) return(NULL)
  
  df <- data.frame(
    date    = as.POSIXct(index(z)),      # keep time: "2016-09-01 06:00:00"
    hs      = as.numeric(core[, "Hobs"])/100,
    swe_obs = as.numeric(core[, "SWEobs"]),
    row.names = NULL
  )
  
  # Optional: if you want to drop rows where hs is NA
  # df <- df[!is.na(df$hs), , drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  df
}

# -------- Iterate stations and write CSVs --------
n_written <- 0L

for (st in names(obs_list)) {
  z  <- obs_list[[st]]
  df <- extract_station(z)
  if (is.null(df)) next

  safe_st <- gsub("[^A-Za-z0-9]", "_", st)
  
  # "<safe_station>_hs_swe_obs.csv"
  out_file <- file.path(out_dir, sprintf("%s_hs_swe_obs.csv", safe_st))
  
  # na = "" → empty fields for NA, like in your example
  write.csv(df, out_file, row.names = FALSE, na = "", quote = FALSE)
  n_written <- n_written + 1L
}

message("Done. Wrote ", n_written, " CSV file(s) to: ", out_dir)