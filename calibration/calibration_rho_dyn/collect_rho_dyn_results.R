# ============================================================
# collect_rho_dyn_results.R
# ------------------------------------------------------------
# Collects every calibration result of the dynamic-rho.max
# deltaSNOW variant (dyn_rho_max = TRUE, 10 parameters) into a
# single summary table.
#
# Searches recursively under dyn_rho_max_res/ for all .rds files
# whose name starts with "opt_results", i.e. over all subsets
# and both optimisers:
#
#   dyn_rho_max_res/<subset>/de_res/opt_results_DE__<weights>.rds
#   dyn_rho_max_res/<subset>/nm_res/opt_results__<weights>.rds
#
# Output columns (in order):
#   subset | dataset | phase | algorithm
#   | w_SWE_NRMSE | w_RHO_NRMSE | w_SWE_NBIAS | w_RHO_NBIAS | w_SWE_KGE | w_RHO_KGE
#   | sigma | mu | rho_h | rho_l | rho_null | eta_null | k | tau | c_ov | k_ov
#   | best_value | iterations | convergence | dyn_rho_max
#   | source_ctime | source_path
#
# A control run with the nixmass defaults of the dynamic model (Winkler et
# al. 2021: Win21 stations, L-BFGS-B) is always appended as one extra row
# (subset "win21", phase "default", algorithm "L-BFGS-B", weights NA), so the
# uncalibrated reference is validated and ranked with the calibrated sets.
#
# The table is the input of run_full_validation_rho_dyn.R.
#
# Usage:
#   Rscript calibration/calibration_rho_dyn/collect_rho_dyn_results.R
# ============================================================

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

# ============================================================
# USER SETTINGS — change paths here
# ============================================================

SEARCH_DIR <- file.path(ROOT, "calibration/calibration_rho_dyn/dyn_rho_max_res")

OUT_CSV    <- file.path(SEARCH_DIR, "all_summaries_dyn_rho_max.csv")
OUT_RDS    <- file.path(SEARCH_DIR, "all_summaries_dyn_rho_max.rds")
# ============================================================

SEARCH_DIR <- normalizePath(SEARCH_DIR, mustWork = TRUE)

message("Searching under : ", SEARCH_DIR)

# -----------------------------------------------------------
# 1) Find all relevant files
# -----------------------------------------------------------
files <- list.files(SEARCH_DIR, pattern = "^opt_results.*\\.rds$",
                    recursive = TRUE, full.names = TRUE)
files <- normalizePath(files, mustWork = FALSE)
files <- files[!grepl("opt_results_summary", files, fixed = TRUE)]

message("Files found: ", length(files))
if (length(files) == 0) message("No calibration results found — only the control run will be written.")

# -----------------------------------------------------------
# 2) Helper functions
# -----------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# dyn_rho_max_res/<subset>/<de_res|nm_res>/<file>.rds
infer_subset <- function(path) basename(dirname(dirname(path)))

infer_dataset <- function(subset) {
  if (grepl("^sp_", subset)) return("SNOWPACK")
  if (grepl("win21", subset, ignore.case = TRUE)) return("Win21")
  NA_character_
}

WEIGHT_NAMES <- c("SWE_NRMSE", "RHO_NRMSE", "SWE_NBIAS", "RHO_NBIAS", "SWE_KGE", "RHO_KGE")

# Same phase table as optimisation_output/helpers/collect_opt_results.R
PHASE_LOOKUP <- data.frame(
  SWE_NRMSE = c(1.0,  0.0,  0.7,  0.5,  0.3,  0.6,  0.7,  0.3,  0.6,  0.7,  0.3,  0.40, 0.80, 0.10, 0.25, 0.0,  0.0,  0.5,  0.0,  0.0 ),
  RHO_NRMSE = c(0.0,  1.0,  0.3,  0.5,  0.7,  0.2,  0.0,  0.5,  0.2,  0.1,  0.5,  0.40, 0.10, 0.80, 0.25, 0.0,  0.5,  0.0,  0.0,  0.0 ),
  SWE_NBIAS = c(0.0,  0.0,  0.0,  0.0,  0.0,  0.2,  0.3,  0.2,  0.0,  0.0,  0.0,  0.10, 0.05, 0.05, 0.25, 0.50, 0.0,  0.0,  0.0,  0.0 ),
  RHO_NBIAS = c(0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.2,  0.2,  0.2,  0.10, 0.05, 0.05, 0.25, 0.50, 0.0,  0.0,  0.0,  0.0 ),
  SWE_KGE   = c(0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  1.0,  0.0 ),
  RHO_KGE   = c(0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  1.0 ),
  label     = c(
    "1A", "1B",
    "2A", "2B", "2C",
    "3A", "3B", "3C",
    "4A", "4B", "4C",
    "5A", "5B", "5C", "5D", "5E",
    "6A", "6B", "6C", "6D"
  ),
  stringsAsFactors = FALSE
)

lookup_phase <- function(w) {
  key_cols <- c("SWE_NRMSE", "RHO_NRMSE", "SWE_NBIAS", "RHO_NBIAS", "SWE_KGE", "RHO_KGE")
  diffs <- abs(sweep(as.matrix(PHASE_LOOKUP[, key_cols]), 2, unname(w[key_cols]), "-"))
  match_rows <- which(rowSums(diffs) < 1e-6)
  if (length(match_rows) == 1) PHASE_LOOKUP$label[match_rows] else NA_character_
}

# The 10 parameters of the dynamic-rho.max model (underscore notation, as in
# the static summary tables).  sigma/mu/rho_h/rho_l replace rho_max.
DSNOW_PAR_NAMES <- c("sigma", "mu", "rho_h", "rho_l",
                     "rho_null", "eta_null", "k", "tau", "c_ov", "k_ov")

normalise_par_name <- function(x) {
  x <- gsub("\\.", "_", x)
  x <- gsub("rho_0$|rho0$",    "rho_null", x)
  x <- gsub("eta_0$|eta0$",    "eta_null", x)
  x <- gsub("^cov$",           "c_ov",     x)
  x <- gsub("^kov$",           "k_ov",     x)
  x
}

# -----------------------------------------------------------
# 3) Loop over all files
# -----------------------------------------------------------
rows <- vector("list", length(files))

for (i in seq_along(files)) {
  path <- files[i]
  fi   <- file.info(path)

  obj <- tryCatch(readRDS(path),
                  error = function(e) { message("ERR ", path, ": ", e$message); NULL })

  w <- setNames(rep(0, length(WEIGHT_NAMES)), WEIGHT_NAMES)
  if (!is.null(obj$weights) && length(obj$weights) > 0) {
    wv  <- obj$weights
    nms <- toupper(names(wv) %||% character(0))
    nms <- gsub("W_|_WEIGHT", "", nms)
    for (j in seq_along(wv)) {
      nm <- nms[j]
      if (nm %in% WEIGHT_NAMES) w[nm] <- as.numeric(wv[j])
    }
  } else {
    stem <- sub("\\.rds$", "", basename(path))
    for (nm in WEIGHT_NAMES) {
      rx <- paste0(nm, "_([0-9p]+)(?:__|$)")
      m  <- regmatches(stem, regexec(rx, stem))[[1]]
      if (length(m) >= 2) w[nm] <- as.numeric(sub("p", ".", m[2], fixed = TRUE))
    }
  }

  best_value  <- as.numeric(obj$best_value %||% NA)
  iterations  <- NA_integer_
  convergence <- NA_character_
  algorithm   <- NA_character_

  opt <- obj$opt
  if (inherits(opt, c("optimx", "data.frame"))) {
    iterations  <- as.integer(opt$fevals[1]    %||% NA)
    convergence <- as.character(opt$convcode[1] %||% NA)
    algorithm   <- rownames(opt)[1] %||%
                   tryCatch(attr(opt, "details")[[1]][1], error = function(e) NA_character_)
  } else if (inherits(opt, "DEoptim")) {
    iterations  <- as.integer(opt$optim$iter    %||% NA)
    convergence <- as.character(opt$optim$nfeval %||% NA)
    algorithm   <- "DE"
    if (is.na(best_value))
      best_value <- as.numeric(opt$optim$bestval %||% NA)
  } else if (is.list(opt) && all(c("par", "value") %in% names(opt))) {
    iterations  <- as.integer(opt$counts[1]    %||% NA)
    convergence <- as.character(opt$convergence %||% NA)
    algorithm   <- "optim"
  }

  par_out <- setNames(as.list(rep(NA_real_, length(DSNOW_PAR_NAMES))), DSNOW_PAR_NAMES)
  bp_src <- if (!is.null(obj$best_par) && length(obj$best_par) > 0) {
    obj$best_par
  } else if (inherits(opt, "DEoptim") && !is.null(opt$optim$bestmem)) {
    opt$optim$bestmem
  } else {
    NULL
  }

  if (!is.null(bp_src)) {
    bp  <- bp_src
    nms <- normalise_par_name(names(bp) %||% character(0))
    if (length(nms) == length(bp)) {
      for (j in seq_along(bp)) {
        nm <- nms[j]
        if (nm %in% DSNOW_PAR_NAMES) par_out[[nm]] <- as.numeric(bp[j])
      }
    } else {
      warning("Unnamed best_par in ", basename(path), " — parameters left NA")
    }
  }

  subset <- as.character(obj$subset %||% infer_subset(path))

  rows[[i]] <- c(
    list(
      subset    = subset,
      dataset   = infer_dataset(subset),
      phase     = lookup_phase(w),
      algorithm = algorithm %||% NA_character_
    ),
    as.list(w),
    par_out,
    list(
      best_value   = best_value,
      iterations   = iterations,
      convergence  = convergence,
      dyn_rho_max  = as.logical(obj$dyn_rho_max %||% NA),
      source_ctime = format(fi$ctime, "%Y-%m-%d %H:%M"),
      # stored relative to the project root, so the table stays portable
      source_path  = sub(paste0("^", ROOT, "/?"), "", path)
    )
  )

  message(sprintf("[%2d/%2d] %-10s %-55s best=%.4g  algo=%s",
                  i, length(files), subset, basename(path),
                  best_value %||% NA, algorithm %||% "?"))
}

# -----------------------------------------------------------
# 3b) Control run: nixmass defaults of the dynamic model
#     (swe.delta.snow, dyn_rho_max = TRUE; Winkler et al. 2021,
#     calibrated on the Win21 stations with L-BFGS-B).
#     Always part of the table, so every summary and the
#     validation include the uncalibrated reference.
# -----------------------------------------------------------
CONTROL_PAR <- c(sigma = 0.03, mu = 80, rho_h = 600, rho_l = 380,
                 rho_null = 80.73706, eta_null = 8543502, k = 0.029297,
                 tau = 0.02356521, c_ov = 0.0005170964, k_ov = 0.3782312)
stopifnot(identical(names(CONTROL_PAR), DSNOW_PAR_NAMES))

rows[[length(rows) + 1]] <- c(
  list(
    subset    = "win21",
    dataset   = "Win21",
    phase     = "default",
    algorithm = "L-BFGS-B"
  ),
  as.list(setNames(rep(NA_real_, length(WEIGHT_NAMES)), WEIGHT_NAMES)),   # no objective weights
  as.list(CONTROL_PAR),
  list(
    best_value   = NA_real_,
    iterations   = NA_integer_,
    convergence  = NA_character_,
    dyn_rho_max  = TRUE,
    source_ctime = NA_character_,
    source_path  = "nixmass::swe.delta.snow defaults (dyn_rho_max = TRUE), Winkler et al. (2021)"
  )
)
message(sprintf("[control] %-10s %-55s %s", "win21", "nixmass defaults (dyn_rho_max = TRUE)", "algo=L-BFGS-B"))

# -----------------------------------------------------------
# 4) Assemble data.frame
# -----------------------------------------------------------
df <- do.call(rbind, lapply(rows, function(r) {
  as.data.frame(r, stringsAsFactors = FALSE)
}))

num_cols <- c(WEIGHT_NAMES, DSNOW_PAR_NAMES, "best_value", "iterations")
for (col in intersect(num_cols, names(df))) {
  df[[col]] <- as.numeric(unlist(df[[col]]))
}

names(df) <- sub("^(SWE_|RHO_)", "w_\\1", names(df))

col_order <- c(
  "subset", "dataset", "phase", "algorithm",
  paste0("w_", WEIGHT_NAMES),
  DSNOW_PAR_NAMES,
  "best_value", "iterations", "convergence", "dyn_rho_max",
  "source_ctime", "source_path"
)
col_order <- intersect(col_order, names(df))
df <- df[, col_order, drop = FALSE]

df <- df[order(df$subset, df$algorithm,
               df$w_SWE_NRMSE, df$w_RHO_NRMSE, df$w_SWE_NBIAS), ]

# -----------------------------------------------------------
# 5) Save + print summary
# -----------------------------------------------------------
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
saveRDS(df,   OUT_RDS)

message("\n=== Summary ===")
print_cols <- c("subset", "phase", "algorithm",
                DSNOW_PAR_NAMES,
                "best_value", "source_ctime")
print(df[, intersect(print_cols, names(df))], row.names = FALSE)

message("\nCSV: ", normalizePath(OUT_CSV, mustWork = FALSE))
message("RDS: ", normalizePath(OUT_RDS,  mustWork = FALSE))
