# ============================================================
# collect_opt_results.R
# ------------------------------------------------------------
# Searches recursively for all .rds files whose name starts
# with "opt_results", extracts key information, and writes
# a single summary table.
#
# Output columns (in order):
#   dataset | phase | algorithm
#   | w_SWE_NRMSE | w_RHO_NRMSE | w_SWE_NBIAS | w_RHO_NBIAS | w_SWE_KGE | w_RHO_KGE
#   | rho_max | rho_null | eta_null | k | tau | c_ov | k_ov
#   | best_value | iterations | convergence
#   | source_ctime | source_path
#
# Usage:
#   Rscript calibration/collect_opt_results.R
# ============================================================

# ============================================================
# USER SETTINGS — change paths here
# ============================================================
SEARCH_DIR <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/AA_opt_out/Win21/data"



OUT_CSV    <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/AA_opt_out/Win21/res/opt_results_summary.csv"
OUT_RDS    <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/AA_opt_out/Win21/res/opt_results_summary.rds"
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
if (length(files) == 0) stop("No matching .rds files found.")

# -----------------------------------------------------------
# 2) Helper functions
# -----------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

infer_dataset <- function(path) {
  if (grepl("calibration_Win21",    path)) return("Win21")
  if (grepl("calibration_SNOWPACK", path)) return("SNOWPACK")
  NA_character_
}

WEIGHT_NAMES <- c("SWE_NRMSE", "RHO_NRMSE", "SWE_NBIAS", "RHO_NBIAS", "SWE_KGE", "RHO_KGE")

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

DSNOW_PAR_NAMES <- c("rho_max", "rho_null", "eta_null", "k", "tau", "c_ov", "k_ov")

normalise_par_name <- function(x) {
  x <- gsub("\\.", "_", x)
  x <- gsub("rho_0$|rho0$",    "rho_null", x)
  x <- gsub("rhomax$",         "rho_max",  x)
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
      for (j in seq_along(bp)) {
        if (j <= length(DSNOW_PAR_NAMES)) par_out[[DSNOW_PAR_NAMES[j]]] <- as.numeric(bp[j])
      }
    }
  }

  rows[[i]] <- c(
    list(
      dataset   = infer_dataset(path),
      phase     = lookup_phase(w),
      algorithm = algorithm %||% NA_character_
    ),
    as.list(w),
    par_out,
    list(
      best_value   = best_value,
      iterations   = iterations,
      convergence  = convergence,
      source_ctime = format(fi$ctime, "%Y-%m-%d %H:%M"),
      source_path  = path
    )
  )

  message(sprintf("[%2d/%2d] %-55s best=%.4g  algo=%s",
                  i, length(files), basename(path),
                  best_value %||% NA, algorithm %||% "?"))
}

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
  "dataset", "phase", "algorithm",
  paste0("w_", WEIGHT_NAMES),
  DSNOW_PAR_NAMES,
  "best_value", "iterations", "convergence",
  "source_ctime", "source_path"
)
col_order <- intersect(col_order, names(df))
df <- df[, col_order, drop = FALSE]

df <- df[order(df$dataset, df$algorithm,
               df$w_SWE_NRMSE, df$w_RHO_NRMSE, df$w_SWE_NBIAS), ]

# -----------------------------------------------------------
# 5) Save + print summary
# -----------------------------------------------------------
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
saveRDS(df,   OUT_RDS)

message("\n=== Summary ===")
print_cols <- c("dataset", "phase", "algorithm",
                paste0("w_", WEIGHT_NAMES),
                DSNOW_PAR_NAMES,
                "best_value", "source_ctime")
print(df[, intersect(print_cols, names(df))], row.names = FALSE)

message("\nCSV: ", normalizePath(OUT_CSV, mustWork = FALSE))
message("RDS: ", normalizePath(OUT_RDS,  mustWork = FALSE))