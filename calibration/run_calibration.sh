#!/usr/bin/env bash
# Usage: ./run_calibration.sh <w_swe_nrmse> <w_rho_nrmse> <w_swe_bias> <w_rho_bias> [w_kge_swe] [w_kge_rho]
#
# Runs all active calibration scripts sequentially with the given weights.
# Results are saved automatically by each R script (tagged by weight combination).
#
# Args:
#   W1  w_swe_nrmse  — weight for NRMSE of SWE               (required)
#   W2  w_rho_nrmse  — weight for NRMSE of bulk density       (required)
#   W3  w_swe_bias   — weight for NBIAS of SWE                (required)
#   W4  w_rho_bias   — weight for NBIAS of bulk density       (required)
#   W5  w_kge_swe    — weight for (1-KGE) of SWE             (optional, default 0.0)
#   W6  w_kge_rho    — weight for (1-KGE) of bulk density    (optional, default 0.0)
#
# Examples:
#   ./run_calibration.sh 1.0 0.0 0.0 0.0          # SWE-only baseline
#   ./run_calibration.sh 0.3 0.7 0.0 0.0          # density-dominant
#   ./run_calibration.sh 0.2 0.6 0.1 0.1          # with bias terms
#   ./run_calibration.sh 0.2 0.5 0.1 0.1 0.1 0.0  # with SWE KGE

set -e

W1=${1:?Error: provide w_swe_nrmse (e.g. 1.0)}
W2=${2:?Error: provide w_rho_nrmse (e.g. 0.0)}
W3=${3:?Error: provide w_swe_bias  (e.g. 0.0)}
W4=${4:?Error: provide w_rho_bias  (e.g. 0.0)}
W5=${5:-0.0}   # w_kge_swe — optional
W6=${6:-0.0}   # w_kge_rho — optional

BASE="$(cd "$(dirname "$0")" && pwd)"

echo "========================================================"
echo "  DeltaSnow calibration"
echo "  SWE_NRMSE=$W1  RHO_NRMSE=$W2  SWE_BIAS=$W3  RHO_BIAS=$W4"
echo "  KGE_SWE=$W5    KGE_RHO=$W6"
echo "  SNOWPACK subset: ${DSNOW_SUBSET:-win21}  (override with DSNOW_SUBSET=...)"
echo "========================================================"

# echo ""
# echo "[SNOWPACK] Nelder-Mead"
# Rscript "$BASE/calibration_snowpack/dsnow_parameter_optimization_nm.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

# echo ""
# echo "[SNOWPACK] Differential Evolution"
# Rscript "$BASE/calibration_snowpack/dsnow_parameter_optimization_de.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

# echo ""
# echo "[Win21] Nelder-Mead"
# Rscript "$BASE/calibration_win21/dsnow_parameter_optimization_nm.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

# echo ""
# echo "[Win21] Differential Evolution"
# Rscript "$BASE/calibration_win21/dsnow_parameter_optimization_de.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

echo ""
echo "[SP] Nelder-Mead"
Rscript "$BASE/calibration_rho_dyn/dsnow_parameter_optimization_nm_dsnow2.0.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

echo ""
echo "========================================================"
echo "  All calibrations complete."
echo "  Results from the steps enabled above are under:"
echo "    SNOWPACK: optimisation_output/${DSNOW_SUBSET:-win21}/data/"
echo "    Win21   : optimisation_output/win21/data/"
echo "  in R_opt_logs/ (Nelder-Mead) and R_opt_logs_DE/ (Differential Evolution),"
echo "  named opt_results__<weights>.rds / opt_results_DE__<weights>.rds"
echo "========================================================"
