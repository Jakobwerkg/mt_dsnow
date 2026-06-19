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
#   ./run_calibration.sh 0.2 0.5 0.1 0.1 0.1 0.0  # with SWE KGE (DE only)

set -e

W1=${1:?Error: provide w_swe_nrmse (e.g. 1.0)}
W2=${2:?Error: provide w_rho_nrmse (e.g. 0.0)}
W3=${3:?Error: provide w_swe_bias  (e.g. 0.0)}
W4=${4:?Error: provide w_rho_bias  (e.g. 0.0)}
W5=${5:-0.0}   # w_kge_swe — optional, DE only
W6=${6:-0.0}   # w_kge_rho — optional, DE only

BASE="$(cd "$(dirname "$0")" && pwd)"

echo "========================================================"
echo "  DeltaSnow calibration with Rain Gauge data and Dynamic RhoMax"
echo "  SWE_NRMSE=$W1  RHO_NRMSE=$W2  SWE_BIAS=$W3  RHO_BIAS=$W4"
echo "  KGE_SWE=$W5    KGE_RHO=$W6    (NM + DE)"
echo "========================================================"

echo ""
echo "[1/2] SNOWPACK — Nelder-Mead"
Rscript "$BASE/calibration_SNOWPACK/dsnow_parameter_optimization.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

echo ""
echo "[2/2] SNOWPACK — Differential Evolution"
Rscript "$BASE/calibration_SNOWPACK/dsnow_parameter_optimization_DE.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

# echo ""
# echo "[3/4] Win21 — Nelder-Mead"
# Rscript "$BASE/calibration_Win21/WIN21_dsnow_paprameter_optimization.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

# echo ""
# echo "[4/4] Win21 — Differential Evolution"
# Rscript "$BASE/calibration_Win21/WIN21_dsnow_parameter_opitmization_DE.R" "$W1" "$W2" "$W3" "$W4" "$W5" "$W6"

echo ""
echo "========================================================"
echo "  All calibrations complete."
echo "  Results saved to:"
echo "    calibration_SNOWPACK/data/R_opt_logs/"
echo "    calibration_SNOWPACK/data/R_opt_logs_DE/"
echo "========================================================"