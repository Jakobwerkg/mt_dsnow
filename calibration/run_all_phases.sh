#!/usr/bin/env bash
# Runs all calibration phases sequentially.
# Each phase calls run_calibration.sh with a specific weight combination.
#
# Weight arguments (positional):
#   1  w_swe_nrmse   weight for NRMSE of SWE
#   2  w_rho_nrmse   weight for NRMSE of bulk density
#   3  w_swe_bias    weight for NBIAS of SWE
#   4  w_rho_bias    weight for NBIAS of bulk density
#   5  w_kge_swe     weight for (1-KGE) of SWE             (optional, default 0.0; NM + DE)
#   6  w_kge_rho     weight for (1-KGE) of bulk density    (optional, default 0.0; NM + DE)
#
# Usage: ./run_all_phases.sh
# Estimated runtime: several hours

set -e

BASE="$(cd "$(dirname "$0")" && pwd)"
RUN="$BASE/run_calibration.sh"

# args: label  w_swe_nrmse  w_rho_nrmse  w_swe_bias  w_rho_bias  [w_kge_swe]  [w_kge_rho]
run_phase() {
  local label=$1 w1=$2 w2=$3 w3=$4 w4=$5 w5=${6:-0.0} w6=${7:-0.0}
  echo ""
  echo "###################################################"
  echo "#  $label"
  echo "#  SWE_NRMSE=$w1  RHO_NRMSE=$w2  SWE_BIAS=$w3  RHO_BIAS=$w4"
  echo "#  KGE_SWE=$w5    KGE_RHO=$w6    (NM + DE)"
  echo "###################################################"
  "$RUN" "$w1" "$w2" "$w3" "$w4" "$w5" "$w6"
}

# # ------------------------------------------------------------
# # PHASE 1 — NRMSE only, no bias
# # ------------------------------------------------------------
run_phase "Phase 1A — SWE-only baseline"           1.0 0.0 0.0 0.0
run_phase "Phase 1B — Density-only extreme"        0.0 1.0 0.0 0.0

# ------------------------------------------------------------
# PHASE 2 — Balanced NRMSE sweep (no bias)
# ------------------------------------------------------------
run_phase "Phase 2A — SWE-dominant   (0.7/0.3)"   0.7 0.3 0.0 0.0
run_phase "Phase 2B — Equal weight   (0.5/0.5)"   0.5 0.5 0.0 0.0
run_phase "Phase 2C — Density-dom.   (0.3/0.7)"   0.3 0.7 0.0 0.0

# ------------------------------------------------------------
# PHASE 3 — SWE bias penalty (NRMSE + SWE bias, no RHO bias)
# ------------------------------------------------------------
run_phase "Phase 3A — Balanced NRMSE + SWE bias   (0.6/0.2/0.2/0.0)" 0.6 0.2 0.2 0.0
run_phase "Phase 3B — SWE-dominant + SWE bias     (0.7/0.0/0.3/0.0)" 0.7 0.0 0.3 0.0
run_phase "Phase 3C — Density-dominant + SWE bias (0.3/0.5/0.2/0.0)" 0.3 0.5 0.2 0.0

# ------------------------------------------------------------
# PHASE 4 — RHO bias penalty (NRMSE + RHO bias, no SWE bias)
# ------------------------------------------------------------
run_phase "Phase 4A — Balanced NRMSE + RHO bias   (0.6/0.2/0.0/0.2)" 0.6 0.2 0.0 0.2
run_phase "Phase 4B — SWE-dominant + RHO bias     (0.7/0.1/0.0/0.2)" 0.7 0.1 0.0 0.2
run_phase "Phase 4C — Density-dominant + RHO bias (0.3/0.5/0.0/0.2)" 0.3 0.5 0.0 0.2

# ------------------------------------------------------------
# PHASE 5 — Both bias terms (NRMSE + SWE bias + RHO bias)
# ------------------------------------------------------------
run_phase "Phase 5A — Balanced NRMSE + both bias   (0.4/0.4/0.1/0.1)"   0.4  0.4  0.10 0.10
run_phase "Phase 5B — SWE-dominant + both bias     (0.8/0.1/0.05/0.05)" 0.8  0.1  0.05 0.05
run_phase "Phase 5C — Density-dominant + both bias (0.1/0.8/0.05/0.05)" 0.1  0.8  0.05 0.05
run_phase "Phase 5D — Bias-focused                 (0.25/0.25/0.25/0.25)" 0.25 0.25 0.25 0.25
run_phase "Phase 5E — Bias-only                    (0.0/0.0/0.5/0.5)"  0.0  0.0  0.50 0.50

# ------------------------------------------------------------
# PHASE 6 — KGE (one run per variable, symmetric substitution)
# ------------------------------------------------------------
run_phase "Phase 6A — KGE SWE + NRMSE RHO balanced (0.0/0.5/0.0/0.0/0.5/0.0)" 0.0 0.5 0.0 0.0 0.5 0.0
run_phase "Phase 6B — NRMSE SWE + KGE RHO balanced (0.5/0.0/0.0/0.0/0.0/0.5)" 0.5 0.0 0.0 0.0 0.0 0.5
run_phase "Phase 6C — KGE SWE                      (0.0/0.0/0.0/0.0/1.0/0.0)" 0.0 0.0 0.0 0.0 1.0 0.0
run_phase "Phase 6D — KGE RHO                      (0.0/0.0/0.0/0.0/0.0/1.0)" 0.0 0.0 0.0 0.0 0.0 1.0


# echo ""
# echo "###################################################"
# echo "#  All phases complete."
# echo "###################################################"