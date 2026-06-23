@echo off
REM Usage: run_calibration.bat <w_swe_nrmse> <w_rho_nrmse> <w_swe_bias> <w_rho_bias> [w_kge_swe] [w_kge_rho]
REM
REM Runs all active calibration scripts sequentially with the given weights.
REM Results are saved automatically by each R script (tagged by weight combination).
REM
REM Args:
REM   W1  w_swe_nrmse  -- weight for NRMSE of SWE               (required)
REM   W2  w_rho_nrmse  -- weight for NRMSE of bulk density       (required)
REM   W3  w_swe_bias   -- weight for NBIAS of SWE                (required)
REM   W4  w_rho_bias   -- weight for NBIAS of bulk density       (required)
REM   W5  w_kge_swe    -- weight for (1-KGE) of SWE             (optional, default 0.0)
REM   W6  w_kge_rho    -- weight for (1-KGE) of bulk density    (optional, default 0.0)
REM
REM Examples:
REM   run_calibration.bat 1.0 0.0 0.0 0.0          # SWE-only baseline
REM   run_calibration.bat 0.3 0.7 0.0 0.0          # density-dominant
REM   run_calibration.bat 0.2 0.6 0.1 0.1          # with bias terms
REM   run_calibration.bat 0.2 0.5 0.1 0.1 0.1 0.0  # with SWE KGE (DE only)

setlocal

if "%~1"=="" (
    echo Error: provide w_swe_nrmse ^(e.g. 1.0^)
    exit /b 1
)
if "%~2"=="" (
    echo Error: provide w_rho_nrmse ^(e.g. 0.0^)
    exit /b 1
)
if "%~3"=="" (
    echo Error: provide w_swe_bias  ^(e.g. 0.0^)
    exit /b 1
)
if "%~4"=="" (
    echo Error: provide w_rho_bias  ^(e.g. 0.0^)
    exit /b 1
)

set W1=%~1
set W2=%~2
set W3=%~3
set W4=%~4
set W5=0.0
set W6=0.0
if not "%~5"=="" set W5=%~5
if not "%~6"=="" set W6=%~6

set BASE=%~dp0
if "%BASE:~-1%"=="\" set BASE=%BASE:~0,-1%

echo ========================================================
echo   DeltaSnow calibration -- Win21 -- Dynamic RhoMax
echo   SWE_NRMSE=%W1%  RHO_NRMSE=%W2%  SWE_BIAS=%W3%  RHO_BIAS=%W4%
echo   KGE_SWE=%W5%    KGE_RHO=%W6%    (NM + DE)
echo ========================================================

echo.
echo [1/2] Win21 -- Nelder-Mead
Rscript "%BASE%\calibration_Win21\WIN21_dsnow_paprameter_optimization.R" %W1% %W2% %W3% %W4% %W5% %W6%
if errorlevel 1 exit /b 1

echo.
echo [2/2] Win21 -- Differential Evolution
Rscript "%BASE%\calibration_Win21\WIN21_dsnow_parameter_opitmization_DE.R" %W1% %W2% %W3% %W4% %W5% %W6%
if errorlevel 1 exit /b 1

echo.
echo ========================================================
echo   All calibrations complete.
echo   Results saved to:
echo     calibration_Win21\data\R_opt_logs\
echo     calibration_Win21\data\R_opt_logs_DE\
echo ========================================================

endlocal
