@echo off
REM Runs all calibration phases sequentially.
REM Each phase calls run_calibration.bat with a specific weight combination.
REM
REM Weight arguments (positional):
REM   1  w_swe_nrmse   weight for NRMSE of SWE
REM   2  w_rho_nrmse   weight for NRMSE of bulk density
REM   3  w_swe_bias    weight for NBIAS of SWE
REM   4  w_rho_bias    weight for NBIAS of bulk density
REM   5  w_kge_swe     weight for (1-KGE) of SWE             (optional, default 0.0; NM + DE)
REM   6  w_kge_rho     weight for (1-KGE) of bulk density    (optional, default 0.0; NM + DE)
REM
REM Usage: run_all_phases.bat
REM Estimated runtime: several hours

setlocal

set BASE=%~dp0
if "%BASE:~-1%"=="\" set BASE=%BASE:~0,-1%
set RUN=%BASE%\run_calibration.bat

REM args: label  w_swe_nrmse  w_rho_nrmse  w_swe_bias  w_rho_bias  [w_kge_swe]  [w_kge_rho]

REM ------------------------------------------------------------
REM PHASE 1 -- NRMSE only, no bias
REM ------------------------------------------------------------
call :run_phase "Phase 1A -- SWE-only baseline"           1.0 0.0 0.0 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 1B -- Density-only extreme"        0.0 1.0 0.0 0.0
if errorlevel 1 exit /b 1

REM ------------------------------------------------------------
REM PHASE 2 -- Balanced NRMSE sweep (no bias)
REM ------------------------------------------------------------
call :run_phase "Phase 2A -- SWE-dominant   (0.7/0.3)"   0.7 0.3 0.0 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 2B -- Equal weight   (0.5/0.5)"   0.5 0.5 0.0 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 2C -- Density-dom.   (0.3/0.7)"   0.3 0.7 0.0 0.0
if errorlevel 1 exit /b 1

REM ------------------------------------------------------------
REM PHASE 3 -- SWE bias penalty (NRMSE + SWE bias, no RHO bias)
REM ------------------------------------------------------------
call :run_phase "Phase 3A -- Balanced NRMSE + SWE bias   (0.6/0.2/0.2/0.0)" 0.6 0.2 0.2 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 3B -- SWE-dominant + SWE bias     (0.7/0.0/0.3/0.0)" 0.7 0.0 0.3 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 3C -- Density-dominant + SWE bias (0.3/0.5/0.2/0.0)" 0.3 0.5 0.2 0.0
if errorlevel 1 exit /b 1

REM ------------------------------------------------------------
REM PHASE 4 -- RHO bias penalty (NRMSE + RHO bias, no SWE bias)
REM ------------------------------------------------------------
call :run_phase "Phase 4A -- Balanced NRMSE + RHO bias   (0.6/0.2/0.0/0.2)" 0.6 0.2 0.0 0.2
if errorlevel 1 exit /b 1
call :run_phase "Phase 4B -- SWE-dominant + RHO bias     (0.7/0.1/0.0/0.2)" 0.7 0.1 0.0 0.2
if errorlevel 1 exit /b 1
call :run_phase "Phase 4C -- Density-dominant + RHO bias (0.3/0.5/0.0/0.2)" 0.3 0.5 0.0 0.2
if errorlevel 1 exit /b 1

REM ------------------------------------------------------------
REM PHASE 5 -- Both bias terms (NRMSE + SWE bias + RHO bias)
REM ------------------------------------------------------------
call :run_phase "Phase 5A -- Balanced NRMSE + both bias   (0.4/0.4/0.1/0.1)"   0.4  0.4  0.10 0.10
if errorlevel 1 exit /b 1
call :run_phase "Phase 5B -- SWE-dominant + both bias     (0.8/0.1/0.05/0.05)" 0.8  0.1  0.05 0.05
if errorlevel 1 exit /b 1
call :run_phase "Phase 5C -- Density-dominant + both bias (0.1/0.8/0.05/0.05)" 0.1  0.8  0.05 0.05
if errorlevel 1 exit /b 1
call :run_phase "Phase 5D -- Bias-focused                 (0.25/0.25/0.25/0.25)" 0.25 0.25 0.25 0.25
if errorlevel 1 exit /b 1
call :run_phase "Phase 5E -- Bias-only                    (0.0/0.0/0.5/0.5)"  0.0  0.0  0.50 0.50
if errorlevel 1 exit /b 1

REM ------------------------------------------------------------
REM PHASE 6 -- KGE (one run per variable, symmetric substitution)
REM ------------------------------------------------------------
call :run_phase "Phase 6A -- KGE SWE + NRMSE RHO balanced (0.0/0.5/0.0/0.0/0.5/0.0)" 0.0 0.5 0.0 0.0 0.5 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 6B -- NRMSE SWE + KGE RHO balanced (0.5/0.0/0.0/0.0/0.0/0.5)" 0.5 0.0 0.0 0.0 0.0 0.5
if errorlevel 1 exit /b 1
call :run_phase "Phase 6C -- KGE SWE                      (0.0/0.0/0.0/0.0/1.0/0.0)" 0.0 0.0 0.0 0.0 1.0 0.0
if errorlevel 1 exit /b 1
call :run_phase "Phase 6D -- KGE RHO                      (0.0/0.0/0.0/0.0/0.0/1.0)" 0.0 0.0 0.0 0.0 0.0 1.0
if errorlevel 1 exit /b 1

goto :eof

:run_phase
REM %~1=label  %~2=w1  %~3=w2  %~4=w3  %~5=w4  %~6=w5(optional)  %~7=w6(optional)
set _LABEL=%~1
set _W1=%~2
set _W2=%~3
set _W3=%~4
set _W4=%~5
set _W5=0.0
set _W6=0.0
if not "%~6"=="" set _W5=%~6
if not "%~7"=="" set _W6=%~7

echo.
echo ###################################################
echo #  %_LABEL%
echo #  SWE_NRMSE=%_W1%  RHO_NRMSE=%_W2%  SWE_BIAS=%_W3%  RHO_BIAS=%_W4%
echo #  KGE_SWE=%_W5%    KGE_RHO=%_W6%    (NM + DE)
echo ###################################################
call "%RUN%" %_W1% %_W2% %_W3% %_W4% %_W5% %_W6%
exit /b %errorlevel%
