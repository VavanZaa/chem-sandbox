@echo off
setlocal
set "ROOT=%~dp0"
set "TOOLS=%ROOT%.tools"
set "MICROMAMBA=%TOOLS%\micromamba.exe"
set "MAMBA_ROOT_PREFIX=%TOOLS%\mamba-root"
set "ENV_PREFIX=%TOOLS%\reaktoro-env"

if not exist "%TOOLS%" mkdir "%TOOLS%"

if not exist "%MICROMAMBA%" (
  echo Downloading the official micromamba win-64 package...
  curl.exe -L --fail --retry 3 -o "%TOOLS%\micromamba-win-64.tar.bz2" https://micro.mamba.pm/api/micromamba/win-64/latest
  if errorlevel 1 exit /b 1
  if not exist "%TOOLS%\micromamba-package" mkdir "%TOOLS%\micromamba-package"
  tar -xf "%TOOLS%\micromamba-win-64.tar.bz2" -C "%TOOLS%\micromamba-package"
  if errorlevel 1 exit /b 1
  copy /Y "%TOOLS%\micromamba-package\Library\bin\micromamba.exe" "%MICROMAMBA%" >nul
  if errorlevel 1 exit /b 1
)

if not exist "%ENV_PREFIX%\python.exe" (
  echo Creating the isolated Reaktoro 2.13 environment...
  "%MICROMAMBA%" create -y -p "%ENV_PREFIX%" -f "%ROOT%runtimes\reaktoro\environment.yml"
  if errorlevel 1 exit /b 1
) else (
  echo Reaktoro environment already exists; preserving the installed lock state.
)

echo Verifying a real electrolyte equilibrium calculation...
pushd "%ROOT%runtimes\reaktoro"
"%ENV_PREFIX%\python.exe" verify_runtime.py
set "VERIFY_ERROR=%ERRORLEVEL%"
popd
if not "%VERIFY_ERROR%"=="0" exit /b %VERIFY_ERROR%

echo Reaktoro runtime is ready. start.cmd will launch it automatically.
endlocal
