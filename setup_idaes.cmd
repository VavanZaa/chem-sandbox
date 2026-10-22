@echo off
setlocal
set "ROOT=%~dp0"
set "IDAES=%ROOT%.venv\Scripts\idaes.exe"
set "EXTENSIONS=%ROOT%.tools\idaes-extensions"

if not exist "%IDAES%" (
  echo IDAES is not installed in the project environment. Run setup.ps1 first.
  exit /b 1
)

if exist "%EXTENSIONS%\ipopt.exe" (
  echo Project-local IDAES/IPOPT extensions already exist; preserving them.
) else (
  echo Downloading the official IDAES Windows solver and property-library bundle...
  pushd "%ROOT%"
  "%IDAES%" get-extensions --to=.tools\idaes-extensions
  set "INSTALL_ERROR=%ERRORLEVEL%"
  popd
  if not "%INSTALL_ERROR%"=="0" exit /b %INSTALL_ERROR%
)

if not exist "%EXTENSIONS%\ipopt.exe" (
  echo IDAES extension verification failed: ipopt.exe is missing.
  exit /b 1
)

echo IDAES/IPOPT process cross-check runtime is ready.
endlocal
