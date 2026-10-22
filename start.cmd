@echo off
setlocal
set "ROOT=%~dp0"
set "PYTHON=%ROOT%.venv\Scripts\python.exe"
set "NODE_ROOT=%ROOT%.tools\node-v22.22.2-win-x64"
set "REAKTORO_PYTHON=%ROOT%.tools\reaktoro-env\python.exe"

if not exist "%PYTHON%" (
  echo Python environment is missing. Run the project setup first.
  pause
  exit /b 1
)

if not exist "%NODE_ROOT%\npm.cmd" (
  echo Project Node runtime is missing. Run the project setup first.
  pause
  exit /b 1
)

if exist "%REAKTORO_PYTHON%" (
  set "CHEM_SANDBOX_REAKTORO_URL=http://127.0.0.1:8011"
  start "Chem Sandbox Electrolytes" /min /d "%ROOT%runtimes\reaktoro" "%REAKTORO_PYTHON%" -m uvicorn service:app --host 127.0.0.1 --port 8011
) else (
  echo Optional Reaktoro runtime is missing. Run setup_reaktoro.cmd to enable electrolyte equilibrium.
)

start "Chem Sandbox API" /min /d "%ROOT%backend" "%PYTHON%" -m uvicorn app.main:app --host 127.0.0.1 --port 8000
start "Chem Sandbox UI" /min /d "%ROOT%frontend" "%NODE_ROOT%\node.exe" "node_modules\vite\bin\vite.js" --host 127.0.0.1 --port 5173

echo Chem Sandbox is starting...
echo App:      http://localhost:5173
echo API docs: http://127.0.0.1:8000/docs
ping 127.0.0.1 -n 4 >nul
if /i not "%CHEM_SANDBOX_NO_BROWSER%"=="1" start "" "http://localhost:5173"
endlocal
