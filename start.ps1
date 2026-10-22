$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$NodeRoot = Join-Path $ProjectRoot ".tools\node-v22.22.2-win-x64"
$PythonExe = Join-Path $ProjectRoot ".venv\Scripts\python.exe"

if (-not (Test-Path $PythonExe)) {
    throw "Python environment is missing. Run .\setup.ps1 first."
}
if (-not (Test-Path (Join-Path $NodeRoot "npm.cmd"))) {
    throw "Project Node runtime is missing. Run .\setup.ps1 first."
}

$env:PATH = "$NodeRoot;$env:PATH"
Start-Process -FilePath $PythonExe -ArgumentList "-m", "uvicorn", "app.main:app", "--reload", "--host", "127.0.0.1", "--port", "8000" -WorkingDirectory (Join-Path $ProjectRoot "backend") -WindowStyle Hidden
Start-Process -FilePath (Join-Path $NodeRoot "npm.cmd") -ArgumentList "run", "dev" -WorkingDirectory (Join-Path $ProjectRoot "frontend") -WindowStyle Hidden

Write-Host "Chem Sandbox is starting..." -ForegroundColor Green
Write-Host "App:      http://localhost:5173"
Write-Host "API docs: http://127.0.0.1:8000/docs"
Start-Sleep -Seconds 3
Start-Process "http://localhost:5173"

