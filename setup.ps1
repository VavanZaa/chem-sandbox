param([switch]$CpuOnly)
$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$NodeVersion = "22.22.2"
$NodeFolder = "node-v$NodeVersion-win-x64"
$NodeRoot = Join-Path $ProjectRoot ".tools\$NodeFolder"
$NodeArchive = Join-Path $ProjectRoot ".tools\$NodeFolder.zip"
$NodeHash = "7c93e9d92bf68c07182b471aa187e35ee6cd08ef0f24ab060dfff605fcc1c57c"

if (-not (Test-Path (Join-Path $NodeRoot "node.exe"))) {
    New-Item -ItemType Directory -Force -Path (Join-Path $ProjectRoot ".tools") | Out-Null
    Invoke-WebRequest "https://nodejs.org/dist/v$NodeVersion/$NodeFolder.zip" -OutFile $NodeArchive
    $ActualHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $NodeArchive).Hash.ToLower()
    if ($ActualHash -ne $NodeHash) { throw "Node.js archive checksum mismatch" }
    Expand-Archive -LiteralPath $NodeArchive -DestinationPath (Join-Path $ProjectRoot ".tools") -Force
    Remove-Item -LiteralPath $NodeArchive -Force
}

if (-not (Test-Path (Join-Path $ProjectRoot ".venv\Scripts\python.exe"))) {
    python -m venv (Join-Path $ProjectRoot ".venv")
}
$PythonExe = Join-Path $ProjectRoot ".venv\Scripts\python.exe"
& $PythonExe -m pip install --upgrade pip "setuptools>=77,<82" wheel
if ($LASTEXITCODE -ne 0) { throw "Python packaging tools installation failed" }
$HasNvidia = $null -ne (Get-Command nvidia-smi -ErrorAction SilentlyContinue)
if ($HasNvidia -and -not $CpuOnly) {
    Write-Host "Installing CUDA 12.8 PyTorch for NVIDIA GPU inference (large download)..."
    & $PythonExe -m pip install --index-url https://download.pytorch.org/whl/cu128 "torch==2.11.0+cu128"
} else {
    & $PythonExe -m pip install --index-url https://download.pytorch.org/whl/cpu "torch==2.11.0+cpu"
}
if ($LASTEXITCODE -ne 0) { throw "PyTorch installation failed" }
& $PythonExe -m pip install -e "$ProjectRoot[test]"
if ($LASTEXITCODE -ne 0) { throw "Project dependencies installation failed" }

$env:PATH = "$NodeRoot;$env:PATH"
Push-Location (Join-Path $ProjectRoot "frontend")
try { & (Join-Path $NodeRoot "npm.cmd") install } finally { Pop-Location }

Write-Host "Dependencies installed. Verifying model checkpoints..." -ForegroundColor Cyan
$env:PYTHONPATH = Join-Path $ProjectRoot "backend"
if ($HasNvidia -and -not $CpuOnly) {
    & $PythonExe (Join-Path $ProjectRoot "backend\scripts\verify_engines.py") cuda
} else {
    & $PythonExe (Join-Path $ProjectRoot "backend\scripts\verify_engines.py") cpu mace_mp ani2x
}
if ($LASTEXITCODE -ne 0) { throw "Model checkpoint verification failed" }
& (Join-Path $ProjectRoot "setup_reaktoro.cmd")
if ($LASTEXITCODE -ne 0) { throw "Isolated Reaktoro runtime setup failed" }
& (Join-Path $ProjectRoot "setup_idaes.cmd")
if ($LASTEXITCODE -ne 0) { throw "Project-local IDAES solver setup failed" }
Write-Host "Setup complete. Run .\start.cmd (no global execution-policy change required)." -ForegroundColor Green
