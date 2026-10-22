@echo off
setlocal
cd /d "%~dp0"
echo Close the Chem Sandbox API before installing GPU dependencies.
echo The installer downloads from PyTorch, verifies SHA-256, and changes only this project's environment.
".venv\Scripts\python.exe" "backend\scripts\install_gpu.py"
if errorlevel 1 exit /b 1
echo GPU setup verified. You can run start.cmd again.
