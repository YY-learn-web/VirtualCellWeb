@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PYTHON_EXE=D:\anaconda3\envs\cell-prj-env\python.exe"

if not exist "%PYTHON_EXE%" (
  set "PYTHON_EXE=python"
)

cd /d "%SCRIPT_DIR%"
echo Starting Virtual Cell API with %PYTHON_EXE%
"%PYTHON_EXE%" main.py
