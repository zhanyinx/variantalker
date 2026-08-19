@echo off
REM MAFigate Windows Launcher
REM Finds Python, sets up venv, installs deps, and launches Streamlit

title MAFigate
setlocal enabledelayedexpansion

set APP_NAME=MAFigate
set APP_DIR=%~dp0streamlit_app
set VENV_DIR=%USERPROFILE%\.mafigate\venv
set LOG_DIR=%USERPROFILE%\.mafigate\logs
set REQUIREMENTS=%APP_DIR%\requirements.txt
set PORT=8501

if not exist "%LOG_DIR%" mkdir "%LOG_DIR%"

REM --- Find Python 3 ---
set PYTHON=
where python >nul 2>&1
if %errorlevel% equ 0 (
    for /f "tokens=*" %%i in ('python -c "import sys; print(sys.version_info.major)" 2^>nul') do (
        if "%%i"=="3" set PYTHON=python
    )
)
if not defined PYTHON (
    where python3 >nul 2>&1
    if %errorlevel% equ 0 set PYTHON=python3
)
if not defined PYTHON (
    REM Check common install locations
    if exist "C:\Python311\python.exe" set PYTHON=C:\Python311\python.exe
    if exist "C:\Python310\python.exe" set PYTHON=C:\Python310\python.exe
    if exist "%LOCALAPPDATA%\Programs\Python\Python311\python.exe" set PYTHON=%LOCALAPPDATA%\Programs\Python\Python311\python.exe
    if exist "%LOCALAPPDATA%\Programs\Python\Python310\python.exe" set PYTHON=%LOCALAPPDATA%\Programs\Python\Python310\python.exe
)
if not defined PYTHON (
    echo.
    echo ERROR: Python 3.9 or later is required but was not found.
    echo.
    echo Please install Python from: https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)

echo Using Python: %PYTHON%

REM --- Create virtual environment ---
if not exist "%VENV_DIR%\Scripts\python.exe" (
    echo Setting up %APP_NAME% environment ^(first run, may take a minute^)...
    %PYTHON% -m venv "%VENV_DIR%"
    if %errorlevel% neq 0 (
        echo ERROR: Failed to create virtual environment.
        pause
        exit /b 1
    )
)

set VENV_PIP=%VENV_DIR%\Scripts\pip.exe
set VENV_STREAMLIT=%VENV_DIR%\Scripts\streamlit.exe

REM --- Install dependencies ---
set MARKER=%VENV_DIR%\.deps_installed
if not exist "%MARKER%" (
    echo Installing dependencies...
    "%VENV_PIP%" install -q -r "%REQUIREMENTS%"
    if %errorlevel% neq 0 (
        echo ERROR: Failed to install dependencies.
        pause
        exit /b 1
    )
    echo. > "%MARKER%"
)

REM --- Find free port ---
:find_port
netstat -an | findstr ":%PORT% " | findstr "LISTENING" >nul 2>&1
if %errorlevel% equ 0 (
    set /a PORT+=1
    if %PORT% gtr 8510 (
        echo ERROR: No free port found.
        pause
        exit /b 1
    )
    goto find_port
)

REM --- Launch Streamlit ---
echo.
echo Starting %APP_NAME% on http://127.0.0.1:%PORT%
echo Press Ctrl+C to stop.
echo.

set STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

REM Open browser after a short delay
start "" cmd /c "timeout /t 5 /nobreak >nul && start http://127.0.0.1:%PORT%"

"%VENV_STREAMLIT%" run "%APP_DIR%\MAFigate.py" ^
    --server.port=%PORT% ^
    --server.address=127.0.0.1 ^
    --server.headless=true

pause
