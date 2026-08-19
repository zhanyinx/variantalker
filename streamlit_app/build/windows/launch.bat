@echo off
REM MAFigate Windows Launcher
REM Runs the bundled interpreter, sets up a venv, installs deps, and launches Streamlit.

title MAFigate
setlocal enabledelayedexpansion

set APP_NAME=MAFigate
set APP_DIR=%~dp0streamlit_app
set PYTHON=%~dp0python\python.exe
set VENV_DIR=%USERPROFILE%\.mafigate\venv
set LOG_DIR=%USERPROFILE%\.mafigate\logs
set REQUIREMENTS=%APP_DIR%\requirements.txt
set PORT=8501

if not exist "%LOG_DIR%" mkdir "%LOG_DIR%"

REM --- The bundled interpreter ---
REM The installer ships a portable CPython beside this script, and it is the only one
REM MAFigate runs: there is no search of the environment and no list of likely install
REM directories to guess from (issue #261). Whatever else is on this machine is irrelevant
REM to the app, which is the point - the version it runs on is the version it was tested on.
REM `!PYTHON!` and not `%PYTHON%`: cmd expands `%VAR%` when it parses the whole block, so an
REM unquoted install path containing a parenthesis - the wizard accepts `C:\Program Files
REM (x86)\MAFigate` - would close this block early and garble the one message a user in this
REM state gets. Delayed expansion happens after parsing, which is why it is enabled above.
if not exist "%PYTHON%" (
    echo.
    echo ERROR: MAFigate's bundled Python is missing from this installation:
    echo   !PYTHON!
    echo.
    echo Reinstall MAFigate to restore it.
    echo.
    pause
    exit /b 1
)

echo Using Python: %PYTHON%

REM --- Create virtual environment ---
REM Rebuilt when it is missing or when it cannot run. The second case is not theoretical:
REM an installation upgraded from a version of MAFigate that used a Python of the user's
REM own left a venv behind pointing at that interpreter, and once it is gone the venv is a
REM directory full of files that cannot start. Reusing one would leave the app running on
REM an interpreter this installer never shipped.
set NEED_VENV=
if not exist "%VENV_DIR%\Scripts\python.exe" set NEED_VENV=1
if not defined NEED_VENV (
    "%VENV_DIR%\Scripts\python.exe" -c "import sys" >nul 2>&1
    if errorlevel 1 (
        echo Existing environment cannot run, recreating it...
        rmdir /s /q "%VENV_DIR%"
        set NEED_VENV=1
    )
)

if defined NEED_VENV (
    echo Setting up !APP_NAME! environment ^(first run, may take a minute^)...
    "%PYTHON%" -m venv "%VENV_DIR%"
    if errorlevel 1 (
        echo ERROR: Failed to create virtual environment.
        pause
        exit /b 1
    )
    REM A fresh environment has nothing installed, whatever an older marker claims.
    if exist "%VENV_DIR%\.deps_installed" del /q "%VENV_DIR%\.deps_installed"
)

set VENV_PIP=%VENV_DIR%\Scripts\pip.exe
set VENV_STREAMLIT=%VENV_DIR%\Scripts\streamlit.exe

REM --- Install dependencies ---
set MARKER=%VENV_DIR%\.deps_installed
if not exist "%MARKER%" (
    echo Installing dependencies...
    "%VENV_PIP%" install -q -r "%REQUIREMENTS%"
    if errorlevel 1 (
        echo ERROR: Failed to install dependencies.
        echo An internet connection is required on first launch.
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
