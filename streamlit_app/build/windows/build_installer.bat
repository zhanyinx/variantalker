@echo off
REM Build MAFigate Windows Installer
REM
REM Prerequisites:
REM   1. Install Inno Setup from https://jrsoftware.org/isinfo.php
REM   2. Add Inno Setup to your PATH, or set ISCC below
REM
REM Output: output\MAFigate-1.0.0-Windows-Setup.exe

setlocal

set ISCC="C:\Program Files (x86)\Inno Setup 6\ISCC.exe"

REM Check for Inno Setup
if not exist %ISCC% (
    where iscc >nul 2>&1
    if %errorlevel% equ 0 (
        set ISCC=iscc
    ) else (
        echo ERROR: Inno Setup not found.
        echo Install from: https://jrsoftware.org/isinfo.php
        echo Then re-run this script.
        pause
        exit /b 1
    )
)

echo === Building MAFigate Windows Installer ===
echo.

REM Create output directory
if not exist output mkdir output

REM Run Inno Setup Compiler
%ISCC% installer.iss

if %errorlevel% equ 0 (
    echo.
    echo === Build complete ===
    echo Output: output\MAFigate-1.0.0-Windows-Setup.exe
) else (
    echo.
    echo === Build FAILED ===
    pause
    exit /b 1
)

pause
