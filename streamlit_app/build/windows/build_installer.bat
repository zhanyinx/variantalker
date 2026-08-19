@echo off
REM Build MAFigate Windows Installer
REM
REM The .exe bundles a portable CPython, so users need NO prerequisites of their own
REM (issue #261). The release is pinned in build\python_release.env, which the macOS build
REM reads too, so the two installers cannot drift to different interpreters.
REM
REM x64 only. ARM Windows runs this build under emulation; a native ARM variant is not
REM built; adding one means shipping an interpreter nobody has launched.
REM
REM Prerequisites:
REM   1. Install Inno Setup from https://jrsoftware.org/isinfo.php
REM   2. Add Inno Setup to your PATH, or set ISCC below
REM   3. curl and tar, both shipped with Windows 10 1803 and later
REM   4. An internet connection on the first build (the download is cached afterwards)
REM   5. Python 3 on PATH - only to read the version, not to build anything
REM
REM Output: output\MAFigate-<version>-Windows-Setup.exe, where <version> is
REM config\constants.py's APP_VERSION and is never typed anywhere else (issue #260).
REM
REM Written with `goto` rather than parenthesised `if (...)` blocks, deliberately. cmd.exe
REM expands `%VAR%` when it parses a whole block, which breaks this script in two ways that
REM are invisible until they happen: `%errorlevel%` inside a block is the value from *before*
REM the block ran, and an unquoted path containing a parenthesis - `C:\Users\me\repo (old)` -
REM closes the block early. One `:fail` label also means the exit path is written once.

setlocal

pushd "%~dp0"

set ISCC="C:\Program Files (x86)\Inno Setup 6\ISCC.exe"

REM Check for Inno Setup
if exist %ISCC% goto have_iscc
where iscc >nul 2>&1
if errorlevel 1 goto no_iscc
set ISCC=iscc
goto have_iscc

:no_iscc
echo ERROR: Inno Setup not found.
echo Install from: https://jrsoftware.org/isinfo.php
echo Then re-run this script.
goto fail

:have_iscc

REM --- Read the pinned CPython release ---
REM One place, shared with build\mac\build_dmg.sh. eol=# skips that file's comments, and
REM the two `if` guards mean a line this loop does not recognise is ignored rather than
REM assigned somewhere unexpected. usebackq is what lets the quoted path be read.
set PYTHON_PIN=%~dp0..\python_release.env
if not exist "%PYTHON_PIN%" goto no_pin

set PYTHON_VERSION=
set PYTHON_TAG=
for /f "usebackq eol=# tokens=1,2 delims==" %%a in ("%PYTHON_PIN%") do (
    if /i "%%a"=="PYTHON_VERSION" set PYTHON_VERSION=%%b
    if /i "%%a"=="PYTHON_TAG" set PYTHON_TAG=%%b
)

if not defined PYTHON_VERSION goto bad_pin
if not defined PYTHON_TAG goto bad_pin
goto pin_read

:no_pin
echo ERROR: build\python_release.env is missing.
echo        It is where the bundled CPython is pinned, for this build and the macOS one.
goto fail

:bad_pin
echo ERROR: build\python_release.env declares no PYTHON_VERSION or no PYTHON_TAG.
echo        Expected two KEY=value lines, unquoted, with no spaces around the `=`.
goto fail

:pin_read

REM No `-shared` in the triple, though issue #261 asked for one. That component only
REM appears on this project's `full` archives; every Windows `install_only` build is the
REM shared one, and `x86_64-pc-windows-msvc-shared-install_only_stripped.tar.gz` is not a
REM published asset of any release. Checked against the pinned release's asset list, where
REM the three Windows install_only_stripped archives are aarch64, i686 and this one.
set PYTHON_TRIPLE=x86_64-pc-windows-msvc-install_only_stripped
set PYTHON_ARCHIVE=cpython-%PYTHON_VERSION%+%PYTHON_TAG%-%PYTHON_TRIPLE%.tar.gz
set PYTHON_BASE_URL=https://github.com/astral-sh/python-build-standalone/releases/download/%PYTHON_TAG%
set PYTHON_CACHE_DIR=%~dp0.python_cache
set PYTHON_ARCHIVE_PATH=%PYTHON_CACHE_DIR%\%PYTHON_ARCHIVE%

REM --- Derive the app version (issue #260) ---
REM
REM The number lives in config\constants.py and reaches this installer only through
REM build\version.py, which writes the version.iss that installer.iss includes. It is
REM regenerated on every build on purpose: an include left over from an earlier checkout
REM would stamp this installer with that checkout's version, silently, which is the
REM failure mode of every generated file that is written once.
REM
REM Nothing here supplies a fallback version. Without this step installer.iss stops with
REM an #error, which is the intended outcome - an installer with no version in Programs
REM and Features, or several releases sharing one filename, is worse than a failed build.
REM
REM Written on `goto` and `if errorlevel 1` like the rest of this file, not the
REM parenthesised `if %errorlevel% neq 0 (...)` it was first drafted as: the note at the
REM top applies to this section too, and tests\test_windows_installer.py enforces it.
set PYTHON=python
%PYTHON% --version >nul 2>&1
if errorlevel 1 set PYTHON=py -3

%PYTHON% ..\version.py --write-iss
if errorlevel 1 goto no_version

REM `for /f` swallows a failure - it simply loops zero times - so the variable is checked
REM rather than trusted: unchecked, a failure here printed "MAFigate--Windows-Setup.exe"
REM and looked like a successful build.
set MAFIGATE_VERSION=
for /f %%v in ('%PYTHON% ..\version.py') do set MAFIGATE_VERSION=%%v
if not defined MAFIGATE_VERSION goto no_version
goto version_read

:no_version
echo ERROR: could not read APP_VERSION from config\constants.py.
echo        Python 3 must be on PATH. This build reads the version; it never guesses one.
goto fail

:version_read

REM --- Stamp the build (issue #263) ---
REM
REM Writes the gitignored config\build_stamp.py, which the About dialog reports; see
REM config\build_identity.py for why the channel and the commit are there at all. Two things
REM about it are local to this script: it runs *before* ISCC, which is what puts the stamp
REM inside the installer (installer.iss ships config\* with recursesubdirs, so nothing here
REM names the stamp file), and its failure is checked rather than trusted, because an .exe
REM whose About dialog calls itself a source checkout is a silent loss of the only build
REM record there is. A build machine with no `git` is not a failure - build_stamp.py says so.
%PYTHON% ..\build_stamp.py --channel windows-installer
if errorlevel 1 goto no_stamp
goto stamped

:no_stamp
echo ERROR: could not write the build stamp (build\build_stamp.py).
echo        Without it this installer would report itself as a source checkout.
goto fail

:stamped

echo === Building MAFigate Windows Installer %MAFIGATE_VERSION% ===
echo     Python: %PYTHON_VERSION% %PYTHON_TRIPLE%
echo.

REM --- Download the interpreter, cached between builds ---
if not exist "%PYTHON_CACHE_DIR%" mkdir "%PYTHON_CACHE_DIR%"
if exist "%PYTHON_ARCHIVE_PATH%" goto have_archive

echo   Downloading Python %PYTHON_VERSION%...
REM --fail so an error page is not cached as if it were a tarball.
curl -L --fail --progress-bar -o "%PYTHON_ARCHIVE_PATH%" "%PYTHON_BASE_URL%/%PYTHON_ARCHIVE%"
if errorlevel 1 goto download_failed
goto extract

:have_archive
echo   Using cached Python %PYTHON_VERSION%
goto extract

:download_failed
REM The partial file goes, so the next run downloads rather than "using the cache".
if exist "%PYTHON_ARCHIVE_PATH%" del /q "%PYTHON_ARCHIVE_PATH%"
echo ERROR: failed to download %PYTHON_ARCHIVE%
goto fail

:extract

REM --- Stage the interpreter beside installer.iss ---
REM The archive's own root directory is `python`, which is what installer.iss ships as
REM {app}\python and what launch.bat then runs. Rebuilt from scratch every time: a stale
REM tree left over from an earlier pin would be shipped without anything noticing.
REM
REM The --exclude list is build_dmg.sh's, component for component, because an installer that
REM shipped CPython's test suite, idlelib and tkinter would be the same release packaged
REM differently on the two platforms - 73 MB unpacked instead of 56, none of the difference
REM reachable from this app. tests\test_windows_installer.py holds the two lists together, so
REM a component trimmed from one build cannot be left in the other.
REM
REM Both tars are bsdtar (Windows 10 1803+ ships it, macOS is bsdtar too), so the patterns
REM behave the same; these were checked against this very archive before being written here.
if exist "python" rmdir /s /q "python"
tar -xzf "%PYTHON_ARCHIVE_PATH%" ^
    --exclude=*/idlelib --exclude=*/idlelib/* ^
    --exclude=*/idle_test --exclude=*/idle_test/* ^
    --exclude=*/tkinter --exclude=*/tkinter/* ^
    --exclude=*/turtledemo --exclude=*/turtledemo/* ^
    --exclude=*/test --exclude=*/test/* ^
    --exclude=*/tests --exclude=*/tests/* ^
    --exclude=*/__pycache__ --exclude=*/__pycache__/* ^
    --exclude=*.pyc
if errorlevel 1 goto extract_failed
if not exist "python\python.exe" goto wrong_layout
goto compile

:extract_failed
echo ERROR: failed to extract %PYTHON_ARCHIVE%
goto fail

:wrong_layout
echo ERROR: %PYTHON_ARCHIVE% did not contain python\python.exe.
echo        The release layout changed; check build\python_release.env.
goto fail

:compile

REM Create output directory
if not exist output mkdir output

REM Run Inno Setup Compiler
%ISCC% installer.iss
if errorlevel 1 goto build_failed

echo.
echo === Build complete ===
echo Output: output\MAFigate-%MAFIGATE_VERSION%-Windows-Setup.exe
echo Bundled Python: %PYTHON_VERSION% - the installed app needs no Python of its own.
popd
pause
exit /b 0

:build_failed
echo.
echo === Build FAILED ===

:fail
popd
pause
exit /b 1
