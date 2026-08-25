#!/bin/bash
# MAFigate macOS Launcher
# This script is the executable inside the .app bundle.
# It uses the bundled Python (no system Python required), sets up a venv,
# and launches the Streamlit app.

set -e

APP_NAME="MAFigate"
BUNDLE_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
RESOURCES_DIR="${BUNDLE_DIR}/Contents/Resources"
APP_DIR="${RESOURCES_DIR}/streamlit_app"
PYTHON_DIR="${RESOURCES_DIR}/python"
VENV_DIR="${HOME}/.mafigate/venv"
BASE_PYTHON_DIR="${HOME}/.mafigate/python"
LOG_DIR="${HOME}/.mafigate/logs"
REQUIREMENTS="${APP_DIR}/requirements.txt"
PORT=8501

# --- Capture file path from macOS "Open With" ---
# When a user right-clicks a .maf file and opens with MAFigate,
# macOS passes the file path as the first argument.
OPEN_FILE=""
if [ -n "$1" ] && [ -f "$1" ]; then
    OPEN_FILE="$1"
fi

mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/mafigate_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

show_error() {
    osascript -e "display dialog \"$1\" with title \"${APP_NAME} Error\" buttons {\"OK\"} default button \"OK\" with icon stop"
    exit 1
}

show_notification() {
    osascript -e "display notification \"$1\" with title \"${APP_NAME}\""
}

show_progress() {
    osascript -e "display dialog \"$1\" with title \"${APP_NAME}\" buttons {\"OK\"} default button \"OK\" giving up after 2"  2>/dev/null &
}

# --- Find Python 3 (bundled first, then system) ---
find_python() {
    local machine_arch
    machine_arch="$(uname -m)"

    # 1. Check for bundled Python (single-arch build)
    if [ -x "${PYTHON_DIR}/bin/python3" ]; then
        echo "${PYTHON_DIR}/bin/python3"
        return 0
    fi

    # 2. Check for bundled Python (universal build — pick matching arch)
    local bundled_arch
    if [ "${machine_arch}" = "arm64" ]; then
        bundled_arch="aarch64"
    else
        bundled_arch="x86_64"
    fi

    if [ -x "${PYTHON_DIR}/${bundled_arch}/python/bin/python3" ]; then
        echo "${PYTHON_DIR}/${bundled_arch}/python/bin/python3"
        return 0
    fi

    # 3. Universal build — try the other architecture (runs via Rosetta on Apple Silicon)
    if [ "${bundled_arch}" = "aarch64" ] && [ -x "${PYTHON_DIR}/x86_64/python/bin/python3" ]; then
        echo "${PYTHON_DIR}/x86_64/python/bin/python3"
        return 0
    fi
    if [ "${bundled_arch}" = "x86_64" ] && [ -x "${PYTHON_DIR}/aarch64/python/bin/python3" ]; then
        echo "${PYTHON_DIR}/aarch64/python/bin/python3"
        return 0
    fi

    # 4. Fallback to system Python
    for candidate in \
        "/opt/homebrew/bin/python3" \
        "/usr/local/bin/python3" \
        "/usr/bin/python3" \
        "$(which python3 2>/dev/null)"; do
        if [ -x "${candidate}" ] 2>/dev/null; then
            version=$("${candidate}" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null)
            major=$(echo "${version}" | cut -d. -f1)
            minor=$(echo "${version}" | cut -d. -f2)
            if [ "${major}" = "3" ] && [ "${minor}" -ge 9 ]; then
                echo "${candidate}"
                return 0
            fi
        fi
    done

    return 1
}

# --- Keep a stable copy of the bundled interpreter ---
# The .app's own path is ephemeral on macOS: Gatekeeper translocation, running
# from the DMG, and a drag to /Applications all move the bundle between
# launches, and a venv built straight from the bundle pins that path in its
# symlinks and in pyvenv.cfg — which is how v1.0.0's second launch died. So
# the bundled runtime is copied once to a path only this user's home controls,
# and the venv is built from the copy. The copy is staged and moved into place
# atomically, health-probed on every launch, and replaced from the bundle
# (always present) whenever it is broken or belongs to a different build.
ensure_stable_python() {
    local bundle_python="$1"
    local source_stamp="$2"
    local bundle_root="${bundle_python%/bin/python3}"
    local stamp_file="${BASE_PYTHON_DIR}/.source_stamp"

    if [ -x "${BASE_PYTHON_DIR}/bin/python3" ] \
        && "${BASE_PYTHON_DIR}/bin/python3" -c "import sys" >/dev/null 2>&1 \
        && [ "$(cat "${stamp_file}" 2>/dev/null)" = "${source_stamp}" ]; then
        return 0
    fi

    log "Copying the bundled Python runtime to ${BASE_PYTHON_DIR}..."
    show_notification "Setting up the ${APP_NAME} runtime..."
    # A leftover .tmp from an interrupted copy is deleted before the next one.
    rm -rf "${BASE_PYTHON_DIR}" "${BASE_PYTHON_DIR}.tmp"
    mkdir -p "${HOME}/.mafigate"
    if ! cp -R "${bundle_root}" "${BASE_PYTHON_DIR}.tmp" >> "${LOG_FILE}" 2>&1; then
        rm -rf "${BASE_PYTHON_DIR}.tmp"
        show_error "Failed to copy the Python runtime.\n\nCheck the log file at:\n${LOG_FILE}"
    fi
    echo "${source_stamp}" > "${BASE_PYTHON_DIR}.tmp/.source_stamp"
    mv "${BASE_PYTHON_DIR}.tmp" "${BASE_PYTHON_DIR}"
}

# --- Set up virtual environment ---
# Recreate the venv when the interpreter's version OR architecture changes,
# when the base interpreter it was built from moved, or when the existing venv
# is broken. The architecture check matters because a venv left behind by a
# different machine/arch (e.g. a synced home directory) points at an
# interpreter that won't run here and would crash on launch.
#
# Every branch that decides to recreate clears the stale venv first: creating
# a venv over a stale directory skips existing symlinks — dangling ones
# included — and then dies executing them (the v1.0.0 second-launch failure).
#
# Reads ${PYTHON} (set by the caller from find_python) and the path constants
# defined at the top of this file.
ensure_venv() {
    PYTHON_STAMP_FILE="${VENV_DIR}/.python_stamp"
    local interpreter_stamp
    interpreter_stamp=$("${PYTHON}" -c 'import sys, platform; print(f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}-{platform.machine()}")')

    # Build the venv from the stable copy whenever the interpreter came out of
    # the bundle. A system Python does not move with the app and is used as-is.
    local base_python="${PYTHON}"
    case "${PYTHON}" in
        "${BUNDLE_DIR}"/*)
            ensure_stable_python "${PYTHON}" "${interpreter_stamp}"
            base_python="${BASE_PYTHON_DIR}/bin/python3"
            ;;
    esac

    # The stamp carries the base interpreter's path as well as its version and
    # architecture, so every venv built before the stable copy existed
    # (v1.0.0 stamps: version-arch only) mismatches and rebuilds whole.
    CURRENT_STAMP="${interpreter_stamp}:${base_python}"

    NEED_NEW_VENV=false
    if [ ! -x "${VENV_DIR}/bin/python" ]; then
        rm -rf "${VENV_DIR}"
        NEED_NEW_VENV=true
    elif [ ! -f "${PYTHON_STAMP_FILE}" ] || [ "$(cat "${PYTHON_STAMP_FILE}" 2>/dev/null)" != "${CURRENT_STAMP}" ]; then
        log "Python interpreter changed, recreating virtual environment..."
        rm -rf "${VENV_DIR}"
        NEED_NEW_VENV=true
    elif ! "${VENV_DIR}/bin/python" -c "import sys" >/dev/null 2>&1; then
        log "Existing virtual environment is broken, recreating..."
        rm -rf "${VENV_DIR}"
        NEED_NEW_VENV=true
    fi

    if [ "${NEED_NEW_VENV}" = true ]; then
        log "Creating virtual environment (first launch, please wait)..."
        show_notification "Setting up ${APP_NAME} (first launch)..."
        if ! "${base_python}" -m venv "${VENV_DIR}" >> "${LOG_FILE}" 2>&1; then
            show_error "Failed to create the Python environment.\n\nCheck the log file at:\n${LOG_FILE}"
        fi
        echo "${CURRENT_STAMP}" > "${PYTHON_STAMP_FILE}"
        # Fresh environment -> force a dependency (re)install below.
        rm -f "${VENV_DIR}/.deps_installed"
    fi
}

# --- Library seam ---
# Sourcing this file with MAFIGATE_LAUNCH_LIB set loads the definitions above
# and stops here, so the launcher's venv logic can be exercised directly
# without booting the app. Everything below this line only runs when the
# bundle executes the script for real.
if [ -n "${MAFIGATE_LAUNCH_LIB:-}" ]; then
    return 0 2>/dev/null || exit 0
fi

PYTHON=$(find_python) || show_error "Python was not found in the application bundle or on the system.\n\nThis should not happen. Please reinstall MAFigate or contact support."

log "Using Python: ${PYTHON}"
log "Python version: $("${PYTHON}" --version 2>&1)"

ensure_venv

VENV_PYTHON="${VENV_DIR}/bin/python"
VENV_PIP="${VENV_DIR}/bin/pip"
VENV_STREAMLIT="${VENV_DIR}/bin/streamlit"

# --- Install dependencies if needed ---
# Each install is guarded with `if ! ...` so a failure (e.g. no internet on
# first launch) shows a dialog instead of silently aborting under `set -e`.
MARKER="${VENV_DIR}/.deps_installed"
if [ ! -f "${MARKER}" ] || [ "${REQUIREMENTS}" -nt "${MARKER}" ]; then
    log "Installing dependencies..."
    show_notification "Installing ${APP_NAME} dependencies (first launch)..."
    if ! "${VENV_PIP}" install --upgrade pip >> "${LOG_FILE}" 2>&1; then
        show_error "Failed to set up the Python environment (pip upgrade).\n\nAn internet connection is required on first launch.\n\nCheck the log file at:\n${LOG_FILE}"
    fi
    if ! "${VENV_PIP}" install -q -r "${REQUIREMENTS}" >> "${LOG_FILE}" 2>&1; then
        show_error "Failed to install dependencies.\n\nAn internet connection is required on first launch.\n\nCheck the log file at:\n${LOG_FILE}"
    fi
    touch "${MARKER}"
fi

# --- Find free port ---
find_free_port() {
    local port=$1
    while [ $port -lt $(($1 + 10)) ]; do
        if ! lsof -i :${port} > /dev/null 2>&1; then
            echo $port
            return 0
        fi
        port=$((port + 1))
    done
    echo $1
}

PORT=$(find_free_port ${PORT})
URL="http://127.0.0.1:${PORT}"

# --- Launch Streamlit ---
log "Starting ${APP_NAME} on ${URL}"
export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
export STREAMLIT_SERVER_HEADLESS=true

# Pass the file path to the Streamlit app if opened via "Open With"
if [ -n "${OPEN_FILE}" ]; then
    export MAFIGATE_OPEN_FILE="${OPEN_FILE}"
    log "Opening file: ${OPEN_FILE}"
fi

"${VENV_STREAMLIT}" run "${APP_DIR}/MAFigate.py" \
    --server.port="${PORT}" \
    --server.address=127.0.0.1 \
    --server.headless=true \
    >> "${LOG_FILE}" 2>&1 &

SERVER_PID=$!

# Clean up server on exit
trap "kill ${SERVER_PID} 2>/dev/null; exit 0" INT TERM HUP

# --- Wait for server and open browser ---
log "Waiting for server to start..."
for i in $(seq 1 60); do
    if curl -s "http://127.0.0.1:${PORT}/_stcore/health" > /dev/null 2>&1; then
        log "Server is ready!"
        open "${URL}"
        break
    fi
    # Check if process died
    if ! kill -0 ${SERVER_PID} 2>/dev/null; then
        show_error "${APP_NAME} failed to start.\n\nCheck the log file at:\n${LOG_FILE}"
    fi
    sleep 0.5
done

# --- Keep running until server exits ---
wait ${SERVER_PID}
