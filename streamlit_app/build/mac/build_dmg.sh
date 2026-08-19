#!/bin/bash
# Build MAFigate.dmg for macOS distribution
#
# Usage: ./build_dmg.sh [VERSION] [--arch ARCH]
#
#   VERSION   App version (default: 1.0.0)
#   --arch    Target architecture: arm64, x86_64, or universal (default: universal)
#
# Examples:
#   ./build_dmg.sh 1.0.0                  # Build universal (Intel + Apple Silicon)
#   ./build_dmg.sh 1.0.0 --arch arm64     # Build for Apple Silicon only
#   ./build_dmg.sh 1.0.0 --arch x86_64    # Build for Intel only
#   ./build_dmg.sh 1.0.0 --arch universal # Build both (larger DMG)
#
# The DMG bundles a portable Python 3.11 — users need ZERO prerequisites.
#
# Prerequisites (build machine only):
#   - macOS with Xcode Command Line Tools
#   - curl (pre-installed on macOS)
#   - Optional: create-dmg (brew install create-dmg) for prettier DMGs
#
# Output: MAFigate-<version>-macOS-<arch>.dmg

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}"
PROJECT_DIR="${SCRIPT_DIR}/../.."
VERSION="${1:-1.0.0}"

# --- Parse architecture argument ---
shift 2>/dev/null || true
TARGET_ARCH=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arch) TARGET_ARCH="$2"; shift 2 ;;
        *) shift ;;
    esac
done

# Default to a universal build so the DMG runs on both Apple Silicon and Intel
# Macs. (A single-arch build crashes on the other architecture when shared.)
if [ -z "${TARGET_ARCH}" ]; then
    TARGET_ARCH="universal"
fi

# --- Python build-standalone configuration ---
PYTHON_VERSION="3.11.15"
PYTHON_TAG="20260303"
PYTHON_BASE_URL="https://github.com/astral-sh/python-build-standalone/releases/download/${PYTHON_TAG}"
PYTHON_CACHE_DIR="${BUILD_DIR}/.python_cache"

mkdir -p "${PYTHON_CACHE_DIR}"

download_python() {
    local arch="$1"  # aarch64 or x86_64
    local filename="cpython-${PYTHON_VERSION}+${PYTHON_TAG}-${arch}-apple-darwin-install_only_stripped.tar.gz"
    local url="${PYTHON_BASE_URL}/${filename}"
    local cache_path="${PYTHON_CACHE_DIR}/${filename}"

    if [ -f "${cache_path}" ]; then
        echo "  Using cached Python for ${arch}" >&2
    else
        echo "  Downloading Python ${PYTHON_VERSION} for ${arch}..." >&2
        curl -L --progress-bar -o "${cache_path}" "${url}"
    fi
    echo "${cache_path}"
}

extract_python() {
    local tarball="$1"
    local dest="$2"
    mkdir -p "${dest}"
    # Skip the components we trim anyway (test suite, idle, tkinter, …). This
    # also avoids writing files that on-access endpoint-security scanners block
    # on managed Macs — notably idlelib/idle.bat, a Windows batch file — which
    # would otherwise abort extraction with "Operation not permitted".
    tar xzf "${tarball}" -C "${dest}" \
        --exclude '*/idlelib' --exclude '*/idlelib/*' \
        --exclude '*/idle_test' --exclude '*/idle_test/*' \
        --exclude '*/tkinter' --exclude '*/tkinter/*' \
        --exclude '*/turtledemo' --exclude '*/turtledemo/*' \
        --exclude '*/test' --exclude '*/test/*' \
        --exclude '*/tests' --exclude '*/tests/*' \
        --exclude '*/__pycache__' --exclude '*/__pycache__/*' \
        --exclude '*.pyc'
}

# --- Determine architectures to bundle ---
declare -a ARCHS
if [ "${TARGET_ARCH}" = "universal" ]; then
    ARCHS=("aarch64" "x86_64")
    ARCH_LABEL="universal"
else
    if [ "${TARGET_ARCH}" = "arm64" ]; then
        ARCHS=("aarch64")
        ARCH_LABEL="arm64"
    else
        ARCHS=("x86_64")
        ARCH_LABEL="x86_64"
    fi
fi

DMG_NAME="MAFigate-${VERSION}-macOS-${ARCH_LABEL}"
STAGING_DIR="${BUILD_DIR}/staging"
APP_BUNDLE="${STAGING_DIR}/MAFigate.app"

echo "=== Building ${DMG_NAME}.dmg ==="
echo "    Python: ${PYTHON_VERSION}"
echo "    Architecture(s): ${ARCHS[*]}"
echo ""

# --- Clean previous build ---
rm -rf "${STAGING_DIR}"
rm -f "${BUILD_DIR}/${DMG_NAME}.dmg"

# --- Create .app bundle ---
echo "Creating application bundle..."
mkdir -p "${APP_BUNDLE}/Contents/MacOS"
mkdir -p "${APP_BUNDLE}/Contents/Resources/streamlit_app"

# Copy app metadata
cp "${BUILD_DIR}/MAFigate.app/Contents/Info.plist" "${APP_BUNDLE}/Contents/"
cp "${BUILD_DIR}/MAFigate.app/Contents/MacOS/launch.sh" "${APP_BUNDLE}/Contents/MacOS/"
chmod +x "${APP_BUNDLE}/Contents/MacOS/launch.sh"

# --- Build the native launcher ---
# A tiny AppKit wrapper runs a real macOS event loop and supervises launch.sh,
# so the app no longer shows as "not responding" and quits/cleans up cleanly.
# If swiftc is unavailable we fall back to launch.sh as the bundle executable
# (the template default), which still works but exhibits the old quit behaviour.
LAUNCHER_SRC="${BUILD_DIR}/MAFigateLauncher.swift"
LAUNCHER_BIN="${APP_BUNDLE}/Contents/MacOS/MAFigateLauncher"

build_native_launcher() {
    [ -f "${LAUNCHER_SRC}" ] || return 1
    command -v swiftc >/dev/null 2>&1 || return 1

    local -a slices=()
    for arch in "${ARCHS[@]}"; do
        local target_arch min
        if [ "${arch}" = "aarch64" ]; then
            target_arch="arm64"; min="11.0"
        else
            target_arch="x86_64"; min="10.15"
        fi
        local slice="${LAUNCHER_BIN}.${target_arch}"
        if swiftc -O -target "${target_arch}-apple-macos${min}" -o "${slice}" "${LAUNCHER_SRC}" 2>/dev/null; then
            slices+=("${slice}")
        else
            echo "  WARNING: swiftc failed to build the ${target_arch} launcher slice." >&2
        fi
    done

    [ "${#slices[@]}" -gt 0 ] || return 1
    if [ "${#slices[@]}" -eq 1 ]; then
        mv "${slices[0]}" "${LAUNCHER_BIN}"
    else
        lipo -create "${slices[@]}" -o "${LAUNCHER_BIN}" || return 1
        rm -f "${slices[@]}"
    fi
    chmod +x "${LAUNCHER_BIN}"
    return 0
}

if build_native_launcher; then
    echo "Built native launcher (responsive app; clean quit)."
    plutil -replace CFBundleExecutable -string "MAFigateLauncher" "${APP_BUNDLE}/Contents/Info.plist"
else
    echo "  WARNING: native launcher not built (swiftc unavailable?);"
    echo "           falling back to launch.sh — the app may show as 'not responding' until quit."
fi

# Copy icon if it exists
if [ -f "${BUILD_DIR}/MAFigate.app/Contents/Resources/icon.icns" ]; then
    cp "${BUILD_DIR}/MAFigate.app/Contents/Resources/icon.icns" "${APP_BUNDLE}/Contents/Resources/"
fi

# --- Bundle Python ---
echo "Bundling Python ${PYTHON_VERSION}..."
PYTHON_DEST="${APP_BUNDLE}/Contents/Resources/python"

for arch in "${ARCHS[@]}"; do
    tarball=$(download_python "${arch}")
    if [ "${#ARCHS[@]}" -eq 1 ]; then
        # Single architecture: extract directly
        extract_python "${tarball}" "${PYTHON_DEST}_tmp"
        mv "${PYTHON_DEST}_tmp/python" "${PYTHON_DEST}"
        rm -rf "${PYTHON_DEST}_tmp"
    else
        # Universal: extract to arch-specific subdirectory
        extract_python "${tarball}" "${PYTHON_DEST}_tmp_${arch}"
        mkdir -p "${PYTHON_DEST}/${arch}"
        mv "${PYTHON_DEST}_tmp_${arch}/python" "${PYTHON_DEST}/${arch}/"
        rm -rf "${PYTHON_DEST}_tmp_${arch}"
    fi
done

# Trim Python to reduce size — remove test suite, idle, tkinter, etc.
echo "Trimming bundled Python..."
find "${PYTHON_DEST}" -type d -name "test" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "tests" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "idle_test" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "idlelib" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "tkinter" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "turtledemo" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find "${PYTHON_DEST}" -name "*.pyc" -delete 2>/dev/null || true

# --- Copy the Streamlit app source ---
echo "Copying application files..."
DEST="${APP_BUNDLE}/Contents/Resources/streamlit_app"

# `vendor` holds the copies of the pipeline's filter code (see vendor/README.md).
# Without it the packaged app has no filter code at all. `make build-mac` gates on
# `check-vendor`, so a drifted copy cannot reach this point.
for item in MAFigate.py requirements.txt components config filters page_modules utils vendor; do
    if [ -e "${PROJECT_DIR}/${item}" ]; then
        cp -R "${PROJECT_DIR}/${item}" "${DEST}/"
    fi
done

if [ -d "${PROJECT_DIR}/resources" ]; then
    cp -R "${PROJECT_DIR}/resources" "${DEST}/"
fi

# Clean up app source
find "${DEST}" -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
find "${DEST}" -name "*.pyc" -delete 2>/dev/null || true
find "${DEST}" -name ".DS_Store" -delete 2>/dev/null || true
rm -rf "${DEST}/tests" 2>/dev/null || true
rm -rf "${DEST}/build" 2>/dev/null || true
# The drift guard lives in tests/ (removed above) and needs bin/, which the bundle
# does not carry — it is skipif-guarded on bin/ so it would skip rather than fail.
# vendor/_sync.py is the developer-side repair tool and needs bin/ too; strip it so
# the shipped vendor package is nothing but the pipeline code it exists to hold.
rm -f "${DEST}/vendor/_sync.py" 2>/dev/null || true

# --- Ad-hoc code signing ---
# This is NOT a Developer ID signature and does NOT satisfy Gatekeeper for
# distribution (recipients still need to bypass quarantine — see the note at the
# end). It does give every nested binary a valid ad-hoc signature, which avoids
# some "killed on launch" failures after the app is copied to another machine.
# For a truly shareable build, codesign with a Developer ID certificate and
# notarize + staple the app and the DMG.
if command -v codesign >/dev/null 2>&1; then
    echo "Ad-hoc code signing the app bundle..."
    codesign --force --deep --sign - "${APP_BUNDLE}" >/dev/null 2>&1 \
        || echo "  WARNING: ad-hoc codesign failed; continuing with an unsigned bundle."
fi

# --- Create DMG ---
echo "Creating DMG..."

if command -v create-dmg &> /dev/null; then
    create-dmg \
        --volname "${DMG_NAME}" \
        --volicon "${APP_BUNDLE}/Contents/Resources/icon.icns" 2>/dev/null \
        --window-pos 200 120 \
        --window-size 600 400 \
        --icon-size 100 \
        --icon "MAFigate.app" 150 190 \
        --app-drop-link 450 190 \
        --hide-extension "MAFigate.app" \
        "${BUILD_DIR}/${DMG_NAME}.dmg" \
        "${STAGING_DIR}/" \
        || {
            echo "create-dmg failed, falling back to hdiutil..."
            hdiutil create -volname "${DMG_NAME}" \
                -srcfolder "${STAGING_DIR}" \
                -ov -format UDZO \
                "${BUILD_DIR}/${DMG_NAME}.dmg"
        }
else
    ln -sf /Applications "${STAGING_DIR}/Applications"
    hdiutil create -volname "${DMG_NAME}" \
        -srcfolder "${STAGING_DIR}" \
        -ov -format UDZO \
        "${BUILD_DIR}/${DMG_NAME}.dmg"
fi

# --- Cleanup ---
rm -rf "${STAGING_DIR}"

# --- Report ---
DMG_SIZE=$(du -sh "${BUILD_DIR}/${DMG_NAME}.dmg" | cut -f1)

echo ""
echo "=== Build complete ==="
echo "Output: ${BUILD_DIR}/${DMG_NAME}.dmg (${DMG_SIZE})"
echo ""
echo "To install:"
echo "  1. Open ${DMG_NAME}.dmg"
echo "  2. Drag MAFigate.app to Applications"
echo "  3. Double-click MAFigate — no Python installation required!"
echo ""
echo "First launch will install pip dependencies (~1 min, needs internet)."
echo "Subsequent launches start in seconds."
echo ""
echo "IMPORTANT — sharing this DMG:"
echo "  This build is ad-hoc signed only (not Developer ID-signed / notarized),"
echo "  so macOS Gatekeeper will quarantine it on other machines. Recipients must"
echo "  either right-click MAFigate.app -> Open (then confirm) on first launch, or"
echo "  run once after copying to Applications:"
echo "      xattr -dr com.apple.quarantine /Applications/MAFigate.app"
echo "  For frictionless sharing, sign with a Developer ID certificate and notarize."
