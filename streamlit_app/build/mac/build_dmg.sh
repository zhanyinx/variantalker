#!/bin/bash
# Build MAFigate.dmg for macOS distribution
#
# Usage: ./build_dmg.sh [--arch ARCH]
#
#   --arch    Target architecture: arm64, x86_64, or universal (default: universal)
#
# Examples:
#   ./build_dmg.sh                  # Build universal (Intel + Apple Silicon)
#   ./build_dmg.sh --arch arm64     # Build for Apple Silicon only
#   ./build_dmg.sh --arch x86_64    # Build for Intel only
#   ./build_dmg.sh --arch universal # Build both (larger DMG)
#
# THERE IS NO VERSION ARGUMENT (issue #260). The version comes from
# config/constants.py's APP_VERSION, read by build/version.py, and it is the same number
# the app's About dialog shows, the Windows installer stamps, and the release tag carries.
# The argument used to default to a literal a whole major version below what the app
# reported; a build could also be handed any number at all from the command line, which is
# two ways to ship a DMG whose filename disagrees with the app inside it.
#
# The DMG bundles a portable Python 3.11 — users need ZERO prerequisites.
#
# Prerequisites (build machine only):
#   - macOS with Xcode Command Line Tools
#   - curl (pre-installed on macOS)
#   - python3, to read the version (macOS ships it with the CLT)
#   - Optional: create-dmg (brew install create-dmg) for prettier DMGs
#
# Output: MAFigate-<version>-macOS-<arch>.dmg

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}"
PROJECT_DIR="${SCRIPT_DIR}/../.."
PYTHON="${PYTHON:-python3}"

# --- Read the version from the app (issue #260) ---
#
# Split from the assignment on purpose. `VERSION="$(...)"` on its own is checked by `set -e`
# — a failing substitution aborts the script — which makes a following `[ -z "$VERSION" ]`
# unreachable, and an unreachable check is decoration. This shape makes both cases real: a
# non-zero exit, and the rc-0-but-silent case a future version.py could introduce.
if ! VERSION="$("${PYTHON}" "${PROJECT_DIR}/build/version.py")" || [ -z "${VERSION}" ]; then
    echo "❌ could not read APP_VERSION from config/constants.py." >&2
    echo "   ${PYTHON} must be a Python 3 and build/version.py must be present." >&2
    exit 1
fi

# --- Parse arguments, and refuse anything that is not one ---
#
# `./build_dmg.sh <version> --arch arm64` is how this script was called until #260, and it
# is still written that way in older tickets and in anyone's shell history. The loop used to
# end in `*) shift ;;` — silently dropping whatever it did not recognise — so a leftover
# version produced a DMG labelled with a number the caller never chose and said nothing.
#
# Refusing every unrecognised token rather than only a leading one, because the version was
# not always leading: with the flag first, the trailing version hit that same catch-all.
TARGET_ARCH=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arch)
            TARGET_ARCH="$2"; shift 2 ;;
        *)
            echo "❌ build_dmg.sh does not understand '$1'." >&2
            echo "   It no longer takes a VERSION argument: the version comes from" >&2
            echo "   config/constants.py's APP_VERSION (this build is ${VERSION})." >&2
            echo "   Usage: ./build_dmg.sh [--arch arm64|x86_64|universal]" >&2
            exit 1 ;;
    esac
done

# Default to a universal build so the DMG runs on both Apple Silicon and Intel
# Macs. (A single-arch build crashes on the other architecture when shared.)
if [ -z "${TARGET_ARCH}" ]; then
    TARGET_ARCH="universal"
fi

# --- Python build-standalone configuration ---
# The release itself is pinned in build/python_release.env, which the Windows build reads
# too (issue #261). Typed here as well, the two could be edited apart and ship one version
# number over two different interpreters.
PYTHON_PIN="${SCRIPT_DIR}/../python_release.env"
if [ ! -f "${PYTHON_PIN}" ]; then
    echo "ERROR: ${PYTHON_PIN} is missing — it is where the bundled CPython is pinned." >&2
    exit 1
fi
# shellcheck source=../python_release.env
source "${PYTHON_PIN}"
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

# --- Stamp the version into the bundle (issue #260) ---
# The template carries a placeholder in both keys rather than a number, so the version the
# Finder shows cannot drift from the one the app's About dialog shows. Both keys, because
# they answer different questions: CFBundleShortVersionString is what a human reads in Get
# Info, CFBundleVersion is what the OS compares between builds.
#
# Then checked, not assumed. `plutil -replace` on a key the template has stopped carrying
# would fail, `set -e` would stop the build, and this grep would never run — but a template
# that grew a *third* placeholder, or a key renamed rather than removed, would leave
# `__APP_VERSION__` in a shipped bundle where it reads as the app's version.
INFO_PLIST="${APP_BUNDLE}/Contents/Info.plist"
plutil -replace CFBundleShortVersionString -string "${VERSION}" "${INFO_PLIST}"
plutil -replace CFBundleVersion -string "${VERSION}" "${INFO_PLIST}"
if grep -qF '__APP_VERSION__' "${INFO_PLIST}"; then
    echo "❌ Info.plist still holds an unreplaced __APP_VERSION__ placeholder." >&2
    grep -nF '__APP_VERSION__' "${INFO_PLIST}" >&2
    exit 1
fi

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

# --- Stamp the build (issue #263) ---
#
# Writes the gitignored config/build_stamp.py, which the About dialog reports; see
# config/build_identity.py for why the channel and the commit are there at all. Two things
# about it are local to this script: it runs *before* the copy list below, which is what puts
# the stamp inside the bundle (that list copies `config` whole, so nothing here names the
# stamp file), and its failure is checked rather than trusted, because a .dmg whose About
# dialog calls itself a source checkout is a silent loss of the only build record there is.
# A build machine with no `git` is not a failure — build_stamp.py says what it did.
if ! "${PYTHON}" "${PROJECT_DIR}/build/build_stamp.py" --channel macos-dmg; then
    echo "❌ could not write the build stamp (build/build_stamp.py)." >&2
    echo "   Without it this .dmg would report itself as a source checkout." >&2
    exit 1
fi

# --- Copy the Streamlit app source ---
echo "Copying application files..."
DEST="${APP_BUNDLE}/Contents/Resources/streamlit_app"

# `vendor` holds the copies of the pipeline's filter code (see vendor/README.md).
# Without it the packaged app has no filter code at all. `make build-mac` gates on
# `check-vendor`, so a drifted copy cannot reach this point.
#
# `.streamlit` is the config that turns Streamlit's usage reporting off, and it has to
# travel *beside MAFigate.py* — that is how Streamlit finds it, since this bundle's
# launcher runs from whatever working directory Finder hands it rather than from the app
# directory (issue #259). A copy list is a deny-by-default list, so a dotted directory
# nobody named here is simply absent from the .dmg; tests/test_telemetry_config.py stands
# this list up and asks Streamlit what it resolves from the result.
for item in MAFigate.py requirements.txt .streamlit components config filters page_modules utils vendor; do
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

# --- Stage the first-open note beside the app ---
# Gatekeeper's alert fires when the app is *opened*, after this window is already on screen,
# so the mounted DMG still reaches a macOS recipient in time. See build/BUILD_INSTRUCTIONS.md
# for why Windows has no equivalent surface.
#
# Copied as .txt rather than kept as .md: a Mac with no developer tooling has no default
# application for .md, and a note that will not open when double-clicked is not a note.
#
# One wording, one home — build/OPENING_MAFIGATE.md. Do not restate any of it in this
# script; what that replaced was two contradictory `xattr` incantations in two files, and
# tests/test_unsigned_artifact_copy.py now fails on a second copy.
cp "${PROJECT_DIR}/build/OPENING_MAFIGATE.md" "${STAGING_DIR}/Opening MAFigate.txt"

# --- Ad-hoc code signing ---
# This is NOT a Developer ID signature and does NOT satisfy Gatekeeper for
# distribution (recipients still need to bypass quarantine — see the note staged
# above). It does give every nested binary a valid ad-hoc signature, which avoids
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
        --window-size 600 440 \
        --icon-size 100 \
        --icon "MAFigate.app" 150 190 \
        --app-drop-link 450 190 \
        --icon "Opening MAFigate.txt" 300 340 \
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
echo "     Blocked on first open? That is expected; the DMG carries the steps."
echo ""
echo "First launch will install pip dependencies (~1 min, needs internet)."
echo "Subsequent launches start in seconds."
echo ""
echo "IMPORTANT — sharing this DMG:"
echo "  This build is ad-hoc signed only (not Developer ID-signed / notarized), so"
echo "  macOS Gatekeeper will quarantine it on every other machine. What recipients"
echo "  should do about that is written once, in build/OPENING_MAFIGATE.md, and rides"
echo "  in this DMG as \"Opening MAFigate.txt\" beside the app. Send them to that,"
echo "  and quote the same file on the release page for Windows recipients — do not"
echo "  paste the steps into a third place."
echo "  For frictionless sharing, sign with a Developer ID certificate and notarize."
