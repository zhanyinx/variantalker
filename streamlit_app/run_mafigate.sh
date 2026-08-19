#!/bin/bash
#
# Run script for MAFigate.
#
# Launches through `$PYTHON -m streamlit` rather than the bare `streamlit` on PATH, and
# refuses to start if a declared dependency is missing from that interpreter.
#
# Issue #162: the old script checked `command -v streamlit` and launched. That check can
# pass while the app is broken, because `streamlit` and `pip` are separate PATH entries
# that need not belong to the same python — on the machine where #162 was measured they
# did not. The app booted in an interpreter with no `plotly`, every Summary chart quietly
# fell back to `st.bar_chart`, and the reader was left with an x-axis reading
# `{"left": 0.0982, "right": 0.14}`. A launcher that starts an app whose charts are
# broken, silently, is the same failure in a different place.
#
# `$PYTHON -m streamlit` cannot make that mistake: the interpreter that is checked is the
# interpreter that runs. It is the same variable, with the same default, that setup.sh
# installs into.
#
# Issue #258: by default that variable now names the virtual environment setup.sh builds in
# this checkout. Which environment, and when it is used, is mafigate_python.sh's to say —
# sourced by both scripts, so the one this launches into is by construction the one that was
# filled, rather than two paths that happen to match today.
#
#     MAFIGATE_PYTHON=/path/to/python ./run_mafigate.sh
#
# The refusal has an escape hatch, because a launcher nobody can override is its own kind
# of trap — a half-working app still beats no app when someone is mid-analysis:
#
#     MAFIGATE_SKIP_DEPS_CHECK=1 ./run_mafigate.sh

set -euo pipefail

# Run from the app directory whatever the caller's working directory is, so MAFigate.py
# and check_dependencies.py are found.
cd "$(dirname "$0")"

PYTHON="${MAFIGATE_PYTHON:-python3}"

# One file decides where the virtual environment is and when it is used. Sourced, not run.
. ./mafigate_python.sh

echo "🚀 Starting MAFigate..."

if mafigate_use_venv; then
    echo "🧪 Using the virtual environment in $MAFIGATE_VENV/"
elif [ -z "${MAFIGATE_PYTHON:-}" ]; then
    # No environment yet, and none named: fall back to the bare default and let the
    # dependency check below be the one that refuses. It has the better message — it can
    # say what is missing — and a launcher that dies on an absent directory would be
    # unhelpful to anyone who installed the requirements some other way.
    echo "⚠️  No virtual environment in $MAFIGATE_VENV/ — falling back to '$PYTHON'."
    echo "   ./setup.sh builds one."
fi

if ! command -v "$PYTHON" > /dev/null 2>&1; then
    echo "❌ '$PYTHON' not found."
    echo "   Set MAFIGATE_PYTHON to the python that should run MAFigate."
    exit 1
fi

echo "🐍 Interpreter: $("$PYTHON" -c 'import sys; print(sys.executable)')"

if [ "${MAFIGATE_SKIP_DEPS_CHECK:-0}" = "1" ]; then
    echo "⚠️  Skipping the dependency check (MAFIGATE_SKIP_DEPS_CHECK=1)."
elif ! "$PYTHON" check_dependencies.py --quiet; then
    echo ""
    echo "   Refusing to launch: MAFigate would run with parts of itself missing."
    echo "   Run ./setup.sh to install them, or set MAFIGATE_SKIP_DEPS_CHECK=1 to"
    echo "   launch anyway."
    exit 1
fi

# Run the application
#
# Loopback, like both packaged launchers (`build/mac/.../launch.sh`,
# `build/windows/launch.bat`). This line read `0.0.0.0` until now, which publishes MAFigate on
# every interface — so anyone who could route to this machine on port 8501 could drive it, and
# open patient data through it. That was never a decision anyone made; it was the developer
# script disagreeing with the two launchers a user actually gets. Kept as one flag rather than
# the machinery that briefly stood here, on the dev's instruction (issue #182).
echo "🌐 Starting Streamlit application..."
exec "$PYTHON" -m streamlit run MAFigate.py --server.port 8501 --server.address 127.0.0.1
