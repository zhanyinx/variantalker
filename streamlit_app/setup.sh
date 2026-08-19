#!/bin/bash
#
# Setup script for MAFigate.
#
# Installs into ONE interpreter, names it out loud, and then checks that what
# requirements.txt asked for actually arrived there.
#
# All three halves are the lesson of issue #162. `pip` and `streamlit` are separate
# entries on PATH and need not belong to the same python — on the machine where #162 was
# measured, `pip` was miniconda's and `streamlit` was Homebrew's. So this script's old
# bare `pip install -r requirements.txt` filled an interpreter the app never ran in, and
# nothing failed: `plotly` was simply absent wherever the app looked, and all six Summary
# charts drew a worse chart in silence for as long as that lasted.
#
# `$PYTHON -m pip` cannot disagree with itself the way two PATH entries can. And
# run_mafigate.sh reads the same variable with the same default, so installing and
# launching land in the same place by construction rather than by luck:
#
#     MAFIGATE_PYTHON=/path/to/python ./setup.sh
#
# No virtual environment is created here. That is a choice about how the developer's
# machine is arranged rather than a bug to fix, and it is left to whoever arranges it —
# point MAFIGATE_PYTHON at a venv's python and this script installs into that venv.

set -euo pipefail

# Run from the app directory whatever the caller's working directory is, so
# requirements.txt and check_dependencies.py are found.
cd "$(dirname "$0")"

PYTHON="${MAFIGATE_PYTHON:-python3}"

echo "🚀 Setting up MAFigate environment..."

if ! command -v "$PYTHON" > /dev/null 2>&1; then
    echo "❌ '$PYTHON' not found."
    echo "   Set MAFIGATE_PYTHON to the python that should run MAFigate."
    exit 1
fi

echo "🐍 Installing into: $("$PYTHON" -c 'import sys; print(sys.executable)')"
echo "📦 Installing dependencies..."
"$PYTHON" -m pip install -r requirements.txt

echo ""
echo "🔍 Checking that what requirements.txt asked for arrived..."
"$PYTHON" check_dependencies.py

echo ""
echo "✅ Setup complete!"
echo ""
echo "To run the application:"
echo "./run_mafigate.sh"
echo ""
echo "Or, to reach the same interpreter by hand:"
echo "$PYTHON -m streamlit run MAFigate.py"
