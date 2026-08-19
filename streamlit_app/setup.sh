#!/bin/bash
#
# Setup script for MAFigate.
#
# Builds ONE interpreter, installs into it, names it out loud, and then checks that what
# requirements.txt asked for actually arrived there.
#
# Installing into one interpreter, naming it, and checking it are all the lesson of issue
# #162. `pip` and `streamlit` are separate entries on PATH and need not belong to the same
# python — on the machine where #162 was
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
# Issue #258: that interpreter is a virtual environment this script builds, in the
# checkout, rather than whichever python3 happens to be first on PATH. Installing MAFigate
# must not be able to disturb the interpreter a collaborator's other work depends on — and
# the desktop installers, from this same codebase, already build an isolated environment,
# so both routes now promise the same thing about the machine.
#
# Where it lives and how it is found is mafigate_python.sh's to say, sourced below and by
# run_mafigate.sh, so that the environment built here is by construction the environment
# launched from. MAFIGATE_PYTHON still short-circuits it entirely.

set -euo pipefail

# Run from the app directory whatever the caller's working directory is, so
# requirements.txt and check_dependencies.py are found.
cd "$(dirname "$0")"

PYTHON="${MAFIGATE_PYTHON:-python3}"

# One file decides where the virtual environment is and when it is used. Sourced, not run.
. ./mafigate_python.sh

echo "🚀 Setting up MAFigate environment..."

if ! command -v "$PYTHON" > /dev/null 2>&1; then
    echo "❌ '$PYTHON' not found."
    echo "   Set MAFIGATE_PYTHON to the python that should run MAFigate."
    exit 1
fi

if [ -n "${MAFIGATE_PYTHON:-}" ]; then
    echo "🔧 MAFIGATE_PYTHON is set, so no virtual environment is built or used."
else
    if mafigate_venv_python > /dev/null; then
        echo "♻️  Reusing the virtual environment in $MAFIGATE_VENV/"
    else
        echo "🧪 Building a virtual environment in $MAFIGATE_VENV/ (your python is untouched)..."
        mafigate_create_venv
    fi
    # One error path for both branches, rather than trusting the branch above to have
    # produced something usable: an environment whose interpreter cannot be run is exactly
    # what a python upgrade leaves behind, and the pip install below would be the one to
    # discover it.
    if ! mafigate_use_venv; then
        echo "❌ No interpreter can be run in $MAFIGATE_VENV/."
        echo "   Remove the directory and run ./setup.sh again, or set MAFIGATE_PYTHON."
        exit 1
    fi
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
echo "Or, to reach the same interpreter by hand (from this directory):"
echo "$PYTHON -m streamlit run MAFigate.py"
