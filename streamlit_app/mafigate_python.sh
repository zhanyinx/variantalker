# Which python MAFigate installs into and runs in.
#
# Sourced by setup.sh and run_mafigate.sh, both of which read
# `PYTHON="${MAFIGATE_PYTHON:-python3}"` first and then let this file have the last word.
# Not executable and not meant to be run: it defines where the virtual environment lives
# and the rule for using it, in one place, so the script that builds the environment and
# the script that launches into it cannot name different directories.
#
# That "one place" is the whole point. Issue #162 was two commands — `pip` and `streamlit`
# — that were each resolved separately and turned out to disagree, so every dependency
# installed into an interpreter the app never ran in and nothing failed. A venv path
# written out twice is the same shape of mistake waiting to be made again.
#
# Issue #258: setup.sh builds this environment instead of installing into whatever python3
# is first on PATH. A collaborator's first run of MAFigate must not be able to break the
# interpreter their other work depends on, and the desktop installers already build an
# isolated environment of their own — so both channels now make one identical promise about
# the machine rather than two different ones.
#
# The environment lives in the checkout, deliberately, and not in the `~/.mafigate` the
# installers use: that path is shared, and two checkouts pointing at it would fight over
# one environment (issue #229 left the choice to here and named this the safer one).
#
# MAFIGATE_PYTHON short-circuits all of it. Name an interpreter and both scripts use it
# untouched, with no environment built and none consulted — a conda environment, a specific
# version, a venv you keep somewhere else. The venv is a default, not a cage.

# Sourced by `sh` as well as by bash — the Makefile reads both facts below out of this file
# rather than restating them — so nothing here uses a bashism, and the two internal variables
# are prefixed rather than declared `local`, which is not POSIX.

# Relative to the app directory, which both scripts cd into before sourcing this file.
#
# Deliberately an assignment and not a `${MAFIGATE_VENV:-...}` knob, despite reading like one
# beside MAFIGATE_PYTHON: *where* the environment goes is the decision this file records, and
# an override could put it back in the shared path two checkouts would fight over. Anyone who
# wants their environment elsewhere names the interpreter instead, which is MAFIGATE_PYTHON.
MAFIGATE_VENV=".venv"

# Where the interpreter sits inside it: `bin` on POSIX, `Scripts` on Windows, which is what
# the same `python3 -m venv` writes when the bash running these scripts is Git Bash.
#
# Tested with -x rather than -e, because a venv outlives the python that built it — a
# Homebrew or apt upgrade is enough — and what is left behind is a symlink pointing at
# nothing. Presence would hand the launcher an interpreter that cannot start; executability
# follows the link and says no, which is what makes re-running setup.sh the repair.
mafigate_venv_python() {
    for _mafigate_candidate in \
        "$MAFIGATE_VENV/bin/python" \
        "$MAFIGATE_VENV/Scripts/python.exe"
    do
        if [ -x "$_mafigate_candidate" ]; then
            printf '%s\n' "$_mafigate_candidate"
            return 0
        fi
    done
    return 1
}

# Point $PYTHON at the environment. Returns non-zero — leaving $PYTHON alone — whenever there
# is nothing to point it at: either the caller named an interpreter, or no environment exists
# yet. It does not distinguish those two, because a caller that cares has to test
# MAFIGATE_PYTHON anyway to know which message to print, and both scripts do.
mafigate_use_venv() {
    if [ -n "${MAFIGATE_PYTHON:-}" ]; then
        return 1
    fi
    _mafigate_venv_python=$(mafigate_venv_python) || return 1
    PYTHON="$_mafigate_venv_python"
}

# Build it, with $PYTHON — through `-m venv` for the same reason everything else here is a
# `-m`: PATH gets no vote in which python this is. Idempotent by construction: run against
# an existing directory, `-m venv` repairs the interpreter links and leaves the installed
# packages alone, which is why setup.sh can be re-run rather than requiring a rm -rf.
mafigate_create_venv() {
    if "$PYTHON" -m venv "$MAFIGATE_VENV"; then
        return 0
    fi
    echo "❌ '$PYTHON -m venv $MAFIGATE_VENV' failed." >&2
    echo "   Some distributions ship venv separately — on Debian and Ubuntu it is the" >&2
    echo "   python3-venv package. Or point MAFIGATE_PYTHON at a python that has it," >&2
    echo "   or at an environment you manage yourself, and nothing is built here." >&2
    return 1
}
