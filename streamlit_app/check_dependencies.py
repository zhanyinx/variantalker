#!/usr/bin/env python3
"""Report whether *this* interpreter has what MAFigate needs to run.

Standalone and stdlib-only, so it can run before anything is installed and on a bare
machine — the same reason ``vendor/_sync.py --check`` is stdlib-only:

    python3 check_dependencies.py            # name the interpreter, list every requirement
    python3 check_dependencies.py --quiet     # say nothing unless something is missing

Exits 1 if any declared runtime dependency is absent. ``setup.sh`` runs it after
installing and ``run_mafigate.sh`` before launching.

**Why this exists.** ``requirements.txt`` declaring a package is not the same as the
package being installed, and "installed" is not a property of a machine — it is a
property of one interpreter. Issue #162 was exactly that gap. On the machine where it
was measured, ``pip`` was miniconda's and ``streamlit`` was Homebrew's: two separate
entries on ``PATH``, belonging to two different pythons. So ``setup.sh``'s bare
``pip install -r requirements.txt`` filled an interpreter the app never ran in, and
nothing failed. ``plotly`` was simply absent everywhere the app looked, every one of the
six Summary charts fell back to ``st.bar_chart``, and the only symptom that reached a
reader was an x-axis reading ``{"left": 0.0982, "right": 0.14}``.

That is why this check names the interpreter it is reporting on, out loud, rather than
just saying yes or no. The interpreter *is* the answer.

**Presence, deliberately not versions.** Every line in ``requirements.txt`` now pins one
exact version (``plotly==6.9.0``, issue #256) — a far stronger claim than the floors the
file used to carry. This check still reports only whether a distribution is installed at
all, which is the gap that was costing the reader a chart, and it stays that way for two
reasons. It runs *before* anything is installed, on a bare machine, where the interesting
answer is "nothing is here" rather than "the wrong thing is here". And comparing versions
means parsing specifiers and ordering releases, which is a packaging library's job and not
something to hand-roll in a stdlib-only script that has to run everywhere this one does.

The trade, said out loud: an interpreter holding some version other than the pinned one
reads as satisfied here. What narrows that is the pin itself — every channel installs from
this file, so a mismatch means either a hand-installed package or a venv older than the
last pin change. Absence is the failure that actually happened (#162); a mismatch has not.

**Distributions, not imports.** The check asks ``importlib.metadata`` what is installed
rather than importing anything, so it needs no table mapping ``PyYAML`` to ``yaml`` and
``streamlit-aggrid`` to ``st_aggrid`` — a table that would be one more thing able to
drift away from ``requirements.txt``. The trade is that a package installed but broken
reads as present here; absence was the failure that happened.
"""

import argparse
import sys
from importlib import metadata
from pathlib import Path

#: Read from the file beside this script, not from the working directory, so both
#: ``./setup.sh`` and ``make -C streamlit_app run`` see the same requirements.
REQUIREMENTS = Path(__file__).resolve().parent / "requirements.txt"

#: The characters that end a distribution name in a requirements line — version
#: specifiers, extras, and environment markers all start with one of these.
_NAME_ENDS = "<>=!~;[, \t"


def declared_requirements(path=REQUIREMENTS):
    """The distribution names ``requirements.txt`` declares for running the app.

    A ``#`` and everything after it is dropped before the line is read, which is what
    keeps the development extras out without naming them: ``pytest``, ``black`` and the
    rest are written there as commented ``pip install`` hints, not as requirements, and
    the app runs without them. It also handles this file's continuation comments, whose
    lines are whitespace followed by ``#`` and so reduce to nothing.
    """
    names = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        name = line
        for index, char in enumerate(line):
            if char in _NAME_ENDS:
                name = line[:index]
                break
        if name:
            names.append(name)
    return names


def missing_requirements(names):
    """Those of ``names`` that are not installed for the running interpreter."""
    absent = []
    for name in names:
        try:
            metadata.version(name)
        except metadata.PackageNotFoundError:
            absent.append(name)
    return absent


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Report whether this interpreter has what MAFigate needs to run."
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="print nothing unless a dependency is missing",
    )
    args = parser.parse_args(argv)

    declared = declared_requirements()
    absent = missing_requirements(declared)

    if not args.quiet:
        print(f"🐍 Interpreter: {sys.executable}")
        for name in declared:
            try:
                print(f"   ✅ {name} {metadata.version(name)}")
            except metadata.PackageNotFoundError:
                print(f"   ❌ {name} — not installed for this interpreter")

    if absent:
        print("")
        print(f"❌ {len(absent)} of {len(declared)} declared dependencies are missing from")
        print(f"   {sys.executable}:")
        print(f"      {', '.join(absent)}")
        print("")
        print("   Install them into *that* interpreter — not whichever pip comes first on")
        print("   PATH, which need not belong to it:")
        print(f"      {sys.executable} -m pip install -r {REQUIREMENTS}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
