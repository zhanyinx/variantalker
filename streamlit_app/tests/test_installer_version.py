"""Where a MAFigate version number may be written, and where it must be derived.

Issue #260. The app reported **2.0.0**, the Windows installer hardcoded **1.0.0** in two
fields, the DMG script defaulted to **1.0.0**, the macOS bundle's ``Info.plist`` carried
**1.0.0** twice, and the Makefile held a fourth copy — five files, two numbers, nothing
reconciling them. Moment two is about to put that number on an artifact a clinician
downloads and on a release tag, so this module makes :data:`config.constants.APP_VERSION`
the only place it is written.

**The failure this guards is a human typing the number twice**, so the assertions are not
"the installers say what the app says" — two literals agreeing today is exactly the state
the tree was already in, and it drifted anyway. Each installer input is instead required
to carry *no* version literal at all, and to name the derivation that fills it in.

Proving the parse
-----------------
Every field assertion below runs against a directive this module has to find first, and
each test names its anchor — ``AppName=MAFigate`` in the Inno script, ``CFBundleName`` in
the plist, the ``DMG_NAME`` line in the shell script. Prior art is the vendor-drift suite,
which fails when it cannot parse the DMG build script's copy list rather than passing over
an empty one: a guard that silently reads nothing reports success on a tree that has lost
the thing it guards. The same reason gives the prose sweep its two vacuity checks — a
claim of *zero* literals is satisfied by a sweep that opened no files.

Proving the derivation
----------------------
"Derives from the constant" is asserted by **moving the constant**: ``build/version.py``
is run against a throwaway tree whose ``APP_VERSION`` is ``9.9.9`` and has to say
``9.9.9``. A tool that had quietly gone back to printing a literal would agree with the
real constant and disagree with this one. Nothing here compares the tool's output to the
number in this repo and calls that derivation.

What this module cannot check
-----------------------------
The Inno Setup preprocessor half. ``installer.iss`` is compiled by ISCC on Windows, which
is not available here, so the assertions about it are over its text: that both
version-bearing fields expand ``{#AppVersion}``, that the script defines that symbol
nowhere itself, and that it stops with ``#error`` when nothing supplied it. That last one
is the important one and it is still only asserted as text — the compile that would prove
it belongs to the release workflow (#264), which builds on Windows.

Out of scope, deliberately: the comments in ``page_modules/help.py`` and
``tests/test_app_identity.py`` that quote the version-bearing copy issue #71 deleted. They
are records of a past state rather than instructions to a builder, and a sweep that forced
them to move would be asking the tree to forget why it looks the way it does. (The
neighbouring passages in ``config/param_migration.py`` and ``page_modules/param_store.py``
*were* rewritten by #260, because they argued from the app version never having moved and
#260 moved it — a stale record is fine, a stale argument is not.)
"""

from __future__ import annotations

import json
import plistlib
import re
import subprocess
import sys
from pathlib import Path

import pytest

from config.constants import APP_VERSION

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

BUILD_DIR = STREAMLIT_APP / "build"
VERSION_TOOL = BUILD_DIR / "version.py"
BUILD_DOCS = BUILD_DIR / "BUILD_INSTRUCTIONS.md"
DMG_SCRIPT = BUILD_DIR / "mac" / "build_dmg.sh"
INFO_PLIST = BUILD_DIR / "mac" / "MAFigate.app" / "Contents" / "Info.plist"
ISS_SCRIPT = BUILD_DIR / "windows" / "installer.iss"
ISS_BAT = BUILD_DIR / "windows" / "build_installer.bat"
MAKEFILE = STREAMLIT_APP / "Makefile"

#: The workflow that builds both installers on a tag (#264). Named here because it is the
#: newest place a version could be typed, and because the sweep below is the only guard that
#: would notice if one were: see :func:`_build_and_docs_files`.
RELEASE_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "mafigate-release.yml"

#: The release tag namespace, and not a cosmetic choice: the public repo's ``v1.x`` tags
#: are the *Nextflow pipeline's*, on their own cadence, so a MAFigate release sharing that
#: namespace would put a second product in the feed pipeline users watch (#229, #260).
#:
#: Written out again here rather than imported from ``build/version.py``. That is the
#: deliberate kind of duplication: a test that imports the value it is checking asserts
#: only that the module equals itself, which is how a renamed namespace would stay green.
TAG_PREFIX = "mafigate-v"

#: The generated Inno Setup include. Written at build time from ``APP_VERSION`` and never
#: committed — a tracked copy would be the fifth literal this module exists to remove.
GENERATED_ISS = BUILD_DIR / "windows" / "version.iss"

#: A dotted version number. Applied only to fields and to MAFigate-adjacent prose, never
#: to a whole build file: ``installer.iss`` names a Python requirement and the DMG script
#: names macOS deployment targets, and neither is a MAFigate version.
VERSION_LITERAL = re.compile(r"\d+\.\d+(?:\.\d+)?")

#: A MAFigate version written out in prose or in an artifact filename — ``MAFigate v2.0.0``,
#: ``MAFigate-1.0.0-Windows-Setup``, ``mafigate-v1.0.0``. Matched by its *shape* rather
#: than against the two numbers the tree happens to hold, because a sweep pinned to today's
#: literals goes blind the moment one of them is bumped — the hole ``test_app_identity``
#: found in issue #71's version guard.
VERSION_IN_PROSE = re.compile(r"MAFigate[ -]v?(\d+\.\d+(?:\.\d+)?)", re.IGNORECASE)

#: A release tag written out. Held to the derived tag by
#: :func:`test_every_release_tag_named_in_prose_is_the_derived_one` rather than banned, so
#: the sweep below skips it — the app README's "coming with the first `mafigate-v1.0.0`
#: release" is copy worth keeping, and pinning it to ``APP_VERSION`` turns it red on a bump
#: instead of letting it rot. Banned *and* pinned would be two guards contradicting each
#: other, and the one that fired first would decide.
TAG_SHAPED = re.compile(r"mafigate-v\d+\.\d+(?:\.\d+)?", re.IGNORECASE)


def _tracked_files() -> list[str]:
    """Every file git tracks, as repo-relative paths.

    Read through git rather than :meth:`Path.rglob`, for the reason
    ``test_delivery_channels_copy`` gives: ``.claude/worktrees`` holds full checkouts of
    this same tree, so a filesystem walk answers questions about a branch this test is not
    running on.

    This is the fourth near-copy of that helper — ``test_public_export``,
    ``test_delivery_channels_copy`` and ``test_public_repo_name`` each carry their own. Not
    shared, on the same ground those three were not: hoisting it into ``conftest`` would
    edit four guard modules for one refactor, and these modules are worked by parallel
    branches. Left as duplication with the duplication named.
    """
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    tracked = [line for line in completed.stdout.splitlines() if line]
    if not tracked:
        raise AssertionError(
            "git ls-files listed nothing — this guard would pass over any tree at all"
        )
    return tracked


def _version_tool(*args: str, app_root: Path | None = None) -> str:
    """``build/version.py``'s stdout, run through this interpreter."""
    command = [sys.executable, str(VERSION_TOOL), *args]
    if app_root is not None:
        command += ["--app-root", str(app_root)]
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 0, (
        f"build/version.py {' '.join(args)} exited {completed.returncode}:\n"
        f"{completed.stderr}"
    )
    return completed.stdout.strip()


def _moved_constant(tmp_path: Path, version: str = "9.9.9") -> Path:
    """A throwaway app root whose ``APP_VERSION`` is *version*, and nothing else.

    The tool has to read a real ``config/constants.py``, so the fake carries the same one
    name in the same place. Everything the app puts around it is irrelevant to a tool that
    is reading one assignment.
    """
    app_root = tmp_path / "streamlit_app"
    (app_root / "config").mkdir(parents=True)
    (app_root / "config" / "constants.py").write_text(
        f'APP_NAME = "MAFigate"\nAPP_VERSION = "{version}"\n', encoding="utf-8"
    )
    return app_root


# ---------------------------------------------------------------------------
# The tool: one reader of the constant, and it really reads it
# ---------------------------------------------------------------------------


def test_the_version_tool_is_present():
    """Every other test in this module runs it, so its absence must say so once."""
    assert VERSION_TOOL.is_file(), (
        f"{VERSION_TOOL} is missing. It is the only thing that turns APP_VERSION into an "
        "installer's version, so without it both installers are back to a literal."
    )


def test_the_version_tool_reports_the_app_version(tmp_path):
    """And reports it by *reading* it, shown by moving the constant out from under it.

    The first assertion alone would be satisfied by a tool that printed ``2.0.0`` from a
    literal of its own — which is the drift this ticket removes, reintroduced one layer
    down. The second is what makes it a derivation: a tree whose constant says ``9.9.9``
    has to produce ``9.9.9``.
    """
    assert _version_tool() == APP_VERSION

    moved = _moved_constant(tmp_path)
    assert _version_tool(app_root=moved) == "9.9.9", (
        "build/version.py did not follow APP_VERSION to a tree that says 9.9.9, so it is "
        "not reading the constant — it is repeating a number of its own."
    )


#: Each way ``config/constants.py`` can fail to yield a version, and what it looks like.
#: Every one of these must stop a build: an installer stamped with an empty version names no
#: release, and several releases sharing one filename overwrite each other.
BROKEN_CONSTANTS = (
    ("absent", None),
    ("no assignment", 'APP_NAME = "MAFigate"\n'),
    ("not a literal", "APP_VERSION = read_version()\n"),
    ("empty", 'APP_VERSION = ""\n'),
    ("not a version", 'APP_VERSION = "TBD"\n'),
)


@pytest.mark.parametrize(
    "description,source", BROKEN_CONSTANTS, ids=[case[0] for case in BROKEN_CONSTANTS]
)
def test_the_version_tool_refuses_rather_than_guesses(tmp_path, description, source):
    """It exits non-zero and says which file it could not read.

    These branches had no test until the review of #260 pointed out that only the success
    paths were exercised — which is the state in which a refusal quietly becomes a
    ``return ""`` and every caller stamps an empty version instead of failing.
    """
    app_root = tmp_path / "streamlit_app"
    if source is not None:
        (app_root / "config").mkdir(parents=True)
        (app_root / "config" / "constants.py").write_text(source, encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(VERSION_TOOL), "--app-root", str(app_root)],
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0, (
        f"version.py accepted a {description} APP_VERSION and printed "
        f"{completed.stdout.strip()!r}, so a build would be stamped with it"
    )
    assert not completed.stdout.strip(), (
        f"version.py printed {completed.stdout.strip()!r} on the {description} case; a "
        "caller reading stdout would take that as the version"
    )
    assert "APP_VERSION" in completed.stderr


def test_the_release_tag_derives_and_stays_out_of_the_pipeline_namespace(tmp_path):
    """``mafigate-v<version>``, from the same one reader.

    The namespace is asserted as *not* the pipeline's rather than merely as "starts with
    mafigate-v", because the collision is the harm: the public repo's ``v1.x`` line is a
    different product's, released on its own cadence, and a pipeline bugfix landing in the
    same feed would read as a new MAFigate.
    """
    tag = _version_tool("--tag")

    assert tag == f"{TAG_PREFIX}{APP_VERSION}"
    assert not re.fullmatch(r"v\d+.*", tag), (
        f"the release tag {tag!r} is in the Nextflow pipeline's own v1.x namespace"
    )

    moved = _moved_constant(tmp_path)
    assert _version_tool("--tag", app_root=moved) == f"{TAG_PREFIX}9.9.9", (
        "the tag does not follow APP_VERSION, so cutting a release would need a human to "
        "type the number a second time"
    )


# ---------------------------------------------------------------------------
# The Windows installer
# ---------------------------------------------------------------------------


def _iss_setup_directives(text: str) -> dict[str, str]:
    """The ``key=value`` directives in the Inno Setup script's ``[Setup]`` section.

    Section-scoped on purpose: ``[Files]`` and ``[Icons]`` carry their own ``Name:``/
    ``Source:`` syntax, and reading the whole file as ``key=value`` pairs would let a
    directive from any of them answer for a ``[Setup]`` field.
    """
    directives: dict[str, str] = {}
    in_setup = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_setup = stripped.lower() == "[setup]"
            continue
        if not in_setup or not stripped or stripped.startswith(";"):
            continue
        key, separator, value = stripped.partition("=")
        if separator:
            directives[key.strip()] = value.strip()
    return directives


#: The two ``[Setup]`` fields that decide what the installed program calls itself and what
#: the downloaded file is named. Both were literals (#260).
ISS_VERSION_FIELDS = ("AppVersion", "OutputBaseFilename")


def test_the_inno_setup_section_is_still_readable():
    """The anchor for every Windows assertion below.

    ``AppName`` is the canary rather than a field count: it is the one directive whose
    value this module knows independently, so a section that stopped parsing cannot look
    like a section with nothing in it.
    """
    directives = _iss_setup_directives(ISS_SCRIPT.read_text(encoding="utf-8"))

    assert directives.get("AppName") == "MAFigate", (
        "could not parse the [Setup] section out of installer.iss — the script's shape "
        f"changed and this module is no longer reading it. Parsed: {sorted(directives)}"
    )
    for field in ISS_VERSION_FIELDS:
        assert field in directives, (
            f"installer.iss has no {field} directive. Both version-bearing fields must "
            "exist for the assertions about them to mean anything."
        )


@pytest.mark.parametrize("field", ISS_VERSION_FIELDS)
def test_the_inno_version_fields_carry_no_literal(field):
    """Neither the program's version nor the download's filename is typed here."""
    value = _iss_setup_directives(ISS_SCRIPT.read_text(encoding="utf-8"))[field]

    found = VERSION_LITERAL.search(value)
    assert not found, (
        f"installer.iss writes a version literal into {field}: {value!r}. That number is "
        "a copy of APP_VERSION and it has already drifted once — expand {#AppVersion} "
        "instead."
    )


@pytest.mark.parametrize("field", ISS_VERSION_FIELDS)
def test_the_inno_version_fields_expand_the_preprocessor_symbol(field):
    """Absence of a literal is not presence of a derivation.

    A field could satisfy the test above by naming no version at all — dropping
    ``AppVersion`` would leave Programs and Features with a blank column, and dropping it
    from ``OutputBaseFilename`` would ship every release under one filename.
    """
    value = _iss_setup_directives(ISS_SCRIPT.read_text(encoding="utf-8"))[field]

    assert "{#AppVersion}" in value, (
        f"installer.iss's {field} is {value!r}, which names no version at all. It must "
        "expand {#AppVersion}, the symbol build/version.py supplies."
    )


def test_the_inno_script_defines_the_symbol_nowhere_itself():
    """The symbol arrives from outside, or the compile stops.

    A ``#define AppVersion "1.0.0"`` in this file would satisfy every assertion above
    while being the same hardcoded literal one line higher up. The include it reads
    instead is generated and untracked (see :data:`GENERATED_ISS`).
    """
    text = ISS_SCRIPT.read_text(encoding="utf-8")

    own_defines = [
        line.strip()
        for line in text.splitlines()
        if re.match(r"\s*#\s*define\s+AppVersion\b", line)
    ]
    assert not own_defines, (
        f"installer.iss defines AppVersion itself: {own_defines}. The whole point of the "
        "symbol is that the number comes from config/constants.py."
    )


def test_the_inno_script_refuses_to_compile_with_no_version_supplied():
    """A missing derivation must stop the build, not ship an unversioned installer.

    Asserted as text because ISCC runs on Windows and this suite does not. That is a real
    limit and it is why the assertion is on the ``#error`` rather than on a compile: the
    guarantee is only as good as the release workflow that compiles it (#264).

    That workflow now exists, and it closes this gap on the one machine that can: before
    building, its Windows job compiles ``installer.iss`` with no ``version.iss`` present and
    requires ISCC to stop **for this reason**, matching the message below rather than merely
    exiting non-zero — which the missing ``python\\`` tree would also do.
    ``tests/test_release_workflow.py`` holds it to that.
    """
    text = ISS_SCRIPT.read_text(encoding="utf-8")

    assert re.search(r"#\s*ifndef\s+AppVersion", text), (
        "installer.iss does not check whether AppVersion was supplied, so compiling it "
        "without the generated include would substitute an empty version silently."
    )
    assert re.search(r"#\s*error", text), (
        "installer.iss has no #error for the undefined case. Inno Setup expands an "
        "unknown symbol rather than failing, so the check needs teeth."
    )


def test_the_batch_file_generates_the_include_before_it_compiles():
    """Order is the whole guarantee, so it is asserted rather than assumed.

    Generating the include *after* ISCC ran would compile whatever the previous build left
    behind — a stale version, silently, which is the failure mode of every generated file.
    Both line numbers have to be found, so a batch file that stopped doing either fails
    here instead of passing over a file whose shape changed.
    """
    lines = ISS_BAT.read_text(encoding="utf-8").splitlines()

    generate = [
        index
        for index, line in enumerate(lines)
        if "version.py" in line and "--write-iss" in line
    ]
    compile_calls = [
        index
        for index, line in enumerate(lines)
        if re.match(r"\s*%ISCC%\s", line)
    ]

    assert generate, (
        "build_installer.bat never runs build/version.py --write-iss, so the version the "
        "installer is stamped with is whatever version.iss happens to hold."
    )
    assert compile_calls, (
        "could not find the ISCC invocation in build_installer.bat — this test is no "
        "longer reading the file it claims to order."
    )
    assert min(generate) < min(compile_calls), (
        f"build_installer.bat compiles at line {min(compile_calls) + 1} and derives the "
        f"version at line {min(generate) + 1}. The stale include would be the one shipped."
    )


def test_the_batch_file_sets_no_version_of_its_own():
    """It prints the version in two messages, and must not be a place one is typed.

    The batch file had no literal assertion at all until this test: ``set
    MAFIGATE_VERSION=1.0.0`` — exactly the fallback its own comment promises does not exist
    — satisfied every other guard in this module. Swept over all its ``set`` variables
    rather than one name, since a fallback would be introduced under whatever name was
    convenient.
    """
    variables = _batch_sets(ISS_BAT.read_text(encoding="utf-8"))

    assert "ISCC" in variables, (
        "could not parse the batch file's `set` lines — this test is no longer reading it. "
        f"Parsed: {sorted(variables)}"
    )
    offenders = {
        name: values
        for name, values in variables.items()
        if any(VERSION_LITERAL.search(value) for value in values)
    }
    assert not offenders, (
        f"build_installer.bat types a version: {offenders}. It must read every version it "
        "prints from build/version.py, or its messages will outlive the number in them."
    )

    derived = [
        line
        for line in ISS_BAT.read_text(encoding="utf-8").splitlines()
        if "MAFIGATE_VERSION" in line and "version.py" in line
    ]
    assert derived, (
        "build_installer.bat never reads MAFIGATE_VERSION from build/version.py, so the "
        "version in its own output messages comes from nowhere."
    )


def test_the_generated_include_derives_the_version(tmp_path):
    """What ``--write-iss`` writes, and that it too follows a moved constant."""
    destination = tmp_path / "version.iss"
    _version_tool("--write-iss", str(destination))

    defined = re.search(
        r'#\s*define\s+AppVersion\s+"([^"]+)"', destination.read_text(encoding="utf-8")
    )
    assert defined, (
        "build/version.py --write-iss wrote no `#define AppVersion \"…\"`, so the include "
        f"installer.iss reads defines nothing:\n{destination.read_text()}"
    )
    assert defined.group(1) == APP_VERSION

    moved_destination = tmp_path / "moved.iss"
    _version_tool("--write-iss", str(moved_destination), app_root=_moved_constant(tmp_path))
    assert '"9.9.9"' in moved_destination.read_text(encoding="utf-8")


def test_the_generated_include_is_not_committed():
    """Tracked, it would be the copy this module exists to remove.

    Two claims, because either alone is weak: git must not be tracking it *today*, and the
    ignore rule must exist so tomorrow's ``git add -A`` cannot commit it by accident.
    """
    relative = GENERATED_ISS.relative_to(REPO_ROOT).as_posix()

    assert relative not in _tracked_files(), (
        f"{relative} is committed. A generated version file in git is a second place the "
        "number is written, and it will be the stale one."
    )

    ignored = subprocess.run(
        ["git", "check-ignore", "-q", relative],
        cwd=REPO_ROOT,
        capture_output=True,
    )
    assert ignored.returncode == 0, (
        f"{relative} is not ignored, so a `git add -A` after a Windows build commits the "
        "generated version. Add it to .gitignore."
    )


# ---------------------------------------------------------------------------
# The macOS installer
# ---------------------------------------------------------------------------


def _shell_assignments(text: str, name: str) -> list[str]:
    """Every right-hand side assigned to *name* in a shell script.

    Not anchored to the start of the line, because ``VERSION`` is assigned inside an ``if``
    — the shape that keeps the emptiness check reachable under ``set -e``. The preceding
    character must still be a boundary, or ``PYTHON_VERSION="3.11.15"`` would be read as an
    assignment to ``VERSION`` and hand this module's literal check a Python release to fail.
    """
    pattern = re.compile(rf"(?:^|[\s!])({re.escape(name)}=(.*))$")
    return [
        match.group(2).strip()
        for match in (pattern.search(line.rstrip()) for line in text.splitlines())
        if match
    ]


#: A make variable assignment, over all four operators and either spacing. ``VERSION ?= x``,
#: ``VERSION?=x``, ``VERSION := x`` and ``VERSION += x`` are one declaration to make and were
#: four different strings to the matcher this replaced, which recognised two of them — so a
#: reintroduced ``VERSION?=1.0.0`` would have passed the guard that exists to catch it.
MAKE_ASSIGNMENT = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:\?|:|\+)?=\s*(.*)$")


def _make_assignments(text: str) -> dict[str, list[str]]:
    """Every top-level make variable, as ``name -> values assigned to it``.

    Recipe lines are skipped: they begin with a tab and are shell, where ``VERSION=1.0.0``
    would be a local of that one command rather than a variable this file declares.
    """
    assignments: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line or line.startswith(("\t", "#", " ")):
            continue
        match = MAKE_ASSIGNMENT.match(line)
        if match:
            assignments.setdefault(match.group(1), []).append(match.group(2).strip())
    return assignments


def _batch_sets(text: str) -> dict[str, list[str]]:
    """Every ``set NAME=VALUE`` in a batch file, as ``name -> values``."""
    sets: dict[str, list[str]] = {}
    pattern = re.compile(r"^\s*set\s+([A-Za-z_][A-Za-z0-9_]*)=(.*)$", re.IGNORECASE)
    for line in text.splitlines():
        match = pattern.match(line)
        if match:
            sets.setdefault(match.group(1), []).append(match.group(2).strip())
    return sets


def test_the_dmg_script_still_names_its_version_the_same_way():
    """The anchor for the macOS assertions: the two lines they read.

    ``DMG_NAME`` is the canary — it is where the version reaches the downloaded file's
    name — and it must interpolate ``VERSION`` rather than compose the name some other
    way, or the ``VERSION`` assertions below are about a variable nothing uses.
    """
    text = DMG_SCRIPT.read_text(encoding="utf-8")

    assert _shell_assignments(text, "VERSION"), (
        "could not find a VERSION assignment in build_dmg.sh — the script's shape changed "
        "and this module is no longer reading it."
    )
    dmg_names = _shell_assignments(text, "DMG_NAME")
    assert dmg_names, "could not find the DMG_NAME assignment in build_dmg.sh"
    assert all("${VERSION}" in value for value in dmg_names), (
        f"build_dmg.sh's DMG_NAME does not interpolate VERSION: {dmg_names}. The version "
        "assertions below would be about a variable the filename ignores."
    )


def test_the_dmg_script_carries_no_version_literal():
    """Its ``VERSION`` is read from the app, not defaulted to a number.

    The old shape was ``VERSION="${1:-1.0.0}"`` — a literal *and* an invitation to type a
    different one on the command line, which is two ways to ship a mislabelled DMG.
    """
    text = DMG_SCRIPT.read_text(encoding="utf-8")

    for value in _shell_assignments(text, "VERSION"):
        assert not VERSION_LITERAL.search(value), (
            f"build_dmg.sh assigns a version literal: VERSION={value}"
        )
        assert "version.py" in value, (
            f"build_dmg.sh's VERSION={value} does not read the app's constant. It must "
            "come from build/version.py, the one reader."
        )


#: Every way a stale caller can still put a version on the command line. The second is the
#: one a text assertion missed: the refusal originally inspected ``$1`` alone, so a version
#: *after* the flag went on being swallowed by the argument loop's catch-all.
STALE_INVOCATIONS = (
    ["1.0.0"],
    ["1.0.0", "--arch", "arm64"],
    ["--arch", "arm64", "1.0.0"],
)


@pytest.mark.parametrize("arguments", STALE_INVOCATIONS, ids=lambda a: " ".join(a))
def test_the_dmg_script_refuses_a_version_argument(arguments):
    """Run, not read. Removing the argument is not enough — ignoring one silently is a bug.

    Asserted by executing the script, because it is bash and this machine has bash. The
    text assertion this replaces (``"no longer takes a VERSION argument" in text``) passed
    on a script whose ``exit 1`` had been deleted and whose ``echo`` remained — a guard on
    the wording rather than on the refusal.

    Bounded by a timeout, and the bound is not a formality: mutation-testing this guard
    found that a refusal branch which loses its ``exit`` also loses its ``shift``, so the
    argument loop spins forever and an unbounded ``subprocess.run`` hangs the whole suite
    instead of failing it. A timeout is also how the *other* regression is caught cheaply —
    a script that stopped refusing would start a real DMG build here, downloading a Python
    and taking minutes. Neither outcome is "it refused", so both count as failure.
    """
    try:
        completed = subprocess.run(
            ["bash", str(DMG_SCRIPT), *arguments],
            capture_output=True,
            text=True,
            cwd=STREAMLIT_APP,
            timeout=45,
        )
    except subprocess.TimeoutExpired:
        raise AssertionError(
            f"build_dmg.sh neither refused nor finished within 45s given {arguments}. It is "
            "either looping on the argument or has started a real build — both mean the "
            "refusal is gone."
        ) from None

    assert completed.returncode != 0, (
        f"build_dmg.sh accepted {arguments}. The argument loop drops what it does not "
        "recognise, so this build would be labelled with a version nobody chose:\n"
        f"{completed.stdout[-2000:]}"
    )
    assert "VERSION argument" in completed.stderr, (
        "build_dmg.sh failed but did not say why. The caller is following a stale example "
        f"and needs to be told where the version comes from now:\n{completed.stderr}"
    )


def test_the_app_bundle_template_carries_no_version_literal():
    """``Info.plist`` is an installer input too, and it held the number twice.

    This is the copy the ticket did not count: the template is copied into the bundle
    verbatim, so its ``1.0.0`` was what macOS showed in Get Info while the app's own About
    dialog said something else.
    """
    plist = plistlib.loads(INFO_PLIST.read_bytes())

    assert plist.get("CFBundleName") == "MAFigate", (
        "could not read the app bundle template — CFBundleName is "
        f"{plist.get('CFBundleName')!r}, so this test is not reading the plist it claims to."
    )
    for key in ("CFBundleVersion", "CFBundleShortVersionString"):
        assert key in plist, f"Info.plist has no {key}; the build cannot replace what is absent"
        assert not VERSION_LITERAL.search(str(plist[key])), (
            f"Info.plist hardcodes {key}={plist[key]!r}. The template must carry a "
            "placeholder the build replaces from APP_VERSION."
        )


@pytest.mark.parametrize("key", ("CFBundleVersion", "CFBundleShortVersionString"))
def test_the_dmg_script_stamps_both_plist_versions(key):
    """Both keys, from ``VERSION``, in the bundle the DMG ships.

    Both rather than one: ``CFBundleShortVersionString`` is what Finder shows and
    ``CFBundleVersion`` is what the OS compares between builds, and the template has no
    number left for either to fall back on.
    """
    text = DMG_SCRIPT.read_text(encoding="utf-8")

    stamped = [
        line.strip()
        for line in text.splitlines()
        if "plutil" in line and key in line and "${VERSION}" in line
    ]
    assert stamped, (
        f"build_dmg.sh never writes {key} from ${{VERSION}}, so the shipped bundle keeps "
        "the template's placeholder and Get Info shows no version at all."
    )


# ---------------------------------------------------------------------------
# The Makefile, the docs, and the prose
# ---------------------------------------------------------------------------


def test_the_makefile_holds_no_version_and_passes_none(tmp_path):
    """It held the fourth copy, as a ``VERSION ?=`` override the build targets passed on.

    Also asserted: the recipe passes ``build_dmg.sh`` no positional argument. A Makefile
    that kept passing ``$(VERSION)`` to a script which now refuses positionals would fail
    the build rather than mislabel it — but it would fail it on every machine, so it is
    worth catching here.
    """
    text = MAKEFILE.read_text(encoding="utf-8")
    assignments = _make_assignments(text)

    assert "PYTHON" in assignments, (
        "could not parse the Makefile's variables — this test is no longer reading the "
        f"file it claims to. Parsed: {sorted(assignments)}"
    )
    # Every variable, not just one called VERSION: the copy could come back under any name,
    # and nothing else this Makefile declares is a dotted number (PYTHON is an interpreter,
    # COV_MIN an integer), so a blanket sweep costs nothing and closes the renaming hole.
    offenders = {
        name: values
        for name, values in assignments.items()
        if any(VERSION_LITERAL.search(value) for value in values)
    }
    assert not offenders, (
        f"the Makefile holds a version literal: {offenders}. The version has one home, "
        "config/constants.py, and one reader, build/version.py."
    )

    invocations = [
        line.strip()
        for line in text.splitlines()
        if "build_dmg.sh" in line and not line.strip().startswith("#")
    ]
    assert invocations, (
        "the Makefile no longer invokes build_dmg.sh — this test is not reading the "
        "recipe it claims to."
    )
    for line in invocations:
        # The recipe is a shell fragment inside an `if`, so what follows the script name is
        # `; \` — make's line continuation, not an argument. Stripped rather than tolerated
        # by a looser assertion, so a real positional still fails.
        tail = line.split("build_dmg.sh", 1)[1].strip().rstrip("\\").strip()
        arguments = tail.lstrip(";").strip()
        assert not arguments or arguments.startswith("--"), (
            f"the Makefile passes {arguments!r} to build_dmg.sh, which takes no version "
            "argument any more"
        )


#: The files that tell someone how to build MAFigate or what they are downloading. Swept
#: by prefix rather than named one by one, so a file added under ``build/`` — or a second
#: workflow — is covered without anyone remembering to add it here.
#:
#: ``.github/workflows/`` joined the list with the release workflow (#264), and it closes a
#: hole that workflow's own guard cannot: ``test_release_workflow`` proves the version is
#: derived by running the workflow's derive step and comparing what it prints to
#: ``APP_VERSION``, and a step that had gone back to a *literal* would agree with the constant
#: today and only disagree after a bump. This sweep is the half that notices the literal. The
#: two existing workflows come along and cost nothing: neither spells a MAFigate version, and
#: the dotted numbers they do carry — a Python series, an action's major — are not one.
def _build_and_docs_files() -> list[Path]:
    prefixes = (
        "streamlit_app/build/",
        "streamlit_app/Makefile",
        "streamlit_app/README.md",
        ".github/workflows/",
    )
    return [
        REPO_ROOT / relative
        for relative in _tracked_files()
        if relative.startswith(prefixes)
    ]


def _prose_versions(text: str) -> list[str]:
    """Every MAFigate version spelled out, minus the tag-shaped ones.

    Tags are removed before matching rather than filtered after, so ``mafigate-v1.0.0``
    cannot leave a bare ``1.0.0`` behind for the sweep to find. See :data:`TAG_SHAPED` for
    why they are pinned instead of banned.
    """
    return VERSION_IN_PROSE.findall(TAG_SHAPED.sub("", text))


def test_the_prose_sweep_reads_the_files_it_claims_to():
    """Vacuity check one: a claim of zero literals is met by opening nothing.

    Two anchors, because the sweep now spans two trees. The release workflow is named as
    well as the build document, since it is the file whose version the *other* guard can
    only check against today's constant — dropping it from these prefixes would leave a
    hardcoded number in the workflow invisible to both.
    """
    swept = _build_and_docs_files()

    assert BUILD_DOCS in swept, (
        f"the sweep did not pick up {BUILD_DOCS.name}, the document that tells a builder "
        f"what to type. It read {[path.name for path in swept]}."
    )
    assert RELEASE_WORKFLOW in swept, (
        f"the sweep did not pick up {RELEASE_WORKFLOW.name}, which names the artifacts a "
        f"clinician downloads. It read {[path.name for path in swept]}."
    )
    assert len(swept) >= 5, f"only {len(swept)} build files swept: {swept}"


def test_the_prose_sweep_catches_a_planted_version():
    """Vacuity check two: the patterns still recognise what they are looking for.

    Both shapes the tree actually wrote — an artifact filename and a spoken version — and
    the tag, which the sweep must *not* claim because a different test owns it. Asserting
    that third case here is what keeps the exemption honest: if ``TAG_SHAPED`` stopped
    matching, tags would silently fall back into this sweep and the pinning test would look
    like the only guard while both were half-applied.
    """
    for planted in (
        "Run MAFigate-1.0.0-Windows-Setup.exe",
        "MAFigate v2.0.0: Advanced MAF File Analysis",
    ):
        assert _prose_versions(planted), f"the sweep no longer notices {planted!r}"

    tagged = "coming with the first mafigate-v1.0.0 release"
    assert TAG_SHAPED.search(tagged), f"the tag pattern no longer notices {tagged!r}"
    assert not _prose_versions(tagged), (
        "a release tag is being counted as loose prose as well as pinned to the "
        "derivation; the two guards now contradict each other"
    )


def test_no_build_file_spells_a_mafigate_version_out():
    """The number is not in the build instructions, the batch echoes, or the comments.

    Prose is where the copies hid in plain sight — ``./build_dmg.sh 1.0.0`` as a worked
    example, ``Output: MAFigate-1.0.0-Windows-Setup.exe`` echoed twice. None of it is read
    by a build, all of it is read by a human, and a human copying a stale example is how
    the mismatch reaches an artifact.
    """
    offenders = {}
    for path in _build_and_docs_files():
        hits = _prose_versions(path.read_text(encoding="utf-8", errors="replace"))
        if hits:
            offenders[path.relative_to(REPO_ROOT).as_posix()] = sorted(set(hits))

    assert not offenders, (
        f"build files spell a MAFigate version out in prose: {offenders}. Write "
        "`MAFigate-<version>-…` and say where the version comes from."
    )


#: What opens a *full-line* comment in each kind of file the sweep reads, by suffix
#: (files with no suffix — the Makefile, ``.gitignore`` — are ``#`` files too, matched
#: by name below). A version literal on such a line cannot reach an artifact: no build
#: reads it, and no worked example points a human at it. Markdown and the plist are
#: deliberately absent — ``#`` in markdown is a heading, and prose read by humans is
#: exactly where a stale copy does its damage — so in those files every line stays swept.
FULL_LINE_COMMENT_MARKERS = {
    ".yml": ("#",),
    ".yaml": ("#",),
    ".sh": ("#",),
    ".py": ("#",),
    ".env": ("#",),
    ".swift": ("//",),
    ".iss": (";",),
    ".bat": ("rem ", "::"),
}


def _is_full_line_comment(path: Path, line: str) -> bool:
    """Whether *line* is a comment from its first character, in *path*'s own syntax."""
    if path.name in {"Makefile", ".gitignore"}:
        markers = ("#",)
    else:
        markers = FULL_LINE_COMMENT_MARKERS.get(path.suffix, ())
    stripped = line.lstrip().lower()
    return any(stripped.startswith(marker) for marker in markers)


def _todays_version_lines(path: Path, text: str) -> list[str]:
    """The lines of *text* that repeat today's ``APP_VERSION`` outside a comment.

    *path* is consulted for its name alone — which comment syntax applies — so the
    vacuity check below can hold the exemption's boundary on planted lines without
    writing files.
    """
    return [
        line.strip()
        for line in TAG_SHAPED.sub("", text).splitlines()
        if re.search(rf"(?<!\d){re.escape(APP_VERSION)}(?!\d)", line)
        and not _is_full_line_comment(path, line)
    ]


def test_no_build_file_repeats_todays_version():
    """The companion sweep, and it catches what the shape-based one cannot.

    :data:`VERSION_IN_PROSE` finds a version only next to the word *MAFigate*, so a bare
    ``1.0.0`` in a build file — a worked example, a message — is invisible to it. This one
    looks for today's ``APP_VERSION`` anywhere in those files instead.

    The two are deliberately complementary rather than redundant. A sweep pinned to the
    current number goes blind the moment the number is bumped, which is why it is not the
    only one; a sweep pinned to a shape cannot see a number written on its own, which is why
    it is not enough. Between them, a *fresh* copy fails today and a *stale* copy fails
    whenever it is next read.

    Nothing else in these files is a dotted number equal to the app version — the bundled
    Python is 3.11.x, the macOS deployment target 10.15 — so this needs no allow-list of
    *paths*, and that is what it must never grow: an exempted path is a copy with
    permission. What issue #422 did exempt is a kind of **line**: a full-line comment, in
    the file's own comment syntax. The launcher fix left the tree remembering which release
    died — ``v1.0.0`` in ``launch.sh``'s and ``launcher-contract.yml``'s comments, naming
    the failure those files exist to prevent — and rewriting that to ``<version>`` would
    delete the one thing the sentence is for. It is the ``build/RELEASES.md`` reasoning one
    level down: a record of history is not an instruction to a builder, and a comment is
    the one shape of line no build and no worked example reads. Everything else still
    fails: a literal in code, in markdown prose, or beside code on the same line — held so
    by the planted sweep below.
    """
    offenders = {}
    for path in _build_and_docs_files():
        hits = _todays_version_lines(
            path, path.read_text(encoding="utf-8", errors="replace")
        )
        if hits:
            offenders[path.relative_to(REPO_ROOT).as_posix()] = hits

    assert not offenders, (
        f"a build file repeats {APP_VERSION}, which is today's APP_VERSION: {offenders}. "
        "Derive it — `make version`, build/version.py — or write `<version>`. A record of "
        "a past release may stay only as a full-line comment."
    )


#: The exemption's boundary, driven from both sides. Left: lines the sweep must still
#: fail on — assignments, values, a make variable, markdown (which has no comment to be
#: in), and the mixed case, code with a trailing comment, because code with a comment on
#: it is code. Right: the one exempt shape, a full-line comment, in every comment syntax
#: the swept tree actually contains.
PLANTED_VERSION_LINES = (
    ("still_caught", "planted.sh", 'VERSION="{v}"'),
    ("still_caught", "planted.sh", 'build_dmg "{v}"  # a worked example'),
    ("still_caught", "planted.yml", "    version: {v}"),
    ("still_caught", "Makefile", "VERSION ?= {v}"),
    ("still_caught", "planted.md", "# MAFigate {v} release notes"),
    ("still_caught", "planted.md", "run the {v} installer"),
    ("still_caught", "planted.bat", "set MAFIGATE_VERSION={v}"),
    ("exempt", "planted.sh", "# which is how v{v}'s second launch died"),
    ("exempt", "planted.yml", "  # v{v} shipped without this"),
    ("exempt", "planted.py", "# the v{v} incident"),
    ("exempt", "planted.bat", "REM the v{v} installer did this"),
    ("exempt", "planted.iss", "; the v{v} installer did this"),
    ("exempt", "planted.swift", "// the v{v} launcher did this"),
    ("exempt", "Makefile", "# the v{v} recipe was wrong"),
)


@pytest.mark.parametrize(
    "expectation,name,template",
    PLANTED_VERSION_LINES,
    ids=[f"{case[0]}-{case[1]}-{index}" for index, case in enumerate(PLANTED_VERSION_LINES)],
)
def test_the_version_sweep_exempts_comments_and_nothing_else(expectation, name, template):
    """Vacuity check for the sweep above: relaxed for comments, and only for comments.

    A guard relaxed until it passes is this repository's recurring failure mode, so the
    comment exemption is not taken on trust: every planted *code* line has to go on
    failing, and only the full-line comment shape may be passed over. Driven through the
    same helper the real sweep calls, so the two cannot drift apart.
    """
    line = template.format(v=APP_VERSION)
    hits = _todays_version_lines(Path(name), line)

    if expectation == "still_caught":
        assert hits, (
            f"the sweep no longer notices {line!r} in {name} — the comment exemption "
            "has swallowed a line a build or a reader can still be pointed at"
        )
    else:
        assert not hits, (
            f"{line!r} in {name} is a full-line comment, the one shape that cannot "
            "reach an artifact; the exemption is not being applied"
        )


def test_the_shipped_examples_carry_the_derived_version():
    """The parameter examples stamp an ``app_version``, so it is held to the constant.

    They are two documents the app ships as documentation of what it writes, and they were
    stamped ``2.0.0`` while the app was being renumbered — which would have made them the
    next stale copies, silently, on the next bump. ``tests/test_param_migration.py`` owns
    everything else about these files (their schema stamp and their agreement with the
    parameter contract); this holds the one field #260 is about.

    Their frozen ``timestamp`` is deliberately *not* held to anything: it records when the
    example was written, has no single source to derive from, and nothing reads it.
    """
    import yaml

    for name in ("example_parameters.json", "example_parameters.yaml"):
        path = STREAMLIT_APP / name
        text = path.read_text(encoding="utf-8")
        document = json.loads(text) if name.endswith(".json") else yaml.safe_load(text)

        assert "app_version" in document, (
            f"{name} no longer carries an app_version, so this guard is watching nothing — "
            "and the app stamps one into every document it writes."
        )
        assert str(document["app_version"]) == APP_VERSION, (
            f"{name} says it was written by MAFigate {document['app_version']} while the "
            f"app is {APP_VERSION}. It is shipped as an example of the app's own output."
        )


def test_the_entry_point_docstring_names_no_version():
    """Issue #71 left the version in the About dialog alone; the docstring kept two copies.

    ``MAFigate.py``'s module docstring opened with ``MAFigate v2.0.0`` and closed with
    ``Version: 2.0.0`` — not rendered anywhere, so #71's surface sweep could not see them,
    and stale for exactly that reason.
    """
    source = (STREAMLIT_APP / "MAFigate.py").read_text(encoding="utf-8")
    docstring = source.split('"""')[1]

    hits = _prose_versions(docstring) + re.findall(r"(?m)^Version:\s*\S+", docstring)
    assert not hits, (
        f"MAFigate.py's module docstring names a version: {hits}. The About dialog is "
        "where that number renders, from APP_VERSION."
    )


def test_every_release_tag_named_in_prose_is_the_derived_one():
    """A tag promised in the README is a copy, so it is held to the derivation.

    The public README tells a reader the installers are "coming with the first
    `mafigate-v1.0.0` release", which is useful copy and a version literal at once. Rather
    than delete the promise, this pins it: bump ``APP_VERSION`` and the sentence goes red
    until someone rewrites it. The private notes are exempt — a ticket body records what
    was true when it was charted, and rewriting history to match the present is how the
    map stops being evidence.

    ``build/RELEASES.md`` is exempt for that same reason and no other (#264): it is the record
    of which releases have actually been *published*, so its entries are history and a bump
    cannot make them wrong. What replaces the pinning there is
    ``test_delivery_channels_copy.TheLedgerIsReadableAtAll``, which refuses an entry naming a
    version ahead of ``APP_VERSION`` — the one direction in which a ledger line can be a
    fiction rather than a memory.
    """
    tag = f"{TAG_PREFIX}{APP_VERSION}"
    named = {}
    for relative in _tracked_files():
        if relative.startswith(("docs/wayfinder/", "streamlit_app/tests/")):
            continue
        if relative == "streamlit_app/build/RELEASES.md":
            continue
        path = REPO_ROOT / relative
        if path.suffix not in {".md", ".py", ".sh", ".bat", ".iss", ".yml", ".toml", ".txt"}:
            continue
        found = {
            hit
            for hit in re.findall(r"mafigate-v\d+\.\d+(?:\.\d+)?", path.read_text(
                encoding="utf-8", errors="replace"
            ))
            if hit != tag
        }
        if found:
            named[relative] = sorted(found)

    assert not named, (
        f"prose names a release tag that is not {tag}, the one APP_VERSION derives: "
        f"{named}"
    )
