"""Drift guard for the vendored pipeline filter code.

``streamlit_app/vendor/`` holds copies of the pipeline's filter code so the app
reaches the *same* PASS/NOPASS verdict as ``bin/filter_variants.py`` instead of
re-deriving it. A copy is only trustworthy while it is provably still a copy, so
this module is the proof.

Three mechanisms, matching the three vendoring units:

* ``pipeline_utils.py`` — whole-file sha256 against ``bin/utils.py``. Intolerant of
  even cosmetic edits, which is acceptable because the fix is one command.
* ``pipeline_filters.py``'s five filter functions — per-symbol ``ast.dump()`` equality
  against the matching ``FunctionDef`` in ``bin/filter_variants.py``. Tolerates
  comment, whitespace and formatting churn in ``bin/``; catches anything semantic.
* ``pipeline_filters.py``'s ``compute_keep`` — the output-column computation, which is
  not a function in ``bin/`` but eight statements inside ``main()``. Five are vendored
  and compared statement by statement; three are guarded from here without being
  vendored. See ``PINNED_COLUMN_STATEMENTS``.

All three skip when ``bin/`` is absent: the guard needs the pipeline checkout, and the
packaged .dmg/.exe deliberately ships neither ``bin/`` nor this test.

**Nothing in this file may be intolerant of comment-only churn in ``bin/``.** The one
exception is ``pipeline_utils.py``'s sha256, which predates this rule and covers a file
copied whole. Everywhere else the comparison is on syntax trees, deliberately: ``bin/``
is frozen for the parity effort but not embalmed, ``.github/workflows/vendor-drift.yml``
runs this whole file on every push, and a guard that turns CI red over a comment is one
that gets switched off. A check that would be red on any case in ``TOLERATED_EDITS``
does not belong here — including a tempting byte-for-byte comparison of the vendored
column statements against ``main()``. That property is real, and it is how ``_sync.py``
lifts them, but policing it continuously costs more than it is worth.

Nothing here imports the vendored modules, pandas, or streamlit — the comparison is
on source text and syntax trees only, so this file is runnable with a bare pytest.
"""

import ast
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent
BIN_DIR = REPO_ROOT / "bin"
VENDOR_DIR = STREAMLIT_APP / "vendor"

# Keep in step with FILTER_SYMBOLS in streamlit_app/vendor/_sync.py.
VENDORED_FILTER_SYMBOLS = [
    "common_filters",
    "has_element_from_list",
    "check_civic_column_exists",
    "somatic_filters",
    "germline_filters",
]

#: The sixth symbol: the output-column computation, lifted from `main()`.
VENDORED_COLUMN_SYMBOL = "compute_keep"

#: Everything vendored. Keep in step with VENDORED_SYMBOLS in _sync.py.
VENDORED_SYMBOLS = VENDORED_FILTER_SYMBOLS + [VENDORED_COLUMN_SYMBOL]

#: The statement `compute_keep` must open with, in place of the pipeline's
#: `keep = KEEP`. The pipeline aliases its module-level column list and then mutates
#: it; the app keeps the module loaded across runs, where that corrupts every call
#: after the first. `tests/test_vendor_compute_keep.py` reproduces both symptoms.
#: Keep in step with LOCAL_COPY_STATEMENT in streamlit_app/vendor/_sync.py.
LOCAL_COPY_STATEMENT = "keep = list(KEEP)"

#: The three statements of the column slice that are guarded here but NOT vendored,
#: in the order they must appear in `main()`. Pinned as source text because there is
#: no vendored copy to compare them against. Keep in step with the pinned entries of
#: COLUMN_SLICE in streamlit_app/vendor/_sync.py.
#:
#: * `keep = KEEP` is the aliasing statement `compute_keep` replaces. If the pipeline
#:   ever changes how it seeds the list, the replacement stops being equivalent and
#:   this is the only thing that would notice.
#: * The two reselections are what apply `keep` to the outgoing frames — what makes it
#:   *the pipeline's output columns* rather than a list the pipeline happens to build.
#:   Deleting either one leaves that frame unnarrowed, and the app would go on
#:   reproducing a column set the pipeline no longer emits.
PINNED_COLUMN_STATEMENTS = [
    "keep = KEEP",
    "out_nopass = out_nopass[keep]",
    "out = out[keep]",
]

#: How many statements of `main()` the column guard covers: the five `compute_keep`
#: vendors plus the three pinned above. Keep in step with the length of COLUMN_SLICE
#: in streamlit_app/vendor/_sync.py — stated here as a plain number on purpose, so
#: that a slice which grew or shrank on that side has to be acknowledged on this one.
COLUMN_SLICE_SIZE = 8

RESYNC_HINT = "Fix with: python streamlit_app/vendor/_sync.py"

needs_pipeline = pytest.mark.skipif(
    not BIN_DIR.is_dir(),
    reason="pipeline bin/ not present (packaged app or partial checkout)",
)

#: CancerVar's clinical-marker table — the third vendored unit, and the only one whose
#: source is not under `bin/`. Guarded by whole-file sha256, like `pipeline_utils.py`:
#: it is data, so it has no symbols to compare and no cosmetic edits to tolerate.
#: Keep in step with MARKERS_SRC / MARKERS_DEST in streamlit_app/vendor/_sync.py.
MARKERS_SRC = REPO_ROOT / "resources" / "CancerVar" / "cancervardb" / "cancervar.out.txt"
MARKERS_DEST = VENDOR_DIR / "cancervar_markers.txt"

#: Its own skip mark, separate from `needs_pipeline`: the marker table lives under
#: `resources/CancerVar/`, so a checkout can have `bin/` without it or the reverse.
#: `_sync.py`'s GUARDED_SOURCES makes the same distinction — see
#: `test_sync_check_agrees_with_this_guard`.
needs_cancervar_db = pytest.mark.skipif(
    not MARKERS_SRC.is_file(),
    reason="resources/CancerVar/cancervardb not present (packaged app or partial checkout)",
)


def _function_defs(path: Path) -> dict:
    """Map top-level function name -> ast node for a Python file."""
    tree = ast.parse(path.read_text())
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


@needs_pipeline
def test_vendor_directory_is_present():
    """The vendored package must exist, or the app has no filter code to run."""
    assert (VENDOR_DIR / "__init__.py").is_file(), f"{VENDOR_DIR} is not a package"
    assert (VENDOR_DIR / "pipeline_utils.py").is_file()
    assert (VENDOR_DIR / "pipeline_filters.py").is_file()


@needs_pipeline
def test_pipeline_utils_is_byte_identical_to_bin_utils():
    """`vendor/pipeline_utils.py` is `bin/utils.py`, to the byte."""
    origin = BIN_DIR / "utils.py"
    vendored = VENDOR_DIR / "pipeline_utils.py"

    origin_hash = hashlib.sha256(origin.read_bytes()).hexdigest()
    vendored_hash = hashlib.sha256(vendored.read_bytes()).hexdigest()

    assert vendored_hash == origin_hash, (
        f"vendor/pipeline_utils.py has drifted from bin/utils.py\n"
        f"  bin/utils.py             sha256 {origin_hash}\n"
        f"  vendor/pipeline_utils.py sha256 {vendored_hash}\n"
        f"{RESYNC_HINT}"
    )


@needs_cancervar_db
def test_cancervar_markers_is_byte_identical_to_the_pipeline_copy():
    """`vendor/cancervar_markers.txt` is CancerVar's marker table, to the byte.

    Not merely tidiness. A MAF's `Therap_list` / `Diag_list` / `Prog_list` cell holds
    **0-based line offsets** into this file — `CancerVar.py:193` reads it with
    `list(csv.reader(...))`, header included, and `:1091` indexes straight into that
    list. So a copy of a different vintage does not fail to resolve: every index after
    an inserted row resolves to the *neighbouring* marker, and the app names a drug the
    file never associated with the variant. There is no cosmetic edit to a data file,
    which is why this is a whole-file hash and not a parse-and-compare.

    Issue #198 measured the two in step: 112,328 indices resolved across 109 real
    somatic MAFs with 0 out of range and 0 whose `Evidence_type` disagreed with the
    criterion that cited them.
    """
    origin_hash = hashlib.sha256(MARKERS_SRC.read_bytes()).hexdigest()
    vendored_hash = hashlib.sha256(MARKERS_DEST.read_bytes()).hexdigest()

    assert MARKERS_DEST.is_file(), (
        f"{MARKERS_DEST.relative_to(REPO_ROOT)} is missing, so the app cannot name any "
        f"of CancerVar's therapeutic, diagnostic or prognostic markers.\n{RESYNC_HINT}"
    )
    assert vendored_hash == origin_hash, (
        f"vendor/cancervar_markers.txt has drifted from the pipeline's copy. MAF marker\n"
        f"indices are line offsets into this file, so the app would name the WRONG drug.\n"
        f"  {MARKERS_SRC.relative_to(REPO_ROOT)} sha256 {origin_hash}\n"
        f"  vendor/cancervar_markers.txt          sha256 {vendored_hash}\n"
        f"{RESYNC_HINT}"
    )


@needs_cancervar_db
def test_cancervar_markers_still_has_its_header_at_row_zero():
    """Row 0 must be the header, because that is what makes the indices 0-based.

    `CancerVar.py:193` does `cancervar_d = list(csv.reader(fh, delimiter="\\t"))` over
    the whole file and never skips a header, so `cancervar_d[0]` *is* the header line
    and `cancervar_d[8989]` is the 8,990th line of the file. A copy stripped of its
    header would shift every index by one while still hashing to something self-
    consistent, so this pins the property the offsets depend on rather than the bytes.

    The sha256 above would also catch a stripped header — but only while the pipeline's
    copy still has one. This says *which* property the app relies on, so that a future
    upstream file without a header fails here, naming the reason, instead of silently
    re-basing every index.
    """
    with MARKERS_DEST.open() as handle:
        first = handle.readline().rstrip("\n").split("\t")

    assert first[9:11] == ["Evidence_type", "Evidence_level"], (
        "vendor/cancervar_markers.txt no longer starts with CancerVar's header row.\n"
        "The app resolves marker indices as 0-based offsets that COUNT that header "
        "(CancerVar.py:193 reads the file with no skip), and it reads column 9 to check "
        f"a marker's Evidence_type and column 10 for its level.\nFirst row: {first}"
    )


@needs_pipeline
@pytest.mark.parametrize("symbol", VENDORED_FILTER_SYMBOLS)
def test_vendored_filter_symbol_matches_pipeline(symbol):
    """Each vendored filter function is semantically identical to the pipeline's."""
    origin_defs = _function_defs(BIN_DIR / "filter_variants.py")
    vendored_defs = _function_defs(VENDOR_DIR / "pipeline_filters.py")

    assert symbol in origin_defs, f"{symbol} no longer exists in bin/filter_variants.py"
    assert symbol in vendored_defs, f"{symbol} is missing from vendor/pipeline_filters.py"

    assert ast.dump(vendored_defs[symbol]) == ast.dump(origin_defs[symbol]), (
        f"vendor/pipeline_filters.py:{symbol} has drifted from "
        f"bin/filter_variants.py:{symbol}\n{RESYNC_HINT}"
    )


@needs_pipeline
def test_no_extra_functions_vendored_silently():
    """The vendored filter module holds exactly the declared symbols, no more.

    An unguarded extra function would be app code masquerading as pipeline code —
    exactly the failure mode `main_utils.KEEP` already demonstrated by drifting to
    39 entries against the pipeline's 45.
    """
    vendored = set(_function_defs(VENDOR_DIR / "pipeline_filters.py"))
    assert vendored == set(VENDORED_SYMBOLS), (
        "vendor/pipeline_filters.py holds functions that no guard covers: "
        f"{sorted(vendored - set(VENDORED_SYMBOLS))}\n"
        "Add them to VENDORED_SYMBOLS (and to _sync.py) or take them out."
    )


# ---------------------------------------------------------------------------
# The column slice
# ---------------------------------------------------------------------------
#
# `_sync.py` locates these eight statements by source-line ANCHOR. This module finds
# them a different way on purpose — by searching `main()` for a statement whose AST
# matches, taking the statements to look for from the vendored file itself and from
# `PINNED_COLUMN_STATEMENTS`. Two mechanisms that must reach the same verdict is the
# same arrangement the rest of this file already runs on, and
# `test_sync_check_agrees_with_this_guard` is what pins them together.


def _statements(body: list):
    """Every statement under `body`, outer before inner, in source order.

    A twin of `_sync.py`'s `_walk_statements`, not an import of it — for the same
    reason `_function_defs` above is a twin of the one in `_sync.py`. If the two guards
    shared their statement enumerator, a walker that missed a nesting form would hide
    the same statements from both of them, and the second opinion would be worth
    nothing. What differs between the two guards is how they pick the eight statements
    out of this stream: `_sync.py` matches source-line anchors, this module matches
    ASTs.
    """
    for statement in body:
        yield statement
        for field in ("body", "orelse", "finalbody"):
            yield from _statements(getattr(statement, field, None) or [])
        for handler in getattr(statement, "handlers", None) or []:
            yield from _statements(handler.body)


def _main_statements() -> list:
    """Every statement inside `bin/filter_variants.py:main()`, at any depth."""
    tree = ast.parse((BIN_DIR / "filter_variants.py").read_text())
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "main":
            return list(_statements(node.body))
    raise AssertionError("bin/filter_variants.py no longer defines main()")


def _compute_keep_body() -> list:
    """`compute_keep`'s statements, docstring dropped."""
    definitions = _function_defs(VENDOR_DIR / "pipeline_filters.py")
    assert VENDORED_COLUMN_SYMBOL in definitions, (
        f"vendor/pipeline_filters.py does not define {VENDORED_COLUMN_SYMBOL}, so the "
        f"app has no output-column computation.\n{RESYNC_HINT}"
    )
    body = definitions[VENDORED_COLUMN_SYMBOL].body
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
        and isinstance(body[0].value.value, str)
    ):
        body = body[1:]
    return body


def _vendored_column_statements() -> list:
    """The statements `compute_keep` carries verbatim: everything but the fix and the
    return."""
    return _compute_keep_body()[1:-1]


def _locate_in_main(statement, what: str):
    """The one statement of `main()` whose AST equals `statement`; assert if not one.

    `ast.dump` ignores line numbers, comments and formatting, so a comment-only edit
    in `bin/` still finds its match. Anything semantic does not.
    """
    wanted = ast.dump(statement)
    matches = [node for node in _main_statements() if ast.dump(node) == wanted]
    assert len(matches) == 1, (
        f"{what} matches {len(matches)} statements in bin/filter_variants.py:main() "
        f"(expected exactly 1). The pipeline's column computation has changed, so the "
        f"app's copy of it is no longer the pipeline's.\n{RESYNC_HINT}\n"
        f"Statement:\n{ast.unparse(statement) if hasattr(ast, 'unparse') else wanted}"
    )
    return matches[0]


@needs_pipeline
def test_compute_keep_takes_a_local_copy_before_the_verbatim_statements():
    """`compute_keep` must open with the fix, and the fix must not be a vendored edit.

    This is the one place the copy is deliberately unfaithful, so it is asserted
    rather than compared: the first statement is the app's, every statement after it
    is the pipeline's, and `bin/` keeps the aliasing form it has always had (pinned
    below as `PINNED_COLUMN_STATEMENTS[0]`).
    """
    body = _compute_keep_body()
    assert body, f"{VENDORED_COLUMN_SYMBOL} has no body beyond its docstring"

    assert ast.dump(body[0]) == ast.dump(ast.parse(LOCAL_COPY_STATEMENT).body[0]), (
        f"{VENDORED_COLUMN_SYMBOL} must start with `{LOCAL_COPY_STATEMENT}`. The "
        f"pipeline's `keep = KEEP` aliases the module-level list and then mutates it, "
        f"which corrupts every call after the first in a long-running app.\n"
        f"{RESYNC_HINT}"
    )
    assert isinstance(body[-1], ast.Return), (
        f"{VENDORED_COLUMN_SYMBOL} must end by returning the column list"
    )


@needs_pipeline
def test_column_guard_covers_eight_statements_of_main():
    """The guard's coverage is a fixed size, so a shrunken slice cannot pass quietly.

    Without this, deleting a statement from `compute_keep` would leave every remaining
    comparison green and the guard would report a column computation it no longer
    covers.
    """
    covered = len(_vendored_column_statements()) + len(PINNED_COLUMN_STATEMENTS)
    assert covered == COLUMN_SLICE_SIZE, (
        f"the column guard covers {covered} statements of main(), expected "
        f"{COLUMN_SLICE_SIZE}. If the pipeline's column computation genuinely changed "
        f"size, update COLUMN_SLICE_SIZE and _sync.py's COLUMN_SLICE_ANCHORS together."
    )


@needs_pipeline
def test_vendored_column_statements_are_present_in_main():
    """Every statement `compute_keep` vendors is still in `main()`, exactly once.

    This is the drift check proper: a column renamed or mistyped in `bin/` leaves the
    vendored statement matching nothing, and a statement duplicated in `main()` leaves
    it matching twice. Either way the app's column set has stopped being derivable
    from the pipeline's.
    """
    vendored = _vendored_column_statements()
    assert vendored, "compute_keep vendors no statements at all"

    located = [
        _locate_in_main(statement, f"vendored column statement {index}")
        for index, statement in enumerate(vendored)
    ]

    linenos = [node.lineno for node in located]
    assert linenos == sorted(linenos), (
        f"the vendored column statements appear in main() in a different order than "
        f"compute_keep runs them (lines {linenos}), so the app would build the column "
        f"list differently from the pipeline.\n{RESYNC_HINT}"
    )


@needs_pipeline
@pytest.mark.parametrize("source", PINNED_COLUMN_STATEMENTS)
def test_pinned_column_statements_are_present_in_main(source):
    """The three guarded-but-not-vendored statements are still in `main()`, unchanged.

    Nothing else covers these. `keep = KEEP` has no vendored counterpart because
    `compute_keep` replaces it; the two reselections have none because the app does its
    own column selection. Deleting a reselection is what would leave the pipeline
    emitting a frame it no longer narrows, with the app still narrowing it.
    """
    _locate_in_main(ast.parse(source).body[0], f"pinned column statement `{source}`")


@needs_pipeline
def test_pinned_column_statements_run_in_the_declared_order():
    """`keep = KEEP` seeds the list before the two reselections consume it."""
    linenos = [
        _locate_in_main(ast.parse(source).body[0], source).lineno
        for source in PINNED_COLUMN_STATEMENTS
    ]
    assert linenos == sorted(linenos), (
        f"the pinned column statements are out of order in main() (lines {linenos})"
    )


DMG_SCRIPT = STREAMLIT_APP / "build" / "mac" / "build_dmg.sh"
ISS_SCRIPT = STREAMLIT_APP / "build" / "windows" / "installer.iss"

needs_build_scripts = pytest.mark.skipif(
    not (DMG_SCRIPT.is_file() and ISS_SCRIPT.is_file()),
    reason="installer scripts not present (build/ is stripped from the packaged app)",
)


def _dmg_copied_sources(dmg: str) -> list[str]:
    """Every *literal* app-source name `build_dmg.sh` copies into the bundle.

    Two shapes to cover: the fixed ``for item in ...`` list, and the standalone
    ``cp -R "${PROJECT_DIR}/<name>"`` calls that sit outside it.

    The loop's own body is one of those calls — ``cp -R "${PROJECT_DIR}/${item}"`` — so
    captures that are shell expansions rather than literal names are dropped. Keeping
    them would leave this function returning a never-empty list, which would make any
    caller's "did we actually parse the script?" check unfalsifiable.
    """
    names = []
    for line in dmg.splitlines():
        stripped = line.strip()
        if stripped.startswith("for item in ") and "; do" in stripped:
            listing = stripped[len("for item in ") : stripped.index("; do")]
            names.extend(listing.split())
        match = re.search(r'cp -R\s+"\$\{PROJECT_DIR\}/([^"]+)"', stripped)
        if match and "$" not in match.group(1):
            names.append(match.group(1))
    return names


def _iss_source_lines(iss: str) -> list[str]:
    """Every ``Source:`` directive in the Inno Setup script, stripped."""
    return [line.strip() for line in iss.splitlines() if line.strip().startswith("Source:")]


#: Names the macOS copy list is known to carry. Asserted by both installer tests as a
#: parse anchor: if `build_dmg.sh` is restructured so `_dmg_copied_sources` stops
#: recognising its shape, these fail loudly instead of every check passing vacuously
#: over an empty list.
DMG_ANCHOR_SOURCES = {"MAFigate.py", "vendor"}


@needs_build_scripts
def test_both_installers_ship_the_vendor_package():
    """Both installers must carry `vendor/`, or the packaged app has no filter code.

    Neither installer globs the source tree — the macOS script copies a fixed list of
    names and the Inno Setup script names one `Source:` line per directory — so a new
    package is invisible to them until it is added by hand. This test is what makes
    that a build failure rather than a silent shipping of an app that cannot filter.
    """
    copied = _dmg_copied_sources(DMG_SCRIPT.read_text())
    assert "MAFigate.py" in copied, (
        "could not parse the app-source copy list out of build_dmg.sh — the script's "
        f"shape changed and this test is no longer reading it. Parsed: {copied}"
    )
    assert "vendor" in copied, (
        "build/mac/build_dmg.sh copies a fixed list of app directories and 'vendor' is "
        f"not among them, so the .dmg would ship without any filter code:\n{copied}"
    )

    source_lines = _iss_source_lines(ISS_SCRIPT.read_text())
    assert source_lines, "could not find any Source: directive in installer.iss"
    assert any("vendor" in line for line in source_lines), (
        "build/windows/installer.iss has no Source: line for vendor\\*, so the .exe "
        "would install without any filter code"
    )


def _vendor_tree_without_bin(tmp_path: Path) -> Path:
    """Copy `vendor/` into a throwaway tree that has no `bin/` above it.

    `_sync.py` derives `REPO_ROOT` from its own location (`parents[1]` of the vendor
    directory), so planting the package at `<tmp>/streamlit_app/vendor/` gives it a
    repo root of `<tmp>` — which has no `bin/`. That is a faithful stand-in for the
    packaged app and for a checkout of the app subtree alone.
    """
    dest = tmp_path / "streamlit_app" / "vendor"
    dest.parent.mkdir(parents=True)
    shutil.copytree(VENDOR_DIR, dest)
    assert not (tmp_path / "bin").exists(), "stand-in tree must not have bin/"
    return dest


def test_sync_check_skips_cleanly_without_bin(tmp_path):
    """The build-gate form must skip, not fail, where `bin/` is absent.

    `--check` is the form the installer build gates, the pre-commit hook and CI all
    run, and it is reachable from trees that carry no pipeline: the app subtree on its
    own, or a bundle being rebuilt from one. Exiting non-zero there would block the
    build over a comparison that is merely impossible to make, not failed.

    Deliberately *not* `@needs_pipeline` — a skip-if-absent mark would remove the only
    test of the absent-`bin/` path exactly where that path is live.
    """
    vendor = _vendor_tree_without_bin(tmp_path)

    result = subprocess.run(
        [sys.executable, str(vendor / "_sync.py"), "--check"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        "`_sync.py --check` failed in a tree with no bin/. The build gates and the "
        "pre-commit hook run this command, so a non-zero exit here breaks builds from "
        f"the app subtree.\n{result.stdout}{result.stderr}"
    )
    assert "SKIP" in result.stdout, (
        "the skip must announce itself, or a build log reads it as a check that ran "
        f"and passed:\n{result.stdout}"
    )
    assert "in sync" not in result.stdout, (
        "a skipped check must not claim the copies were verified in sync:\n"
        f"{result.stdout}"
    )


def test_sync_repair_still_fails_without_bin(tmp_path):
    """The *repair* form must keep failing where `bin/` is absent.

    Asymmetric on purpose. `--check` can honestly report "cannot compare"; a rewrite
    has nothing to copy from, so silently succeeding would leave the caller believing
    the copies were refreshed against a pipeline that was never there.
    """
    vendor = _vendor_tree_without_bin(tmp_path)

    result = subprocess.run(
        [sys.executable, str(vendor / "_sync.py")],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0, (
        "`_sync.py` (rewrite mode) exited 0 with no bin/ to copy from, so it reported "
        f"success for work it could not have done.\n{result.stdout}{result.stderr}"
    )


@needs_build_scripts
def test_windows_installer_does_not_exclude_the_marker_table():
    """The Inno Setup `vendor\\*` line excludes `_sync.py`; it must not exclude the table.

    `test_both_installers_ship_the_vendor_package` proves `vendor/` ships, but that line
    carries an `Excludes:` clause, and a 1.2 MB data file is exactly the sort of thing
    someone trims to shrink an installer. Losing it would not break the app — the marker
    disclosure would simply never appear in the .exe, which is the failure mode issue
    #198 rejected when it chose vendoring over reading `resources/` in place.
    """
    vendor_lines = [
        line for line in _iss_source_lines(ISS_SCRIPT.read_text()) if "vendor" in line
    ]
    assert vendor_lines, "installer.iss has no Source: line for vendor\\*"
    for line in vendor_lines:
        assert MARKERS_DEST.name not in line, (
            f"installer.iss excludes {MARKERS_DEST.name} from the .exe, so the packaged "
            f"app could not name any CancerVar marker:\n  {line}"
        )


@needs_build_scripts
def test_neither_installer_ships_the_test_suite_or_fixtures():
    """Neither installer may carry `tests/` — the suite or the parity fixtures.

    Three separate reasons, any one of which is sufficient:

    * The fixtures include `germline_reference.maf` and `somatic_reference.maf`, cut
      from the real GERSOM reference run. Real variant calls have no business in a
      binary handed to end users.
    * `tests/fixtures/parity/` is tens of thousands of lines the app never reads at
      runtime, in an installer whose whole point is to be self-contained and small.
    * The drift guard needs `bin/`, which no bundle carries, so shipping it would put
      a permanently-skipping test suite in front of users.

    Neither installer globs the source tree, so today this holds by construction — the
    macOS script copies a fixed list of names and Inno Setup names one `Source:` line
    per directory. That is exactly why it needs a test: the way it would break is
    someone adding `tests` to a list, and nothing else would notice.
    """
    dmg = DMG_SCRIPT.read_text()
    copied = _dmg_copied_sources(dmg)
    assert DMG_ANCHOR_SOURCES <= set(copied), (
        "could not parse the app-source copy list out of build_dmg.sh — the script's "
        "shape changed, so the offender check below would pass over an empty list "
        f"rather than because nothing offends. Parsed: {copied}"
    )

    offenders = [n for n in copied if n == "tests" or n.startswith("tests/") or "fixtures" in n]
    assert not offenders, (
        f"build/mac/build_dmg.sh copies {offenders} into the .dmg, which puts the test "
        "suite and the parity fixtures — including real reference variant calls — into "
        "a bundle handed to end users."
    )

    # Belt and braces in the script itself, and worth pinning: it is what protects the
    # bundle if `tests/` ever arrives inside one of the directories that IS copied.
    assert 'rm -rf "${DEST}/tests"' in dmg, (
        "build_dmg.sh no longer strips ${DEST}/tests after copying. Restore it: the "
        "copy list is an allowlist, but this is what catches a tests/ directory that "
        "arrives nested inside a copied package."
    )

    source_lines = _iss_source_lines(ISS_SCRIPT.read_text())
    assert source_lines, "could not find any Source: directive in installer.iss"
    iss_offenders = [ln for ln in source_lines if re.search(r"[\\/]tests[\\/*]|fixtures", ln)]
    assert not iss_offenders, (
        "build/windows/installer.iss has Source: line(s) that would install the test "
        f"suite or the parity fixtures:\n{iss_offenders}"
    )


# ---------------------------------------------------------------------------
# Proving the column guard actually fires
# ---------------------------------------------------------------------------
#
# A guard that has never been seen to fail is indistinguishable from a guard that
# cannot fail. These tests plant a throwaway repo, edit `bin/` inside it, and check
# that both implementations of the rule — `_sync.py --check` and the checks above —
# reach the verdict the edit deserves.
#
# The green cases matter as much as the red ones: a guard that fires on a reformatted
# comment gets switched off by the third person it inconveniences.


def _repo_with_bin(tmp_path: Path) -> Path:
    """A throwaway repo carrying both sides of the copy, safe to mutate.

    `_sync.py` derives its repo root from its own location, so `bin/` beside
    `streamlit_app/vendor/` is all it needs to run against this tree instead of the
    real one.
    """
    shutil.copytree(BIN_DIR, tmp_path / "bin")
    vendor = tmp_path / "streamlit_app" / "vendor"
    vendor.parent.mkdir(parents=True)
    shutil.copytree(VENDOR_DIR, vendor)
    return tmp_path


def _edit_pipeline(repo: Path, old: str, new: str) -> None:
    """Replace `old` with `new` in the throwaway `bin/filter_variants.py`.

    Asserts the target text was there exactly once. A mutation test whose mutation
    silently failed to apply is a test that proves nothing.
    """
    path = repo / "bin" / "filter_variants.py"
    source = path.read_text()
    assert source.count(old) == 1, (
        f"expected to find {old!r} exactly once in bin/filter_variants.py, found "
        f"{source.count(old)}. This mutation test no longer edits what it means to."
    )
    path.write_text(source.replace(old, new))


def _sync_check(repo: Path) -> subprocess.CompletedProcess:
    """Run the stdlib-only build gate against the throwaway repo."""
    return subprocess.run(
        [sys.executable, str(repo / "streamlit_app" / "vendor" / "_sync.py"), "--check"],
        capture_output=True,
        text=True,
    )


#: The checks above, as callables, so the mutation tests can run them against a
#: throwaway tree. Every one of them reads `BIN_DIR` / `VENDOR_DIR` when called rather
#: than at import, which is what makes redirecting them possible.
#:
#: Every column check in this module belongs here. If you add one that cannot go in —
#: because it would go red on one of `TOLERATED_EDITS` — that is the signal to not add
#: it at all. See the note on comment tolerance in this module's docstring.
COLUMN_GUARD_CHECKS = (
    test_compute_keep_takes_a_local_copy_before_the_verbatim_statements,
    test_column_guard_covers_eight_statements_of_main,
    test_vendored_column_statements_are_present_in_main,
    test_pinned_column_statements_run_in_the_declared_order,
    test_no_extra_functions_vendored_silently,
)


def _pytest_guard_failures(monkeypatch, repo: Path) -> list:
    """Run this module's own guard against `repo`; return what it complains about."""
    monkeypatch.setattr(sys.modules[__name__], "BIN_DIR", repo / "bin")
    monkeypatch.setattr(sys.modules[__name__], "VENDOR_DIR", repo / "streamlit_app" / "vendor")

    failures = []
    for check in COLUMN_GUARD_CHECKS:
        try:
            check()
        except AssertionError as error:
            failures.append(f"{check.__name__}: {error}")
    for source in PINNED_COLUMN_STATEMENTS:
        try:
            test_pinned_column_statements_are_present_in_main(source)
        except AssertionError as error:
            failures.append(f"pinned {source!r}: {error}")
    return failures


@needs_pipeline
def test_unmutated_throwaway_repo_is_green(tmp_path, monkeypatch):
    """The mutation harness must start from green, or every red below is meaningless."""
    repo = _repo_with_bin(tmp_path)

    result = _sync_check(repo)
    assert result.returncode == 0, (
        f"the copied-but-unedited tree already reports drift, so the mutation tests "
        f"below prove nothing.\n{result.stdout}{result.stderr}"
    )
    assert not _pytest_guard_failures(monkeypatch, repo)


#: Edits to `bin/` that must turn the guard red, and why each one matters.
BREAKING_EDITS = {
    "renamed_column": (
        'keep.append("RENOVO_Class")',
        'keep.append("RENOVO_Classification")',
    ),
    "mistyped_column": (
        'keep.append("gnomAD_genome_AF")',
        'keep.append("gnomAD_genome_A")',
    ),
    "reseeded_column_list": (
        "    keep = KEEP\n",
        "    keep = KEEP.copy()\n",
    ),
    "altered_nopass_reselection": (
        "        out_nopass = out_nopass[keep]\n",
        "        out_nopass = out_nopass.reindex(columns=keep)\n",
    ),
    "deleted_pass_reselection": (
        "        out = out[keep]\n",
        "",
    ),
}


@needs_pipeline
@pytest.mark.parametrize("edit", sorted(BREAKING_EDITS))
def test_column_guard_fires_on(edit, tmp_path, monkeypatch):
    """Each of these edits to `bin/` must be reported by both forms of the guard.

    `renamed_column` and `mistyped_column` are the two ways a column name goes wrong:
    the pipeline deliberately renames one, or someone fat-fingers one. Both leave the
    app emitting a column the pipeline no longer does, and neither changes anything
    the app would crash on.

    The last three are the statements the guard covers *without* vendoring, which is
    the coverage that would otherwise be missing. `reseeded_column_list` is the
    pipeline changing how it seeds the list, which would make `compute_keep`'s local
    copy no longer an equivalent substitution. The other two are the reselections that
    apply the list: alter or delete one and the pipeline stops narrowing that frame,
    while the app goes on narrowing it to a column set the pipeline no longer emits.
    """
    old, new = BREAKING_EDITS[edit]
    repo = _repo_with_bin(tmp_path)
    _edit_pipeline(repo, old, new)

    result = _sync_check(repo)
    assert result.returncode != 0, (
        f"`_sync.py --check` stayed green after {edit} in bin/filter_variants.py, so "
        f"the build gate would ship an app whose column set no longer matches the "
        f"pipeline's.\n{result.stdout}{result.stderr}"
    )

    failures = _pytest_guard_failures(monkeypatch, repo)
    assert failures, (
        f"the pytest guard stayed green after {edit} in bin/filter_variants.py, while "
        f"`_sync.py --check` went red. The two implementations of the rule disagree."
    )


#: Edits to `bin/` that must NOT turn the guard red.
TOLERATED_EDITS = {
    "comment_inside_a_vendored_statement": (
        '        keep.append("RENOVO_Class")\n',
        '        # RENOVO columns are germline-only\n        keep.append("RENOVO_Class")\n',
    ),
    "comment_above_a_pinned_statement": (
        "        out = out[keep]\n",
        "        # narrow to the reported columns\n        out = out[keep]\n",
    ),
    "reformatted_vendored_statement": (
        '    if "gnomAD_exome_AF_raw" in out.columns.values:\n'
        '        keep.append("gnomAD_exome_AF_raw")\n',
        '    if "gnomAD_exome_AF_raw" in out.columns.values:\n'
        "        keep.append(\n"
        '            "gnomAD_exome_AF_raw"\n'
        "        )\n",
    ),
    "edit_outside_the_slice": (
        'f"filtered.{args.output}.nopass.tsv"',
        'f"filtered.{args.output}.nopass.txt"',
    ),
}


@needs_pipeline
@pytest.mark.parametrize("edit", sorted(TOLERATED_EDITS))
def test_column_guard_tolerates(edit, tmp_path, monkeypatch):
    """Churn that changes no behaviour must leave the guard green.

    Comments and re-wrapping are invisible to `ast.dump`, which is the whole reason
    this half of the guard is an AST comparison rather than a hash. `bin/` is frozen
    for the parity effort but not embalmed, and a guard that fires on a reflowed line
    is one that gets routed around.

    `edit_outside_the_slice` is the other half of the claim: the guard covers eight
    statements of `main()`, not `main()` itself, so an unrelated change to the same
    function is none of its business.
    """
    old, new = TOLERATED_EDITS[edit]
    repo = _repo_with_bin(tmp_path)
    _edit_pipeline(repo, old, new)

    result = _sync_check(repo)
    assert result.returncode == 0, (
        f"`_sync.py --check` went red after {edit} in bin/filter_variants.py, which "
        f"changes no behaviour. A guard that cries wolf over formatting gets switched "
        f"off.\n{result.stdout}{result.stderr}"
    )
    assert not _pytest_guard_failures(monkeypatch, repo), (
        f"the pytest guard went red after {edit} while `_sync.py --check` stayed "
        f"green. The two implementations of the rule disagree."
    )


@needs_pipeline
def test_guard_fires_on_an_undeclared_symbol_in_the_vendored_package(tmp_path, monkeypatch):
    """A function added to `vendor/` without being declared must be reported.

    This is how app code gets in among pipeline code: something plausible is written
    into the vendored module, no declaration is added, and from then on it is the only
    thing in there that nothing compares against `bin/`.
    """
    repo = _repo_with_bin(tmp_path)
    vendored = repo / "streamlit_app" / "vendor" / "pipeline_filters.py"
    vendored.write_text(
        vendored.read_text() + '\n\ndef tidy_column_names(keep):\n    return sorted(keep)\n'
    )

    result = _sync_check(repo)
    assert result.returncode != 0, (
        f"`_sync.py --check` stayed green with an undeclared function in "
        f"vendor/pipeline_filters.py.\n{result.stdout}{result.stderr}"
    )
    assert "UNGUARDED" in result.stdout, (
        f"the finding must name the problem as an unguarded symbol:\n{result.stdout}"
    )

    failures = _pytest_guard_failures(monkeypatch, repo)
    assert any("no_extra_functions" in failure for failure in failures), (
        f"the pytest guard did not report the undeclared symbol. Failures: {failures}"
    )


@needs_pipeline
def test_guard_fires_when_a_vendored_column_statement_is_deleted(tmp_path, monkeypatch):
    """Dropping a statement from `compute_keep` must not shrink the guard quietly.

    Every remaining statement would still match `main()` perfectly. Only the count
    check notices, which is why it exists.
    """
    repo = _repo_with_bin(tmp_path)
    vendored = repo / "streamlit_app" / "vendor" / "pipeline_filters.py"
    vendored.write_text(
        vendored.read_text().replace(
            '    if "gnomAD_genome_AF_raw" in out.columns.values:\n'
            '        keep.append("gnomAD_genome_AF_raw")\n'
            '    elif "gnomAD_genome_AF" in out.columns.values:\n'
            '        keep.append("gnomAD_genome_AF")\n',
            "",
        )
    )

    result = _sync_check(repo)
    assert result.returncode != 0, (
        f"`_sync.py --check` stayed green after a vendored column statement was "
        f"deleted.\n{result.stdout}{result.stderr}"
    )
    assert _pytest_guard_failures(monkeypatch, repo), (
        "the pytest guard stayed green after a vendored column statement was deleted"
    )


@needs_pipeline
def test_sync_check_agrees_with_this_guard():
    """`_sync.py --check` must reach the same verdict as the tests above.

    The build gate and the pre-commit hook call `_sync.py --check` rather than pytest,
    because it is stdlib-only and needs neither pytest nor pandas installed. That means
    the rule is implemented twice, so this test pins the two together: if they ever
    diverge, one of them is lying about whether the app still matches the pipeline.
    """
    result = subprocess.run(
        [sys.executable, str(VENDOR_DIR / "_sync.py"), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        "`_sync.py --check` reports drift. If the tests above are green, the two "
        "implementations of the drift rule disagree and one of them is lying; if they "
        f"are red too, fix the drift they name.\n{result.stdout}{result.stderr}"
    )
