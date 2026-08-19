#!/usr/bin/env python3
"""Regenerate the vendored copies of the pipeline's filter code.

This is the *fix* command for a red drift guard (see
``streamlit_app/tests/test_vendor_drift.py``). It never edits ``bin/`` — the
pipeline is frozen — it only refreshes the copies under this directory.

    python streamlit_app/vendor/_sync.py          # rewrite the vendored files
    python streamlit_app/vendor/_sync.py --check  # report drift, write nothing

Three units, matching the three guard mechanisms:

* ``pipeline_utils.py`` is ``bin/utils.py`` byte-for-byte. No header is added,
  because the guard is a whole-file sha256.
* ``pipeline_filters.py`` is a generated header plus the verbatim source of the
  six vendored symbols lifted out of ``bin/filter_variants.py``. That file cannot
  be copied whole: its line 10 imports ``database.database_utils``, which pulls
  in ``psycopg``, a Postgres driver the packaged .dmg/.exe will never carry.
* ``cancervar_markers.txt`` is
  ``resources/CancerVar/cancervardb/cancervar.out.txt`` byte-for-byte — data, not
  code, and guarded by a whole-file sha256 for the same reason ``pipeline_utils.py``
  is. See ``MARKERS_SRC``.

Five of those six symbols are plain functions. The sixth, ``compute_keep``, is not a
function in ``bin/`` at all — it is a slice of statements out of ``main()``, wrapped
in a ``def`` here so the app can call it. See ``COLUMN_SLICE_ANCHORS``.

The units come from **two independent source trees** (``GUARDED_SOURCES``), so absence
is per source: a checkout with ``bin/`` but no ``resources/CancerVar`` compares the
first two units and skips the third, rather than skipping everything.
"""

import argparse
import ast
import hashlib
import sys
from pathlib import Path

VENDOR_DIR = Path(__file__).resolve().parent
REPO_ROOT = VENDOR_DIR.parents[1]
BIN_DIR = REPO_ROOT / "bin"

# Vendored whole-file, guarded by sha256.
UTILS_SRC = BIN_DIR / "utils.py"
UTILS_DEST = VENDOR_DIR / "pipeline_utils.py"

# Vendored per-symbol, guarded by ast.dump() equality.
FILTERS_SRC = BIN_DIR / "filter_variants.py"
FILTERS_DEST = VENDOR_DIR / "pipeline_filters.py"

# CancerVar's clinical-marker table, vendored whole-file and guarded by sha256 — the
# `pipeline_utils.py` treatment, because this is one opaque blob with no symbols to
# compare. It is *data*, and it does not come from `bin/`: it is the file CancerVar's
# own `config.ini` names as `cancervar_markers`, which `CancerVar.py:193` reads with
# `list(csv.reader(...))` — header included, so row 0 is the header and an index in a
# MAF's `Therap_list` / `Diag_list` / `Prog_list` cell is a plain 0-based offset into
# that list.
#
# Vendored rather than read in place for the reason README.md gives for the filter
# code: the packaged .dmg/.exe carries no pipeline, so `resources/CancerVar/` is not
# there to read. Without a copy here the marker detail would exist only in a developer
# checkout — see issue #198.
#
# The sha256 is what makes the indices safe to resolve. They are offsets into *this
# file's* line order, so a copy of a different vintage does not fail to resolve, it
# resolves to the wrong drug — which is worse than naming nothing.
CANCERVAR_DB_DIR = REPO_ROOT / "resources" / "CancerVar" / "cancervardb"
MARKERS_SRC = CANCERVAR_DB_DIR / "cancervar.out.txt"
MARKERS_DEST = VENDOR_DIR / "cancervar_markers.txt"

# The five function symbols the app needs, in source order. Keep in step with
# VENDORED_FILTER_SYMBOLS in tests/test_vendor_drift.py.
FILTER_SYMBOLS = [
    "common_filters",
    "has_element_from_list",
    "check_civic_column_exists",
    "somatic_filters",
    "germline_filters",
]

# The sixth symbol, lifted from main() rather than copied from a def. See below.
COLUMN_SYMBOL = "compute_keep"

# Everything vendored, and therefore everything the "no unguarded symbol" check
# tolerates. Keep in step with VENDORED_SYMBOLS in tests/test_vendor_drift.py.
VENDORED_SYMBOLS = FILTER_SYMBOLS + [COLUMN_SYMBOL]

# --------------------------------------------------------------------------------
# The column slice
# --------------------------------------------------------------------------------
#
# The pipeline's output column list is not a function. It is eight statements strewn
# through `main()`: one that seeds `keep`, five that reshape it, and two that apply it
# to a frame on the way out.
#
# The eight are located by ANCHOR rather than by line number, so that unrelated edits
# to `main()` cannot silently shift the slice onto different statements. Each anchor is
# a tuple of line prefixes that must match the statement's own source lines in order;
# one line is enough for seven of them, and the germline branch needs two because
# `main()` tests `args.sample_type == "germline"` twice.
#
# A missing or duplicated anchor is a hard error, never a silent empty slice: the shape
# of `main()` changed, and a human has to look.
#
# Each entry carries its own disposition, so the vendored set, the pinned set and the
# size of the slice all derive from this one table rather than from integers that would
# have to be updated in step with it. `PINNED` means guarded but not vendored, with the
# source text that statement must keep; `VENDORED` means `compute_keep` carries it
# verbatim and the guard compares it against `bin/`.
VENDORED = None

COLUMN_SLICE = (
    # The aliasing statement `compute_keep` replaces — see LOCAL_COPY_STATEMENT. It has
    # no vendored counterpart, but it is watched anyway: if the pipeline ever changes
    # how it seeds the list, the replacement stops being equivalent and nothing else
    # would notice.
    (("keep = KEEP",), "keep = KEEP"),
    # The five that build the list. `compute_keep` is these, verbatim.
    (('if args.sample_type == "germline":', 'keep.remove("Tumor_Sample_Barcode")'), VENDORED),
    (("civic_columns_exist = any(",), VENDORED),
    (("if args.skip_civic or not civic_columns_exist:",), VENDORED),
    (('if "gnomAD_exome_AF_raw" in out.columns.values:',), VENDORED),
    (('if "gnomAD_genome_AF_raw" in out.columns.values:',), VENDORED),
    # The two reselections that apply `keep` to the outgoing frames. They are what make
    # the computed list *the pipeline's output columns* rather than a list the pipeline
    # happens to build. Delete either one and the pipeline stops narrowing that frame —
    # the app would go on reproducing a column set the pipeline no longer emits, with
    # nothing else in the guard covering the loss. Not vendored: the app does its own
    # column selection.
    (("out_nopass = out_nopass[keep]",), "out_nopass = out_nopass[keep]"),
    (("out = out[keep]",), "out = out[keep]"),
)

#: Anchors only, in slice order.
COLUMN_SLICE_ANCHORS = tuple(anchor for anchor, _ in COLUMN_SLICE)

#: Positions `compute_keep` carries verbatim. Not assumed contiguous.
VENDORED_SLICE_INDICES = tuple(
    index
    for index, (_, disposition) in enumerate(COLUMN_SLICE)
    if disposition is VENDORED
)

#: Positions guarded but NOT vendored, mapped to the source text they must keep. Pinned
#: as text because there is no vendored counterpart to compare them against. Keep in
#: step with PINNED_COLUMN_STATEMENTS in tests/test_vendor_drift.py.
UNVENDORED_SLICE_STATEMENTS = {
    index: disposition
    for index, (_, disposition) in enumerate(COLUMN_SLICE)
    if disposition is not VENDORED
}

#: The statement `compute_keep` substitutes for slice position 0, and the whole point of
#: this symbol. `main()` writes `keep = KEEP`, aliasing the module-level list and then
#: mutating it in place; harmless in a process that exits, wrong in an app that keeps
#: the module loaded. The copy is taken here rather than by editing the vendored
#: statements, so every statement the guard compares is still an unedited copy. Keep in
#: step with LOCAL_COPY_STATEMENT in tests/test_vendor_drift.py.
LOCAL_COPY_STATEMENT = "keep = list(KEEP)"

# Not part of the vendored surface — free to differ from bin/, and deliberately
# does not import psycopg.
FILTERS_HEADER = '''"""Vendored verbatim from ``bin/filter_variants.py`` — DO NOT EDIT BY HAND.

Regenerate with ``python streamlit_app/vendor/_sync.py``; the drift guard in
``streamlit_app/tests/test_vendor_drift.py`` compares this file to ``bin/`` by
``ast.dump()`` equality, so any hand-edit here that changes behaviour turns the
suite red.

The first five functions are copies of the pipeline's own, compared symbol for
symbol. ``compute_keep`` is not: the pipeline has no such function, only eight
statements inside ``main()``. Five of them are lifted into it verbatim and
compared statement by statement; its **first** statement is deliberately not the
pipeline's, and its docstring says why.

Only the import header is local: ``bin/filter_variants.py`` cannot be imported by
the app, because it pulls in ``psycopg`` via ``database.database_utils``.
"""

import os

import pandas as pd

from .pipeline_utils import CLINVAR_PATHO, KEEP, has_clinvar_term
'''

#: Wrapper for the column slice. Everything outside the `{statements}` hole is local:
#: `main()` has no such function, and the first statement is deliberately *not* the
#: pipeline's — see `LOCAL_COPY_STATEMENT`.
COLUMN_TEMPLATE = '''def compute_keep(args, out):
    """The pipeline's output column list — `keep` in `bin/filter_variants.py:main()`.

    The statements after the first are lifted verbatim out of `main()`, so they read
    the pipeline's own names. Supply both as call-site shims: `args` is anything with
    `sample_type` and `skip_civic` attributes, and `out` is the frame whose `.columns`
    decide the CIViC and gnomAD branches. Returns the ordered column list.

    THE FIRST STATEMENT IS NOT THE PIPELINE'S. `main()` writes `keep = KEEP`, which
    aliases the module-level list from `pipeline_utils` and then mutates it in place.
    A pipeline process runs `main()` once and exits, so it never sees the consequence;
    the app keeps this module loaded for the whole session, where a germline call
    leaves the shared list three entries short and the next somatic call silently
    returns germline columns. See `tests/test_vendor_compute_keep.py`.

    The copy is taken here, outside the verbatim block, rather than by editing the
    vendored statements — so the drift guard still compares unedited copies.
    """
    {local_copy}

{statements}

    return keep
'''


def extract_symbols(source: str, names: list) -> dict:
    """Return {name: verbatim source segment} for top-level functions in `names`."""
    tree = ast.parse(source)
    found = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names:
            found[node.name] = ast.get_source_segment(source, node)
    missing = [n for n in names if n not in found]
    if missing:
        raise SystemExit(f"ERROR: symbols not found in {FILTERS_SRC}: {', '.join(missing)}")
    return found


class ColumnSliceError(Exception):
    """`main()` no longer has the shape the column slice is located by."""


def _walk_statements(body: list):
    """Every statement under `body`, outer before inner, in source order."""
    for statement in body:
        yield statement
        for field in ("body", "orelse", "finalbody"):
            yield from _walk_statements(getattr(statement, field, None) or [])
        for handler in getattr(statement, "handlers", None) or []:
            yield from _walk_statements(handler.body)


def _source_segment(source: str, node) -> str:
    """The node's verbatim source, including the indentation of its own first line.

    `ast.get_source_segment(..., padded=True)` pads the first line only for nodes that
    span several lines; for a single-line node it ignores `padded` entirely. Re-adding
    the indentation here is what lets a lifted statement drop straight into a function
    body at the same depth and still be byte-identical to the pipeline's.
    """
    segment = ast.get_source_segment(source, node, padded=True)
    if segment is None:
        raise ColumnSliceError(f"no source segment for the statement at line {node.lineno}")
    if node.end_lineno == node.lineno:
        segment = " " * node.col_offset + segment
    return segment


def _matches_anchor(source: str, node, anchor: tuple) -> bool:
    """True when the node's first source lines start with each line of `anchor`."""
    lines = [line.strip() for line in _source_segment(source, node).splitlines()]
    if len(lines) < len(anchor):
        return False
    return all(line.startswith(want) for line, want in zip(lines, anchor))


def column_slice(source: str) -> list:
    """The eight statements of `main()` that build and apply the output column list.

    Raises `ColumnSliceError` rather than returning a short list: an anchor that finds
    no statement, or more than one, means `main()` was reshaped and the guard no longer
    knows what it is guarding. Failing there is the point — a silently empty slice
    would make every downstream comparison pass vacuously.
    """
    tree = ast.parse(source)
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "main":
            statements = list(_walk_statements(node.body))
            break
    else:
        raise ColumnSliceError(f"{FILTERS_SRC.name} no longer defines main()")

    picked = []
    for index, anchor in enumerate(COLUMN_SLICE_ANCHORS):
        matches = [s for s in statements if _matches_anchor(source, s, anchor)]
        if len(matches) != 1:
            raise ColumnSliceError(
                f"column slice statement {index} matched {len(matches)} statements in "
                f"main() (expected exactly 1); anchor: {' / '.join(anchor)}"
            )
        picked.append(matches[0])

    linenos = [node.lineno for node in picked]
    if linenos != sorted(linenos):
        raise ColumnSliceError(
            f"the column slice statements are out of source order in main(): {linenos}"
        )
    return picked


def render_compute_keep(source: str) -> str:
    """Wrap the vendored statements of the column slice in a callable function."""
    # One blank line between statements, as in `main()`. Blank lines belong to no
    # statement, so this changes nothing the guard compares — it is for the reader.
    located = column_slice(source)
    statements = "\n\n".join(
        _source_segment(source, located[index]) for index in VENDORED_SLICE_INDICES
    )
    # rstrip so this joins with the function segments, which carry no trailing newline.
    return COLUMN_TEMPLATE.format(
        local_copy=LOCAL_COPY_STATEMENT, statements=statements
    ).rstrip("\n")


def render_filters() -> str:
    source = FILTERS_SRC.read_text()
    segments = extract_symbols(source, FILTER_SYMBOLS)
    parts = [segments[name] for name in FILTER_SYMBOLS] + [render_compute_keep(source)]
    return f"{FILTERS_HEADER}\n\n" + "\n\n\n".join(parts) + "\n"


def _function_defs(path: Path) -> dict:
    """Map top-level function name -> ast node for a Python file."""
    tree = ast.parse(path.read_text())
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _without_docstring(body: list) -> list:
    """`body` with a leading docstring removed, if it has one."""
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
        and isinstance(body[0].value.value, str)
    ):
        return body[1:]
    return body


def _compute_keep_findings(vendored_defs: dict, origin_slice: list) -> list:
    """Compare the vendored `compute_keep` against the column slice it was lifted from.

    Its body must be exactly: the docstring, the local copy that replaces the
    pipeline's aliasing statement, the vendored statements, and `return keep`. The
    length is checked before anything is compared, so a body of an unexpected shape
    reports that rather than comparing misaligned statements against each other.
    """
    if COLUMN_SYMBOL not in vendored_defs:
        return [f"MISSING {FILTERS_DEST.relative_to(REPO_ROOT)}: {COLUMN_SYMBOL}"]

    where = f"{FILTERS_DEST.relative_to(REPO_ROOT)}: {COLUMN_SYMBOL}"
    body = _without_docstring(vendored_defs[COLUMN_SYMBOL].body)

    wanted = [origin_slice[index] for index in VENDORED_SLICE_INDICES]
    if len(body) != len(wanted) + 2:
        return [
            f"DRIFT {where} has {len(body)} statements after its docstring; expected "
            f"{len(wanted) + 2} (the local copy, {len(wanted)} vendored, and the return)"
        ]

    findings = []
    if ast.dump(body[0]) != ast.dump(ast.parse(LOCAL_COPY_STATEMENT).body[0]):
        findings.append(
            f"DRIFT {where} no longer starts with `{LOCAL_COPY_STATEMENT}`. That local "
            f"copy is the fix for the pipeline's list aliasing; without it a germline "
            f"call corrupts the column set of every call that follows it."
        )
    if not isinstance(body[-1], ast.Return):
        findings.append(f"DRIFT {where} no longer ends with a return statement")

    for index, vendored, origin in zip(VENDORED_SLICE_INDICES, body[1:-1], wanted):
        if ast.dump(vendored) != ast.dump(origin):
            findings.append(
                f"DRIFT {where}: column slice statement {index} "
                f"(bin/filter_variants.py line {origin.lineno})"
            )
    return findings


def markers_findings() -> list:
    """Compare the vendored CancerVar marker table against the pipeline's copy.

    A whole-file sha256, for the same reason ``pipeline_utils.py`` gets one: there are
    no symbols to compare, only bytes. Here the intolerance is not merely acceptable
    but the point — the MAF cells hold *line offsets* into this file, so a copy that
    differs by one inserted row resolves every later index to the wrong marker and
    names the wrong drug. There is no cosmetic edit to a data file.
    """
    if not MARKERS_DEST.exists():
        return [f"MISSING {MARKERS_DEST.relative_to(REPO_ROOT)}"]
    want = hashlib.sha256(MARKERS_SRC.read_bytes()).hexdigest()
    got = hashlib.sha256(MARKERS_DEST.read_bytes()).hexdigest()
    if want == got:
        return []
    return [
        f"DRIFT {MARKERS_DEST.relative_to(REPO_ROOT)}: sha256 {got[:16]} != "
        f"{MARKERS_SRC.relative_to(REPO_ROOT)} {want[:16]}. The app resolves MAF "
        f"marker indices as line offsets into this file, so a mismatch names the "
        f"wrong drug rather than failing to name one."
    ]


#: The independent source trees the vendored copies are compared against, as
#: ``(label, path)``. Two, not three: ``pipeline_utils.py`` and ``pipeline_filters.py``
#: both come out of ``bin/`` and stand or fall together, while the marker table comes
#: from ``resources/CancerVar/`` and can be absent on its own.
#:
#: A table rather than two ``if``s so that :data:`SOURCE_COUNT` derives from it. A
#: hand-written count is the thing that goes stale when a third source is added, and a
#: count that is too high would make the all-skipped branch unreachable — leaving
#: ``--check`` claiming "in sync" in a tree where nothing was compared at all.
GUARDED_SOURCES = (
    ("pipeline filter code", BIN_DIR),
    ("CancerVar marker table", MARKERS_SRC),
)

#: How many independent sources exist. Derived — never written down.
SOURCE_COUNT = len(GUARDED_SOURCES)


def skipped_units() -> list:
    """Sources that are absent, so no comparison against them is possible.

    Reported separately from findings, and never silently: the README's rule is that a
    skip must be loud enough that a build log cannot read it as a check that ran and
    passed. Per source, because the two are independent — the app subtree alone has
    neither, a pipeline checkout has both, and either can be missing on its own.
    """
    return [
        f"SKIP {label}: {path} not present (packaged app or partial checkout)"
        for label, path in GUARDED_SOURCES
        if not path.exists()
    ]


def drift_report() -> list:
    """Return a list of human-readable drift findings; empty means in sync.

    This deliberately applies the *same* rule as
    ``streamlit_app/tests/test_vendor_drift.py`` — whole-file sha256 for
    ``pipeline_utils.py`` and ``cancervar_markers.txt``, per-symbol ``ast.dump()``
    equality for ``pipeline_filters.py``, and the eight-statement column slice — so
    that the dependency-free build gate and the pytest guard cannot disagree about what
    counts as drift. A comment-only edit in ``bin/`` is drift to neither.
    ``tests/test_vendor_drift.py`` asserts the two agree.

    Each unit is compared only where its own source tree is present; what was skipped
    is :func:`skipped_units`' business, so that an impossible comparison never reads
    here as a passing one.
    """
    findings = []
    if BIN_DIR.is_dir():
        findings.extend(_pipeline_findings())
    if MARKERS_SRC.exists():
        findings.extend(markers_findings())
    return findings


def _pipeline_findings() -> list:
    """The two units vendored out of ``bin/``: ``pipeline_utils`` and ``pipeline_filters``."""
    findings = []

    if not UTILS_DEST.exists():
        findings.append(f"MISSING {UTILS_DEST.relative_to(REPO_ROOT)}")
    else:
        want = hashlib.sha256(UTILS_SRC.read_bytes()).hexdigest()
        got = hashlib.sha256(UTILS_DEST.read_bytes()).hexdigest()
        if want != got:
            findings.append(
                f"DRIFT {UTILS_DEST.relative_to(REPO_ROOT)}: "
                f"sha256 {got[:16]} != bin/utils.py {want[:16]}"
            )

    # The eight-statement column slice, checked whether or not the vendored file is
    # readable: three of the eight have no vendored counterpart, so their only
    # description lives here.
    filters_source = FILTERS_SRC.read_text()
    origin_slice = None
    try:
        origin_slice = column_slice(filters_source)
    except ColumnSliceError as error:
        findings.append(f"SHAPE bin/filter_variants.py: {error}")
    else:
        for index, pinned in UNVENDORED_SLICE_STATEMENTS.items():
            if ast.dump(origin_slice[index]) != ast.dump(ast.parse(pinned).body[0]):
                findings.append(
                    f"DRIFT bin/filter_variants.py line {origin_slice[index].lineno}: "
                    f"column slice statement {index} no longer reads `{pinned}`. It is "
                    f"guarded but not vendored, so nothing else would notice."
                )

    if not FILTERS_DEST.exists():
        findings.append(f"MISSING {FILTERS_DEST.relative_to(REPO_ROOT)}")
    else:
        origin = _function_defs(FILTERS_SRC)
        vendored = _function_defs(FILTERS_DEST)
        for name in FILTER_SYMBOLS:
            if name not in origin:
                findings.append(f"GONE bin/filter_variants.py no longer defines {name}")
            elif name not in vendored:
                findings.append(f"MISSING {FILTERS_DEST.relative_to(REPO_ROOT)}: {name}")
            elif ast.dump(vendored[name]) != ast.dump(origin[name]):
                findings.append(f"DRIFT {FILTERS_DEST.relative_to(REPO_ROOT)}: {name}")
        if origin_slice is not None:
            findings.extend(_compute_keep_findings(vendored, origin_slice))
        for extra in sorted(set(vendored) - set(VENDORED_SYMBOLS)):
            findings.append(f"UNGUARDED {FILTERS_DEST.relative_to(REPO_ROOT)}: {extra}")

    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="report drift and exit non-zero instead of rewriting the vendored files",
    )
    args = parser.parse_args()

    # Asymmetric on purpose. --check is the form the installer build gates, the
    # pre-commit hook and CI all run, and it is reachable from a tree that carries no
    # pipeline — the app subtree on its own, or a bundle rebuilt from one. There the
    # comparison is impossible to make rather than failed, so blocking the build would
    # be wrong; it skips, and says so loudly enough that a build log cannot read the
    # skip as a check that ran and passed. This mirrors the `needs_pipeline` skipif in
    # tests/test_vendor_drift.py, so both forms of the guard behave the same way.
    #
    # The rewrite form below has no such out: with no bin/ to copy from there is
    # nothing it could do, and exiting 0 would report success for work it never did.
    skips = skipped_units()

    if args.check:
        for skip in skips:
            print(skip)
        findings = drift_report()
        for finding in findings:
            print(finding)
        if findings:
            print(
                f"\n{len(findings)} finding(s). The app's filter results can no longer be "
                f"trusted to match the pipeline's.\nFix: python streamlit_app/vendor/_sync.py"
            )
            return 1
        # Every unit skipped: nothing was compared, so say that and nothing stronger.
        # The wording is asserted by `test_sync_check_skips_cleanly_without_bin`, which
        # requires the phrase "in sync" to be absent here — a build log must not be able
        # to read an impossible comparison as a passing one.
        if len(skips) == SOURCE_COUNT:
            print("The vendored copies were not compared against the pipeline.")
            return 0
        print(
            f"vendored copies are in sync "
            f"({SOURCE_COUNT - len(skips)} of {SOURCE_COUNT} sources compared)"
        )
        return 0

    # Rewrite mode has no skip: with no source to copy from there is nothing it could
    # do, and exiting 0 would report success for work it never did. Both source trees
    # are required, so a repair either refreshes every unit or refuses.
    if skips:
        for skip in skips:
            print(skip.replace("SKIP", "ERROR", 1))
        raise SystemExit(
            "ERROR: cannot repair from a tree that is missing a source — run this from "
            "a full pipeline checkout."
        )

    planned = {
        UTILS_DEST: UTILS_SRC.read_bytes(),
        FILTERS_DEST: render_filters().encode(),
        # Bytes, not text: this is the file whose sha256 is the guard, and a round trip
        # through `read_text`/`write_text` would re-encode line endings on Windows and
        # turn the guard red for a copy that is otherwise correct.
        MARKERS_DEST: MARKERS_SRC.read_bytes(),
    }
    written = []
    for dest, wanted in planned.items():
        if (dest.read_bytes() if dest.exists() else None) != wanted:
            dest.write_bytes(wanted)
            written.append(dest)
            print(f"wrote {dest.relative_to(REPO_ROOT)}")
    if not written:
        print("vendored copies already in sync — nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
