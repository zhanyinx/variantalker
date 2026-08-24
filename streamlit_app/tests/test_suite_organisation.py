"""The suite's own organisation, asserted rather than described.

``tests/README.md`` states the governing rule — **the harness owns every behaviour with a
pipeline counterpart, the unit suite owns every behaviour without one**, plus a ``bin/``-free
net for the claims that have a counterpart but must hold where the pipeline is absent — and
then names every module against it in one table.

A rule written down is a rule until the first hurried afternoon. This file makes the table
load-bearing:

* every test module on disk has a row, and every row names a module that exists — so a new
  file cannot land without someone deciding which instrument it belongs to, and a deleted
  one cannot leave a row behind that reads as coverage;
* the instrument is one of the four the rule allows;
* the table's ``needs bin/`` column agrees with what the module actually does, read off its
  own ``skipif`` reasons;
* the one exemption to the first of those — the modules the public export strips, which
  this table names and a published clone does not have — is itself checked rather than
  taken on trust, wherever there is a tree that can check it.

That last bullet is issue #348, and the failure it fixes was in the *published* repository
rather than this one: a row for a stripped module read to this file as a deleted one, so a
clone of what we publish failed two of its own tests. The fix is not a softer first bullet.
A row for a deleted file really is coverage that is not there, and that assertion keeps its
teeth for every module but the named few — which have to go on proving they are stripped.

That last one is the one with teeth, and issue #24 is why. A module-level
``skipif(not PIPELINE_AVAILABLE)`` silently removed the entire parity suite from every
checkout without ``bin/`` — including the packaged app — and the suite went on reporting
green. The failure was not the skip; it was that nothing anywhere said the skip existed. A
``net`` or ``unit`` module that quietly acquires a pipeline gate is that failure returning,
and it now has to get past this file first.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
README = TESTS_DIR / "README.md"
REPO_ROOT = TESTS_DIR.parents[1]
DENY_LIST = REPO_ROOT / ".publicignore"

#: The four instruments.
INSTRUMENTS = {"harness", "net", "unit", "guard"}

#: The modules this table names that a **published** clone does not have, **by name**.
#:
#: This table travels and the one-way public export does not carry everything the table
#: names, so in an exported tree a row here points at a file that is not on disk. Before
#: issue #348 that made the published repository's own ``make test`` red — two failures,
#: 2942 passing — which is the first thing a contributor checks and the worst possible
#: first impression. The guard was *right*: a row for an absent module really is coverage
#: that is not there. So the answer is not to soften it into silence but to teach it the
#: one shape of absence that is deliberate, and to keep it strict about every other.
#:
#: Named rather than derived — the same argument as
#: :data:`MODULES_ALLOWED_TO_NEED_THE_PIPELINE` above it, and for the same reason: an
#: exemption a module could infer from how it is filed is an exemption every module can
#: take. Adding a name here is a deliberate act and should be argued for in the commit
#: that does it.
#:
#: An exemption granted on trust is the vacuous guard this suite keeps finding, so this
#: one is **checked wherever it can be**:
#:
#: * in this repository, where ``.publicignore`` is present,
#:   :func:`test_every_declared_exemption_is_real` requires each name below to be matched
#:   by a real deny-list rule *and* to exist on disk. So the exemption goes red the moment
#:   it stops being true in either direction — the deny-list stops naming the module, or
#:   the module is simply deleted and the row left behind, which is the failure the whole
#:   check exists for;
#: * in an exported tree, where ``.publicignore`` is itself stripped, the absence cannot
#:   be re-derived from anything that travels, so exactly these names may be missing and
#:   nothing else can be.
#:
#: One name. ``test_public_export.py`` is the export's own proof: it performs real exports
#: from a private checkout, reading the export script and ``.publicignore``, both
#: of which the deny-list also strips. It cannot travel, and its own docstring says so.
MODULES_STRIPPED_FROM_THE_PUBLIC_TREE = {
    "test_public_export.py",
}

#: The modules allowed to need the pipeline, **by name**. Named rather than derived from
#: the instrument, so that no module can exempt itself by how it is filed — the first
#: version of this file said "harness and guard may need bin/", which quietly exempted
#: this module too, since it files itself as a guard. Seven names, each of which has to
#: compare something against the pipeline checkout to do its job at all:
#:
#: * the three parity modules run ``bin/filter_variants.py`` in a subprocess;
#: * ``test_vendor_drift.py`` compares ``vendor/`` to ``bin/``;
#: * ``test_param_contract.py`` compares the contract to ``nextflow.config``;
#: * ``test_chromosome_rule_contract.py`` compares ``config/chromosome_spelling.py``'s
#:   rule to the one ``bin/add_guidelines_escat_to_funcotator.py`` applies (issue #149).
#:   The fifth name, added deliberately: the app's chromosome rule is a *copy* of the
#:   pipeline's, exactly as the contract is a copy of ``nextflow.config``'s defaults, and
#:   a copy with no guard on it diverges silently — here, into MAFigate spelling a
#:   chromosome differently from the report the file came out of. The claim is "the app
#:   agrees with the pipeline", so it cannot be checked without the pipeline. Its twenty
#:   sibling assertions about the app's own behaviour stay in
#:   ``test_chromosome_spelling.py``, which needs no ``bin/`` at all.
#:
#: * ``test_predictor_cutoff_contract.py`` compares the five predictor thresholds
#:   ``components/predictor_context.py`` *prints on screen* to the assignments
#:   ``resources/CancerVar/CancerVar.py`` and ``resources/InterVar/Intervar.py`` make
#:   (issue #190). The sixth name, and the same argument as the fifth with one addition:
#:   the panel does not merely agree with the number, it **states** it to a clinician —
#:   *"damaging — below CancerVar's 0.05 cutoff"* — so a drifted copy is not a quiet
#:   behavioural divergence but a false sentence on a clinical screen. Map #184 keeps
#:   ``resources/`` read-only for this effort, which makes the app's copy the only side
#:   that can move, and this the only thing that would notice. ``resources/`` rather than
#:   ``bin/``, but the same dependency: it ships with the pipeline and not with the app.
#:   Its sibling assertions about the app's own behaviour stay in
#:   ``test_predictor_context.py``, which needs no pipeline at all.
#:
#: * ``parity/test_mutation_coverage.py`` measures what the parity fixtures still *catch*,
#:   by re-injecting each of the map's divergences into the app and asking whether the
#:   pipeline still agrees (issue #242). The seventh name, and the first added rather than
#:   removed since this list was written, so it owes the larger argument: a ``bin/``-free
#:   version of it would have to compare the app against a recorded answer, which is a mock
#:   of the exact side the measurement exists to consult — and the failure it guards
#:   against is a fixture set that agrees with everything because it exercises nothing,
#:   which no recording can detect. What it *removed* from the suite,
#:   ``test_attribution_coverage.py``, needed no pipeline, and the half of it that did not
#:   invert lives on in ``parity/test_baseline_integrity.py``, which still does not.
#:
#: Everything else must hold in a checkout with no pipeline in it. Adding an eighth name is
#: a deliberate act, and should be argued for in the commit that does it — as the seventh
#: is, above.
MODULES_ALLOWED_TO_NEED_THE_PIPELINE = {
    "parity/test_parity.py",
    "parity/test_absent_columns.py",
    "parity/test_mutation_coverage.py",
    "test_vendor_drift.py",
    "test_param_contract.py",
    "test_chromosome_rule_contract.py",
    "test_predictor_cutoff_contract.py",
}

#: What marks a gate as a *pipeline* gate rather than any other precondition. Matched
#: against the reason string, which is the one part of a gate that is written for a human
#: and so the one part that is reliably descriptive. ``vendor/ not present`` and
#: ``installer scripts not present`` are gates too, and deliberately do not match: they are
#: about this app's own files, and they hold in a checkout that has no pipeline.
PIPELINE_TOKENS = ("pipeline", "bin/", "nextflow.config")


def _modules_on_disk() -> set[str]:
    """Every module pytest would collect, as a path relative to ``tests/``."""
    return {
        path.relative_to(TESTS_DIR).as_posix()
        for path in TESTS_DIR.rglob("test_*.py")
        if "__pycache__" not in path.parts
    }


def _table() -> dict[str, tuple[str, str]]:
    """The README's table: ``module -> (instrument, needs_bin)``.

    Parsed between the header row and the blank line that ends the table, so a future
    table elsewhere in the file cannot be read as this one.
    """
    rows: dict[str, tuple[str, str]] = {}
    lines = README.read_text().splitlines()
    try:
        start = next(
            i for i, line in enumerate(lines) if line.startswith("| module | instrument |")
        )
    except StopIteration:  # pragma: no cover - the assertion below reports it
        return rows

    for line in lines[start + 2 :]:
        if not line.startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        rows[cells[0].strip("`")] = (cells[1], cells[2])
    return rows


def _literal_text(node: ast.AST) -> str:
    """The literal part of a string node — plain or f-string.

    An f-string is read for the text around its holes, because that is where the human
    part of the message lives: ``f"{column} is a filter input the pipeline does not
    emit"`` says "pipeline" whatever ``column`` turns out to be.
    """
    if isinstance(node, ast.Constant):
        return str(node.value)
    if isinstance(node, ast.JoinedStr):
        return " ".join(
            str(part.value) for part in node.values if isinstance(part, ast.Constant)
        )
    return ""


def _skip_reasons(path: Path) -> list[str]:
    """Every reason a module can skip for, statically.

    **Both** spellings, and the second is not a nicety. This module was first written
    reading only ``pytest.mark.skipif(..., reason=...)``, which left it blind to
    ``pytest.skip("…")`` inside a fixture — the very pattern the same change introduced
    in ``test_param_contract.py``, on the argument that a gate a call site can forget is
    a gate that will be forgotten. A module that adopted the better pattern would have
    read as needing no pipeline, passed this file vacuously, and skipped in silence:
    issue #24's failure walking back in through the door held open for it.

    Non-literal reasons still cannot be read — a reason built by calling a function, say.
    That is the known edge, and it is why the ``needs bin/`` column exists to be declared
    by a human as well as checked here.
    """
    reasons = []
    for node in ast.walk(ast.parse(path.read_text())):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
        if name not in {"skipif", "skip"}:
            continue
        if name == "skip" and node.args:
            reasons.append(_literal_text(node.args[0]))
        for keyword in node.keywords:
            if keyword.arg == "reason":
                reasons.append(_literal_text(keyword.value))
    return [reason for reason in reasons if reason]


def _is_pipeline_gated(module: str) -> bool:
    return any(
        any(token in reason.lower() for token in PIPELINE_TOKENS)
        for reason in _skip_reasons(TESTS_DIR / module)
    )


TABLE = _table()


def test_the_readme_table_is_readable_at_all():
    """A parse failure here would make every other test in this file vacuously green.

    Asserted against a known row rather than a row *count*, so a half-parsed table cannot
    slip through under a threshold. This row is a good canary because it is the one the
    whole suite is organised around: the harness, which does need the pipeline.
    """
    assert TABLE.get("parity/test_parity.py") == ("harness", "yes"), (
        f"tests/README.md's table did not parse as expected — the row for the harness "
        f"read {TABLE.get('parity/test_parity.py')!r}, and {len(TABLE)} rows were read in "
        "total. The header line or the column order has changed, and every check below is "
        "now reading something other than the table it thinks it is."
    )


def test_every_module_has_a_home_and_every_home_has_a_module():
    """No test file without an instrument; no instrument without a test file.

    With one exemption in the second direction, for the modules the public export
    deliberately does not carry: see :data:`MODULES_STRIPPED_FROM_THE_PUBLIC_TREE`, and
    :func:`test_every_declared_exemption_is_real` for what stops it becoming a hiding
    place. The exemption is subtracted rather than added — a row that is *not* on the
    list and *not* on disk still fails, which is every ordinary deleted file.
    """
    on_disk = _modules_on_disk()
    listed = set(TABLE)

    assert not (on_disk - listed), (
        f"test modules with no row in tests/README.md: {sorted(on_disk - listed)}. Decide "
        "which instrument each belongs to — harness (has a pipeline counterpart and "
        "compares against it live), net (has a counterpart, asserted without bin/), unit "
        "(no counterpart), guard (a copy held against its source) — and add the row."
    )
    missing = listed - on_disk - MODULES_STRIPPED_FROM_THE_PUBLIC_TREE
    assert not missing, (
        f"tests/README.md names modules that do not exist: {sorted(missing)}. A row for a "
        "deleted file reads as coverage that is not there. If the module is absent because "
        "the public export strips it, name it in MODULES_STRIPPED_FROM_THE_PUBLIC_TREE — "
        "which then has to prove the strip is real."
    )


def _deny_list_rules() -> list[str]:
    """``.publicignore``'s rules — its non-comment, non-blank lines."""
    return [
        line.strip()
        for line in DENY_LIST.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def _is_denied(rules: list[str], path: str) -> bool:
    """Whether some rule in ``rules`` covers ``path``, relative to the repository root.

    Deliberately the *narrow* half of ``.publicignore``'s stated grammar — a rule ending
    in ``/`` strips that directory and everything under it, any other rule names one exact
    file — and not an import of the real matcher, because the real matcher lives in
    the export script, which the deny-list also strips. A check that has to hold
    in the exported tree cannot import something the exported tree does not have.

    That makes this a copy, which the house rule in ``README.md`` forbids where the
    original is reachable. It is tolerable only because the copy is doing strictly less:
    it answers "does a rule name this path" and never decides what to strip. The full
    grammar, and the claim that every rule still matches something tracked, stay with
    ``test_public_export.py`` — which reads the real matcher, in the tree that has it.
    """
    return any(
        path == rule.lstrip("/")
        if not rule.endswith("/")
        else path.startswith(rule.lstrip("/"))
        for rule in rules
    )


def test_every_declared_exemption_is_real():
    """An exemption that cannot go red is not an exemption, it is a hole.

    :data:`MODULES_STRIPPED_FROM_THE_PUBLIC_TREE` lets a row name a module that is not on
    disk, which is precisely the thing the sibling check above exists to refuse. So in the
    one tree that can tell — this one, the private checkout, which is also the only tree
    anyone edits — the exemption has to earn itself in both directions:

    * the deny-list really does name the module. A name left here after the rule went away
      would exempt a module nothing strips, and the row would then be false in the public
      tree as well as here;
    * the module really is on disk. This is the direction that matters most, because it is
      the ordinary failure the whole check protects: someone deletes the test, leaves the
      row, and the exemption quietly absorbs it as though it had been stripped.

    In an exported tree neither can be asked — ``.publicignore`` is stripped along with
    the module — so this skips there rather than passing vacuously.
    """
    if not DENY_LIST.exists():
        pytest.skip(
            "no .publicignore in this tree: it is stripped along with the modules it "
            "names, so the strip cannot be re-derived here"
        )

    rules = _deny_list_rules()
    assert rules, (
        f"{DENY_LIST} parsed to no rules at all, so both assertions below would be "
        "quantified over nothing and pass whatever the exemption said."
    )

    for module in sorted(MODULES_STRIPPED_FROM_THE_PUBLIC_TREE):
        tracked_path = f"streamlit_app/tests/{module}"
        assert _is_denied(rules, tracked_path), (
            f"{module} is named in MODULES_STRIPPED_FROM_THE_PUBLIC_TREE, but no "
            f".publicignore rule covers {tracked_path}. Either the rule was renamed — in "
            "which case the module is travelling under a new name and the deny-list is "
            "the thing to fix — or the exemption has outlived it and should go."
        )
        assert (TESTS_DIR / module).exists(), (
            f"{module} is exempted as stripped-by-the-export, but it does not exist in "
            "this tree either. An exemption is not a way to retire a module: delete the "
            "row in tests/README.md, the .publicignore rule and this name together."
        )


@pytest.mark.parametrize("module", sorted(TABLE))
def test_the_instrument_is_one_of_the_four(module):
    instrument, _ = TABLE[module]
    assert instrument in INSTRUMENTS, (
        f"{module} is filed as {instrument!r}, which is not one of {sorted(INSTRUMENTS)}"
    )


@pytest.mark.parametrize("module", sorted(TABLE))
def test_the_needs_bin_column_matches_what_the_module_does(module):
    """The column is a constraint, not a description.

    Read off the module's own ``skipif`` reasons rather than trusted, because the whole
    point of the ``bin/``-free net is that it holds where ``bin/`` does not exist, and a
    net module that gained a pipeline gate would go on *passing* — by skipping — while
    asserting nothing. That is issue #24's failure, and a green suite is exactly how it
    presented the first time.
    """
    _, declared = TABLE[module]
    assert declared in {"yes", "no"}, (
        f"{module}: the 'needs bin/' cell reads {declared!r}, not 'yes' or 'no'"
    )

    # The cell is still checked for legality above, because that is a claim about the
    # table and the table travels. What cannot be checked is the half read off the
    # module's own source, there being no source here to read: in an exported tree this
    # parametrisation used to raise an uncaught FileNotFoundError (issue #348).
    if not (TESTS_DIR / module).exists():
        assert module in MODULES_STRIPPED_FROM_THE_PUBLIC_TREE  # the check above
        pytest.skip(f"{module} is not in this tree — the public export strips it")

    actual = _is_pipeline_gated(module)
    assert actual == (declared == "yes"), (
        f"{module}: tests/README.md says it {'needs' if declared == 'yes' else 'does not need'} "
        f"bin/, but its skip reasons say otherwise ({_skip_reasons(TESTS_DIR / module)}). "
        f"A gate counts as a pipeline gate when its reason mentions one of {PIPELINE_TOKENS}."
    )


@pytest.mark.parametrize("module", sorted(TABLE))
def test_only_the_named_modules_need_the_pipeline(module):
    """The governing rule's boundary, in the direction that can rot silently.

    Seven modules have to reach the pipeline to do their job, and they are named in
    :data:`MODULES_ALLOWED_TO_NEED_THE_PIPELINE`. Every other module must hold in a
    checkout with no pipeline in it — which is what keeps the packaged app's filter path
    asserted, and what ``make test-cov`` then actually runs.

    Keyed on the module's name and not on how it is filed. An instrument-based exemption
    ("harness and guard may need bin/") let a module exempt itself by choosing a label,
    and this file was the first to do it: it files itself as a guard.
    """
    _, declared = TABLE[module]
    if module in MODULES_ALLOWED_TO_NEED_THE_PIPELINE:
        assert declared == "yes", (
            f"{module} is named as one of the modules that needs the pipeline, but the "
            "table says it does not. If the dependency really is gone, take it out of "
            "MODULES_ALLOWED_TO_NEED_THE_PIPELINE — the list is meant to shrink."
        )
        return
    assert declared == "no", (
        f"{module} declares that it needs bin/, but it is not one of the modules "
        f"allowed to: {sorted(MODULES_ALLOWED_TO_NEED_THE_PIPELINE)}. Either the "
        "dependency is accidental and should go, or this is genuinely another module "
        "that compares something against the pipeline — in which case add it there, in a "
        "commit that says why."
    )
