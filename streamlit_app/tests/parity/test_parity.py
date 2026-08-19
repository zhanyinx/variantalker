"""Parity assertions against the recorded baseline.

Two layers, and the distinction matters:

* **Baseline tests** (green today) assert the measured numbers equal ``baseline.json``.
  They pin what the app does *now*, so an unrelated change that moves parity shows up
  as a diff to explain rather than as silence.
* **Parity tests** assert the destination — identical PASS sets, identical output
  columns. Every case is there, so these are plain assertions: no case carries an
  expected-failure marker and there is no longer any machinery to give one.

**Why the ratchet is gone rather than merely unused** (issue #41). While divergences were
open, a diverging case carried ``xfail(strict=True)`` so the suite went red the moment the
case started passing, and whoever fixed it re-recorded ``baseline.json`` and dropped the
marker in the same commit. That ratchet did its job — it is what kept issue #28's run of
divergence fixes honest, and it was proven in both directions on a real one-line fix for
divergence #7. With the last divergence closed it inverts: a marker that can still be
handed out is a marker a *regression* can be absorbed into, by re-recording the baseline
and calling the case "not yet at parity" again. Issue #35 made exactly this argument for
:func:`test_column_parity` when column parity became structural; issue #41 applies it to
the rows. A case that stops passing must now go red as a failure, which is the only
reading left that is true.

Green here is not green-by-skipping, and that claim needs the checks below rather than
this sentence — everything in this module needs the pipeline checkout and **skips
without it**, since the packaged .dmg/.exe ships neither ``bin/`` nor these tests. Two
things hold the line where this module cannot:

* ``tests/parity/test_baseline_integrity.py`` reads the committed baseline with no
  ``bin/`` and asserts that no case diverges — so "nothing is xfailed" is a measured
  property of a checkout without the pipeline, not an absence of evidence;
* ``tests/test_filter_entry_point.py`` asserts the verdicts themselves against the
  pipeline verdicts frozen in ``fixtures/parity/MANIFEST.json``, which is the layer that
  keeps the filter path asserted where this module is skipped.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(STREAMLIT_APP))

from filters.variant_filters import (  # noqa: E402
    FREQUENCY_COLUMNS,
    MAFIGATE_FILTER,
    apply_filters,
    frequency_mask,
)
from tests.parity import harness as H  # noqa: E402
from tests.parity.contract import (  # noqa: E402
    APP_VARIANT_CLASSIFICATIONS,
    CASES,
    KEY_COLUMNS,
    app_params,
)

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not H.PIPELINE_AVAILABLE,
        reason="pipeline bin/ not present (packaged app or partial checkout)",
    ),
]

_BASELINE = H.load_baseline()["cases"]
_results: dict[str, H.ParityResult] = {}


def result_for(name: str) -> H.ParityResult:
    """One pipeline subprocess per case per session, not per assertion."""
    if name not in _results:
        _results[name] = H.compare(H.CASES_BY_NAME[name])
    return _results[name]


def _case_params(predicate=lambda record: True):
    """Parametrise over the cases the baseline records, filtered by ``predicate``.

    It takes no ``xfail_when`` any more, and the omission is the point rather than a
    simplification — see the module docstring. ``predicate`` selects which cases an
    assertion *applies to* (a case where the pipeline raised has no column list to
    compare), which is a different question from whether a case is allowed to fail.
    """
    return [
        pytest.param(case.name, id=case.name)
        for case in CASES
        if predicate(_BASELINE[case.name])
    ]


# ---------------------------------------------------------------------------
# The instrument's own preconditions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", sorted({c.fixture for c in CASES}))
def test_comparison_key_is_unique(fixture):
    """The four-column variant key must identify a row, or every diff is meaningless.

    Issue #18 measured it unique within any single reference MAF and deduplicated on it
    while pooling; this asserts the property actually holds on what was committed.
    ``Tumor_Sample_Barcode`` cannot be added as a tiebreaker — it is ``__UNKNOWN__`` in
    every germline reference row.
    """
    frame = H.read_maf(H.FIXTURE_DIR / fixture)
    keys = H.variant_keys(frame)
    duplicates = {k for k in keys if keys.count(k) > 1}
    assert not duplicates, (
        f"{fixture}: {len(duplicates)} duplicate variant keys on "
        f"{KEY_COLUMNS} — the PASS-set diff cannot be trusted"
    )


#: The harness's loading premise — that its frame is the frame the app loads — used to
#: be asserted here, as the *impossibility* of loading with the app's own loader
#: (``test_app_loader_cannot_reach_parity``, issue #10). Issue #32 made it an identity
#: instead, and moved it to ``test_loader_premise.py``, which carries no
#: ``skipif(not PIPELINE_AVAILABLE)`` because the claim needs no ``bin/``.


def test_harness_is_excluded_from_packaged_builds():
    """Neither installer may ship ``tests/`` — this harness or its 1.8 MiB of fixtures.

    Both already exclude it, and for different reasons, so both are checked: the macOS
    script removes ``tests`` from the staged bundle, while ``installer.iss`` is an
    allowlist of ``Source:`` lines that has never named it. Asserted rather than
    assumed, because the fixtures carry real germline variant calls from 50 patients
    (see the fixture set's privacy note) and must not leave the repository inside an
    installer.
    """
    dmg = STREAMLIT_APP / "build" / "mac" / "build_dmg.sh"
    iss = STREAMLIT_APP / "build" / "windows" / "installer.iss"
    if not (dmg.exists() and iss.exists()):
        pytest.skip("installer scripts not present (build/ is stripped from the app)")

    assert 'rm -rf "${DEST}/tests"' in dmg.read_text(), (
        "build/mac/build_dmg.sh no longer strips tests/ from the bundle — the parity "
        "fixtures would ship inside the .dmg"
    )

    sources = [
        line for line in iss.read_text().splitlines() if line.startswith("Source:")
    ]
    offenders = [line for line in sources if "tests" in line]
    assert not offenders, (
        f"build/windows/installer.iss gained a Source: line covering tests/: {offenders}"
    )


def test_app_vocabulary_unchanged():
    """``contract.APP_VARIANT_CLASSIFICATIONS`` must still mirror the app's own list.

    The contract holds its own copy so the baseline stays interpretable after issue #12
    moves or deletes the app's list; this keeps the copy honest until then.
    """
    from config.vocabularies import VARIANT_CLASSIFICATIONS

    live = [v for v in VARIANT_CLASSIFICATIONS if v != "All"]
    assert live == APP_VARIANT_CLASSIFICATIONS, (
        "config/vocabularies.py:VARIANT_CLASSIFICATIONS changed. If this is issue #12's "
        "fix landing, re-measure the baseline; otherwise the mirror has drifted."
    )


#: The map's divergence #12 — "the app matches CIViC with an exact ``isin`` that cannot
#: match list-repr at all" — was measured not to exist: the app's own
#: ``has_element_from_list`` was a substring test behaviourally identical to the
#: pipeline's, and ``test_civic_matching_already_agrees`` pinned the two implementations
#: against each other value by value. Issue #33 deleted the app's copy, so there is no
#: longer a pair to compare: the app calls the pipeline's function. The claim is now
#: structural rather than measured, which is strictly stronger, and the ``civic_present``
#: and ``civic_skipped`` cases continue to check it end-to-end.


# ---------------------------------------------------------------------------
# The baseline
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", _case_params())
def test_baseline_matches(name):
    """The measured numbers still equal the committed baseline.

    Regenerate deliberately:
    ``python streamlit_app/tests/parity/harness.py --update-baseline``
    """
    measured = result_for(name).to_json()
    expected = _BASELINE[name]
    assert measured == expected, (
        f"{name}: parity measurement moved.\n"
        f"  expected {expected}\n"
        f"  measured {measured}\n"
        "If this is a fix landing, re-record the baseline; if not, it is a regression."
    )


# ---------------------------------------------------------------------------
# The destination
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", _case_params())
def test_row_parity(name):
    """Identical PASS sets, and agreeing failure behaviour where either side fails.

    "Agreeing" rather than "identical" since issue #38: on a MAF whose depth or VAF columns
    cannot be read as numbers the pipeline raises ``TypeError`` and the app **refuses** with
    a typed exception naming every offending column and value. Those are the same
    non-verdict, and requiring the strings to match would have meant the only route to
    parity was reproducing the pipeline's traceback. ``ParityResult.errors_in_parity`` holds
    the rule, including what it still refuses to accept.
    """
    result = result_for(name)
    assert result.errors_in_parity, (
        f"{name}: pipeline error {result.pipeline_error!r} but app error "
        f"{result.app_error!r} (refused columns {result.app_refused_columns!r}) — "
        "one side fails where the other does not"
    )
    assert result.only_pipeline == 0 and result.only_app == 0, (
        f"{name}: {result.only_pipeline} variants PASS only in the pipeline, "
        f"{result.only_app} only in the app.\nAttribution: {result.attribution}"
    )


@pytest.mark.parametrize(
    "name", _case_params(predicate=lambda record: record["app_error"] is None)
)
def test_diagnostics_decomposition_reproduces_the_verdict(name):
    """The app's four-cell account must agree with the PASS set it actually produced.

    The decomposition is what the app tells the user about *why* each row is in the report,
    so a decomposition that did not add up to the verdict would be a confidently wrong
    explanation — the failure mode the 28 deleted diagnostic strings had, where the
    classification count read 68,593 against a true 20,386.

    Two claims, and both matter:

    * it is a **partition** — the four cells sum to the row count, so no row is
      double-counted or dropped;
    * the three passing cells sum to the **union PASS** count, which where the case is in
      row parity is the *pipeline's* PASS count too.

    Asserted as a partition rather than against four recorded numbers, deliberately: fixed
    cell counts would be a second baseline to keep in step, and the property that makes the
    channel trustworthy is this arithmetic, not any particular value. `reference.py` runs
    the stronger per-variant form of the same check over all 1,100 reference measurements.
    """
    result = result_for(name)
    cells = result.app_cells
    assert cells is not None, f"{name}: the app returned no diagnostics"

    assert sum(cells.values()) == result.rows, (
        f"{name}: cells {cells} sum to {sum(cells.values())}, not the {result.rows} rows "
        "handed in — the decomposition is not a partition"
    )

    union_pass = cells["criteria_only"] + cells["rescue_only"] + cells["both"]
    assert union_pass == result.app_pass, (
        f"{name}: the decomposition says {union_pass} rows pass but the verdict column "
        f"says {result.app_pass}"
    )
    if result.in_parity and result.pipeline_error is None:
        assert union_pass == result.pipeline_pass, (
            f"{name}: union PASS {union_pass} against the pipeline's "
            f"{result.pipeline_pass} — note this is the *union* cell, not the criteria "
            f"path, which is {cells['criteria_only'] + cells['both']}"
        )


#: This was the first assertion to be taken off the ratchet, by issue #35: column parity
#: became *structural* — both sides run the same vendored ``compute_keep`` — so a case
#: could not drift back without someone breaking that, and leaving the marker wired up
#: would have meant a future regression could be absorbed by re-recording the baseline.
#: Issue #41 extended the same reasoning to :func:`test_row_parity` and removed the
#: machinery outright, so the ``predicate`` below is now the only thing ``_case_params``
#: does. It still earns its keep: a case where the pipeline raised has no column list.
@pytest.mark.parametrize(
    "name", _case_params(predicate=lambda record: record["pipeline_error"] is None)
)
def test_column_parity(name):
    """The app's column list must open with the pipeline's ``keep``, element for element.

    A prefix, not equality: since issue #35 the app is a deliberate superset, appending
    its own derived columns after everything the pipeline emits. The app's list is
    ``config.columns.resolve_visible_columns``, imported and called; the pipeline's is
    read off a real output file.

    Compared position by position, never by length — see ``columns_in_parity``.

    Since issue #117 the pipeline's list is taken less the columns the app answers itself
    (``ParityResult.expected_prefix``). Every other pipeline column must still be there,
    in the pipeline's order, and the names let off are recorded in the baseline — so
    dropping a further column fails here, and dropping this one deliberately shows up as a
    baseline diff rather than as silence.
    """
    result = result_for(name)
    pipeline, app = result.expected_prefix, result.app_columns
    prefix = app[: len(pipeline)]
    assert prefix == pipeline, (
        f"{name}: the app's columns do not open with the pipeline's.\n"
        f"  first mismatch at index "
        f"{next((i for i, (p, a) in enumerate(zip(pipeline, prefix)) if p != a), len(prefix))}\n"
        f"  pipeline-only: {[c for c in pipeline if c not in app]}\n"
        f"  app-only:      {[c for c in prefix if c not in pipeline]}\n"
        f"  lengths:       pipeline {len(pipeline)}, app {len(app)}"
    )


# ---------------------------------------------------------------------------
# Neutrality of the app-only extras (issue #16)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", sorted({c.fixture for c in CASES}))
def test_frequency_mask_is_tautological_at_1_0(fixture):
    """At ``max_freq_population = 1.0`` the frequency mask is all-True *by algebra*.

    Not by the ``if max_freq_population < 1.0:`` guard — issues #16 and #37 require the
    stronger claim, because the composition they specify has no guard to hide behind: the
    mask now gates the pipeline's whole verdict, so a mask that is merely *usually* all-True
    at 1.0 would cost parity rather than merely widen the criteria path. Allele frequencies
    live in [0, 1], so ``(v <= 1.0) | v.isna()`` is a tautology.

    This calls the shipped mask with the guard bypassed rather than restating the
    expression. An earlier version of this test kept its own copy of both the expression and
    the frequency-column list, which made it a test of the copy.
    """
    frame = H.read_maf(H.FIXTURE_DIR / fixture)
    present = [c for c in FREQUENCY_COLUMNS if c in frame.columns]
    if not present:
        pytest.skip(f"{fixture} carries no frequency columns")

    mask = frequency_mask(frame, 1.0)
    assert bool(mask.all()), (
        f"{fixture}: the frequency mask is not all-True at threshold 1.0 over "
        f"{present} — a value outside [0, 1] breaks issue #16's neutrality argument"
    )


@pytest.mark.parametrize(
    "name",
    [
        "somatic_defaults",
        "germline_defaults",
        "civic_present",
        "somatic_synthetic",
        "germline_synthetic",
    ],
)
def test_frequency_extra_is_neutral_end_to_end(name):
    """Declaring the extra at 1.0 gives the same verdicts as not declaring it at all."""
    case = H.CASES_BY_NAME[name]
    frame = H.load_frame(case)

    declared = H.run_app(case, frame)
    params = app_params(case)
    del params["max_freq_population"]

    labelled, _ = apply_filters(frame, params)
    omitted = dict(zip(H.variant_keys(labelled), labelled[MAFIGATE_FILTER]))

    assert declared.verdicts == omitted, (
        f"{name}: max_freq_population=1.0 changed the verdict set — the app-only "
        "frequency extra is not neutral at its parity value"
    )


# ---------------------------------------------------------------------------
# The refusal (issue #38)
# ---------------------------------------------------------------------------


def test_the_unreadable_fixture_is_refused_against_the_pipelines_raise():
    """``dot_numeric`` is in parity *because* the app refuses where the pipeline raises.

    The one case in the set where neither side produces a verdict, and the case whose
    manifest note used to read "both sides must raise". Both sides no longer raise the same
    thing, on purpose: the pipeline dies with ``TypeError: '>=' not supported between
    instances of 'str' and 'int'``, naming neither the column nor the value, and the app
    refuses with a message naming all three offending columns and the ``.`` in each.

    That is parity, not divergence — the app is neither fussier nor laxer than the pipeline
    about which files can be filtered — and this test says so directly rather than leaving
    it to be inferred from ``errors_in_parity`` returning ``True``.
    """
    result = result_for("dot_numeric")

    assert result.pipeline_error is not None, (
        "dot_numeric no longer makes the pipeline raise, so there is nothing for the "
        "app's refusal to be in parity with — the fixture has changed shape"
    )
    assert result.pipeline_error.startswith("TypeError")
    assert result.app_error == "UnreadableNumericColumns"
    assert result.app_refused_columns == ["t_alt_count", "t_ref_count", "tumor_f"]
    assert result.in_parity


def test_no_other_case_is_refused():
    """The refusal fires on exactly one fixture. Everything else still gets a verdict.

    The direction that matters most: a validator that refused more widely would be fussier
    than the pipeline, and would do it silently — every over-refused case would simply stop
    producing a PASS set, and ``in_parity`` would go on reading ``True`` for the cases that
    still worked.
    """
    refused = [
        r.case for r in (result_for(c.name) for c in H.CASES) if r.app_refused_columns
    ]
    assert refused == ["dot_numeric"]
