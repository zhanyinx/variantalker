"""The filter entry point, asserted with neither ``bin/`` nor the clinical drive.

This is the first layer of the ``bin/``-free net the spec (issue #28) asks for, and it
exists because ``tests/parity/`` is skipped *entirely* on any checkout without ``bin/``
— so before this file, the filter path was unasserted on every packaged build and every
pipeline-less CI job. The two test files this replaces (``test_filters.py``,
``test_integration.py``) were the only other thing touching the filter module, and they
tested three functions that no longer exist as app code.

The oracle is **already in git**: ``tests/fixtures/parity/MANIFEST.json`` records, per
named fixture row, the verdict ``bin/filter_variants.py`` reached at the contract's
default parameters. Those recorded verdicts were produced by the pipeline, so asserting
against them is a real parity assertion that needs no pipeline present.

Three layers, weakest last — deliberately in this order, because the trap in this
codebase is the coincidence:

1. **Per-row verdict equality on the named rows**, addressed *positionally*. This is
   the strong layer. Addressing by position is what turns S1's "every row returned, in
   input order, index preserved" guarantee into a load-bearing assertion rather than a
   claim: if the entry point ever reorders, drops or re-indexes a row, the names stop
   lining up and these tests fail.
2. **Oracle-free contract checks** — every row returned, order and index preserved, both
   verdict columns present, the reason agreeing with the verdict, the four diagnostic
   cells partitioning the frame, and the visible column list opening with the pipeline's
   own (computed from the vendored ``compute_keep``, since there is no ``bin/`` to read
   a header from).
3. **Aggregate PASS counts**, which are the *weakest* layer and are labelled as such
   below, with the measured reason.

Issue #41 closed the net at layer 2's column check and named the weakest layer's reason;
where each behaviour is asserted, and why it is asserted *here* rather than in the
harness, is ``tests/README.md``.

Every assertion here **quotes the cell, not the number**: the somatic reference is
legitimately both 408 and 411 rows depending on whether the criteria path or the union
is meant, and conflating them is how a real discrepancy survived several tickets.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.pipeline_params import pipeline_params  # noqa: E402
from filters.variant_filters import (  # noqa: E402
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    NOPASS,
    PASS,
    REASON_BOTH,
    REASON_CRITERIA,
    REASON_REJECTED,
    REASON_RESCUE,
    apply_filters,
)
from vendor.pipeline_utils import read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"
MANIFEST = json.loads((FIXTURE_DIR / "MANIFEST.json").read_text())

#: The fixtures carrying per-row recorded verdicts, and the key each row's verdict is
#: recorded under. ``somatic_civic.maf`` records two verdicts per row — one per
#: ``skip_civic`` setting — which is why the manifest holds 59 verdicts across 47 named
#: rows. Both prose numbers in the closed tickets are right about different things;
#: this table is derived from the manifest so neither has to be trusted.
VERDICT_FIXTURES = [
    ("somatic_synthetic.maf", "somatic", "filter", {}),
    ("germline_synthetic.maf", "germline", "filter", {}),
    ("somatic_civic.maf", "somatic", "filter_skip_civic_false", {"skip_civic": False}),
    ("somatic_civic.maf", "somatic", "filter_skip_civic_true", {"skip_civic": True}),
]

#: Fixtures with no recorded per-row verdicts, used for the contract layer only.
CONTRACT_FIXTURES = [
    ("somatic_reference.maf", "somatic"),
    ("germline_reference.maf", "germline"),
    ("somatic_synthetic.maf", "somatic"),
    ("germline_synthetic.maf", "germline"),
    ("somatic_civic.maf", "somatic"),
    ("somatic_gnomad_genome.maf", "somatic"),
]


def _params(arm: str, **overrides) -> dict:
    """The contract, as the entry point consumes it.

    Taken from ``config/pipeline_params.py`` rather than written out here. That module is
    guarded against ``nextflow.config`` by ``test_param_contract.py``, so the parameters
    these assertions run under cannot drift from the pipeline's without a test saying so
    — where a copied dict here would drift silently. This is the house rule: where the
    app must agree with the pipeline about a list, derive it or guard it, never copy it.
    """
    params = pipeline_params(arm)
    params.update(overrides)
    return params


def _run(fixture: str, arm: str, **overrides):
    frame = read_maf(str(FIXTURE_DIR / fixture))
    labelled, diagnostics = apply_filters(frame, _params(arm, **overrides))
    return frame, labelled, diagnostics


# ---------------------------------------------------------------------------
# Layer 1 — per-row verdict equality against the recorded pipeline verdicts
# ---------------------------------------------------------------------------


def _verdict_cases():
    """One parameter per (fixture, skip_civic setting, named row): 59 in total."""
    cases = []
    for fixture, arm, verdict_key, overrides in VERDICT_FIXTURES:
        record = MANIFEST["fixtures"][fixture]
        for position, name in enumerate(record["rows_named"]):
            cases.append(
                pytest.param(
                    fixture,
                    arm,
                    verdict_key,
                    overrides,
                    position,
                    name,
                    id=f"{fixture.removesuffix('.maf')}-{verdict_key}-{name}",
                )
            )
    return cases


@pytest.mark.parametrize(
    "fixture,arm,verdict_key,overrides,position,name", _verdict_cases()
)
def test_recorded_row_verdict_reproduces(
    fixture, arm, verdict_key, overrides, position, name
):
    """The entry point reaches the verdict the pipeline reached, row by row.

    Addressed by *position*, not by any key column: that is what makes the row-order and
    index guarantees load-bearing rather than merely documented. ``name`` comes from the
    manifest's ``rows_named``, which is in file order, so ``position`` is both the
    fixture's row and the returned frame's row — and only stays so while the entry point
    returns every row in input order.
    """
    expected = MANIFEST["fixtures"][fixture]["expected"][name][verdict_key]
    _, labelled, _ = _run(fixture, arm, **overrides)

    actual = labelled.iloc[position][MAFIGATE_FILTER]
    assert actual == expected, (
        f"{fixture} row {position} ({name}): the pipeline recorded {expected} for this "
        f"row at the contract defaults ({verdict_key}), the app returned {actual}. "
        f"Reason column says {labelled.iloc[position][MAFIGATE_REASON]!r}."
    )


def test_every_recorded_verdict_is_covered():
    """All 59 recorded verdicts are asserted, so none can be quietly skipped.

    Two closed tickets quote 47 and 59 for "the recorded per-row verdicts", and the ticket
    warns to read the manifest and trust it rather than either number. Read: both are
    right about different things — 47 named rows, 59 verdicts, because the CIViC fixture
    records two verdicts per row (one per ``skip_civic`` setting). Counted from the
    manifest here so the reconciliation is checked rather than asserted in prose.

    The pair is then pinned deliberately: it is a tripwire, so that growing a fixture
    fails here — where the message says what to do — instead of silently leaving new rows
    unasserted. The last clause is what actually prevents that: it compares the
    parametrised count against the manifest, so the two can never drift.
    """
    recorded = sum(
        1
        for record in MANIFEST["fixtures"].values()
        for row in record.get("expected", {}).values()
        for key in row
        if key.startswith("filter")
    )
    named = sum(
        len(record["expected"])
        for record in MANIFEST["fixtures"].values()
        if "expected" in record
    )
    assert (named, recorded) == (47, 59), (
        f"the fixture manifest now records {recorded} verdicts across {named} named "
        "rows, not 59 across 47. If a fixture grew, extend VERDICT_FIXTURES and this "
        "expectation together."
    )
    assert len(_verdict_cases()) == recorded, (
        f"{len(_verdict_cases())} of {recorded} recorded verdicts are parametrised — "
        "VERDICT_FIXTURES is missing a fixture or a skip_civic setting"
    )


# ---------------------------------------------------------------------------
# Layer 2 — the contract, with no oracle needed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_every_row_is_returned_in_input_order_with_its_index(fixture, arm):
    """S1's central guarantee: nothing is dropped, reordered or re-indexed.

    This is what lets a caller derive the passed and failed frames itself, what makes
    layer 1's positional addressing valid, and what makes the verdict column meaningful
    as a column rather than as a filtered frame.
    """
    frame, labelled, _ = _run(fixture, arm)

    assert len(labelled) == len(frame), (
        f"{fixture}: {len(frame)} rows in, {len(labelled)} out — the entry point must "
        "return every input row and let the caller split them"
    )
    assert list(labelled.index) == list(frame.index), (
        f"{fixture}: the returned index is not the input index in input order"
    )


@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_both_verdict_columns_are_present_and_populated(fixture, arm):
    """Every row carries a PASS/NOPASS verdict and one of the four reasons."""
    _, labelled, _ = _run(fixture, arm)

    assert MAFIGATE_FILTER in labelled.columns
    assert MAFIGATE_REASON in labelled.columns
    assert set(labelled[MAFIGATE_FILTER]) <= {PASS, NOPASS}
    assert set(labelled[MAFIGATE_REASON]) <= {
        REASON_CRITERIA,
        REASON_RESCUE,
        REASON_BOTH,
        REASON_REJECTED,
    }
    assert not labelled[MAFIGATE_FILTER].isna().any()
    assert not labelled[MAFIGATE_REASON].isna().any()


@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_the_reason_cannot_disagree_with_the_verdict(fixture, arm):
    """Both columns come off the same masks, so the pairing is exact by construction.

    Asserted anyway, because "derived from the very masks that made the decision" is a
    claim about the implementation that a later refactor could break while leaving both
    columns individually plausible.
    """
    _, labelled, _ = _run(fixture, arm)

    passed = labelled[MAFIGATE_FILTER] == PASS
    rejected = labelled[MAFIGATE_REASON] == REASON_REJECTED
    assert bool((passed ^ ~rejected).sum() == 0), (
        f"{fixture}: a row is PASS with reason {REASON_REJECTED!r}, or NOPASS with a "
        "passing reason — the verdict and the reason have come apart"
    )


@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_the_four_cells_partition_the_frame(fixture, arm):
    """The decomposition is a partition: the cells sum to the row count.

    Asserted as a partition rather than against fixed numbers, on purpose. Fixed cell
    counts would be a second baseline to maintain, and the property that actually matters
    — that the four cells are a Boolean rearrangement of the masks that made the
    decision, and so cannot overlap or leave a row out — is exactly this.
    """
    _, labelled, diagnostics = _run(fixture, arm)

    cells = diagnostics.cells()
    assert sum(cells.values()) == len(labelled) == diagnostics.rows, (
        f"{fixture}: cells {cells} do not sum to {len(labelled)} rows"
    )
    assert diagnostics.passed == (labelled[MAFIGATE_FILTER] == PASS).sum()

    # And each cell equals the count of its own reason label, which is the other
    # direction of the same claim: the channel and the column agree.
    for reason, count in (
        (REASON_CRITERIA, cells["criteria_only"]),
        (REASON_RESCUE, cells["rescue_only"]),
        (REASON_BOTH, cells["both"]),
        (REASON_REJECTED, cells["rejected"]),
    ):
        assert (labelled[MAFIGATE_REASON] == reason).sum() == count, (
            f"{fixture}: diagnostics says {count} rows are {reason!r} but the column "
            "says otherwise"
        )


@pytest.mark.parametrize("skip_civic", [False, True])
@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_the_visible_columns_open_with_the_pipeline_list_on_this_frame(
    fixture, arm, skip_civic
):
    """The last of layer 2's contract checks: the column prefix, with no ``bin/``.

    The pipeline's half is not read off a pipeline output file here — there is none
    without ``bin/`` — and it is not read off ``config/columns.py`` either, which would
    be asserting the resolver against itself. It is computed by calling the **vendored**
    ``compute_keep``: the pipeline's own function, held byte-for-byte and guarded against
    ``bin/`` by ``test_vendor_drift.py``. So this is a real claim about the pipeline's
    column list on a checkout that has no pipeline, which is exactly what the net is for.

    Over the *fixtures'* frames rather than a synthetic column list, which is the
    difference from ``test_column_resolver.py``: ``compute_keep`` branches on which
    columns the frame actually carries (the CIViC strip, the gnomAD genome branch), so a
    real MAF exercises branches a hand-written list would have to remember to include.

    Both ``skip_civic`` settings, because the column list is a function of all three of
    arm, flag and frame. Only ``somatic_civic.maf`` carries CIViC columns at all, so it is
    the only fixture where the flag can change the answer — and it is exactly the fixture
    where the strip has to be right, since every other one reaches the same list by the
    columns simply being absent. A sweep fixed at ``False`` would have left the branch
    unexercised on real data and looked complete.

    A prefix, not equality — the app is a deliberate superset since issue #35, appending
    its own extras after everything the pipeline emits.

    Less ``PIPELINE_COLUMNS_THE_APP_REPLACES`` since issue #117, the columns the app
    answers itself and keeps out of the default view. The subtraction is taken from the
    app's own constant, so this still checks the resolver against the *vendored* pipeline
    function and not against itself: what the constant can do is excuse a name, and what
    it cannot do is put a name back in the pipeline's list or change the order of the rest.
    """
    from types import SimpleNamespace

    from config.columns import (
        PIPELINE_COLUMNS_THE_APP_REPLACES,
        resolve_visible_columns,
    )
    from vendor.pipeline_filters import compute_keep

    frame, labelled, _ = _run(fixture, arm)

    # The pipeline's half is computed from the frame the *pipeline* would be handed — the
    # one loaded off disk, before the app appends its two verdict columns. ``compute_keep``
    # branches on which columns are present, so handing it the labelled frame would be
    # asking it about a MAF that does not exist. Harmless today, since it looks for no name
    # the app adds; wrong in principle, and the branch it takes is the whole point here.
    pipeline = [
        col
        for col in compute_keep(
            SimpleNamespace(sample_type=arm, skip_civic=skip_civic), frame
        )
        if col not in PIPELINE_COLUMNS_THE_APP_REPLACES
    ]

    # The app's half sees what the app sees: the labelled frame, verdict columns included.
    app = resolve_visible_columns(
        sample_type=arm, skip_civic=skip_civic, available_columns=list(labelled.columns)
    )

    assert app[: len(pipeline)] == pipeline, (
        f"{fixture} (skip_civic={skip_civic}): the app's visible columns do not open "
        f"with the pipeline's list.\n"
        f"  first mismatch at index "
        f"{next((i for i, (p, a) in enumerate(zip(pipeline, app)) if p != a), None)}\n"
        f"  pipeline-only: {[c for c in pipeline if c not in app]}\n"
        f"  app-only in the prefix: "
        f"{[c for c in app[: len(pipeline)] if c not in pipeline]}"
    )
    assert MAFIGATE_FILTER in app[len(pipeline) :], (
        f"{fixture}: the verdict column is not among the app's extras — the grid would "
        "show a filtered frame with nothing saying why"
    )


@pytest.mark.parametrize("fixture,arm", CONTRACT_FIXTURES)
def test_the_input_frame_is_not_mutated(fixture, arm):
    """The caller's frame is left alone — it is Streamlit session state.

    The deleted implementation added a ``CancerVar`` column to the caller's frame in
    place, which the parity harness had to defend against by copying every fixture it
    loaded.
    """
    frame, labelled, _ = _run(fixture, arm)

    assert MAFIGATE_FILTER not in frame.columns
    assert MAFIGATE_REASON not in frame.columns


#: Which fixtures can tell the two cells apart. They differ exactly where the
#: unconditional pathogenic rescue admits a row the criteria path did not — measured at
#: contract defaults as criteria path 19 against union PASS 20 on the somatic reference,
#: and 27 against 31 on the germline.
#:
#: Those four numbers are recorded here and **not asserted**, deliberately. Fixed cell
#: counts would be a second baseline to keep in step with ``parity/baseline.json``, which
#: is the thing this suite says not to do everywhere else it touches the decomposition.
#: The property that matters is the *inequality*, and that is what the test asserts.
CELL_SEPARATING_FIXTURES = [
    ("somatic_reference.maf", "somatic"),
    ("germline_reference.maf", "germline"),
]

#: And where it cannot: at contract defaults these three have an empty ``rescue_only``
#: cell, so the criteria path and the union coincide and an assertion written on the
#: wrong one still passes. Named rather than left to be discovered.
CELL_BLIND_FIXTURES = [
    ("somatic_synthetic.maf", "somatic"),
    ("germline_synthetic.maf", "germline"),
    ("somatic_civic.maf", "somatic"),
]


@pytest.mark.parametrize("fixture,arm", CELL_SEPARATING_FIXTURES)
def test_the_criteria_path_and_the_union_are_different_numbers_here(fixture, arm):
    """The suite's vocabulary rule is only enforceable while some fixture separates them.

    "Quote the cell, not the number" is the rule — the somatic reference baseline is
    legitimately **408** on the criteria path and **411** on the union, and calling either
    "the somatic baseline" is how a real discrepancy stayed open across three tickets. A
    rule about which cell you mean has no teeth on data where the two cells hold the same
    number: an assertion written on the wrong one passes, and keeps passing until someone
    runs it on the reference.

    This pins the discriminating power itself, and nothing else. If a fixture edit ever
    collapsed these two to one number, every cell-sensitive assertion in the suite would
    quietly become untestable here, and this is where that is noticed. What it does *not*
    do is pin the two counts: those live in ``parity/baseline.json`` already, and a second
    copy of them here would be a second baseline to re-record.
    """
    _, _, diagnostics = _run(fixture, arm)

    assert diagnostics.criteria_path < diagnostics.passed, (
        f"{fixture}: the criteria path and the union PASS cell now hold the same number, "
        "so this fixture can no longer tell an assertion written on one from an "
        "assertion written on the other. Restore a row the pathogenic rescue admits."
    )


@pytest.mark.parametrize("fixture,arm", CELL_BLIND_FIXTURES)
def test_the_cell_blind_fixtures_are_still_the_ones_we_think(fixture, arm):
    """The other half: these three cannot separate the cells, and that is recorded.

    Not a defect — they exist to isolate other behaviour, and the synthetic pair is
    deliberately built without a rescued row. Asserted so the list above stays true: a
    fixture that quietly gained a rescued row would become cell-separating without anyone
    noticing it had, and the useful fixture for a cell-sensitive assertion would still be
    documented as useless.
    """
    _, _, diagnostics = _run(fixture, arm)

    assert diagnostics.rescue_only == 0, (
        f"{fixture} now rescues {diagnostics.rescue_only} row(s), so it separates the "
        "criteria path from the union after all — move it to CELL_SEPARATING_FIXTURES"
    )
    assert diagnostics.criteria_path == diagnostics.passed


def test_the_pathogenic_rescue_is_unconditional_and_the_flag_is_passed_down():
    """``skip_pathogenic`` reaches the vendored code; the app does not branch on it.

    The vendored functions already return an all-``False`` rescue mask when the flag is
    set, so there is nothing for the app to decide. The observable consequence is that
    setting the flag can only ever *remove* rescued rows, never add any — asserted here
    as a subset relation rather than as a count, because counts can coincide.
    """
    frame = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))

    kept, kept_diag = apply_filters(frame, _params("somatic", skip_pathogenic=False))
    skipped, skipped_diag = apply_filters(
        frame, _params("somatic", skip_pathogenic=True)
    )

    assert skipped_diag.rescue_only == 0 and skipped_diag.both == 0, (
        "skip_pathogenic=True still rescued rows — the flag is not reaching the "
        f"vendored functions (cells: {skipped_diag.cells()})"
    )
    # The criteria cell is untouched by the flag: it is the other half of the union.
    assert skipped_diag.criteria_only + skipped_diag.both == (
        kept_diag.criteria_only + kept_diag.both
    ), "skip_pathogenic changed the criteria path, which it must not touch"

    passing = set(skipped.index[skipped[MAFIGATE_FILTER] == PASS])
    assert passing <= set(kept.index[kept[MAFIGATE_FILTER] == PASS]), (
        "skipping the pathogenic rescue admitted a row that keeping it did not"
    )


def test_the_no_gene_filter_sentinel_is_the_one_the_vendored_code_honours():
    """``NO_GENE_FILE`` is a literal copied out of the vendored bodies — so guard it.

    The vendored gene clause is ``if somatic_genes != "null":``, and the app has to hand it
    that exact string to mean "no gene filter". A copied literal is what the house rule
    forbids ("derive it or guard it, never copy it"), and it cannot be derived: in ``bin/``
    it is a bare string inside two function bodies, not a named constant.

    So it is guarded *behaviourally* instead, which is stronger than comparing spellings —
    it asserts the consequence. The mask the vendored function returns for the sentinel must
    be the all-True "no restriction" mask, identical to what a gene list containing every
    symbol in the MAF would produce. If ``bin/`` ever changes the sentinel, this fails.
    """
    from filters.variant_filters import NO_GENE_FILE, _gene_file
    from vendor.pipeline_filters import somatic_filters

    frame = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))
    contract = pipeline_params("somatic")

    def specific_mask(genes_path):
        mask, _ = somatic_filters(
            frame,
            vaf=-1.0,  # so the VAF gate cannot mask the gene clause
            somatic_genes=genes_path,
            cancervar_keep=contract["filter_cancervar"],
            civic_keep=contract["filter_civic"],
            escat_keep=contract["filter_escat"],
            clinvar_keep=contract["filter_clinvar"],
            skip_civic=True,
            skip_pathogenic=True,
        )
        return mask

    with _gene_file([]) as sentinel_path:
        assert sentinel_path == NO_GENE_FILE
        unrestricted = specific_mask(sentinel_path)

    every_symbol = sorted(set(frame["Hugo_Symbol"].dropna().astype(str)))
    with _gene_file(every_symbol) as listed_path:
        listing_everything = specific_mask(listed_path)

    assert unrestricted.equals(listing_everything), (
        f"passing {NO_GENE_FILE!r} no longer means 'no gene filter' to the vendored "
        "code — the sentinel in bin/ has changed and filters/variant_filters.py's "
        "NO_GENE_FILE is stale"
    )


def test_the_gene_adapter_removes_its_temp_file():
    """The adapter cleans up, and does so around the call rather than before it.

    Two failures this guards, both silent: a file left behind on every re-filter click on a
    platform that does not reap its temp directory, and — the worse one — a file removed
    *before* the vendored call, which checks for the file itself and applies no gene filter
    at all when it is missing, widening the report with no warning.
    """
    from filters.variant_filters import _gene_file

    with _gene_file(["BRCA1", "TP53"]) as path:
        assert os.path.exists(path), (
            "the gene file must still exist inside the context manager — the vendored "
            "body does its own existence check, so an early close means no gene filter"
        )
        assert Path(path).read_text().split() == ["BRCA1", "TP53"]

    assert not os.path.exists(path), f"the adapter leaked {path}"


def test_the_germline_arm_has_no_escat_clause():
    """The single largest divergence on real data — 540 rows, 51% — closes by routing.

    The app used to OR ``ESCAT.isin(filter_escat)`` into the *germline* guideline block.
    The pipeline's germline guidelines are ``InterVar | ClinVar | RENOVO``, with no ESCAT
    term at all. Handing the germline arm an ESCAT keep-list must therefore change
    nothing whatsoever: the parameter has no germline meaning.

    The keep-list is taken from the *somatic* contract rather than written out, for the
    same reason :func:`_params` gives: a copied list here could go stale against
    ``nextflow.config`` and quietly turn this into a test that passes an empty or
    unmatched list, which would prove nothing.
    """
    escat_keep = pipeline_params("somatic")["filter_escat"]
    assert escat_keep, "the somatic contract has no ESCAT terms — this test is vacuous"

    frame = read_maf(str(FIXTURE_DIR / "germline_reference.maf"))

    without, _ = apply_filters(frame, _params("germline"))
    with_escat, _ = apply_filters(frame, _params("germline", filter_escat=escat_keep))

    assert list(without[MAFIGATE_FILTER]) == list(with_escat[MAFIGATE_FILTER]), (
        "an ESCAT keep-list changed the germline verdicts — the clause the pipeline "
        "does not have is still in the germline guideline block"
    )


# ---------------------------------------------------------------------------
# Layer 3 — aggregates. DELIBERATELY THE WEAKEST LAYER.
# ---------------------------------------------------------------------------
#
# Read this before adding to it. A PASS *total* can match while the sets differ, and both
# halves of that sentence were measured in this repository rather than imagined:
#
#   * `somatic_depth_500` sat at **pipeline 10, app 10** in the pre-#33 baseline — equal
#     totals — while one row passed only in the pipeline and a *different* row passed only
#     in the app. A totals-only check called that case parity. It was two divergences
#     pointing opposite ways, which is the same shape as the PIK3CA row that hid
#     divergence #6 for the whole effort. The historical table in
#     `tests/parity/README.md` is where those figures live.
#   * Column parity was once measured 40-against-40 while the contents differed by the
#     substitution -`variantalker_naive` +`gnomAD_exome_AF`.
#
# So these tests are a coarse smoke check that the union is not wildly off, and nothing
# more. Layer 1 is where per-row correctness is actually established — if you are
# tempted to assert a new number here, assert a row there instead.


@pytest.mark.parametrize(
    "fixture,arm,verdict_key,overrides,expected_union_pass",
    [
        ("somatic_synthetic.maf", "somatic", "filter", {}, 7),
        ("germline_synthetic.maf", "germline", "filter", {}, 7),
        ("somatic_civic.maf", "somatic", "filter_skip_civic_false", {}, 8),
        ("somatic_civic.maf", "somatic", "filter_skip_civic_true", {"skip_civic": True}, 0),
    ],
)
def test_union_pass_total_agrees_with_the_recorded_verdicts(
    fixture, arm, verdict_key, overrides, expected_union_pass
):
    """The **union PASS** cell, quoted as that cell and derived from the manifest.

    The number is named ``expected_union_pass`` and not ``expected_pass`` because the
    criteria path and the union are different cells with different totals, and every
    assertion in this suite has to say which one it means.
    """
    recorded = sum(
        1
        for row in MANIFEST["fixtures"][fixture]["expected"].values()
        if row[verdict_key] == PASS
    )
    assert recorded == expected_union_pass, (
        f"{fixture}: the manifest now records {recorded} union-PASS rows for "
        f"{verdict_key}, not {expected_union_pass} — update this case deliberately"
    )

    _, labelled, diagnostics = _run(fixture, arm, **overrides)
    assert diagnostics.passed == expected_union_pass
    assert (labelled[MAFIGATE_FILTER] == PASS).sum() == expected_union_pass
