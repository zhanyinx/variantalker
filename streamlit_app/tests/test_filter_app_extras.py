"""Behaviour the pipeline has no counterpart for, so the harness cannot judge it.

The governing rule for this suite: **the parity harness owns every behaviour with a
pipeline counterpart; the unit suite owns every behaviour without one.**
``test_filter_entry_point.py`` is the third instrument — pipeline-matching behaviour
asserted against verdicts frozen in git, needing no ``bin/``. This file is the second: the
parts of the filter entry point that are the app's own and that no comparison with
``bin/filter_variants.py`` could ever check.

Three groups, and each exists because deleting ``test_filters.py`` would otherwise have
left it unasserted:

1. **The population-frequency layer.** The pipeline has no frequency filter at any value,
   so this is app-owned entirely. It is also the *only* filter logic left in the module,
   and it is the regime every live preset uses — ``test_param_contract.py`` asserts
   ``max_freq_population < 1.0`` for all four of SOFT and CLINICAL. The parity harness runs
   it only at 1.0, where it is deliberately inert, so without this file the app's one piece
   of real filter logic would be exercised nowhere.
2. **The transitional parameter fallbacks.** ``_Settings`` accepts the older
   ``keep_pathogenic``, ``vaf_threshold_germline`` and ``somatic_genes``/``germline_genes``
   spellings. The app opens on the contract since #36, but the SOFT and CLINICAL presets
   and every saved parameter file still carry these names, so the fallbacks hold the
   running app together until #40 migrates them. They have no pipeline counterpart at all
   — the pipeline never saw these names.
3. **What happens when the arm does not match the MAF.** Not a parity question either: the
   pipeline would also fail, but the app is what a clinician is sitting in front of.

Frequency assertions quote the **union**, not the criteria path, and that changed with
#37: the layer used to be joined into the criteria mask and is now composed on top of the
pipeline's verdict. Quoting the wrong cell is how a real discrepancy survived several
tickets, so every assertion here names the cell it means.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.pipeline_params import pipeline_params  # noqa: E402
from filters import variant_filters  # noqa: E402
from filters.variant_filters import (  # noqa: E402
    FREQUENCY_COLUMNS,
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    PASS,
    REASON_BOTH,
    REASON_REJECTED,
    REASON_RESCUE,
    apply_filters,
    frequency_mask,
    pathogenic_exemption,
)
from page_modules import data_loading  # noqa: E402
from tests.fakes import note_texts  # noqa: E402
from vendor.pipeline_utils import CLINVAR_PATHO, has_clinvar_term, read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"


@pytest.fixture(scope="module")
def somatic():
    return read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))


@pytest.fixture(scope="module")
def germline():
    return read_maf(str(FIXTURE_DIR / "germline_reference.maf"))


def params(arm: str, **overrides) -> dict:
    """The contract, with overrides. Derived, never copied — see the house rule."""
    merged = pipeline_params(arm)
    merged.update(overrides)
    return merged


def pipeline_verdict(maf: pd.DataFrame, arm: str, **overrides) -> pd.Series:
    """The verdict with the frequency layer inert, i.e. the pipeline's own PASS set.

    The layer is neutral at 1.0 — asserted end-to-end below and algebraically in
    ``tests/parity/test_parity.py`` — so running the entry point there is the honest way to
    obtain the mask the layer is supposed to compose *on top of*, without reaching into the
    vendored calls and rebuilding the gene-file plumbing in the test.
    """
    labelled, _ = apply_filters(maf, params(arm, max_freq_population=1.0, **overrides))
    return labelled[MAFIGATE_FILTER] == PASS


#: The two arms, as ``(arm, fixture)`` — the fixture is named for its arm, so one parameter
#: does both jobs and ``request.getfixturevalue`` resolves it.
BOTH_ARMS = pytest.mark.parametrize("arm", ["somatic", "germline"])


# ---------------------------------------------------------------------------
# 1. The population-frequency layer (app-owned; no pipeline counterpart)
# ---------------------------------------------------------------------------


def test_the_frequency_layer_is_inert_at_its_parity_value(somatic):
    """1.0 must change nothing, or having the layer available costs parity.

    The end-to-end half of the neutrality claim. The algebraic half — that the comparison
    is a tautology over [0, 1] with the ``< 1.0`` guard bypassed — lives in
    ``tests/parity/test_parity.py``, because it is about the mask expression rather than
    about the entry point.
    """
    wide, wide_diag = apply_filters(somatic, params("somatic", max_freq_population=1.0))
    omitted_params = params("somatic")
    del omitted_params["max_freq_population"]
    omitted, omitted_diag = apply_filters(somatic, omitted_params)

    assert list(wide[MAFIGATE_FILTER]) == list(omitted[MAFIGATE_FILTER])
    assert wide_diag.cells() == omitted_diag.cells()


def test_a_real_frequency_threshold_removes_rows_from_the_union(somatic):
    """Below 1.0 the layer does something, and it does it to the whole verdict.

    Asserted as a strict subset rather than against a count: the point is the direction —
    tightening the threshold can only ever remove rows, never add any — and a count here
    would be a second baseline to maintain.

    The union is the cell that moves now. Before #37 the mask was joined into the criteria
    mask, so a common variant kept by pathogenic rescue stayed in the report; it is now
    composed on top of the pipeline's verdict, so the rescue path is masked too — with the
    ClinVar exemption as its one carve-out.
    """
    loose, loose_diag = apply_filters(somatic, params("somatic", max_freq_population=1.0))
    tight, tight_diag = apply_filters(
        somatic, params("somatic", max_freq_population=0.01)
    )

    assert tight_diag.passed < loose_diag.passed, (
        "a 0.01 threshold removed nothing from the union on a fixture where 44 of 82 "
        f"rows exceed it (union PASS {tight_diag.passed} vs {loose_diag.passed})"
    )

    # Every row the tight threshold still admits was admitted by the loose one.
    assert set(tight.index[tight[MAFIGATE_FILTER] == PASS]) <= set(
        loose.index[loose[MAFIGATE_FILTER] == PASS]
    )


@BOTH_ARMS
def test_the_layer_composes_on_top_of_the_pipeline_verdict(arm, request):
    """The composition, asserted as the equation #37 specifies.

    ``final_pass = pipeline_pass & (frequency acceptable | ClinVar-pathogenic)``.

    This is the whole point of the ticket, and it is asserted against the equation rather
    than against a row count because the two rejected alternatives are both expressible as
    counts that look reasonable. Joining the mask into the vendored common filters leaves
    115 of the 136 disputed reference rows in the report — the union re-admits them through
    pathogenic rescue — while masking the union with no exemption drops a genuinely
    pathogenic low-penetrance allele along with the polymorphisms. Only the equation tells
    those three apart. The reference measurement is in ``docs/wayfinder/issue-37/``.
    """
    maf = request.getfixturevalue(arm)
    threshold = 0.01

    passed = pipeline_verdict(maf, arm)
    expected = passed & (frequency_mask(maf, threshold) | pathogenic_exemption(maf))

    tight, _ = apply_filters(maf, params(arm, max_freq_population=threshold))
    assert list(tight[MAFIGATE_FILTER] == PASS) == list(expected)


def test_the_exemption_is_the_vendored_clinvar_call_and_nothing_else(somatic):
    """The exemption is the pipeline's own pathogenic test, reused not restated.

    Two claims, because either one alone can hold while the pair drifts: the mask equals
    what the vendored function returns over the vendored constant, and the constant the
    module holds *is* the vendored object rather than a copy of its contents.

    The column name is spelled out rather than imported from
    ``variant_filters.CLINVAR_SIGNIFICANCE``, deliberately and against the house rule about
    copies: importing it would make this a test of consistency with itself. The next test
    pins that constant against the vendored source instead, which is the direction that can
    actually be wrong.
    """
    expected = somatic["ClinVar_VCF_CLNSIG"].apply(
        lambda value: has_clinvar_term(value, CLINVAR_PATHO)
    )
    assert list(pathogenic_exemption(somatic)) == list(expected)
    assert bool(expected.any()), "fixture carries no ClinVar-pathogenic rows"

    assert variant_filters.CLINVAR_PATHO is CLINVAR_PATHO, (
        "the module holds its own copy of the pathogenic list — that is the second list "
        "#37 exists to prevent"
    )


def test_the_exemption_reads_a_column_the_vendored_filters_read():
    """``CLINVAR_SIGNIFICANCE`` must name a column the vendored code actually subscripts.

    The one identifier #37 copies rather than imports. The copy is safe because a drift
    cannot be silent — the vendored calls read the same column first, so a rename raises
    there before the exemption runs — but "cannot be silent" is a claim about the vendored
    source, so it is checked against the vendored source rather than asserted in a comment.

    By AST over the subscripts, not by substring: the point is that the vendored code reads
    *this column*, which a mention in one of its docstrings would not establish.
    """
    vendored = (STREAMLIT_APP / "vendor" / "pipeline_filters.py").read_text()

    subscripted = {
        node.slice.value
        for node in ast.walk(ast.parse(vendored))
        if isinstance(node, ast.Subscript)
        and isinstance(node.slice, ast.Constant)
        and isinstance(node.slice.value, str)
    }
    assert variant_filters.CLINVAR_SIGNIFICANCE in subscripted, (
        f"the vendored filters no longer read {variant_filters.CLINVAR_SIGNIFICANCE!r} — "
        "the pathogenic exemption is filtering on a column the pipeline does not, which is "
        f"the drift this constant's comment says is unreachable. Read: {sorted(subscripted)}"
    )

    # The gene column, guarded here rather than in a second near-identical test. Its drift
    # *would* be silent, which is why it earns a guard more than the one above: the app tests
    # `GENE_SYMBOL in maf.columns` before deciding whether to write a gene file at all, so a
    # renamed column reads as "this MAF has no gene symbols", drops the restriction, and
    # widens the report — with a warning that says so, which is exactly the wrong warning.
    assert variant_filters.GENE_SYMBOL in subscripted, (
        f"the vendored filters no longer read {variant_filters.GENE_SYMBOL!r} — the app's "
        "gene adapter is testing for a column the pipeline does not match on, and would "
        f"silently stop applying gene restrictions. Read: {sorted(subscripted)}"
    )


def test_the_module_holds_no_second_pathogenic_list():
    """No list, tuple or set literal in the filter module names a ClinVar term.

    The AST rather than the text, so prose about ``Pathogenic/Pathogenic,_low_penetrance``
    in a docstring — which the module has, because the reference's witness variant is worth
    naming — cannot fail this, while a re-typed keep-list cannot pass it. The same shape of
    guard as ``tests/test_vendor_drift.py``, applied to a value instead of to code.
    """
    source = (STREAMLIT_APP / "filters" / "variant_filters.py").read_text()
    terms = set(CLINVAR_PATHO)

    restated = [
        element.value
        for node in ast.walk(ast.parse(source))
        if isinstance(node, (ast.List, ast.Tuple, ast.Set))
        for element in node.elts
        if isinstance(element, ast.Constant) and element.value in terms
    ]
    assert not restated, (
        f"filters/variant_filters.py restates ClinVar term(s) {restated} in a literal "
        "collection instead of importing CLINVAR_PATHO from the vendored module"
    )


@BOTH_ARMS
def test_the_exemption_keeps_the_pathogenic_allele_and_drops_the_polymorphisms(arm, request):
    """At a moved threshold, the exemption splits exactly along clinical significance.

    Both fixtures carry the reference's witness: ``SERPINA1`` chr14:94847262
    ``Pathogenic/Pathogenic,_low_penetrance|other`` at a population frequency of 2.3%, the
    only ClinVar-pathogenic call in the reference above 0.01. Without the exemption a 0.01
    threshold discards it — a genuinely pathogenic low-penetrance allele dropped for being
    common, which is the failure mode this carve-out exists for.

    The rows on the other side of the split are what makes it a real test rather than a
    tautology: on the somatic fixture five common variants lose their PASS, four of them
    ClinVar ``Benign`` and one ``drug_response``, all surviving only on a
    potential-significance CancerVar tier.

    The witness is *derived* — the pathogenic rows above the threshold — rather than looked
    up by position, so the test still means what it says if the fixtures are regenerated.
    """
    maf = request.getfixturevalue(arm)
    threshold = 0.01

    unfiltered, _ = apply_filters(maf, params(arm, max_freq_population=1.0))
    passed = unfiltered[MAFIGATE_FILTER] == PASS
    common = ~frequency_mask(maf, threshold)
    exempt = pathogenic_exemption(maf)

    kept = passed & common & exempt
    dropped = passed & common & ~exempt
    assert int(kept.sum()) >= 1, "fixture has no common pathogenic row to exempt"
    assert int(dropped.sum()) >= 1, "fixture has no common non-pathogenic row to drop"

    # The witness is a *low-penetrance* call: its ClinVar string is not itself a member of
    # CLINVAR_PATHO, and only matches once split on the separators the vendored function
    # knows. That is the property that makes it the hard case — an exemption written as an
    # exact-match `isin` would drop it while looking correct on every ordinary Pathogenic
    # row — so it is asserted rather than only named in prose.
    significances = set(maf.loc[kept, "ClinVar_VCF_CLNSIG"].dropna())
    assert any(value not in set(CLINVAR_PATHO) for value in significances), (
        "the exempted rows are all exact CLINVAR_PATHO matches, so this fixture no longer "
        f"witnesses the split-only low-penetrance case: {significances}"
    )

    # And the dropped rows must include one the *pathogenic rescue* was carrying, or the
    # test cannot tell the new composition from the old one: before #37 the mask was joined
    # into the criteria path, where the union re-admitted exactly these rows.
    rescued = unfiltered[MAFIGATE_REASON].isin([REASON_RESCUE, REASON_BOTH])
    assert int((dropped & rescued).sum()) >= 1, (
        "no dropped row came from the pathogenic union, so this fixture cannot witness "
        "the mask reaching the rescue path"
    )

    tight, _ = apply_filters(maf, params(arm, max_freq_population=threshold))
    verdicts = tight[MAFIGATE_FILTER]

    lost = maf.loc[kept & (verdicts != PASS), ["Hugo_Symbol", "ClinVar_VCF_CLNSIG"]]
    assert set(verdicts[kept]) == {PASS}, (
        "the frequency layer dropped a ClinVar-pathogenic variant for being common: "
        f"{lost.to_dict('records')}"
    )
    assert PASS not in set(verdicts[dropped]), (
        "a common, non-pathogenic variant kept its PASS — the exemption is wider than "
        "ClinVar-pathogenic, or the mask never reached the rescue path"
    )


def test_the_exemption_survives_skip_pathogenic(somatic):
    """Turning off the pipeline's rescue must not turn off the app's exemption.

    Two different questions that both say "pathogenic". ``skip_pathogenic`` empties the
    *pipeline's* rescue mask — which variants the pipeline reports — while the exemption
    bounds what the **app's own** frequency filter may throw away. A user asking for a
    stricter report is not asking the frequency cut-off to start overruling ClinVar, so a
    pathogenic variant that reaches the report on its own criteria keeps its exemption.

    Worth pinning because the natural mis-implementation — reusing the rescue mask as the
    exemption — passes every other test in this file and fails only here.
    """
    threshold = 0.01
    strict = params("somatic", skip_pathogenic=True)

    unfiltered, _ = apply_filters(somatic, {**strict, "max_freq_population": 1.0})
    criteria_pass = unfiltered[MAFIGATE_FILTER] == PASS
    witness = criteria_pass & ~frequency_mask(somatic, threshold) & pathogenic_exemption(somatic)
    assert int(witness.sum()) >= 1, (
        "no common pathogenic row survives on criteria alone here, so this fixture cannot "
        "witness the interaction"
    )

    tight, _ = apply_filters(somatic, {**strict, "max_freq_population": threshold})
    assert set(tight.loc[witness, MAFIGATE_FILTER]) == {PASS}, (
        "with skip_pathogenic the exemption stopped applying — the frequency layer is "
        "reading the pipeline's rescue mask instead of ClinVar"
    )


def test_the_warning_reports_the_removals_and_names_the_exemption(somatic):
    """The count is unreadable unless the message says pathogenic rows were spared.

    A user reading "5 variants removed" on a fixture where 6 of the passing rows are above
    the threshold has been told something that looks wrong. The exemption is why, so the
    warning states it and gives its own count.
    """
    threshold = 0.01
    passed = pipeline_verdict(somatic, "somatic")
    common = ~frequency_mask(somatic, threshold)
    exempt = pathogenic_exemption(somatic)
    removed = int((passed & common & ~exempt).sum())
    exempted = int((passed & common & exempt).sum())
    assert removed and exempted, "fixture cannot witness both halves of the message"

    _, diagnostics = apply_filters(somatic, params("somatic", max_freq_population=threshold))
    notes = [t for t in note_texts(diagnostics) if "frequency" in t.lower()]
    assert len(notes) == 1, f"expected one frequency note, got {note_texts(diagnostics)}"
    note = notes[0]

    assert f"{removed} variant" in note, f"the removed count is not {removed}: {note}"
    assert "exempt" in note.lower(), f"the exemption is not mentioned: {note}"
    assert f"{exempted} " in note, f"the exempted count is not {exempted}: {note}"


def test_the_frequency_layer_never_contradicts_the_reason_column(somatic):
    """A row the layer removes must read as rejected, not as having met the criteria.

    This is the composition's one real hazard, and #37 sharpened it: the mask now gates the
    union, and a mask applied to the union *after* the four cells were derived would produce
    rows marked NOPASS whose reason still claimed they passed on criteria — a report that
    contradicts itself, which is precisely what this seam exists to prevent. It is instead
    applied to both deciding masks first, which is the same mask on the verdict and leaves
    the cells a partition.
    """
    tight, diagnostics = apply_filters(
        somatic, params("somatic", max_freq_population=0.01)
    )

    passed = tight[MAFIGATE_FILTER] == PASS
    rejected = tight[MAFIGATE_REASON] == REASON_REJECTED
    assert int((passed & rejected).sum()) == 0
    assert int((~passed & ~rejected).sum()) == 0
    assert sum(diagnostics.cells().values()) == len(tight)


def test_a_missing_frequency_passes_rather_than_dropping_the_row(somatic):
    """Missing means *not seen in the panel*, not *unmeasurable*.

    The opposite of the rule for depth and VAF, deliberately: there a missing value means
    we cannot tell whether the call is sound, so the pipeline drops the row. Missing
    frequency rates run 13–35% per column on the reference, and dropping them would invert
    the filter and discard rare pathogenic variants precisely *because* they are rare.

    Asserted on the fixture's 25 rows that carry no frequency at all: at a threshold so
    tight that nothing measurable survives, they must all still be judged on their other
    criteria rather than removed by this layer.

    Since #122 this is a claim about rows **no** column measured, not about rows with *a*
    blank, and it is the whole reason the mask carries its ``~measured`` fallback: with the
    blank rescue gone, this is the only thing standing between an unmeasured variant and
    the floor. The rows with a blank *and* a populated value are the next test's.
    """
    present = [c for c in FREQUENCY_COLUMNS if c in somatic.columns]
    assert present, "fixture carries no frequency columns — this test is vacuous"

    numeric = pd.concat(
        [pd.to_numeric(somatic[c], errors="coerce") for c in present], axis=1
    )
    all_missing = numeric.isna().all(axis=1)
    assert bool(all_missing.any()), "fixture has no all-missing frequency rows"

    unfiltered, _ = apply_filters(somatic, params("somatic", max_freq_population=1.0))
    floored, _ = apply_filters(
        somatic, params("somatic", max_freq_population=0.0000001)
    )

    # For the all-missing rows the layer is a no-op, so their verdicts are untouched.
    assert list(floored.loc[all_missing, MAFIGATE_FILTER]) == list(
        unfiltered.loc[all_missing, MAFIGATE_FILTER]
    ), "a row with no frequency data was removed by the frequency layer"


def _frequency_frame(**columns) -> pd.DataFrame:
    """A frame carrying only the named frequency columns.

    ``frequency_mask`` reads nothing else, so anything more would be columns whose values
    cannot affect the answer. Written out rather than taken from a fixture because the
    cases below are *combinations* — a blank beside a common value, a blank beside a rare
    one, the same column absent instead of empty — and pinning each to whichever reference
    row happens to hold it today makes the test a hostage to the fixture rather than a
    statement of the rule.
    """
    return pd.DataFrame(columns)


def test_a_blank_does_not_rescue_a_variant_the_populated_columns_call_common():
    """Issue #122's change, stated at the smallest scale that can express it.

    Before #122 the mask read ``(v <= threshold) | v.isna()``, which counts a blank as
    zero, so one empty column carried a row whatever the populated ones said. Measured
    through the whole filter over 171 real MAFs that was holding 4,167 rows across 534
    distinct variants in Broad reports — ``ATM`` at gnomAD 1.0, ``PTEN`` at 1000G 1.0 —
    each of them in the report only because some other column was empty.

    Row 2 is the control that keeps the assertion from being satisfiable by a mask that
    simply drops everything with a blank in it.
    """
    frame = _frequency_frame(
        gnomAD_exome_AF=[float("nan"), 0.44, float("nan")],
        Freq_1000g2015aug_all=[0.44, float("nan"), float("nan")],
    )

    assert list(frequency_mask(frame, 0.01)) == [False, False, True], (
        "a blank in one column still rescues a variant another column calls common"
    )


def test_a_blank_still_cannot_sink_a_variant_another_column_calls_rare():
    """The half of the missing-value rule #122 did **not** change.

    #28 and #37 chose it because a panel that never saw a variant is not evidence the
    variant is common, and missing rates run 13–35% per column. That reasoning survives
    #122 untouched: a blank contributes nothing, so it cannot sink a row either, and the
    minimum over the columns that *did* measure the variant still decides.

    Row 1 is the case that separates this rule from "every column must agree" — the option
    #122 rejected, under which a single disagreeing panel vetoes. It matters because
    ``FREQUENCY_COLUMNS`` includes the pre-QC ``_raw`` columns, which #126 measured as
    calling a variant common that the adjusted column does not on 1,052 germline and 1,549
    somatic rows at 1%: under a veto every one of those would be dropped. ``MAX_AF``, a
    maximum across subpopulations, was the other name this argument cited, and #126 removed
    it from the list — a veto on it would have turned a 1% cut into "no subpopulation
    above 1%", and it was never a panel's opinion in the first place.
    """
    frame = _frequency_frame(
        gnomAD_exome_AF=[float("nan"), 0.44],
        Freq_1000g2015aug_all=[0.001, 0.001],
    )

    assert list(frequency_mask(frame, 0.01)) == [True, True], (
        "a blank, or a single disagreeing panel, sank a variant another column calls rare"
    )


def test_an_absent_column_and_an_empty_one_reach_the_same_verdict():
    """The asymmetry #122 dissolved rather than fixed.

    Under the old rule a column present-but-empty could rescue a row while the same column
    *missing* could not, because a missing column is not in the comprehension at all and so
    has no blank to contribute. Nothing ever argued for that, and it is what decided
    ``HMGCR``/``NOS3`` one way and ``FUT2`` the other. Now neither contributes, so the two
    agree without a second mechanism to hold them together — which is why this is asserted
    as an equality between the two frames rather than as a value for each.
    """
    empty = _frequency_frame(
        gnomAD_exome_AF=[float("nan")], Freq_1000g2015aug_all=[0.44]
    )
    absent = _frequency_frame(Freq_1000g2015aug_all=[0.44])

    assert list(frequency_mask(empty, 0.01)) == list(frequency_mask(absent, 0.01)), (
        "an empty column and an absent one still disagree about the same variant"
    )
    assert list(frequency_mask(empty, 0.01)) == [False]


def test_the_maf_blank_spelling_is_read_as_a_blank_and_not_as_a_number():
    """``.`` is MAF's own blank, and it is the exact shape #122 was opened about.

    Every file carrying ``HMGCR`` or ``NOS3`` spells ``Freq_esp6500siv2_all`` as ``.``, so
    a mask that read that token as a value — or that failed to coerce it and treated the
    whole column as unmeasured — would answer this ticket's own case wrongly while passing
    every test written with ``NaN``. Row 0 is ``NOS3``'s shape at its real frequency.
    """
    frame = _frequency_frame(
        Freq_esp6500siv2_all=[".", "."],
        Freq_1000g2015aug_all=["0.765575", "0.001"],
    )

    assert list(frequency_mask(frame, 0.01)) == [False, True]


def test_the_mask_is_all_true_at_1_0_including_rows_no_panel_measured():
    """Neutrality at 1.0 is what the composition rests on, and #122 rebuilt the algebra.

    ``tests/parity/test_parity.py`` asserts this over the fixtures; the case it cannot
    reach is the new ``~measured`` branch, which is a *different* term in the expression
    rather than the same comparison on other data. A mask that returned the populated
    comparison alone would still be all-True on every measured row and would silently drop
    every unmeasured one at 1.0 — costing parity, which is the failure the guard in
    ``apply_filters`` deliberately does not hide.
    """
    frame = _frequency_frame(
        gnomAD_exome_AF=[float("nan"), 1.0, 0.0],
        Freq_1000g2015aug_all=[float("nan"), float("nan"), 1.0],
    )

    assert bool(frequency_mask(frame, 1.0).all())


def test_a_pre_qc_value_cannot_veto_a_variant_the_adjusted_column_calls_rare():
    """Issue #126's decision, in the direction that makes keeping ``_raw`` safe.

    The two ``gnomAD_*_AF_raw`` columns are frequencies over every genotype gnomAD called,
    including the ones its own quality filters reject, so they run higher than the adjusted
    value at the cut: measured over 164 real MAFs, ``_raw`` calls a variant common that the
    adjusted column does not on **1,052 germline and 1,549 somatic rows** at 1%, against 35
    and 32 rows the other way. Under the minimum rule every one of those is kept, which is
    the whole reason the pre-QC column can sit in the list beside the adjusted one.

    Row 0 is that case. Row 1 is the control: with the adjusted column silent, the pre-QC
    value is the only gnomAD evidence there is and it decides — which is what a MAF that has
    been through the pipeline looks like, ``bin/filter_variants.py`` keeping
    ``gnomAD_exome_AF_raw`` in preference to ``gnomAD_exome_AF``. Dropping the ``_raw`` names
    from the list fails this second assertion, and that is the trade #126 measured: 435 Broad
    and 91 Stringent germline rows, on files where nothing else speaks for gnomAD at all.
    """
    frame = _frequency_frame(
        gnomAD_exome_AF=[0.001, float("nan")],
        gnomAD_exome_AF_raw=[0.44, 0.44],
    )

    assert list(frequency_mask(frame, 0.01)) == [True, False], (
        "a pre-QC value either vetoed a variant the adjusted column calls rare, or stopped "
        "deciding a variant no other column measured"
    )


def test_max_af_is_not_a_measurement_this_filter_reads():
    """Issue #126 removed it, and this is the removal stated as behaviour.

    ``MAX_AF`` is the highest frequency a variant reaches in any *one* population, so it is
    not a panel's opinion the way the other seven are — it is already a maximum, and under
    the minimum rule it could only ever rescue a row, never remove one, except in exactly the
    case below where it is the only populated column. #122 refused the max rule partly
    because its presence would have turned a 1% cut into "no subpopulation above 1%".

    It is also absent from **all 164** real MAFs measured, so nothing it could have done had
    ever been exercised, and removing it moved 0 report rows on either arm at either preset.
    A frame carrying it alone must therefore be a frame with nothing to judge by — the
    ``None`` that makes ``apply_filters`` say so rather than filter on it.

    Stated as that behaviour and **not** as ``"MAX_AF" not in FREQUENCY_COLUMNS``, which is a
    restatement of the constant it is checking and would pass over a filter that read the
    column by some other route. Row 1 is below the threshold, so a mask built over this frame
    would be ``[False, True]`` — the failure this asserts against is a real verdict, not an
    empty one.
    """
    frame = _frequency_frame(MAX_AF=[0.44, 0.001])

    assert frequency_mask(frame, 0.01) is None, (
        "MAX_AF is being read as a population-frequency measurement again"
    )


def test_the_load_page_and_the_filter_agree_what_a_frequency_column_is():
    """One list, because they make the same claim about the same file (issue #126).

    ``page_modules/data_loading.py`` warns *No Population Frequency Columns Found* and
    ``frequency_mask`` returns ``None``; before #126 the warning was driven by a hand-written
    second copy of these eight names in ``config/columns.py`` and the mask by this list, so
    the two could drift and one of them be wrong about the file on screen. Removing
    ``MAX_AF`` from one copy and not the other would have done exactly that.

    Driven through :func:`page_modules.data_loading.validate_required_columns` rather than
    asserted over the source, because what has to hold is that the *warning* follows this
    list. A grep for the constant's name passes on a file that imports it and then reads
    something else, which is the shape the weaker version of this test had; nothing else in
    the suite calls this function at all, so there was no existing guard to lean on.

    ``MAX_AF`` is the case that separates the two lists: it was in the copy and is not in
    this list, so a load page still reading the copy stays silent about a MAF the filter has
    nothing to judge by. The control is a MAF carrying one real frequency column, where the
    sentence must not be said.
    """
    from config.columns import OPTIONAL_COLUMNS

    assert "population_frequency" not in OPTIONAL_COLUMNS, (
        "a second list of the frequency columns is back in config/columns.py; "
        "data_loading.py reads FREQUENCY_COLUMNS, and two lists can disagree"
    )

    def warnings_for(frame):
        said: list[str] = []
        fake = SimpleNamespace(
            warning=lambda text: said.append(str(text)),
            error=lambda text: said.append(str(text)),
            session_state=SimpleNamespace(
                filter_params={"sample_type": "somatic"}, get=lambda *a: None
            ),
        )
        with patch.object(data_loading, "st", fake):
            data_loading.validate_required_columns(frame)
        return said

    complete = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))
    stripped = complete.drop(
        columns=[c for c in FREQUENCY_COLUMNS if c in complete.columns]
    )

    only_max_af = stripped.copy()
    only_max_af["MAX_AF"] = 0.44
    assert any("Population Frequency" in w for w in warnings_for(only_max_af)), (
        "a MAF whose only frequency-ish column is MAX_AF was not told the frequency filter "
        "has nothing to read — the load page and the filter disagree about the file"
    )

    one_real_column = stripped.copy()
    one_real_column[FREQUENCY_COLUMNS[0]] = 0.001
    assert not any("Population Frequency" in w for w in warnings_for(one_real_column)), (
        f"the load page claims there are no frequency columns over a MAF carrying "
        f"{FREQUENCY_COLUMNS[0]}, which frequency_mask does read"
    )


def test_a_maf_with_no_frequency_columns_is_reported_not_filtered(somatic):
    """No columns to judge by is a fact to state, not a reason to drop everything."""
    stripped = somatic.drop(
        columns=[c for c in FREQUENCY_COLUMNS if c in somatic.columns]
    )
    with_layer, diagnostics = apply_filters(
        stripped, params("somatic", max_freq_population=0.01)
    )
    without_layer, _ = apply_filters(stripped, params("somatic", max_freq_population=1.0))

    assert list(with_layer[MAFIGATE_FILTER]) == list(without_layer[MAFIGATE_FILTER])
    assert any("frequency columns" in t for t in note_texts(diagnostics)), (
        f"the user was not told the filter could not be applied: {note_texts(diagnostics)}"
    )


# ---------------------------------------------------------------------------
# 2. The transitional parameter fallbacks (deleted when #40 lands)
# ---------------------------------------------------------------------------


def test_the_legacy_keep_pathogenic_flag_still_inverts_correctly(somatic):
    """``keep_pathogenic=False`` must mean ``skip_pathogenic=True``.

    Getting this polarity wrong is worth +2 rows on each arm, silently — and silently
    *wider*, which is the invisible direction. The live presets all still carry
    ``keep_pathogenic``, so this is the path the running app takes today.
    """
    legacy = params("somatic")
    del legacy["skip_pathogenic"]
    legacy["keep_pathogenic"] = False

    contract = params("somatic", skip_pathogenic=True)

    legacy_out, legacy_diag = apply_filters(somatic, legacy)
    contract_out, contract_diag = apply_filters(somatic, contract)

    assert legacy_diag.cells() == contract_diag.cells()
    assert list(legacy_out[MAFIGATE_FILTER]) == list(contract_out[MAFIGATE_FILTER])
    assert legacy_diag.rescue_only == 0 and legacy_diag.both == 0


def test_the_legacy_germline_vaf_key_wins_on_the_germline_arm(germline):
    """A legacy dict carrying *both* VAF keys must use the germline one.

    ``vaf_threshold`` holds the somatic value in such a dict, so preferring it would apply
    0.01 where 0.2 was meant and widen the germline report. Three sources once disagreed
    here — the config and both presets said 0.2 while a widget fallback said 0.3.
    """
    both_keys = params("germline", vaf_threshold=0.01)
    both_keys["vaf_threshold_germline"] = 0.2

    contract = params("germline", vaf_threshold=0.2)

    from_legacy, legacy_diag = apply_filters(germline, both_keys)
    from_contract, contract_diag = apply_filters(germline, contract)

    assert legacy_diag.cells() == contract_diag.cells(), (
        "the germline arm did not prefer vaf_threshold_germline — it is filtering at the "
        "somatic threshold"
    )
    assert list(from_legacy[MAFIGATE_FILTER]) == list(from_contract[MAFIGATE_FILTER])


@pytest.mark.parametrize(
    "arm,legacy_key,fixture_name",
    [
        ("somatic", "somatic_genes", "somatic"),
        ("germline", "germline_genes", "germline"),
    ],
)
def test_the_legacy_per_arm_gene_string_is_still_read(
    arm, legacy_key, fixture_name, request
):
    """The old comma-separated string form must reach the vendored gene clause.

    The contract's one unprefixed ``filter_genes`` holds a list; the presets and saved
    files hold a per-arm comma string. Both have to arrive at the same gene filter, or the
    running app silently stops filtering by gene — which widens, invisibly.

    Only the *comma* form is asserted, deliberately: making a one-symbol-per-line paste
    work is #34's, and asserting it here would claim a fix this ticket did not make.
    """
    frame = request.getfixturevalue(fixture_name)
    symbols = sorted(set(frame["Hugo_Symbol"].dropna().astype(str)))[:3]
    assert len(symbols) == 3, "fixture has too few gene symbols for this test"

    legacy = params(arm)
    del legacy["filter_genes"]
    legacy[legacy_key] = ",".join(symbols)

    from_legacy, legacy_diag = apply_filters(frame, legacy)
    from_list, list_diag = apply_filters(frame, params(arm, filter_genes=symbols))

    assert legacy_diag.cells() == list_diag.cells()
    assert list(from_legacy[MAFIGATE_FILTER]) == list(from_list[MAFIGATE_FILTER])


def test_a_gene_list_restricts_the_criteria_path_case_insensitively(somatic):
    """Gene matching is case-insensitive, and it got that way by routing.

    The app used to compare symbols verbatim, so a lowercase paste silently emptied the
    report. Nothing in the app upper-cases anything now — the vendored clause upper-cases
    both sides — so this asserts the consequence of the routing rather than a fix of the
    app's own.

    Uses ``skip_pathogenic`` because otherwise the pathogenic rescue re-admits exactly the
    rows the gene list excluded, and the case bug could be present with the totals equal.
    """
    symbols = sorted(set(somatic["Hugo_Symbol"].dropna().astype(str)))[:3]

    def criteria_path(gene_values):
        _, diagnostics = apply_filters(
            somatic,
            params("somatic", filter_genes=gene_values, skip_pathogenic=True),
        )
        return diagnostics.criteria_path

    unrestricted = apply_filters(
        somatic, params("somatic", skip_pathogenic=True)
    )[1].criteria_path

    upper = criteria_path([s.upper() for s in symbols])
    lower = criteria_path([s.lower() for s in symbols])

    assert upper == lower, (
        f"case changed the result ({upper} vs {lower}) — gene matching is still "
        "case-sensitive somewhere"
    )
    assert upper < unrestricted, (
        f"restricting to {symbols} did not narrow the criteria path ({upper} vs "
        f"{unrestricted}) — the gene list never reached the vendored clause"
    )


# ---------------------------------------------------------------------------
# 3. When the arm does not match the MAF
# ---------------------------------------------------------------------------


def test_a_mismatched_arm_reports_loudly_rather_than_filtering_on_nothing(somatic):
    """Germline parameters against a somatic MAF must never produce a quiet report.

    This is a real state the running app reaches: ``MAFigate.py`` seeds ``filter_params``
    from ``~/.mafigate/cached_parameters.json``, so a user whose last session was germline
    opens a somatic MAF on the germline arm. #40 stopped the cache reaching session state;
    #39 turned the refusal this used to assert into a report plus an escalated warning, which
    is the direction this test was written expecting.

    What is being pinned is unchanged across that switch, because it was never the exception:
    the app's own germline filter guarded every source with ``if column in maf.columns``, so
    with all three germline sources absent it skipped the entire guideline block and reported
    a confident, *silent*, wrong cut. That is the failure. A report the user is told to
    distrust is not it — so the assertion is on the telling, not on the raising.
    """
    assert "InterVar" not in somatic.columns, (
        "the somatic fixture grew an InterVar column — pick another absent germline "
        "source for this test"
    )
    labelled, diagnostics = apply_filters(somatic, params("germline"))

    assert len(labelled) == len(somatic)
    assert "InterVar" in diagnostics.filled_columns
    assert "InterVar" in diagnostics.degraded_columns
    assert any(
        "PATHOGENIC RETENTION" in warning and "InterVar" in warning
        for warning in note_texts(diagnostics)
    )


def test_an_unknown_arm_is_rejected_by_name():
    """A typo in ``sample_type`` must not silently select an arm."""
    with pytest.raises(ValueError, match="somatic or germline"):
        apply_filters(pd.DataFrame({"Hugo_Symbol": ["BRCA1"]}), {"sample_type": "both"})
