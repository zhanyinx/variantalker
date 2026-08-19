"""``Clinical_Summary`` against the pipeline's own answer to the same question.

The annotation workflow runs ``bin/generate_clinical_summary.py`` as a process of its own,
downstream of RENOVO and upstream of ``filter_maf``, and writes ``variantalker_naive``.
``compute_keep`` keeps that column on both arms, so it is in every filtered frame — beside
``Clinical_Summary``, which ``components/clinical_summary.py`` derives by re-implementing that
same script.

Two columns on one screen answering one question is only safe if they agree, and until issue
#108 they did not: the app read ``filter_renovo``, a filter *parameter* name, where the pipeline
reads ``RENOVO_Class``, so RENOVO reached the app's column on no row ever rendered while
reaching the pipeline's on every one. That is what this module holds shut, and it can hold it
**exactly** rather than approximately, because the oracle is already in git: the committed
``germline_reference.maf`` carries the pipeline's column for all 94 of its rows.

**Issue #117 then took that column off the screen, and this module is deliberately unchanged
by it.** The report shows one verdict now — ``config.columns.PIPELINE_COLUMNS_THE_APP_REPLACES``
keeps the pipeline's copy out of the default view, on the measurement that the two agree on
every somatic row and on 51 of 52 germline MAFs. What that changes is which columns the *grid*
opens with; the column is still in the frame and still in this fixture, so the comparison below
reads it exactly as before. A check against the pipeline's own answer must not depend on the
app choosing to display that answer — if it did, hiding a column would silently retire the
guard that says the app agrees with it.

**This is a net, not a harness.** It has a pipeline counterpart and must hold where ``bin/`` is
absent — the fixture carries the pipeline's output, so nothing here runs the pipeline.

Five things this module deliberately does *not* do. (It said three while listing four: the
count was not updated when issue #109 added the CIViC divergence. The fifth entry arrived
with #285's third divergence, and stays now that issue #296 has closed that divergence —
what it records is no longer a disagreement but a limit of this fixture, and the limit is
what the test below leans on.)

It does not compare the two columns on the non-pathogenicity classes. Since issue #98 the app
classifies ClinVar's non-standard terms by the institute's term table, giving them **six classes
of their own** below the ladder where the pipeline pools every one of them into two buckets,
``Not_Provided`` and ``Unknown`` — so disagreement there is a decision this map took, recorded at
:data:`CLINICAL_VALUE_MAPPING`, not drift. The row filter below is written so that a row it
excludes cannot pass unnoticed: it asserts the excluded set is *empty* on this fixture, so the
day one appears, this module fails and someone decides.

It also does not see the **second** deliberate divergence, and that is a limit rather than a
choice: issue #109 reads CIViC's ``E`` as ``Uncertain_Significance`` where the pipeline tiers it
*Benign*, and neither committed reference MAF carries a CIViC column, so no fixture here can
witness it either way. It is recorded in the module docstring under test instead. A MAF with
CIViC annotation would make the two columns disagree on rows whose strongest CIViC level is
``E``, and that disagreement is correct.

And it cannot witness ``Pathogenic/Pathogenic,_low_penetrance|other`` either way, which was
the **third** divergence until issue #296 closed it. Both reference MAFs carry that cell —
germline row 56 and somatic row 46 — and issue #285 taught the app to read it as
``Pathogenic`` where ``bin/generate_clinical_summary.py`` split on the first separator only,
mapped the cell to tier 7 (``Unknown``, its *worst* tier) and let ``min()`` defer to whatever
else the row carried: ``Benign`` on these two rows, ``Uncertain_Significance`` on all 8 rows
of the real corpus that hold the cell. **#296 ported the app's rule into ``bin/``**, so the
two now agree here and the fixture's stored ``Pathogenic`` is what both would produce.

What has not changed is that this fixture cannot *hold* the two to it. The constructed set's
verdict column is written by *calling the app*
(``build_parity_fixtures.fill_pipeline_verdict``), so on this row it tracks the app rather
than witnessing it — it read ``Pathogenic`` before ``bin/`` did, and would go on reading it
if ``bin/`` regressed. That is why the agreement is recorded here, and measured against the
real corpus on issue #295, rather than waited for below. MAFs annotated before #296 keep the
old reading; per issue #117 the app's column is the current one and a file's own is an
earlier vintage.

It does not read ``bin/`` to derive the expected tiers. ``test_vendor_drift.py`` owns
app-versus-pipeline *code* equality for the filter, and this function is not vendored — it is a
re-implementation, deliberately divergent on ClinVar. What is checkable is the *output*, which
is what the fixture holds.

And it does not assert the six RENOVO tiers by restating them. It asserts each is load-bearing:
that RENOVO decides the label on rows of the real fixture, and that the two low-confidence
classes either side of RENOVO's 0.5 boundary land on one class.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NamedTuple

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from components import clinical_summary  # noqa: E402
from components.clinical_summary import (  # noqa: E402
    CLINICAL_CIRCLE_SOURCES,
    CLINICAL_HIERARCHY,
    CLINICAL_VALUE_MAPPING,
    NO_CLINICAL_DATA,
    _CLASS_GLYPHS,
    _SUMMARY_LABELS,
    _SUMMARY_SOURCE_COLUMNS,
    generate_clinical_summary,
    pathogenicity_circle_glyphs,
)
from components.clinical_badges import CLINICAL_ROW  # noqa: E402
from components.variant_detail import _POP_FREQ_COLUMNS  # noqa: E402
from config.pipeline_params import pipeline_params  # noqa: E402
from config.vocabularies import RENOVO_OPTIONS  # noqa: E402
from vendor.pipeline_utils import read_maf  # noqa: E402

FIXTURE = STREAMLIT_APP / "tests" / "fixtures" / "parity" / "germline_reference.maf"

#: The pipeline's column, and the app's RENOVO column.
PIPELINE_VERDICT = "variantalker_naive"
RENOVO_COLUMN = "RENOVO_Class"

#: The seven names ``bin/generate_clinical_summary.py``'s ``NUMERIC_TO_CLINICAL`` can write.
#: The five that are points on the ACMG ladder are the ones the two columns must agree on;
#: ``Not_Provided`` and ``Unknown`` are the pipeline's two pooled buckets, which issue #98
#: replaced with six classes from the institute's term table and so may differ.
PIPELINE_LADDER = ("Pathogenic", "Likely_Pathogenic", "Uncertain_Significance", "Likely_Benign", "Benign")
PIPELINE_POOLED = ("Not_Provided", "Unknown")

class RenovoTier(NamedTuple):
    """One of RENOVO's six classes, and everything this module asserts about it.

    A named tuple rather than a bare 5-tuple because most tests here want two of the five
    fields, and ``for value, _lo, hi, _c, _g in ...`` names three things it does not use.
    """

    value: str
    score_from: float
    score_to: float
    clinical_class: str
    circle: str


#: RENOVO's six classes, strongest first: the interval of RENOVO's own ``PL_score`` each
#: occupies, the class the app reads it as, and the circle that class draws.
#:
#: The score intervals are measured over the **211,744 rows carrying both a class and a numeric
#: score**, in the 64 MAFs on this machine carrying ``RENOVO_Class``. The score is in the MAF
#: beside the class — ``bin/add_renovo_to_maf.py`` copies RENOVO's ``PL_score`` in as
#: ``RENOVO_pls`` — which is what makes the confidence prefix a measurable thing rather than a
#: reading of three letters, and what lets
#: :func:`test_the_score_intervals_hold_on_the_committed_fixture` check these numbers instead of
#: comparing them with themselves.
#:
#: **The class names are written out here on purpose, and that is not the transcription this
#: repo keeps deleting.** They are the *subject* of this module, not its vocabulary: a guard
#: that reads the expected class out of :data:`CLINICAL_VALUE_MAPPING` asserts the mapping
#: equals itself, which is #77's ``deepcopy(X) == X`` and #100's legend-compared-with-itself in
#: a third place. The first draft of this file did exactly that, and the mutation that exposed
#: it is the one the committed fixture cannot see: ``IP Pathogenic`` occurs on **no row of
#: germline_reference.maf** — 41 real MAFs and 300 corpus rows carry it — so re-grading that one
#: value back to ``Pathogenic`` passed every derived assertion here. Six literals is the price
#: of the mutation being caught.
#:
#: Ordered strongest first so the *ordering* is checkable as a property as well: the class a
#: value maps to may never be stronger than the class of a value scoring higher.
#: Bounds are given to five decimals and rounded **outward** — lower bounds down, upper bounds
#: up — so each interval contains every score measured in its class. Five rather than four
#: because four is not enough to keep the bins apart: ``HP Benign`` tops out at 0.0091836 and
#: ``IP Benign`` starts at 0.0092063, which both round to 0.0092, and the disjointness check
#: below failed on that the moment it was written. Do not "tidy" these to four.
RENOVO_TIERS = (
    RenovoTier("HP Pathogenic", 0.88928, 1.00000, "Pathogenic", "🔴"),
    RenovoTier("IP Pathogenic", 0.78540, 0.88858, "Likely_Pathogenic", "🟠"),
    RenovoTier("LP Pathogenic", 0.50000, 0.78453, "Uncertain_Significance", "🟡"),
    RenovoTier("LP Benign", 0.23501, 0.49985, "Uncertain_Significance", "🟡"),
    RenovoTier("IP Benign", 0.00920, 0.23494, "Likely_Benign", "🟢"),
    RenovoTier("HP Benign", 0.00000, 0.00919, "Benign", "🔵"),
)

#: Just the values, for parametrising. Derived from the table so the two cannot part company —
#: a parametrised guard whose table has quietly stopped covering a case is a guard that has
#: quietly stopped applying, which issue #103 found had already happened once here.
RENOVO_VALUES = [tier.value for tier in RENOVO_TIERS]


@pytest.fixture(scope="module")
def germline() -> pd.DataFrame:
    return read_maf(str(FIXTURE))


def _labels(frame: pd.DataFrame) -> list[str]:
    """The summary label of every row of ``frame``, in row order."""
    return [generate_clinical_summary(frame.loc[i]) for i in frame.index]


# --- the premise: this fixture can see the defect ----------------------------------------


def test_the_fixture_carries_the_pipelines_answer_and_not_the_parameter_name(germline):
    """The oracle is present, and the name the app used to read is not.

    Both halves matter. Without ``variantalker_naive`` there is nothing to compare against;
    and if a MAF *did* carry a ``filter_renovo`` column, the old code would have worked by
    accident and this module would be asserting something other than what it claims.
    """
    assert PIPELINE_VERDICT in germline.columns
    assert RENOVO_COLUMN in germline.columns
    assert "filter_renovo" not in germline.columns


def test_renovo_decides_the_label_on_enough_rows_to_be_load_bearing(germline):
    """Dropping RENOVO from the read must move the label on a substantial share of rows.

    Without this, :func:`test_the_summary_matches_the_pipelines_own_column` could pass on a
    fixture where RENOVO never mattered — green, and blind to the whole of issue #108. It is
    the same fault as a guard whose parametrising table has quietly stopped covering anything.

    Measured at 23 of 94 rows on this fixture. Asserted as a floor rather than the exact
    number, because the claim is *this instrument has purchase*, and pinning 23 would make a
    fixture edit look like a regression.
    """
    with_renovo = _labels(germline)
    without_renovo = _labels(germline.drop(columns=[RENOVO_COLUMN]))
    moved = sum(a != b for a, b in zip(with_renovo, without_renovo))
    assert moved >= 20, f"RENOVO decides only {moved} of {len(germline)} labels"


# --- the claim ----------------------------------------------------------------------------


def test_the_summary_matches_the_pipelines_own_column(germline):
    """Row for row, the app's label is the pipeline's verdict.

    The exact check issue #108 was resolved by. Before it, 71 of these 94 rows agreed; after
    it, all 94 do — and on the wider corpus, 99,210 of 99,210 comparable rows across the 52
    real germline MAFs whose ``variantalker_naive`` was written by a RENOVO-aware pipeline.

    **"All 94" is 94 of 94 on this fixture and not a claim about real files.** Row 56 carries
    the ClinVar cell the app and ``bin/`` differed on between issues #285 and #296; they agree
    on it now, but this test would have stayed green either way, because the constructed set's
    verdict column is written by calling the app. The module docstring's last entry says so —
    read it before treating a pass here as agreement with the pipeline.
    """
    verdicts = germline[PIPELINE_VERDICT].astype(str).str.strip()

    # No row may fall outside the comparison unnoticed — see the module docstring.
    unexpected = sorted(set(verdicts) - set(PIPELINE_LADDER) - set(PIPELINE_POOLED))
    assert not unexpected, f"{PIPELINE_VERDICT} holds names this module cannot place: {unexpected}"
    pooled = verdicts.isin(PIPELINE_POOLED)
    assert not pooled.any(), (
        f"{pooled.sum()} rows carry one of the pipeline's pooled buckets, where issue #98's "
        "classes may legitimately differ — decide what this module should say about them"
    )

    # The pipeline's five ladder names *are* `CLINICAL_HIERARCHY` names, so the translation to an
    # app label is a lookup in the app's own label table rather than a second table of five emoji
    # strings written out here — which is how the drift this module exists to catch would re-enter
    # through its own oracle.
    mismatches = [
        (i, got, _SUMMARY_LABELS[verdicts[i]])
        for i, got in zip(germline.index, _labels(germline))
        if got != _SUMMARY_LABELS[verdicts[i]]
    ]
    assert not mismatches, f"{len(mismatches)} of {len(germline)} rows disagree: {mismatches[:5]}"


def test_the_summary_reads_renovo_by_the_name_the_pipeline_emits():
    """The column set names ``RENOVO_Class``, and no member of it is a filter parameter.

    The second half is issue #108's fourth question made structural. The defect was invisible
    to every kind of check the app had, because the name was *plausible* and the function
    returns a label either way — nothing raised, nothing warned, and one entry of a five-entry
    table simply never matched. What distinguishes it is the category error: a parameter name
    among column names. That is checkable, and here it is checked.
    """
    assert RENOVO_COLUMN in _SUMMARY_SOURCE_COLUMNS

    parameters = set(pipeline_params("germline")) | set(pipeline_params("somatic"))
    named_parameters = sorted(set(_SUMMARY_SOURCE_COLUMNS) & parameters)
    assert not named_parameters, (
        f"{named_parameters} are filter parameter names, not MAF columns — the shape of the "
        "defect issue #108 fixed"
    )


#: Every table in the app that pairs a MAF column with how to draw it, walked per row.
#:
#: Named rather than discovered, because the *syntax* does not give them away — which is the
#: whole reason this defect survived three times. A sweep for column-reading shapes
#: (``row["X"]``, ``"X" in row``) finds none of them: the column name sits in a collection, and
#: the read happens in a loop over that collection somewhere else entirely. An AST sweep for
#: collections-of-column-names misses ``CLINICAL_ROW`` too, because each of its entries is a
#: 3-tuple holding one column name, one label and one colour, so no single tuple looks like a
#: table of columns. (It was ``variant_detail._CLINICAL_BADGE_CONFIG`` until issue #204 moved the
#: row's membership into :mod:`components.clinical_badges` beside its vocabularies.)
#:
#: Adding a table here is a deliberate act. The cost of forgetting is on the record: three
#: surfaces held ``filter_renovo`` — the circles (issue #95), ``Clinical_Summary`` (#108) and the
#: variant-detail badges (#108, found by review) — and each was invisible to the others.
COLUMN_RENDER_TABLES = {
    "components.clinical_summary._SUMMARY_SOURCE_COLUMNS": lambda: list(_SUMMARY_SOURCE_COLUMNS),
    "components.clinical_summary.CLINICAL_CIRCLE_SOURCES": lambda: [
        entry[2] for entry in CLINICAL_CIRCLE_SOURCES
    ],
    "components.clinical_badges.CLINICAL_ROW": lambda: [
        entry[0] for entry in CLINICAL_ROW
    ],
    "components.variant_detail._POP_FREQ_COLUMNS": lambda: [entry[0] for entry in _POP_FREQ_COLUMNS],
}


@pytest.mark.parametrize("table", sorted(COLUMN_RENDER_TABLES))
def test_no_render_table_names_a_filter_parameter(table):
    """Issue #108's fourth question, asked of every such table instead of one.

    The ticket asked *whether anything else names a parameter where it means a column*, and the
    answer was yes — one more, in the variant-detail badges, where the RENOVO badge had therefore
    never been drawn. So the answer is a guard rather than a sweep I ran once: a filter parameter
    name is a category error wherever it appears in a table of column names, and the sixteen
    parameters are enumerable from the contract.

    Deliberately **not** also asserting that each column appears in a reference fixture. That
    rule was drafted and measured first: it fails on ``CIViC_Evidence_Level`` and
    ``gnomAD_genome_AF``, which neither committed reference MAF carries and which are both
    perfectly real columns — 34 of 103 real somatic MAFs carry the CIViC block. A guard needing
    two standing exceptions on the day it is written is a guard that will acquire a third
    silently, and absence is measured against what the pipeline emits, not against a list of
    names (issue #90).
    """
    parameters = set(pipeline_params("germline")) | set(pipeline_params("somatic"))
    columns = COLUMN_RENDER_TABLES[table]()
    named_parameters = sorted(set(columns) & parameters)
    assert not named_parameters, (
        f"{table} names {named_parameters}, which are filter parameters and not MAF columns; a "
        "table walked against row.index will never match one, and nothing raises"
    )


# --- the six RENOVO tiers -----------------------------------------------------------------


def test_the_control_offers_exactly_the_classes_the_report_can_name():
    """Every RENOVO term the filter offers is nameable, and nothing else is mapped.

    Issue #103's rule — *the control offers exactly what the report can name* — applied to the
    one source it was not applied to. Guarded in both directions, so neither list can grow a
    member the other lacks: a term offered and unnameable would report a variant as carrying an
    unrecognised annotation, and a term nameable and unofferable is a keep-list option missing
    from a vocabulary.
    """
    assert set(RENOVO_OPTIONS) == set(RENOVO_VALUES)
    for value in RENOVO_OPTIONS:
        assert value in CLINICAL_VALUE_MAPPING, f"{value} is offered but cannot be named"


@pytest.mark.parametrize("tier", RENOVO_TIERS, ids=RENOVO_VALUES)
def test_a_renovo_call_on_its_own_names_its_tier(tier):
    """A variant whose only annotation is a RENOVO class gets that class, not a shrug.

    Before issue #108 every one of these six returned ``🔍 No Clinical Data`` — said of a
    variant whose only distinguishing feature is the RENOVO prediction it carries. On the
    germline arm that sentence covered 2,033 of the 218,349 rows in the 64 MAFs carrying
    ``RENOVO_Class``, and now covers 1.

    The expected class is the literal from :data:`RENOVO_TIERS`, so this fails on a re-grade of
    any of the six — including ``IP Pathogenic``, which the committed fixture cannot witness.
    """
    row = pd.Series({RENOVO_COLUMN: tier.value})
    assert generate_clinical_summary(row) != NO_CLINICAL_DATA
    assert CLINICAL_VALUE_MAPPING[tier.value] == tier.clinical_class
    assert generate_clinical_summary(row) == _SUMMARY_LABELS[tier.clinical_class]


def test_a_stronger_renovo_score_never_maps_to_a_weaker_class():
    """The six classes are ordered by RENOVO's own score, and the mapping respects that order.

    A property rather than six equalities, and it holds a claim the equalities do not: that
    the *score* is what orders them. The mapping this fix replaced satisfied it too — what that
    one failed is the boundary claim below.
    """
    ranks = [CLINICAL_HIERARCHY.index(CLINICAL_VALUE_MAPPING[v]) for v in RENOVO_VALUES]
    assert ranks == sorted(ranks), dict(zip(RENOVO_VALUES, ranks))


def test_the_two_low_confidence_classes_land_on_one_class():
    """``LP Pathogenic`` and ``LP Benign`` read the same, because the score says they are.

    They are the low-confidence pair either side of RENOVO's decision boundary:
    ``LP Pathogenic`` bottoms out at 0.50000 and ``LP Benign`` tops out at 0.49985. Reading them
    three ranks apart — 🟠 Likely Pathogenic against 🟢 Likely Benign, which is what the app
    did until issue #108 — puts the ACMG ladder's widest gap across RENOVO's narrowest score
    difference. The pipeline tiers both 3, with the comment "borderline / conservative".
    """
    assert CLINICAL_VALUE_MAPPING["LP Pathogenic"] == CLINICAL_VALUE_MAPPING["LP Benign"]
    assert CLINICAL_VALUE_MAPPING["LP Pathogenic"] == "Uncertain_Significance"


def test_the_score_intervals_hold_on_the_committed_fixture(germline):
    """RENOVO's own score, read off the fixture, puts every class where this table says.

    This is what makes :data:`RENOVO_TIERS`' intervals an assertion rather than a decoration.
    The first draft checked them against each other — ``scores == sorted(scores, reverse=True)``
    over a list read from the same literal table, and a boundary check that subtracted two of its
    own numbers — which is the table compared with itself, the fault this module's own docstring
    warns about and `/code-review` caught anyway.

    ``germline_reference.maf`` carries ``RENOVO_pls``, so the intervals are checkable against
    committed data: every row's score must fall inside its class's interval, and the classes must
    be *separated* by score, which is the claim that HP/IP/LP is a confidence on one number.
    Five of the six classes occur here; ``IP Pathogenic`` does not, which is why its tier is
    pinned as a literal above rather than rested on this test.
    """
    scores = pd.to_numeric(germline["RENOVO_pls"], errors="coerce")
    intervals = {tier.value: (tier.score_from, tier.score_to) for tier in RENOVO_TIERS}

    # The pinned intervals must be six *bins*: ordered strongest first and pairwise disjoint.
    # Checking rows fall inside their interval cannot catch a widened one — a bigger interval
    # only makes that easier — and a mutation stretching `LP Benign` up to 0.6 survived the
    # first version of this test for exactly that reason. This is a coherence invariant of the
    # table rather than a comparison of it with itself: it fails on any overlap.
    bounds = [intervals[value] for value in RENOVO_VALUES]
    for (lo, hi) in bounds:
        assert lo < hi, (lo, hi)
    for (lower_lo, lower_hi), (upper_lo, _upper_hi) in zip(bounds[1:], bounds[:-1]):
        assert lower_hi < upper_lo, (lower_hi, upper_lo)

    seen = set()
    for value, group in scores.groupby(germline[RENOVO_COLUMN]):
        group = group.dropna()
        if value not in intervals or group.empty:
            continue
        seen.add(value)
        lo, hi = intervals[value]
        assert group.min() >= lo, (value, group.min(), lo)
        assert group.max() <= hi, (value, group.max(), hi)
    assert len(seen) >= 5, f"only {sorted(seen)} witnessed — the fixture has stopped covering these"

    # Separated by score, in the order the table lists them: each class's highest observed score
    # is below the next-stronger class's lowest. This is what "six bins of one score" means, and
    # it fails if two classes are ever given overlapping intervals.
    observed = [(value, scores[germline[RENOVO_COLUMN] == value].dropna()) for value in RENOVO_VALUES]
    observed = [(value, group) for value, group in observed if not group.empty]
    for (stronger, above), (weaker, below) in zip(observed, observed[1:]):
        assert below.max() < above.min(), (weaker, below.max(), stronger, above.min())


# --- what the fix cost the other derived column -------------------------------------------


@pytest.mark.parametrize("tier", RENOVO_TIERS, ids=RENOVO_VALUES)
def test_the_renovo_circle_follows_the_same_six_tiers(tier):
    """The ``RN`` glyph is the class's glyph, for every RENOVO value.

    Both derived columns read :data:`CLINICAL_VALUE_MAPPING`'s RENOVO entries, so re-grading
    three of them for ``Clinical_Summary`` re-graded this column too — on 4,485 of the 218,348
    rows carrying a RENOVO value (``LP Benign`` 🟢→🟡, ``LP Pathogenic`` 🟠→🟡, ``IP Pathogenic``
    🔴→🟠). Unlike the summary's reading, the circle has been *live* since issue #95, and nothing
    pinned its glyphs, so that would have been a silent re-grade of a rendered clinical scale —
    the failure issue #100 pinned the ESCAT circles against. Pinned to the glyph *and* to the
    class it comes from, so the two columns cannot be given different RENOVO readings either.
    """
    row = pd.Series({RENOVO_COLUMN: tier.value})
    sources = [("RN", "RENOVO", RENOVO_COLUMN)]
    assert pathogenicity_circle_glyphs(row, sources) == [tier.circle]
    assert _CLASS_GLYPHS[tier.clinical_class] == tier.circle


# --- no source outranks another -----------------------------------------------------------


def test_the_strongest_class_wins_whatever_source_asserts_it():
    """The label is the strongest class any source asserts — there is no source priority.

    ``_SUMMARY_SOURCE_COLUMNS`` used to be a dict mapping each column to a priority 1–5, and
    those numbers decided nothing: they were collected into a local list nothing read, while
    the label came from walking :data:`CLINICAL_HIERARCHY`. Issue #108 asked what priority
    RENOVO should have and the answer was that there is no priority to give it. This makes that
    a checkable claim rather than a deleted comment: a machine-learning prediction of
    ``HP Pathogenic`` outranks a curated ``Benign``, in whichever order the columns arrive.

    That escalation is the pipeline's behaviour too, not a reading this app invented: on the
    corpus rows where RENOVO predicts ``HP Pathogenic`` over a curated benign call,
    ``variantalker_naive`` reads ``Pathogenic``.
    """
    values = {RENOVO_COLUMN: "HP Pathogenic", "ClinVar_VCF_CLNSIG": "Benign", "InterVar": "Benign"}
    assert generate_clinical_summary(pd.Series(values)) == _SUMMARY_LABELS["Pathogenic"]
    assert generate_clinical_summary(pd.Series(dict(reversed(list(values.items()))))) == (
        _SUMMARY_LABELS["Pathogenic"]
    )


def test_the_order_of_the_source_columns_carries_no_weight(germline, monkeypatch):
    """Reading the five columns in reverse gives every row of the fixture the same label.

    The claim :data:`_SUMMARY_SOURCE_COLUMNS` makes in words, asserted. It earns its place by
    the mutation that found it: reversing that tuple was the one mutation of seven that the
    suite did **not** catch, and rightly so — it is semantically equivalent, which is precisely
    what "the priorities decided nothing" means. A guard that cannot fail is worthless, but a
    mutation that cannot break anything is a property worth pinning, and this is the difference.

    It also fails if selection ever goes back to being decided by the position of a source in
    that tuple, which is the shape the deleted priority numbers pretended to have.
    """
    forward = _labels(germline)
    monkeypatch.setattr(
        clinical_summary, "_SUMMARY_SOURCE_COLUMNS", tuple(reversed(_SUMMARY_SOURCE_COLUMNS))
    )
    assert _labels(germline) == forward


def test_adding_a_source_can_only_strengthen_a_label(germline):
    """Reading one more source never weakens a label, on every row of the real fixture.

    The property that bounds this whole change: the hierarchy takes the strongest class, so
    wiring RENOVO in could move a label up the ladder or leave it, never down. It is why the
    fix needed no worry about a curated call being *lost*, and why ``🔍 No Clinical Data`` could
    only shrink. Asserted over the fixture rather than argued in a comment.
    """
    ranks = {label: rank for rank, label in enumerate(_SUMMARY_LABELS[n] for n in CLINICAL_HIERARCHY)}
    for i, before, after in zip(
        germline.index, _labels(germline.drop(columns=[RENOVO_COLUMN])), _labels(germline)
    ):
        if before == NO_CLINICAL_DATA:
            continue
        assert ranks[after] <= ranks[before], (i, before, after)
