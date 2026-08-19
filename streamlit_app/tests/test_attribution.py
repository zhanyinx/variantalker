"""What a per-variant explanation may claim, and what stops it claiming something false.

Issue #147. This surface says *which of your settings this variant falls outside*, which is a
claim about the filter — so the failure mode it has to be held against is not "does it
render" but the one every parameter echo on this map has fallen into: describing a cut that
did not happen. Issue #28 deleted eight on-screen echoes and **every one was wrong**; #137's
recap was false on almost every real MAF until review caught two clauses it reported as
applied when the file's columns had dropped them whole.

So the guards here are the ways that class of defect gets in:

* the isolation does not isolate — a "neutral" value that still excludes rows, or one clause's
  override moving another clause;
* the answer names a setting in words the interface does not use;
* the answer names *one* setting where the variant falls outside several;
* the per-row answer and the whole-report answer disagree about the same row;
* a second copy of a fact the app already states once.

Every assertion below is mutation-checked and the mutation is named, because a green test
over a wired-up module proves nothing on its own (issues #67, #83, #90, #140).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.param_labels import label_of  # noqa: E402
from config.vocabularies import GUIDELINE_SOURCES as VOCAB_SOURCES  # noqa: E402
from filters.attribution import (  # noqa: E402
    CLASSIFICATION_GROUP,
    CLASSIFICATION_SOURCES,
    CLAUSE_COLUMNS,
    CLAUSE_LABEL_IDS,
    CLAUSES,
    IN_REPORT_BECAUSE,
    NEUTRAL,
    NOT_IN_FILE,
    SAYS_NOTHING,
    _all_neutral,
    attribute_report,
    explain_variant,
    label_for,
)
from filters.variant_filters import (  # noqa: E402
    MAFIGATE_FILTER,
    PASS,
    REASON_BOTH,
    REASON_CRITERIA,
    REASON_RESCUE,
    apply_filters,
)

ARMS = ("somatic", "germline")

#: A frame with one obliging row per arm, built so that every clause can be made to fail on
#: its own. The values are chosen against the presets rather than being round numbers: the
#: point of each test below is which clause moves, so the starting row has to *pass*.
_BASE = {
    "Hugo_Symbol": "BRCA1",
    "Chromosome": "chr17",
    "Start_Position": 41276045,
    "Reference_Allele": "A",
    "Tumor_Seq_Allele2": "T",
    "Variant_Classification": "Missense_Mutation",
    "t_alt_count": 300,
    "t_ref_count": 300,
    "tumor_f": 0.45,
    "ClinVar_VCF_CLNSIG": "Pathogenic",
}

_ARM_COLUMNS = {
    "somatic": {"CancerVar": "Tier_I_strong", "ESCAT": "I-A", "CIViC_Evidence_Level": "A"},
    "germline": {"InterVar": "Pathogenic", "RENOVO_Class": "HP Pathogenic"},
}


def a_frame(arm: str, **overrides) -> pd.DataFrame:
    """One row that passes on ``arm``, with ``overrides`` applied."""
    row = dict(_BASE)
    row.update(_ARM_COLUMNS[arm])
    row.update(overrides)
    return pd.DataFrame([row])


def params_for(arm: str, **overrides) -> dict:
    """Parameters that keep :func:`a_frame`'s row, with ``overrides`` applied."""
    base = {
        "sample_type": arm,
        "min_depth": 50,
        "vaf_threshold": 0.05,
        "vaf_threshold_germline": 0.05,
        "filter_variant_classification": ["Silent", "Intron"],
        "filter_clinvar": ["Pathogenic", "Likely_pathogenic"],
        "filter_cancervar": ["Tier_I_strong", "Tier_II_potential"],
        "filter_escat": ["I-A", "I-B"],
        "filter_civic": ["A", "B"],
        "filter_intervar": ["Pathogenic", "Likely pathogenic"],
        "filter_renovo": ["HP Pathogenic"],
        "filter_genes": [],
        "skip_civic": False,
        "skip_pathogenic": False,
        "max_freq_population": 1.0,
    }
    base.update(overrides)
    return base


def passes(frame, params) -> bool:
    labelled, _ = apply_filters(frame, params)
    return bool((labelled[MAFIGATE_FILTER] == PASS).iloc[0])


# ---------------------------------------------------------------------------
# The tables, against the facts the app already states elsewhere
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_source_names_come_from_their_control(arm):
    """Each source is called what its own control calls it, not something new.

    The failure this prevents is a second vocabulary: issue #79 catalogued five places where
    Help had drifted from the tooltip beside the control, and the copy beside the control was
    right every time. So the name here has to be a substring of the label there.

    Mutated by renaming ``RENOVO`` to ``Renovo`` — the plausible drift, since the control says
    ``Keep RENOVO Classifications`` — and this fails.
    """
    for source in CLASSIFICATION_SOURCES[arm]:
        label = label_of(source.label_id)
        assert source.name in label, (
            f"{source.name!r} is not the word {label!r} uses for this source"
        )


@pytest.mark.parametrize("arm", ARMS)
def test_the_classification_sources_are_the_arms_own(arm):
    """This module's source table holds exactly the arm's guideline parameters.

    ``config.vocabularies.GUIDELINE_SOURCES`` already states which keep-lists each arm's
    filter ORs together, and ``tests/test_param_contract.py`` keeps *that* copy honest
    against the vendored signatures. This asserts the two agree, so the table here is a
    widening of one fact rather than a second copy of it — the failure ``vendor/README.md``
    records for a hand-copied ``KEEP`` that drifted to 39 entries against 45.

    Mutated by dropping the ESCAT entry from the somatic tuple, and this fails.
    """
    assert {source.param_key for source in CLASSIFICATION_SOURCES[arm]} == set(
        VOCAB_SOURCES[arm]
    )


def test_every_clause_can_be_named_and_neutralised():
    """No clause is missing a label or a neutral value, on either arm.

    A clause with no label would render as a ``KeyError`` inside a dialog; a clause with no
    neutral value would be silently skipped, and the answer would be missing a closed path
    without saying so.

    Mutated by deleting the ``depth`` entry from ``NEUTRAL``, and this fails.
    """
    assert set(CLAUSES) == set(NEUTRAL) | {"guidelines"}
    assert set(CLAUSES) == set(CLAUSE_LABEL_IDS) | {"guidelines"}
    for arm in ARMS:
        for clause in CLAUSES:
            assert label_for(clause, arm)
    assert label_for("guidelines", "somatic") == CLASSIFICATION_GROUP


def test_the_reason_table_covers_every_way_into_the_report():
    """Every cell of the union that means PASS has something to say about it.

    ``MAFigate_reason`` has four cells and three of them are in the report. A cell missing
    here would draw a bare "In this report" with no account, which is the state
    :func:`components.variant_detail._render_report_standing` falls back to.

    Mutated by removing ``REASON_RESCUE``, and this fails.
    """
    assert set(IN_REPORT_BECAUSE) == {REASON_CRITERIA, REASON_BOTH, REASON_RESCUE}


# ---------------------------------------------------------------------------
# The isolation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_all_neutral_keeps_a_row_every_clause_would_drop(arm):
    """With every setting neutral and retention off, a row failing all of them still passes.

    This is the property the whole mechanism rests on: if a "neutral" value still excluded
    rows, every clause would read as failing and the answer would name all of them.

    The row here is built to fail each clause at once — excluded classification, below depth,
    below VAF, no matching guideline call, and common in every frequency column.

    Mutated by changing ``NEUTRAL["depth"]`` to ``{"min_depth": 1}`` with the row's counts at
    zero, and this fails.
    """
    hostile = {
        "Variant_Classification": "Silent",
        "t_alt_count": 0,
        "t_ref_count": 0,
        "tumor_f": 0.0,
        "ClinVar_VCF_CLNSIG": "Benign",
        "gnomAD_exome_AF": 0.99,
    }
    hostile.update(
        {"CancerVar": "Tier_IV_benign", "ESCAT": "X", "CIViC_Evidence_Level": "E"}
        if arm == "somatic"
        else {"InterVar": "Benign", "RENOVO_Class": "HP Benign"}
    )
    frame = a_frame(arm, **hostile)
    params = params_for(arm, max_freq_population=0.01)

    assert not passes(frame, params), "the hostile row should be rejected as set"
    assert passes(frame, _all_neutral(frame, arm, params)), (
        "every setting neutral and retention off should keep it"
    )


@pytest.mark.parametrize("arm", ARMS)
@pytest.mark.parametrize(
    "clause,overrides",
    [
        ("exclude_classifications", {"Variant_Classification": "Silent"}),
        ("depth", {"t_alt_count": 1, "t_ref_count": 1}),
        ("vaf", {"tumor_f": 0.001}),
    ],
)
def test_one_broken_column_names_exactly_that_clause(arm, clause, overrides):
    """A row that fails one clause and no other is told about that one clause.

    The 5–30% case: exactly one setting, so there is a one-sentence answer and it must be the
    right sentence. ``ClinVar_VCF_CLNSIG`` is left at ``Pathogenic`` here, so pathogenic
    retention would carry the row — which is why each override also turns retention off.

    Mutated by having ``_isolate`` restore the whole parameter dict rather than one clause,
    and every one of these fails with the full set of clauses named.
    """
    frame = a_frame(arm, **overrides)
    params = params_for(arm, skip_pathogenic=True)

    explanation = explain_variant(frame, params)
    assert not explanation.in_report
    assert [outcome.clause for outcome in explanation.failing] == [clause]
    assert explanation.sole_reason is not None
    assert explanation.sole_reason.label == label_for(clause, arm)
    # One clause failing is exactly the case where relaxing it alone brings the row back.
    assert [outcome.clause for outcome in explanation.recoverable] == [clause]


@pytest.mark.parametrize("arm", ARMS)
def test_two_broken_columns_name_both_and_offer_no_single_change(arm):
    """The 70–95% case: several settings, none of them recoverable alone.

    This is the finding that decides the copy. A surface built on the leave-one-out alone
    would say nothing here, and one that named a single clause would be true of that clause
    and misleading about the variant.

    Mutated by having :attr:`VariantExplanation.recoverable` return every failing clause, and
    this fails on the emptiness assertion.
    """
    frame = a_frame(arm, t_alt_count=1, t_ref_count=1, tumor_f=0.001)
    params = params_for(arm, skip_pathogenic=True)

    explanation = explain_variant(frame, params)
    assert {outcome.clause for outcome in explanation.failing} == {"depth", "vaf"}
    assert explanation.sole_reason is None
    assert explanation.recoverable == ()


@pytest.mark.parametrize("arm", ARMS)
def test_a_variant_no_guideline_speaks_for_names_the_block_once(arm):
    """The classification block is one clause with its sources beneath it, not four clauses.

    The block is False only when *every* source is False, so "which source rejected it" has
    no answer. Naming it once and listing the sources is what the dev chose over a five-label
    chain; the sources still appear, so nothing is hidden.

    Mutated by adding one clause per source to ``CLAUSES``, and this fails on the length.
    """
    silent = (
        {"CancerVar": "Tier_IV_benign", "ESCAT": "X", "CIViC_Evidence_Level": "E"}
        if arm == "somatic"
        else {"InterVar": "Benign", "RENOVO_Class": "HP Benign"}
    )
    frame = a_frame(arm, ClinVar_VCF_CLNSIG="Benign", **silent)
    params = params_for(arm, skip_pathogenic=True)

    explanation = explain_variant(frame, params)
    assert [outcome.clause for outcome in explanation.failing] == ["guidelines"]
    outcome = explanation.failing[0]
    assert outcome.label == CLASSIFICATION_GROUP
    assert [name for name, _ in outcome.values] == [
        source.name for source in CLASSIFICATION_SOURCES[arm]
    ]


def test_a_source_the_file_lacks_is_said_to_be_missing_not_blank():
    """An absent column and a blank value are different sentences.

    Measured on real data: 415 somatic rows are admitted by widening ESCAT and not by
    widening CancerVar, because CancerVar says nothing there — so "changing this setting
    would bring it back" is false for a source the file is silent on, and the two kinds of
    silence are not the same thing to a user deciding whether to re-annotate.

    Mutated by having ``shown_value`` return :data:`SAYS_NOTHING` for an absent column too,
    and this fails.
    """
    frame = a_frame("somatic", ClinVar_VCF_CLNSIG="Benign", CancerVar="Tier_IV_benign",
                    ESCAT=".")
    frame = frame.drop(columns=["CIViC_Evidence_Level"])
    params = params_for("somatic", skip_pathogenic=True)

    explanation = explain_variant(frame, params)
    values = dict(explanation.failing[0].values)
    assert values["CIViC"] == NOT_IN_FILE
    assert values["ESCAT"] == SAYS_NOTHING
    assert values["CancerVar"] == "Tier_IV_benign"


def test_a_filled_input_is_named_as_the_files_gap():
    """Where a filter input was absent and filled, the answer says so.

    Otherwise it names a setting the user chose for a decision their file's missing column
    made — the ticket's fifth bullet, and #136's line that a warning states what MAFigate did.

    Mutated by intersecting ``plan.filled`` with the empty set, and this fails.

    Every guideline source has to be silenced for the block to fail at all — the first draft
    of this test dropped ``CancerVar`` alone and the row stayed in the report on its ESCAT
    call, which is the OR working exactly as #136 describes it.
    """
    frame = a_frame(
        "somatic",
        ClinVar_VCF_CLNSIG="Benign",
        ESCAT="X",
        CIViC_Evidence_Level="E",
    ).drop(columns=["CancerVar"])
    params = params_for("somatic", skip_pathogenic=True)

    explanation = explain_variant(frame, params)
    assert not explanation.in_report
    assert [outcome.clause for outcome in explanation.failing] == ["guidelines"]
    assert "CancerVar" in explanation.filled_inputs


# ---------------------------------------------------------------------------
# The variant that *is* in the report
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_a_passing_variant_is_told_why_it_is_in(arm):
    """The dialog is shared by both tabs, so the block must be true of a passing row too.

    Mutated by having ``explain_variant`` return an empty failing set without setting
    ``in_report``, and this fails.
    """
    explanation = explain_variant(a_frame(arm), params_for(arm))
    assert explanation.in_report
    assert explanation.reason in IN_REPORT_BECAUSE
    assert explanation.failing == ()


@pytest.mark.parametrize("arm", ARMS)
def test_a_rescued_variant_says_it_did_not_meet_the_criteria(arm):
    """A row in the report only by pathogenic retention is not told it met the criteria.

    This is the cell nothing in the app reported before #28's decomposition, and it is 19.2%
    of the germline reference report. Saying "it met your filter criteria" of these rows would
    be the report contradicting itself.

    Mutated by mapping ``REASON_RESCUE`` onto the ``REASON_CRITERIA`` wording, and this fails.
    """
    # Below depth, so the criteria path is closed, but ClinVar-pathogenic so retention keeps it.
    frame = a_frame(arm, t_alt_count=1, t_ref_count=1, ClinVar_VCF_CLNSIG="Pathogenic")
    explanation = explain_variant(frame, params_for(arm))
    assert explanation.in_report
    assert explanation.reason == REASON_RESCUE
    assert "did not meet" in IN_REPORT_BECAUSE[explanation.reason]


def test_explain_variant_refuses_more_than_one_row():
    """A silently-averaged answer over several rows would be worse than a refusal.

    Mutated by taking ``.iloc[0]`` instead of raising, and this fails.
    """
    frame = pd.concat([a_frame("somatic"), a_frame("somatic")], ignore_index=True)
    with pytest.raises(ValueError, match="one-row frame"):
        explain_variant(frame, params_for("somatic"))


# ---------------------------------------------------------------------------
# The whole-report answer
# ---------------------------------------------------------------------------


def test_the_report_answer_agrees_with_the_per_row_answer():
    """The same row is not told two different things by the two surfaces.

    They are one mechanism at two scales, and the reason to build them together is that they
    cannot then disagree. This asserts it over a frame holding one row per failure shape.

    Mutated by having :func:`attribute_report` isolate against the *live* parameter dict
    rather than the one it was handed, and this fails.
    """
    rows = [
        a_frame("somatic").iloc[0],
        a_frame("somatic", Variant_Classification="Silent").iloc[0],
        a_frame("somatic", t_alt_count=1, t_ref_count=1).iloc[0],
        a_frame("somatic", tumor_f=0.001, ClinVar_VCF_CLNSIG="Benign",
                CancerVar="Tier_IV_benign", ESCAT="X", CIViC_Evidence_Level="E").iloc[0],
    ]
    frame = pd.DataFrame(rows).reset_index(drop=True)
    params = params_for("somatic", skip_pathogenic=True)

    report = attribute_report(frame, params)
    labelled, _ = apply_filters(frame, dict(params))
    rejected = labelled[MAFIGATE_FILTER] != PASS

    assert report.rows == len(frame)
    assert report.left_out == int(rejected.sum())
    assert report.in_report == len(frame) - report.left_out

    per_label: dict[str, int] = {}
    for index in frame.index[rejected]:
        for outcome in explain_variant(frame.loc[[index]], params).failing:
            per_label[outcome.label] = per_label.get(outcome.label, 0) + 1
    assert dict(report.excluded_by) == per_label


def test_the_report_counts_overlap_rather_than_partition():
    """A row outside two settings is counted by both, which is the framing #147 chose.

    Stated as a test because the alternative is a silent change of meaning: if these ever
    summed to :attr:`left_out` they would have become the partitioning framing, and the
    caption promising they add to more than it would be false.

    Mutated by counting only rows whose failing set has one member, and this fails.
    """
    frame = pd.DataFrame(
        [a_frame("somatic", t_alt_count=1, t_ref_count=1, tumor_f=0.001).iloc[0]]
    ).reset_index(drop=True)
    report = attribute_report(frame, params_for("somatic", skip_pathogenic=True))

    assert report.left_out == 1
    assert sum(count for _, count in report.excluded_by) > report.left_out
    assert {label for label, _ in report.excluded_by} == {
        label_for("depth", "somatic"),
        label_for("vaf", "somatic"),
    }


def test_a_report_that_left_nothing_out_attributes_nothing():
    """No exclusions, no table of zeroes beneath a complete report.

    Mutated by returning a zero entry per clause, and this fails.
    """
    report = attribute_report(a_frame("somatic"), params_for("somatic"))
    assert report.left_out == 0
    assert report.excluded_by == ()


def test_the_report_counts_are_descending():
    """Ordered by how much each setting excludes, because that is the question it answers.

    Mutated by sorting on the label, and this fails.

    The fixture is built so that alphabetical order and descending order **disagree**, which
    the first draft of it did not do: two rows failing only VAF and one failing only the
    classification block put the larger count on ``VAF Threshold`` and the smaller on
    ``Clinical classification filters``, so a label sort inverts them. Sorted by count the
    first draft happened to come out alphabetical too, and the mutation survived.
    """
    frame = pd.DataFrame(
        [
            a_frame("somatic", tumor_f=0.001).iloc[0],
            a_frame("somatic", tumor_f=0.001, Start_Position=41276046).iloc[0],
            a_frame(
                "somatic",
                ClinVar_VCF_CLNSIG="Benign",
                CancerVar="Tier_IV_benign",
                ESCAT="X",
                CIViC_Evidence_Level="E",
            ).iloc[0],
        ]
    ).reset_index(drop=True)
    report = attribute_report(frame, params_for("somatic", skip_pathogenic=True))

    assert report.excluded_by == (
        (label_for("vaf", "somatic"), 2),
        (CLASSIFICATION_GROUP, 1),
    )
    assert CLASSIFICATION_GROUP < label_for("vaf", "somatic"), (
        "the fixture must put the smaller count on the alphabetically earlier label, "
        "or a label sort would pass this"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_the_clauses_read_the_columns_they_claim_to(arm):
    """Every column named beside a clause is one the filter reads for that clause.

    A column shown beside the wrong setting is the quietest way this surface could mislead:
    the number would be real and the sentence about it false.

    Mutated by moving ``tumor_f`` from ``vaf`` to ``depth``, and this fails.
    """
    from filters.absent_columns import REQUIRED_INPUTS
    from filters.variant_filters import GENE_SYMBOL

    # ``Hugo_Symbol`` is deliberately **not** in ``REQUIRED_INPUTS``: it is the one filter
    # input the *app* guards rather than fills, because the vendored body only reads it when
    # a gene file exists and the app is what writes that file (``_gene_symbols``). So it is a
    # column a clause reads without being one absence fills, and it is named from the same
    # constant the filter reads it by rather than spelled again here.
    readable = (
        set(REQUIRED_INPUTS[arm]) | set(CLAUSE_COLUMNS["frequency"]) | {GENE_SYMBOL}
    )
    for clause, columns in CLAUSE_COLUMNS.items():
        for column in columns:
            assert column in readable, f"{column} is not a filter input on {arm}"
    assert CLAUSE_COLUMNS["vaf"] == ("tumor_f",)
    assert set(CLAUSE_COLUMNS["depth"]) == {"t_alt_count", "t_ref_count"}


def test_this_module_does_not_need_streamlit():
    """``filters/`` stays importable without Streamlit, and this is its newest member.

    The parity harness imports this package in an environment that has no Streamlit, which is
    the line #54 drew when it decided ``page_modules/param_store.py`` could *not* move into
    ``config/``. This module is the one most likely to break it: what it produces is rendered
    prose, so reaching for a component helper — ``components.variant_detail.shown``, say — is
    a natural mistake that nothing else here would catch.

    Mutated by adding ``import streamlit`` to ``filters/attribution.py``, and this fails.
    """
    import builtins
    import importlib
    import sys

    real_import = builtins.__import__

    def refuse_streamlit(name, *args, **kwargs):
        if name == "streamlit" or name.startswith("streamlit."):
            raise ImportError("streamlit is not installed in this environment")
        return real_import(name, *args, **kwargs)

    saved = {
        name: module
        for name, module in sys.modules.items()
        if name == "filters.attribution"
    }
    for name in saved:
        del sys.modules[name]
    builtins.__import__ = refuse_streamlit
    try:
        importlib.import_module("filters.attribution")
    finally:
        builtins.__import__ = real_import
        sys.modules.update(saved)
