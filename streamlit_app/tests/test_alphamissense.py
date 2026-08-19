"""AlphaMissense's own section: its column, its band, and the claims it may not make (issue #203).

The section is small; what it must not say is not. Four claims are live here, and each is one the
prototype rendered wrongly before the dev chose against it:

* it reads **``am_pathogenicity``**, not the dbNSFP ``AlphaMissense_score`` this panel's naming
  convention would have picked — they are different annotations;
* it names the band **from the score** rather than echoing ``am_class``, whose spelling 100 of 139
  real files overstate;
* it leads with **provenance**, and that sentence is *asked* of the panel's own wiring rather than
  written into a string that can outlive its truth;
* it claims **no ACMG evidence strength**, which is what keeps a fourth predictor section from
  reading as a fourth PP3/BP4 surface after issue #201 refused RENOVO one.

``tests/test_clinical_badges.py`` holds the other half — that the badge it replaced is drawn
nowhere. ``docs/wayfinder/issue-204/README.md`` records what each guard here was mutated into.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from components.alphamissense import (  # noqa: E402
    AM_COLUMN,
    AM_LIKELY_BENIGN_BELOW,
    AM_LIKELY_PATHOGENIC_FROM,
    _provenance_caption,
    _scale_caption,
    alphamissense_reading,
    read_by_any_classifier,
)
from components.reference_scales import CLASSIFIER_INPUTS, CLINGEN_SVI_TABLE  # noqa: E402
from config.contaminated_columns import PREDICTOR_COLUMNS  # noqa: E402


# ---------------------------------------------------------------------------
# The column it reads
# ---------------------------------------------------------------------------


def test_it_reads_am_pathogenicity_and_not_the_dbnsfp_column():
    """The two are **different annotations**, and the natural wiring would have been wrong.

    Of 55,353 real rows carrying both, only 6.3% agree to 1e-9, 33.9% differ beyond rounding, and
    **346 rows across 84 files fall in different AlphaMissense classes** — one germline variant
    holds 0.3475 in one column and 0.287 in the other, *ambiguous* beside *likely benign*. And
    ``am_pathogenicity`` covers strictly more files (55 germline, 84 somatic against 53 and 51), so
    the dbNSFP spelling would have gone silent on 15,340 variants where the class speaks.
    """
    assert AM_COLUMN == "am_pathogenicity"

    source = (STREAMLIT_APP / "components" / "alphamissense.py").read_text()
    read_literals = {
        node.value
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    reading = alphamissense_reading(
        pd.Series({"AlphaMissense_score": "0.9", "am_class": "likely_pathogenic"}), frozenset()
    )
    assert reading is None, "the section read a column that is not am_pathogenicity"
    assert "AlphaMissense_score" not in {
        column for _name, column, *_rest in CLINGEN_SVI_TABLE
    }, "the dbNSFP column has been wired into the ClinGen table, which issue #203 ruled out"
    assert "AlphaMissense_score" in read_literals, (
        "the module no longer records which column it deliberately does not read"
    )


def test_the_column_it_scores_is_measured_for_contamination():
    """A table that scores a column asks whether the column holds a score (issue #194).

    ``am_pathogenicity`` is clean on every file — maximum chromosome fraction exactly 0.0000,
    measured for this ticket — but a predictor column ANNOVAR did not supply arrives holding the
    **chromosome**, not a ``.``, so a plausible-looking number can be nothing at all. The cheap
    move is to measure rather than to argue an exemption.
    """
    assert AM_COLUMN in PREDICTOR_COLUMNS


def test_a_contaminated_column_is_warned_about_and_not_banded():
    """A chromosome must not be given an AlphaMissense band.

    Asked *before* the value is parsed, which is the whole of one bug next door: on a ``chrX``
    variant the contaminated cell holds ``X``, which ``float()`` refuses — so parsing first drops
    the row through the unreadable branch and the warning is never reached.
    """
    row = pd.Series({AM_COLUMN: "X", "Chromosome": "chrX"})
    reading = alphamissense_reading(row, frozenset({AM_COLUMN}))

    assert reading is not None
    assert reading.side == "chromosome"
    assert "⚠" in reading.value
    assert not reading.readable
    for band in ("likely benign", "likely pathogenic", "ambiguous"):
        assert band not in reading.says


# ---------------------------------------------------------------------------
# The band, named from the score
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "score,expected",
    [
        ("0.3396", "likely benign"),
        ("0.3401", "ambiguous"),
        ("0.5626", "ambiguous"),
        ("0.5640", "likely pathogenic"),
        ("0.0", "likely benign"),
        ("1.0", "likely pathogenic"),
    ],
)
def test_the_band_is_named_at_alphamissenses_own_cutoffs(score, expected):
    """The four values are the corpus's own class boundaries, to four decimal places.

    ``am_class`` is an exact function of ``am_pathogenicity`` at 0.34 / 0.564: ``likely_benign``
    tops out at 0.3396 and ``ambiguous`` starts at 0.3401; ``ambiguous`` tops out at 0.5626 and
    ``likely_pathogenic`` starts at 0.5640. Testing the boundaries rather than the middles is what
    makes this a test of the cutoffs.
    """
    reading = alphamissense_reading(pd.Series({AM_COLUMN: score}), frozenset())
    assert reading.says.startswith(expected)


def test_the_cutoffs_are_the_published_pair():
    assert (AM_LIKELY_BENIGN_BELOW, AM_LIKELY_PATHOGENIC_FROM) == (0.34, 0.564)


def test_the_band_never_echoes_the_files_own_class_spelling():
    """100 of 139 real files spell ``am_class`` in a vocabulary AlphaMissense does not publish.

    They write ``benign``/``ambiguous``/``pathogenic``; AlphaMissense publishes only the
    ``likely_`` form. Echoing the label would reproduce an overstatement, while naming the band
    from the score does not disagree with the annotator — the two vocabularies occupy identical
    score ranges and agree on the class, so this is declining to reproduce a spelling rather than
    contradicting one.
    """
    row = pd.Series({AM_COLUMN: "0.9", "am_class": "pathogenic"})
    reading = alphamissense_reading(row, frozenset())

    assert reading.says.startswith("likely pathogenic")
    assert reading.value == "0.9"
    # The overstated spelling must not travel onto the row by any route.
    assert "pathogenic —" in reading.says
    assert not reading.says.startswith("pathogenic")


def test_a_variant_the_predictor_did_not_score_draws_nothing():
    """The section makes one claim, and there is no such claim to make on an unscored variant."""
    assert alphamissense_reading(pd.Series({"Hugo_Symbol": "TP53"}), frozenset()) is None
    assert alphamissense_reading(pd.Series({AM_COLUMN: "."}), frozenset()) is None
    assert alphamissense_reading(pd.Series({AM_COLUMN: ""}), frozenset()) is None
    assert alphamissense_reading(pd.Series({AM_COLUMN: "not a number"}), frozenset()) is None


# ---------------------------------------------------------------------------
# What the section claims, and what it must not
# ---------------------------------------------------------------------------


def test_the_provenance_claim_is_asked_of_the_wiring_rather_than_asserted():
    """*"Neither classifier read AlphaMissense"* must not be able to outlive its own truth.

    That is the defect issue #203 found one section down, where ``_empty_caption``'s hand-written
    predictor list goes stale. Here the sentence is derived from
    :data:`components.reference_scales.CLASSIFIER_INPUTS`, which
    ``tests/test_reference_scales.py`` already holds to what ``predictor_context`` actually reads —
    so the chain runs from the classifier's real inputs to this caption with no transcription in it.
    """
    assert read_by_any_classifier() == ()

    every_input = {
        column for inputs in CLASSIFIER_INPUTS.values() for _name, column in inputs
    }
    assert AM_COLUMN not in every_input
    assert "AlphaMissense_score" not in every_input

    assert "did not read AlphaMissense" in _provenance_caption("InterVar")
    assert "contributed to the classification above" in _provenance_caption("InterVar")
    assert "contributed to the tier above" in _provenance_caption("CancerVar")


def test_the_provenance_claim_goes_the_other_way_if_a_classifier_ever_reads_it(monkeypatch):
    """The derivation is only worth having if it can produce the opposite sentence.

    A caption computed from a table that can only ever be empty is a hand-written caption with
    extra steps, so the table is given an AlphaMissense input and the sentence has to change.
    """
    import components.alphamissense as module

    monkeypatch.setitem(
        module.CLASSIFIER_INPUTS, "InterVar", (("AlphaMissense", AM_COLUMN),)
    )
    assert read_by_any_classifier() == ("InterVar",)
    assert "InterVar read AlphaMissense" in _provenance_caption("InterVar")
    assert "did not read" not in _provenance_caption("InterVar")


def test_a_row_carrying_no_classifier_is_not_told_about_one():
    """There is no classification above for anything to have contributed to.

    The same split :func:`components.reference_scales._provenance_captions` makes, and the same
    reason: naming a tool that drew nothing on this page would be false.
    """
    caption = _provenance_caption(None)
    assert "InterVar" not in caption and "CancerVar" not in caption
    assert "No guideline classifier scored this variant" in caption


def test_the_section_claims_no_acmg_evidence_strength():
    """The reconciliation issue #203 asked #204 for, in copy rather than by reopening either call.

    Issue #201 refused RENOVO a predictor table partly because all three the panel drew were
    PP3/BP4 surfaces; this is a fourth that is not one. What makes both calls stand is that this
    section names AlphaMissense's own three-class boundaries and disclaims the ACMG strengths
    outright — what RENOVO was refused was a seat among *calibrated* predictors, and this claims no
    calibration at all.
    """
    caption = _scale_caption()
    assert "AlphaMissense's own three-class boundaries" in caption
    assert "not ACMG PP3/BP4 evidence strengths" in caption

    for strength in ("Supporting", "Moderate", "Strong", "Very Strong"):
        assert strength not in caption

    reading = alphamissense_reading(pd.Series({AM_COLUMN: "0.9"}), frozenset())
    for strength in ("PP3", "BP4", "Supporting", "Moderate", "Strong"):
        assert strength not in reading.says


def test_it_is_not_a_row_in_the_clingen_table():
    """A fifth row there makes that section's own copy false, and the prototype rendered it.

    With no thresholds the tally reports *"BP4 (benign): 0/1 predictors"* — asserting AlphaMissense
    *failed* a ClinGen threshold that does not exist for it — and with AlphaMissense's own cutoffs
    wearing ClinGen's strength names it reports *"1/1"* and prints *Supporting*, indistinguishable
    on screen from a calibrated row.
    """
    named = {name for name, *_rest in CLINGEN_SVI_TABLE}
    assert "AlphaMissense" not in named
    assert named == {"REVEL", "BayesDel (no AF)", "VEST4", "MutPred2", "CADD"}


def test_the_section_draws_after_both_evidence_sections_and_above_the_clingen_table():
    """Position, asserted on the panel's own source order — the thing the prototype got wrong.

    Mounted under the guideline badges, the shape read as part of the classification; the shape
    survived and the position did not. Above the ClinGen table so that table stays the panel's last
    classification row, which is the claim its own module makes about itself.
    """
    source = (STREAMLIT_APP / "components" / "variant_detail.py").read_text()
    tree = ast.parse(source)
    panel = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "render_variant_detail_panel"
    )
    positions = {
        node.func.id: node.lineno
        for node in ast.walk(panel)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert positions["render_acmg_evidence"] < positions["render_alphamissense"]
    assert positions["render_cbp_evidence"] < positions["render_alphamissense"]
    assert positions["render_alphamissense"] < positions["render_reference_scale"]
