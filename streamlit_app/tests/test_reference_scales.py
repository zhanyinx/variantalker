"""The ClinGen SVI reference scale says what it is, and does not claim to be the verdict.

Issue #191. The section's whole job is to be read as *separate*, so most of these guards are
about its **copy** and its **position** rather than its arithmetic — the arithmetic moved
unchanged from ``variant_detail`` and ``tests/test_contaminated_columns.py`` still holds it.

The load-bearing one is :func:`test_no_predictor_this_table_scores_is_read_by_either_classifier`.
The label's first sentence is a factual claim about wiring — *"InterVar read none of these
predictors"* — and it is true only while the two column sets stay disjoint. Nothing else in the
app would notice them overlapping.
"""

from __future__ import annotations

import os
import sys
from unittest.mock import MagicMock

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components import reference_scales  # noqa: E402
from components.reference_scales import (  # noqa: E402
    CLASSIFIER_INPUTS,
    CLINGEN_SVI_TABLE,
    render_reference_scale,
)


def _row(**overrides):
    """A germline-looking variant with a REVEL score that reaches PP3 Supporting."""
    cells = {
        "Hugo_Symbol": "BRCA1",
        "Chromosome": "chr17",
        "Variant_Classification": "Missense_Mutation",
        "REVEL_score": 0.70,
    }
    cells.update(overrides)
    return pd.Series(cells)


def _draw(monkeypatch, row, tool="InterVar", contaminated=()):
    """Every string the section draws, plus the expander's label and expanded flag."""
    drawn = []
    opened = {}

    fake = MagicMock()
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))

    def expander(label, *a, **k):
        opened["label"] = label
        opened["expanded"] = k.get("expanded", a[0] if a else None)
        return MagicMock()

    fake.expander.side_effect = expander
    monkeypatch.setattr(reference_scales, "st", fake)

    render_reference_scale(row, tool, frozenset(contaminated))
    return drawn, opened


# ---------------------------------------------------------------------------
# The claim the label makes
# ---------------------------------------------------------------------------


def test_no_predictor_this_table_scores_is_read_by_either_classifier():
    """The label says the classifier read none of these. This is that sentence, as a test.

    Issue #191 measured the two sets as disjoint and built the whole label on it: leading with
    *provenance* — nothing here fed the verdict — is what lets the section separate itself from
    the classification above without claiming either calibration is more correct, which map #184
    puts out of scope as upstream work.

    Read out of :mod:`components.predictor_context` itself rather than transcribed, so a
    predictor added to either criterion fails here instead of quietly falsifying the copy. That
    is not hypothetical: ticket #203 is considering wiring ``AlphaMissense_score`` into a
    predictor table, and #190 has already had to add a column-spelling resolver.
    """
    from components.predictor_context import _cbp10_readings, _pp3_bp4_readings

    empty = pd.Series({}, dtype=object)
    scored = {col for _, col, _, _, _ in CLINGEN_SVI_TABLE if col is not None}

    for tool, readings in (
        ("InterVar", _pp3_bp4_readings(empty, frozenset())),
        ("CancerVar", _cbp10_readings(empty, frozenset())),
    ):
        read = {r.column for r in readings}
        assert not (read & scored), (
            f"{tool} reads {sorted(read & scored)}, which this table also scores — "
            "the label's 'read none of these predictors' is now false"
        )


def test_the_label_names_exactly_the_inputs_the_classifier_actually_reads():
    """The sentence listing the tool's own predictors cannot fall behind the section drawing them.

    A stale list here would be the quieter half of the same defect: the disjointness claim could
    still hold while the label pointed the reader at four inputs out of five.
    """
    from components.predictor_context import _cbp10_readings, _pp3_bp4_readings

    empty = pd.Series({}, dtype=object)
    for tool, readings in (
        ("InterVar", _pp3_bp4_readings(empty, frozenset())),
        ("CancerVar", _cbp10_readings(empty, frozenset())),
    ):
        assert CLASSIFIER_INPUTS[tool] == tuple((r.name, r.column) for r in readings), tool


@pytest.mark.parametrize(
    "tool,verdict",
    [("InterVar", "classification"), ("CancerVar", "tier")],
)
def test_the_label_leads_with_the_fact_that_nothing_here_fed_the_verdict(
    monkeypatch, tool, verdict
):
    """Provenance first — issue #191's decision, and the one framing that takes no side."""
    drawn, _ = _draw(monkeypatch, _row(), tool=tool)
    captions = [t for t in drawn if "<table" not in t]

    assert captions, "the section drew no label at all"
    first = captions[0]
    assert tool in first
    assert "read none of these predictors" in first
    assert verdict in first, first


def test_the_label_names_the_calibration_and_its_source():
    """A reader has to be able to look the thresholds up; the scale is named and dated."""
    captions = reference_scales._provenance_captions("InterVar")
    joined = " ".join(captions)

    assert "ClinGen SVI" in joined
    assert "Pejaver" in joined and "2022" in joined
    assert "PP3" in joined and "BP4" in joined


def test_on_a_somatic_variant_the_label_says_the_calibration_is_germline():
    """Ticket constraint: the copy has to survive being read on a somatic variant.

    ClinGen SVI calibrated PP3/BP4 for ACMG/AMP germline classification and this table draws
    unchanged on somatic files, where the guideline in force is AMP/ASCO/CAP. Said outright, not
    left for the reader to infer from a heading that says ACMG on a page that says AMP/ASCO/CAP
    one section above.
    """
    joined = " ".join(reference_scales._provenance_captions("CancerVar"))

    assert "germline" in joined
    assert "AMP/ASCO/CAP" in joined
    assert "different guideline" in joined


def test_the_germline_label_does_not_carry_the_somatic_guideline_sentence():
    """The pair with the last one: an unconditional sentence would pass it and be noise here."""
    joined = " ".join(reference_scales._provenance_captions("InterVar"))

    assert "AMP/ASCO/CAP" not in joined


def test_with_no_classifier_on_the_row_no_tool_is_named():
    """Naming a tool that drew nothing would be a false sentence, not a harmless one.

    Reached when the row carries neither classifier: 15 of 188 real MAFs are in that state. The
    scale still introduces itself — it is on screen and has to say what it is — but it has no
    verdict to disclaim against.
    """
    joined = " ".join(reference_scales._provenance_captions(None))

    assert "InterVar" not in joined and "CancerVar" not in joined
    assert "ClinGen SVI" in joined


# ---------------------------------------------------------------------------
# Where it sits
# ---------------------------------------------------------------------------


def test_the_section_is_collapsed(monkeypatch):
    """Empty on ~70% of variants, and never the page's answer — so it opens closed."""
    _, opened = _draw(monkeypatch, _row())

    assert opened["expanded"] is False


def test_the_collapsed_label_says_whether_there_is_anything_inside(monkeypatch):
    """A closed expander that cannot say whether it is empty makes the reader open it
    to find out."""
    _, with_score = _draw(monkeypatch, _row())
    _, without = _draw(monkeypatch, _row(REVEL_score="."))

    assert "(1)" in with_score["label"], with_score["label"]
    assert "none for this variant" in without["label"], without["label"]


def test_the_reference_scale_is_drawn_after_both_evidence_sections():
    """Position, asserted on the panel's own source order.

    Issue #191's placement decision is not decoration: a table directly beneath a verdict reads
    as that verdict's evidence whatever its heading says, and these two point opposite ways on
    810 germline and 2,632 somatic real rows. Asserted as *order* rather than mere presence,
    because moving the call back above the evidence sections would leave every copy guard here
    passing.
    """
    import inspect

    from components import variant_detail

    source = inspect.getsource(variant_detail.render_variant_detail_panel)
    positions = {
        name: source.index(name)
        for name in (
            "_render_guideline_classifications",
            "render_acmg_evidence",
            "render_cbp_evidence",
            "render_reference_scale",
        )
    }

    assert positions["_render_guideline_classifications"] < positions["render_acmg_evidence"]
    assert positions["render_acmg_evidence"] < positions["render_reference_scale"]
    assert positions["render_cbp_evidence"] < positions["render_reference_scale"]


def test_the_guideline_section_no_longer_draws_the_clingen_table():
    """The move is a move, not a copy — the old heading is gone from the verdict's function."""
    import inspect

    from components import variant_detail

    source = inspect.getsource(variant_detail._render_guideline_classifications)

    assert "In Silico Predictors" not in source
    assert "REVEL" not in source


def test_the_table_draws_on_a_somatic_row_too(monkeypatch):
    """Not arm-gated: it is neither arm's scale, so the tool only changes the words."""
    drawn, opened = _draw(monkeypatch, _row(), tool="CancerVar")

    assert [t for t in drawn if "<table" in t], "the table did not draw on the somatic arm"
    assert "(1)" in opened["label"]


# ---------------------------------------------------------------------------
# The empty state
# ---------------------------------------------------------------------------


def test_the_empty_state_names_the_variant_rather_than_claiming_no_data(monkeypatch):
    """*"No in silico predictions available"* read as a gap in knowledge. It is usually not one.

    On 99.9% of empty germline rows and 93.6% of empty somatic ones the file carries the columns
    and holds no value on that row — 82-85% of them introns and silent changes. So the caption
    names what the file carries and what this variant is, and stops there.
    """
    row = _row(REVEL_score=".", Variant_Classification="Intron")
    drawn, _ = _draw(monkeypatch, row)
    captions = " ".join(drawn)

    assert "REVEL" in captions
    assert "Intron" in captions
    assert "No in silico predictions available" not in captions


def test_the_empty_state_distinguishes_a_file_that_carries_nothing(monkeypatch):
    """A different absence, and a different sentence: nothing about this variant explains it."""
    row = pd.Series({"Hugo_Symbol": "BRCA1", "Variant_Classification": "Intron"})
    drawn, _ = _draw(monkeypatch, row)
    captions = " ".join(drawn)

    assert "carries none of the ClinGen SVI predictors" in captions
    assert "Intron" not in captions.split("carries none")[1]


def test_the_empty_caption_survives_a_row_with_no_classification(monkeypatch):
    """``Variant_Classification`` is not guaranteed; the sentence drops the clause,
    not the caption."""
    row = _row(REVEL_score=".")
    del row["Variant_Classification"]
    drawn, _ = _draw(monkeypatch, row)
    captions = " ".join(t for t in drawn if "<table" not in t)

    assert "records no value for this variant" in captions


def test_an_empty_section_still_says_what_the_scale_is(monkeypatch):
    """Opening a collapsed section to find an unexplained apology is the worst of both."""
    drawn, _ = _draw(monkeypatch, _row(REVEL_score="."))

    assert any("ClinGen SVI" in t for t in drawn)


# ---------------------------------------------------------------------------
# The arithmetic, which moved unchanged
# ---------------------------------------------------------------------------


def test_a_score_is_placed_on_the_clingen_scale(monkeypatch):
    """REVEL 0.70 is PP3 Supporting (>=0.644) and the table says so."""
    drawn, _ = _draw(monkeypatch, _row())
    table = "\n".join(t for t in drawn if "<table" in t)

    assert "REVEL" in table
    assert "Supporting" in table


def test_the_table_headings_say_clingen_not_acmg_amp(monkeypatch):
    """The strength columns are ClinGen's units, and the heading no longer implies otherwise.

    They used to read *ACMG/AMP Benign (BP4)* / *ACMG/AMP Pathogenic (PP3)*, directly under the
    verdict — the exact pair of words the guideline row above uses for the classification the
    tool actually made. Half the germline and 63% of the somatic rows that reach PP3 here reach
    Moderate or Strong, strengths InterVar's PP3 — only ever *Supporting* — cannot express, so
    the same heading over both was the invitation to read them as one scale.
    """
    drawn, _ = _draw(monkeypatch, _row())
    table = "\n".join(t for t in drawn if "<table" in t)

    assert "ClinGen Pathogenic (PP3)" in table
    assert "ClinGen Benign (BP4)" in table
    assert "ACMG/AMP Pathogenic" not in table
