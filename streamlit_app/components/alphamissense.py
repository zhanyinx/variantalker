"""AlphaMissense: an in-silico missense predictor, in its own section and on its own scale.

Issue #203. AlphaMissense used to draw **twice** — a brown badge in the clinical row and a
red/green/amber one in the guideline row, the same value in two colours inches apart. It is neither
a database opinion nor a guideline verdict, so the dev gave it neither badge and its own one-row
section instead, in issue #190's landed ``Predictor | This file's value | How AlphaMissense read
it`` shape.

Beside the ClinGen SVI table, not inside it
--------------------------------------------
A fifth row in :data:`components.reference_scales.CLINGEN_SVI_TABLE` was possible and was
rejected on what the prototype **rendered** rather than on an argument. That table is calibrated
strictly to ClinGen SVI (Pejaver et al. 2022) and says so in its own label, so an uncalibrated row
makes the section's copy false: with no thresholds the tally reports *"BP4 (benign): 0/1
predictors"* — asserting AlphaMissense *failed* a ClinGen threshold that does not exist for it —
and with AlphaMissense's own cutoffs wearing ClinGen's strength names it reports *"1/1"* and prints
*Supporting*, indistinguishable on screen from a calibrated row.

A published AlphaMissense calibration does exist — Stenton/Nadeau et al. 2024, named in
``reference_scales`` and deliberately untranscribed — so a calibrated row is possible rather than
impossible. Map #199 rules transcribing a second calibration **out of scope**: it is new
capability, not drawing an annotation once.

Not a PP3/BP4 surface, and it says so
--------------------------------------
Issue #201 refused RENOVO a predictor table partly because all three the panel drew were **PP3/BP4**
surfaces — a classifier's own criterion inputs (:mod:`components.predictor_context`) or ClinGen's
calibration (:mod:`components.reference_scales`) — and RENOVO has no such calibration. This section
is a fourth that is not one, and both calls stand: the reconciliation is that this table names
**AlphaMissense's own three-class boundaries** and never an ACMG evidence strength, which is what
the second caption below says outright. What RENOVO was refused was a seat among *calibrated*
predictors; what this has is a section that claims no calibration at all.

It reads ``am_pathogenicity``, and that is not the obvious column
------------------------------------------------------------------
The ClinGen table addresses every row by dbNSFP column name, and ``AlphaMissense_score`` **is** in
the MAFs — 53 germline and 51 somatic files. Wiring it would have been the natural move and would
have been wrong: the two columns are **different annotations**. Of 55,353 rows carrying both, only
6.3% agree to 1e-9, 33.9% differ beyond rounding, and **346 rows across 84 files fall in different
AlphaMissense classes** — one real germline variant holds 0.3475 in one column and 0.287 in the
other, *ambiguous* beside *likely benign*. ``am_pathogenicity`` also covers strictly more files —
55 germline and **84** somatic — so the dbNSFP spelling would have gone silent on 15,340 variants.

The band is named from the score, not echoed from ``am_class``
---------------------------------------------------------------
``am_class`` is an **exact function** of ``am_pathogenicity`` at :data:`AM_LIKELY_BENIGN_BELOW` and
:data:`AM_LIKELY_PATHOGENIC_FROM`, and the score is present on every row the class is — 0
exceptions on either arm — so the label carries no information the score lacks. It carries one
liability: **100 of 139 files** spell it ``benign``/``ambiguous``/``pathogenic`` for a tool that
publishes only the ``likely_`` form. Naming the band from the score is not disagreeing with the
annotator, which map #199 rules out of scope — the two vocabularies occupy identical score ranges
and agree on the class — it is declining to reproduce a spelling the tool does not use.
"""

from typing import Optional

import pandas as pd
import streamlit as st

from config.contaminated_columns import held_value
from config.missing_values import says_nothing

from .predictor_context import _Reading, _format_score
from .predictor_context import _table_html as _context_table_html
from .reference_scales import CLASSIFIER_INPUTS

#: The column the section reads. See the module docstring for why it is not ``AlphaMissense_score``.
AM_COLUMN = "am_pathogenicity"

#: AlphaMissense's own published three-class cutoffs, recovered from the corpus: ``am_class`` is an
#: exact function of ``am_pathogenicity`` at these two numbers — ``likely_benign`` tops out at
#: 0.3396 and ``ambiguous`` starts at 0.3401; ``ambiguous`` tops out at 0.5626 and
#: ``likely_pathogenic`` starts at 0.5640.
#:
#: They are **not** ACMG evidence strengths and the captions say so. See the module docstring.
AM_LIKELY_BENIGN_BELOW = 0.34
AM_LIKELY_PATHOGENIC_FROM = 0.564

#: The heading. Says what kind of thing this is in its own words, because there is no scale name
#: to lead with — unlike the ClinGen section, which is named for its calibration.
_HEADING = "**AlphaMissense (in silico missense predictor):**"


def _says(score: float) -> tuple:
    """What AlphaMissense itself makes of this score, in its own vocabulary.

    Returns:
        tuple: ``(sentence, side)`` — ``side`` is :class:`components.predictor_context._Reading`'s,
        used only to grey an unreadable row. A band is never a damaging/benign *vote*: the rows
        carry no colour, for the reason issue #190 gave for its own two tables.
    """
    if score >= AM_LIKELY_PATHOGENIC_FROM:
        return (
            f"likely pathogenic — at or above its own {AM_LIKELY_PATHOGENIC_FROM} cutoff",
            "neither",
        )
    if score < AM_LIKELY_BENIGN_BELOW:
        return f"likely benign — below its own {AM_LIKELY_BENIGN_BELOW} cutoff", "neither"
    return (
        f"ambiguous — between its own {AM_LIKELY_BENIGN_BELOW} and "
        f"{AM_LIKELY_PATHOGENIC_FROM} cutoffs",
        "neither",
    )


def read_by_any_classifier() -> tuple:
    """Which classifiers read an AlphaMissense column, as the panel's own wiring records it.

    Returns:
        tuple: the tool names whose criterion inputs include an AlphaMissense column — empty on
        this tree, which is what the provenance caption asserts.

    Asked of :data:`components.reference_scales.CLASSIFIER_INPUTS` rather than written into the
    sentence, so *"Neither classifier read AlphaMissense"* cannot outlive its own truth. That is
    the defect issue #203 found in ``reference_scales._empty_caption``, whose hand-written
    predictor list goes stale, and ``tests/test_alphamissense.py`` holds the claim to this
    function.
    """
    columns = {AM_COLUMN, "AlphaMissense_score", "am_class", "AlphaMissense_pred"}
    return tuple(
        tool
        for tool, inputs in CLASSIFIER_INPUTS.items()
        if columns.intersection({column for _name, column in inputs})
    )


def _provenance_caption(tool: Optional[str]) -> str:
    """Why this section is here and what it did not contribute to.

    Provenance-led, which is issue #191's framing and the only one of the true framings that takes
    no side: it is a fact about wiring, checkable, and it blocks the misreading that a predictor
    table beneath a verdict is that verdict's evidence.

    A row carrying neither classifier gets a different sentence, because there is no
    classification above for anything to have contributed to and naming one would be false — the
    same split :func:`components.reference_scales._provenance_captions` makes.
    """
    if read_by_any_classifier():
        # Unreachable on this tree, and it must not be a silent branch: if a classifier ever does
        # read AlphaMissense, the sentence below stops being true and this one is what says so.
        read = " and ".join(read_by_any_classifier())
        return f"{read} read AlphaMissense, so this score is one of its inputs."
    if tool is None:
        return (
            "No guideline classifier scored this variant, so nothing here is compared."
        )
    verdict = "tier" if tool == "CancerVar" else "classification"
    return (
        f"{tool} did not read AlphaMissense. Nothing here contributed to the {verdict} above."
    )


def _scale_caption() -> str:
    """The cutoffs named as AlphaMissense's own, and refused as ACMG strengths.

    This sentence is the reconciliation the module docstring describes: it is what keeps a fourth
    predictor section from reading as a fourth PP3/BP4 surface.
    """
    return (
        f"The cutoffs are AlphaMissense's own three-class boundaries "
        f"({AM_LIKELY_BENIGN_BELOW} / {AM_LIKELY_PATHOGENIC_FROM}), not ACMG PP3/BP4 evidence "
        f"strengths. No published calibration onto those strengths is transcribed here."
    )


def alphamissense_reading(row: pd.Series, untrustworthy: frozenset) -> Optional[_Reading]:
    """This variant's AlphaMissense row, or ``None`` when there is nothing to draw.

    Args:
        row: the variant.
        untrustworthy: the columns this file holds the chromosome in (issue #194).

    Returns:
        A one-row reading, or ``None`` when the column is absent, the cell says nothing, or the
        cell cannot be read as a number.

    ``am_pathogenicity`` joined ``config.contaminated_columns.PREDICTOR_COLUMNS`` with this
    section, so #194's detector measures it at load time. It is **clean on every file** — maximum
    chromosome fraction exactly 0.0000 — but a table that scores a column asks rather than
    concludes: a column ANNOVAR did not supply arrives holding the *chromosome*, not a ``.``, so a
    score can look plausible and be nothing.
    """
    if AM_COLUMN in untrustworthy:
        held = held_value(row, AM_COLUMN)
        if held is None:
            return None
        return _Reading(
            "AlphaMissense",
            AM_COLUMN,
            f"⚠ {held}",
            f"this file's {AM_COLUMN} column holds the chromosome, not a score — so this row "
            "is not a reading",
            "chromosome",
        )

    if AM_COLUMN not in row.index:
        return None
    value = row[AM_COLUMN]
    if says_nothing(value):
        return None
    try:
        score = float(str(value).strip())
    except (TypeError, ValueError):
        return None

    says, side = _says(score)
    return _Reading("AlphaMissense", AM_COLUMN, _format_score(score), says, side)


def render_alphamissense(
    row: pd.Series, tool: Optional[str], untrustworthy: frozenset
) -> None:
    """Draw AlphaMissense's own section for one variant.

    Args:
        row: the variant.
        tool: the classifier whose verdict the panel drew above — ``"InterVar"``, ``"CancerVar"``,
            or ``None`` when the row carries neither. Names the tool in the provenance caption; it
            does **not** gate the section, which is neither arm's scale.
        untrustworthy: as :func:`alphamissense_reading`.

    Draws nothing at all when the file carries no readable ``am_pathogenicity`` for this variant.
    The section makes one claim — *this predictor read this variant this way* — and there is no
    such claim to make on a variant it did not score; the columns the file lacks are already
    reported at load time. That is the prototype's behaviour, and it is what the dev chose from.
    """
    reading = alphamissense_reading(row, untrustworthy)
    if reading is None:
        return

    st.markdown("---")
    st.markdown(_HEADING)
    st.markdown(
        _context_table_html([reading], "This file's value", "How AlphaMissense read it"),
        unsafe_allow_html=True,
    )
    st.caption(_provenance_caption(tool))
    st.caption(_scale_caption())
