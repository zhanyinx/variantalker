"""The predictor scores each classifier actually read, beside the verdict it wrote.

Issue #190. Two criteria on the variant page state a computational verdict and show nothing behind
it: InterVar's **PP3/BP4** on a germline file, CancerVar's **CBP10** on a somatic one. This module
draws the columns each one reads, as *this file* holds them.

It re-applies no thresholds. Map #184's standing decision is that the app reports the tool's own
answer; the verdict here is read from the evidence string, exactly as
:mod:`components.acmg_evidence` and :mod:`components.cbp_evidence` read it. What this module adds is
the *description* of each cell — which side of the tool's own cutoff it fell on — and the cutoffs
are transcribed from the vendored scripts under a drift guard
(``tests/test_predictor_cutoff_contract.py``), not chosen here.

The two arms read different tools, and treat a missing one oppositely
--------------------------------------------------------------------
This is the module's central finding, measured over 323,926 real rows in 167 MAFs (109 somatic /
58 germline), and it is why one shared table would have been wrong:

* **CancerVar (CBP10, ``CancerVar.py:1590``)** counts damaging and benign calls over five
  predictors and grades on the counts — ``dam > 2 -> +1``, ``dam > 4 -> +2``, ``ben > 2 -> -1``.
  An input it cannot read increments ``var``, and ``var`` **cannot change the result**: the branch
  that would use it (``if var > 2``) assigns to a misspelled local (``Prep``), which issue #185
  found. So on the somatic arm a missing predictor is genuinely **inert**.
* **InterVar (PP3/BP4, ``Intervar.py:1540``/``:1847``)** does the opposite: an unreadable column
  *supplies* a point. An unreadable ``GERP++_RS`` or ``dbscSNV_*`` gives BP4 its point outright
  ("absent means there are gaps in the multiple alignment, so not conserved"), and an unreadable
  ``MetaSVM_score`` gives **PP3** its first point unconditionally — the ``except ValueError``
  branch loops over ``["synon", "coding-synon"]`` and no ANNOVAR consequence contains the literal
  ``coding-synon``, so the second iteration always fires. Measured: BP4 fired on 97,223 of 211,553
  germline rows, and on **63,222** of those all three of its points came from columns the file does
  not carry. PP3 fired 2,377 times, **1,278** of them with ``MetaSVM_score`` unreadable.

So a blank cell means "this changed nothing" on one arm and "this is what fired the criterion" on
the other, and each table says which.

The numbers in the file are the numbers the tool branched on — with a measured exception
----------------------------------------------------------------------------------------
Transcribing both classifiers and running them over the MAFs' own cells:

* **germline: 211,553 rows, 0 disagreements** on PP3 and 0 on BP4. Every column InterVar reads is
  on all 58 files.
* **somatic: 112,299 comparable rows, 3,584 disagreements**, and they are entirely a
  column-presence story. Grouped by which of CBP10's five inputs a file carries:

  ===================================  =====  =======  ========
  inputs the file carries              files    rows    differ
  ===================================  =====  =======  ========
  all five                                57   89,783         0
  four (no ``Polyphen2_HDIV_pred``)       37   13,468         0
  two (``SIFT_score``, ``GERP++_RS``)      9    7,803     2,939
  none                                     6    1,245       645
  ===================================  =====  =======  ========

  The 37 four-input files disagree on nothing because of the ``var`` typo above — a missing
  ``Polyphen2_HDIV_pred`` costs CancerVar nothing either. Only the 15 files carrying two inputs or
  none cannot account for what CancerVar recorded, and on those the disagreement is total: every
  row where CBP10 fired differs.

:func:`_beyond_what_this_file_holds` names that case, and it is **arithmetic on the counts, not a
re-derivation**: ``+2`` needs five damaging calls and ``+1`` needs more than two, so a file with
fewer readable inputs than the recorded value needs cannot account for it whatever the values are.
Checked against the transcription over 112,299 rows: it flags **3,584 rows, with 0 false alarms and
0 misses** — exactly the set where the file's own cells do not reproduce ``EVS[9]``.

No substitution, ever
---------------------
``Polyphen2_HDIV_pred`` is absent on 52 of 109 somatic files while ``Polyphen2_HDIV_score`` — the
near-neighbour a substituted row would reach for — is on 101. It is not a substitute, and the
corpus says so twice over: of those 52 files, **35 hold the chromosome in
``Polyphen2_HDIV_score``** (issue #194's failure mode), 8 have no score column either, and only 9
have one present and clean. A row that quietly swapped the score in for the prediction would print
a chromosome where PolyPhen-2's call belongs on two thirds of the files that needed it. The row
says the column is not in the file instead.

Why these rows carry no colour
------------------------------
:mod:`components.cbp_evidence` colours a criterion's score cell, because a CBP score *is* an
evidence strength on CancerVar's scale. A per-predictor call is not: five coloured rows would read
as five votes on pathogenicity and invite exactly the adding-up that #191 exists to prevent, and
#188 already settled that a signed axis is the wrong one for this evidence. The rows are plain, and
only the absent and unreadable ones are muted.
"""

from typing import NamedTuple, Optional

import pandas as pd
import streamlit as st

from config.columns import spelled_in
from config.contaminated_columns import held_value
from config.missing_values import says_nothing

from .acmg_evidence import parse_intervar
from .cbp_evidence import parse_cancervar

# ---------------------------------------------------------------------------
# The cutoffs, transcribed from the vendored classifiers
# ---------------------------------------------------------------------------
#
# Every number below is read out of ``resources/CancerVar/CancerVar.py`` or
# ``resources/InterVar/Intervar.py`` and appears on screen attributed to the tool that owns it.
# ``tests/test_predictor_cutoff_contract.py`` re-reads both scripts and fails if any of them moves,
# which is what makes printing a threshold safe: map #184 rule 2 keeps ``resources/`` read-only, so
# the app can only ever be a copy, and a copy needs a guard or it drifts in silence — the same
# relationship ``tests/test_chromosome_rule_contract.py`` keeps with ``bin/``.
#
# The two tools do **not** agree on how to read the one column they share: CancerVar's CBP10 calls
# ``GERP++_RS`` conserved at ``>= 2`` and InterVar's PP3 at ``> 2``. Two constants, therefore, and
# not one — a shared cutoff would be a claim the scripts do not make.

#: ``CancerVar.py:1590`` — ``sift_cutoff``. Below it is damaging; at or above it is tolerated.
CANCERVAR_SIFT_DAMAGING_BELOW = 0.05

#: ``CancerVar.py:1590`` — ``cutoff_conserv``, applied as ``>=``.
CANCERVAR_GERP_CONSERVED_FROM = 2.0

#: ``Intervar.py:1540`` — ``metasvm_cutoff``, applied as ``> 0`` for PP3 and ``< 0`` for BP4.
INTERVAR_METASVM_CUTOFF = 0.0

#: ``Intervar.py:1540`` — ``cutoff_conserv``, applied as ``>`` for PP3 and ``<=`` for BP4.
INTERVAR_GERP_CONSERVED_ABOVE = 2.0

#: ``Intervar.py:1540`` — ``dbscSNV_cutoff``, applied as ``>`` for PP3 and ``<=`` for BP4.
INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE = 0.6


# ---------------------------------------------------------------------------
# What one row of either table holds
# ---------------------------------------------------------------------------

#: What the value cell says when the column is not in the file. Not an em dash: the panel's em
#: dash means *this cell is empty*, and the whole point of these rows is which of two different
#: things an empty cell is.
_ABSENT = "not in this file"

#: What it says when the column is there and the cell is blank — a variant-level absence rather
#: than a file-level one, which is a different fact and gets different words.
_BLANK = "no value for this variant"


class _Reading(NamedTuple):
    """One predictor's cell, and what the classifier that reads it made of it.

    ``side`` is one of ``"damaging"``, ``"benign"``, ``"neither"``, ``"absent"``, ``"blank"`` or
    ``"chromosome"``. Only the first two are counts either classifier acts on; the rest exist so
    the caption can say *why* a row did not count, which differs per arm.
    """

    name: str
    column: str
    value: str
    says: str
    side: str

    @property
    def readable(self) -> bool:
        """Whether the classifier could have read a value here at all."""
        return self.side in ("damaging", "benign", "neither")


def _numeric(row: pd.Series, column: str) -> Optional[float]:
    """The cell as a float, or ``None``.

    Splits on ``;`` for the same reason :func:`components.reference_scales.parse_score`
    does — the dbNSFP multi-value convention — and issue #186 measured **zero** such cells in
    210,792 rows, as did this module's own pass over 323,926. It is defensive and has never fired.
    """
    spelling = spelled_in(row.index, column)
    if spelling is None:
        return None
    value = row[spelling]
    if says_nothing(value):
        return None
    text = str(value).strip().split(";")[0].strip()
    try:
        return float(text)
    except (ValueError, TypeError):
        return None


def _letter(row: pd.Series, column: str) -> Optional[str]:
    """The cell as the classifier compares it: whole, stripped, and **not** ``;``-split.

    CancerVar tests these columns with ``==`` against a single letter, so a multi-valued cell
    matches nothing at all — neither damaging, nor benign, nor the ``.`` that counts as unknown.
    Splitting here would make the app read a value the tool did not. Measured over 60,417 real
    ``_pred`` cells: every one is a single letter from the tool's own vocabulary, and none carries
    a separator.
    """
    spelling = spelled_in(row.index, column)
    if spelling is None:
        return None
    value = row[spelling]
    if says_nothing(value):
        return None
    return str(value).strip()


def _format_score(number: float) -> str:
    """A score, formatted as the ClinGen table formats one, so the panel reads as one panel."""
    return f"{number:.3g}" if abs(number) < 100 else f"{number:.1f}"


def _absent_or_blank(row: pd.Series, column: str) -> str:
    """Which kind of nothing this is: no column, or a column with nothing in it."""
    return _BLANK if _present(row, column) else _ABSENT


def _present(row: pd.Series, column: str) -> bool:
    """Whether the file has this column at all, under any spelling it might carry."""
    return spelled_in(row.index, column) is not None


def _kind_of_nothing(row: pd.Series, column: str) -> str:
    """``"blank"`` when the column is there and the cell is not, ``"absent"`` when it is not there.

    The two are different facts and the captions turn on which: a file that never carried
    ``Polyphen2_HDIV_pred`` is a file-level gap, and a variant with an empty cell in a column the
    file does carry is a variant-level one. Both are unreadable to the classifier, which is why
    :attr:`_Reading.readable` is false for each.
    """
    return "blank" if _present(row, column) else "absent"


# ---------------------------------------------------------------------------
# CancerVar CBP10 — five predictors, counted
# ---------------------------------------------------------------------------

#: The three letter columns CBP10 reads, and what each letter means *to CancerVar*.
#:
#: Read straight off ``check_PreP``'s string comparisons, and they are exhaustive of what the
#: corpus holds: 60,417 non-missing cells across both arms, every one of ``D``/``P``/``B`` for
#: PolyPhen-2, ``D``/``T`` for FATHMM and ``H``/``M``/``L``/``N`` for MutationAssessor. A value
#: outside its column's table is drawn with :data:`_UNRECOGNISED_LETTER`, because CancerVar's
#: ``==`` would match none of its branches either — such a cell is counted by nothing.
#:
#: ``P`` is spelled out per column deliberately: issue #186 found it means *possibly damaging* in
#: PolyPhen-2 and *polymorphism* in MutationTaster, so one shared letter table would be wrong.
_LETTER_MEANINGS = {
    "Polyphen2_HDIV_pred": {
        "D": ("damaging", "probably damaging (D)"),
        "P": ("damaging", "possibly damaging (P) — which CancerVar counts as damaging"),
        "B": ("benign", "benign (B)"),
    },
    "FATHMM_pred": {
        "D": ("damaging", "damaging (D)"),
        "T": ("benign", "tolerated (T)"),
    },
    "MutationAssessor_pred": {
        "H": ("damaging", "high functional impact (H)"),
        "M": ("damaging", "medium functional impact (M) — which CancerVar counts as damaging"),
        "L": ("benign", "low functional impact (L)"),
        "N": ("benign", "neutral (N)"),
    },
}

#: What a letter outside the tool's own vocabulary is said to be. Never seen in the corpus.
_UNRECOGNISED_LETTER = "not one of the values CancerVar recognises, so it counted for neither side"

#: How many damaging or benign calls each recorded CBP10 value needs, from ``CancerVar.py:1590``.
#: ``+1`` needs ``dam > 2``, ``+2`` needs ``dam > 4``, ``-1`` needs ``ben > 2``. The basis of
#: :func:`_beyond_what_this_file_holds`.
_CBP10_CALLS_NEEDED = {2: 5, 1: 3, -1: 3, 0: 0}

#: Where CBP10 sits in CancerVar's evidence vector: zero-based index 9. Named because the ticket
#: that charted this section called it CBP9 — ``EVS[9]`` is the tenth criterion.
CBP10_INDEX = 9


def _cbp10_readings(row: pd.Series, untrustworthy: frozenset) -> list:
    """CBP10's five inputs, in the order ``check_PreP`` reads them."""
    readings = []

    readings.append(_numeric_reading(
        row,
        untrustworthy,
        name="SIFT",
        column="SIFT_score",
        damaging=lambda score: score < CANCERVAR_SIFT_DAMAGING_BELOW,
        damaging_says=f"damaging — below CancerVar's {CANCERVAR_SIFT_DAMAGING_BELOW} cutoff",
        benign_says=(
            f"tolerated — at or above CancerVar's {CANCERVAR_SIFT_DAMAGING_BELOW} cutoff"
        ),
        absent_says="CancerVar counted it for neither side",
    ))

    for column, name in (
        ("Polyphen2_HDIV_pred", "PolyPhen-2 HDIV"),
        ("FATHMM_pred", "FATHMM"),
        ("MutationAssessor_pred", "MutationAssessor"),
    ):
        readings.append(_letter_reading(row, untrustworthy, name=name, column=column))

    readings.append(_numeric_reading(
        row,
        untrustworthy,
        name="GERP++ RS",
        column="GERP++_RS",
        damaging=lambda score: score >= CANCERVAR_GERP_CONSERVED_FROM,
        damaging_says=(
            f"conserved — at or above CancerVar's {CANCERVAR_GERP_CONSERVED_FROM:g} cutoff"
        ),
        benign_says=f"not conserved — below CancerVar's {CANCERVAR_GERP_CONSERVED_FROM:g} cutoff",
        absent_says="CancerVar counted it for neither side",
    ))
    return readings


def _numeric_reading(
    row: pd.Series,
    untrustworthy: frozenset,
    *,
    name: str,
    column: str,
    damaging,
    damaging_says: str,
    benign_says: str,
    absent_says: str,
) -> _Reading:
    """One numeric predictor's row."""
    if column in untrustworthy:
        held = held_value(row, column)
        if held is not None:
            return _Reading(name, column, f"⚠ {held}", _CHROMOSOME_SAYS.format(column=column),
                            "chromosome")
    score = _numeric(row, column)
    if score is None:
        return _Reading(name, column, _absent_or_blank(row, column), absent_says,
                        _kind_of_nothing(row, column))
    if damaging(score):
        return _Reading(name, column, _format_score(score), damaging_says, "damaging")
    return _Reading(name, column, _format_score(score), benign_says, "benign")


def _letter_reading(
    row: pd.Series, untrustworthy: frozenset, *, name: str, column: str
) -> _Reading:
    """One ``_pred`` predictor's row."""
    if column in untrustworthy:
        held = held_value(row, column)
        if held is not None:
            return _Reading(name, column, f"⚠ {held}", _CHROMOSOME_SAYS.format(column=column),
                            "chromosome")
    letter = _letter(row, column)
    if letter is None:
        return _Reading(
            name,
            column,
            _absent_or_blank(row, column),
            "CancerVar counted it for neither side",
            _kind_of_nothing(row, column),
        )
    meaning = _LETTER_MEANINGS[column].get(letter)
    if meaning is None:
        return _Reading(name, column, letter, _UNRECOGNISED_LETTER, "neither")
    side, says = meaning
    return _Reading(name, column, letter, says, side)


#: The sentence a contaminated cell gets, in both tables. Says what was measured and no more —
#: the same restraint :func:`components.reference_scales._untrustworthy_note` documents.
_CHROMOSOME_SAYS = (
    "this file's {column} column holds the chromosome, not a score — so this row is not a reading"
)


def _beyond_what_this_file_holds(cbp10: int, readable: int) -> bool:
    """Whether the recorded CBP10 needs more calls than this file has readable inputs.

    Arithmetic on the counts, never a re-derivation: ``+2`` requires five damaging calls and
    ``+1`` more than two, so a file with fewer readable inputs than that cannot account for the
    value whatever the values are. Verified against a faithful transcription of ``check_PreP`` over
    112,299 real rows — 3,584 flagged, 0 false alarms, 0 misses.
    """
    return readable < _CBP10_CALLS_NEEDED.get(cbp10, 0)


# ---------------------------------------------------------------------------
# InterVar PP3 / BP4 — four columns, three points
# ---------------------------------------------------------------------------

#: The two ``ExonicFunc``/``Func`` tokens ``check_BP4`` looks for, and the one that vetoes them.
#: Transcribed rather than reworded: the veto is why ``nonsynonymous_SNV`` does not satisfy a test
#: for ``"synon"``.
_SYNONYMOUS_TOKENS = ("synon", "coding-synon")
_SYNONYMOUS_VETO = "nonsynon"

#: The consequence columns, in InterVar's own order. It concatenates ``Func.refGene`` and
#: ``ExonicFunc.refGene``; the MAF carries ``ExonicFunc.refGene`` on all 58 germline files and
#: ``Func.refGene`` on only 4, and reading whichever are present gives the same answer, because no
#: ``Func.refGene`` value contains ``synon``. Empirically: transcribing BP4 this way reproduced it
#: on 211,553 of 211,553 rows.
_CONSEQUENCE_COLUMNS = ("Func.refGene", "ExonicFunc.refGene")


def _looks_synonymous(row: pd.Series) -> bool:
    """``check_BP4``'s consequence test, transcribed."""
    parts = [
        str(row[column]).strip()
        for column in _CONSEQUENCE_COLUMNS
        if column in row.index and not says_nothing(row[column])
    ]
    text = " ".join(parts).lower()
    if _SYNONYMOUS_VETO in text:
        return False
    return any(token in text for token in _SYNONYMOUS_TOKENS)


def _metasvm_reading(row: pd.Series, untrustworthy: frozenset) -> _Reading:
    """MetaSVM, whose *absence* is a reading in its own right.

    An unreadable ``MetaSVM_score`` gives PP3 its first point unconditionally (see the module
    docstring), and gives BP4 its first point only on a synonymous consequence. Both are stated,
    because a reader looking at an empty cell beside a fired criterion is owed the reason the cell
    being empty was the reason.
    """
    column = "MetaSVM_score"
    if column in untrustworthy:
        held = held_value(row, column)
        if held is not None:
            return _Reading("MetaSVM", column, f"⚠ {held}",
                            _CHROMOSOME_SAYS.format(column=column), "chromosome")
    score = _numeric(row, column)
    if score is None:
        says = (
            "InterVar counts an unreadable MetaSVM score as one of PP3's three points"
        )
        if _looks_synonymous(row):
            says += ", and — on a synonymous consequence — as one of BP4's three points too"
        return _Reading("MetaSVM", column, _absent_or_blank(row, column), says,
                        _kind_of_nothing(row, column))
    if score > INTERVAR_METASVM_CUTOFF:
        return _Reading(
            "MetaSVM", column, _format_score(score),
            f"above InterVar's {INTERVAR_METASVM_CUTOFF:g} cutoff — one of PP3's three points",
            "damaging",
        )
    if score < INTERVAR_METASVM_CUTOFF:
        return _Reading(
            "MetaSVM", column, _format_score(score),
            f"below InterVar's {INTERVAR_METASVM_CUTOFF:g} cutoff — one of BP4's three points",
            "benign",
        )
    return _Reading(
        "MetaSVM", column, _format_score(score),
        f"exactly InterVar's {INTERVAR_METASVM_CUTOFF:g} cutoff — a point for neither criterion",
        "neither",
    )


def _gerp_reading(row: pd.Series, untrustworthy: frozenset) -> _Reading:
    """GERP++ as InterVar reads it — ``> 2``, where CancerVar uses ``>= 2``."""
    column = "GERP++_RS"
    if column in untrustworthy:
        held = held_value(row, column)
        if held is not None:
            return _Reading("GERP++ RS", column, f"⚠ {held}",
                            _CHROMOSOME_SAYS.format(column=column), "chromosome")
    score = _numeric(row, column)
    if score is None:
        return _Reading(
            "GERP++ RS", column, _absent_or_blank(row, column),
            "InterVar reads an unreadable score as not conserved — one of BP4's three points",
            _kind_of_nothing(row, column),
        )
    if score > INTERVAR_GERP_CONSERVED_ABOVE:
        return _Reading(
            "GERP++ RS", column, _format_score(score),
            f"conserved — above InterVar's {INTERVAR_GERP_CONSERVED_ABOVE:g} cutoff, one of "
            "PP3's three points",
            "damaging",
        )
    return _Reading(
        "GERP++ RS", column, _format_score(score),
        f"not conserved — at or below InterVar's {INTERVAR_GERP_CONSERVED_ABOVE:g} cutoff, one "
        "of BP4's three points",
        "benign",
    )


def _dbscsnv_reading(
    row: pd.Series, untrustworthy: frozenset, *, name: str, column: str
) -> _Reading:
    """One of the two splice scores. They are **one** point, taken together.

    PP3 takes the point if *either* exceeds 0.6; BP4 takes it only if *neither* does, and takes it
    outright if either column is unreadable. The row states its own column's comparison and the
    caption states the pairing, because a row that claimed the point on its own would be wrong
    half the time.
    """
    if column in untrustworthy:
        held = held_value(row, column)
        if held is not None:
            return _Reading(name, column, f"⚠ {held}", _CHROMOSOME_SAYS.format(column=column),
                            "chromosome")
    score = _numeric(row, column)
    if score is None:
        return _Reading(
            name, column, _absent_or_blank(row, column),
            "InterVar reads an unreadable splice score as not splice-altering — which gives BP4 "
            "its third point",
            _kind_of_nothing(row, column),
        )
    if score > INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE:
        return _Reading(
            name, column, _format_score(score),
            f"splice-altering — above InterVar's {INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE:g} "
            "cutoff",
            "damaging",
        )
    return _Reading(
        name, column, _format_score(score),
        f"not splice-altering — at or below InterVar's "
        f"{INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE:g} cutoff",
        "benign",
    )


def _pp3_bp4_readings(row: pd.Series, untrustworthy: frozenset) -> list:
    """PP3/BP4's four columns, in the order ``check_PP3`` reads them."""
    return [
        _metasvm_reading(row, untrustworthy),
        _gerp_reading(row, untrustworthy),
        _dbscsnv_reading(row, untrustworthy, name="dbscSNV ADA", column="dbscSNV_ADA_SCORE"),
        _dbscsnv_reading(row, untrustworthy, name="dbscSNV RF", column="dbscSNV_RF_SCORE"),
    ]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _table_html(readings: list, value_heading: str, says_heading: str) -> str:
    """``Predictor | <value_heading> | <says_heading>``.

    The same three-column shape as :func:`components.cbp_evidence._cbp_table_html`, with explicit
    dark text for the same reason: the light row backgrounds would render white-on-white in
    Streamlit's dark theme. No colour by side — see the module docstring.
    """
    th = (
        "background-color:#4472C4;color:#fff;padding:6px 10px;"
        "text-align:left;font-weight:600;"
    )
    parts = [
        '<table style="width:100%;border-collapse:collapse;font-size:0.88em;'
        'table-layout:fixed;">',
        '<colgroup><col style="width:22%;"><col style="width:22%;">'
        '<col style="width:56%;"></colgroup>',
        f'<thead><tr><th style="{th}">Predictor</th>'
        f'<th style="{th}">{value_heading}</th>'
        f'<th style="{th}">{says_heading}</th></tr></thead><tbody>',
    ]
    for i, reading in enumerate(readings):
        bg = "#f7f9fc" if i % 2 == 0 else "#ffffff"
        color = "#1a1a1a" if reading.readable else "#999999"
        if reading.side == "chromosome":
            color = "#b8860b"
        cell = f"padding:5px 10px;vertical-align:top;color:{color};"
        parts.append(
            f'<tr style="background-color:{bg};">'
            f'<td style="{cell}font-weight:600;">{reading.name}'
            f'<br><span style="font-weight:400;font-size:0.82em;color:#8a8a8a;">'
            f"{reading.column}</span></td>"
            f'<td style="{cell}">{reading.value}</td>'
            f'<td style="{cell}">{reading.says}</td></tr>'
        )
    parts.append("</tbody></table>")
    return "\n".join(parts)


def _tally(readings: list) -> str:
    """This file's own counts, said as counts and attributed to the file rather than the tool."""
    damaging = sum(1 for r in readings if r.side == "damaging")
    benign = sum(1 for r in readings if r.side == "benign")
    unread = len(readings) - damaging - benign
    parts = [f"**{damaging}** damaging", f"**{benign}** benign"]
    if unread:
        parts.append(f"**{unread}** that CancerVar could not read here")
    return ", ".join(parts)


def render_cbp10_inputs(row: pd.Series, cancervar_result, untrustworthy: frozenset) -> None:
    """The five predictors behind CancerVar's CBP10, as this file holds them.

    Args:
        row: the variant.
        cancervar_result: the file's ``CancerVar and Evidence`` cell — the same value
            :func:`components.cbp_evidence.render_cbp_evidence` is given, so the value named here
            and the criterion drawn one section above cannot come apart.
        untrustworthy: the columns this file holds the chromosome in, from the load-time verdict
            (issue #194). Passed in rather than read from the session, so this module needs no
            session state and the panel keeps the one place that reads it.

    Renders into the current Streamlit container. Draws nothing at all when the evidence string
    does not parse — :func:`components.cbp_evidence.render_cbp_evidence` has already said so
    directly above, and a second sentence saying it again is noise.
    """
    parsed = parse_cancervar(cancervar_result)
    if parsed is None:
        return

    cbp10 = parsed["scores"][CBP10_INDEX]
    readings = _cbp10_readings(row, untrustworthy)
    readable = sum(1 for r in readings if r.readable)

    st.markdown("**The predictors behind CBP10 (CancerVar):**")
    st.markdown(_table_html(readings, "This file's value", "How CancerVar read it"),
                unsafe_allow_html=True)

    st.caption(
        f"CancerVar recorded **CBP10 = {cbp10:+d}** from these five. It grades on the counts — "
        f"more than two damaging is +1, all five is +2, more than two benign is −1 — and the "
        f"cutoffs above are its own, not the guideline's. This file holds {_tally(readings)}."
    )

    if _beyond_what_this_file_holds(cbp10, readable):
        st.caption(
            f"⚠ Only {readable} of the five columns CancerVar reads are in this file, and "
            f"CBP10 = {cbp10:+d} needs {_CBP10_CALLS_NEEDED[cbp10]}. CancerVar saw more than "
            "this file kept, so what is above cannot account for its answer — the answer is "
            "still CancerVar's, read from the evidence string it wrote."
        )


def render_pp3_bp4_inputs(row: pd.Series, intervar_result, untrustworthy: frozenset) -> None:
    """The four columns behind InterVar's PP3 and BP4, as this file holds them.

    Args:
        row: the variant.
        intervar_result: the file's ``InterVar and Evidence`` cell, the same value
            :func:`components.acmg_evidence.render_acmg_evidence` is given.
        untrustworthy: as :func:`render_cbp10_inputs`.

    Drawn on every germline variant, not only where a criterion fired: the four columns are on all
    58 germline files measured, and the numbers explain a criterion that did *not* fire as much as
    one that did. PP3 fires on 1.1% of real rows and BP4 on 46.0%, so gating on them would take the
    section off screen on the majority of variants with nothing to say why.
    """
    parsed = parse_intervar(intervar_result)
    if parsed is None:
        return

    criteria = parsed.get("criteria") or {}
    pp3 = bool(criteria.get("PP3"))
    bp4 = bool(criteria.get("BP4"))
    readings = _pp3_bp4_readings(row, untrustworthy)

    st.markdown("**The predictors behind PP3 and BP4 (InterVar):**")
    st.markdown(_table_html(readings, "This file's value", "How InterVar read it"),
                unsafe_allow_html=True)

    recorded = (
        "**PP3 and BP4**" if pp3 and bp4
        else "**PP3**" if pp3
        else "**BP4**" if bp4
        else "neither PP3 nor BP4"
    )
    st.caption(
        f"InterVar recorded {recorded} for this variant. It counts three points — MetaSVM, "
        f"GERP++, and the two dbscSNV scores together — and takes PP3 on two of them, BP4 only on "
        f"all three. The cutoffs above are InterVar's own."
    )

    from_absence = _points_from_absence(readings, row)
    if bp4 and from_absence:
        st.caption(
            f"⚠ {from_absence} of BP4's three points came from a column this file does not "
            f"carry a readable value in. InterVar reads an unreadable conservation or splice "
            f"score as evidence of no impact, so an empty cell above is part of why BP4 fired "
            f"rather than a gap in the reasoning."
        )
    elif pp3 and not readings[0].readable:
        st.caption(
            "⚠ One of PP3's three points came from this file having no readable MetaSVM score. "
            "InterVar treats an unreadable MetaSVM score as support for a deleterious effect, so "
            "an empty cell above is part of why PP3 fired."
        )


def _points_from_absence(readings: list, row: pd.Series) -> int:
    """How many of BP4's three points an unreadable column supplied.

    The three points are MetaSVM, GERP++, and the dbscSNV pair — so the pair contributes at most
    one, and contributes it whenever *either* of its two columns is unreadable, which is what
    ``check_BP4``'s ``except ValueError`` does.
    """
    metasvm, gerp, ada, rf = readings
    points = 0
    if not metasvm.readable and _looks_synonymous(row):
        points += 1
    if not gerp.readable:
        points += 1
    if not ada.readable or not rf.readable:
        points += 1
    return points
