"""The ClinGen SVI PP3/BP4 reference scale — a separate calibration, and not the page's answer.

REVEL, BayesDel (no AF), VEST4 and CADD, scored against the ClinGen Sequence Variant
Interpretation thresholds (Pejaver et al. 2022, *Am J Hum Genet* 109:2163) for ACMG/AMP PP3 and
BP4. Map #184 keeps this table and issue #191 places and labels it: it lives here, in its own
module and its own collapsed section *below* the classifier's evidence, rather than beneath the
guideline verdict where it used to sit.

Why it is not the verdict's evidence (issue #191)
-------------------------------------------------
The move and the label answer one hazard: this table and the classifier above it both appear to
answer *"do the computational predictors support pathogenicity"*, and they answer it with
**different tools** and **different cutoffs**. Drawn directly under the verdict, it read as the
reasoning behind it, and a clinician could reasonably add the two together.

Three measurements say why that reading is wrong, all over the real corpus — 188 MAFs, 211,634
germline and 116,417 somatic rows (``docs/wayfinder/issue-191/``):

* **Nothing here fed the verdict.** The four columns this table reads are **disjoint** from the
  nine either classifier reads. InterVar's PP3/BP4 read ``MetaSVM_score``, ``GERP++_RS`` and the
  two ``dbscSNV_*`` scores; CancerVar's CBP10 reads ``SIFT_score``, three ``*_pred`` letters and
  ``GERP++_RS``. Not one of REVEL, BayesDel, VEST4 or CADD appears in either list. This is the
  fact the label leads with, because it is checkable and it takes no side: it says the table did
  not contribute, never that the classifier's calibration is worse.
* **They point opposite ways often enough to matter.** Against InterVar's own PP3/BP4: 810 rows
  where this table says benign and PP3 fired, 24 the other way. Against CancerVar's CBP10: 2,632
  rows where this table says benign and CBP10 called it damaging, 59 the other way — **13.0% of
  the somatic rows where both take a side**. A reader adding them up is not resolving a close
  call; they are averaging two different questions.
* **It speaks strengths the classifier cannot.** Of the rows reaching PP3 here, half the germline
  and 63% of the somatic ones reach Moderate or Strong. InterVar's PP3 is only ever *Supporting*,
  so those cells are not a louder version of the criterion above — they are a different scale's
  units, and the label must not let them be read as the same.

It is a **germline** calibration on both arms. ClinGen SVI calibrated PP3/BP4 for ACMG/AMP
germline classification, and this table draws unchanged on somatic files, where the guideline in
force is AMP/ASCO/CAP. The somatic copy says so outright rather than leaving the reader to notice
the heading says ACMG on a page that says AMP/ASCO/CAP one section above.

Empty far more often than not
-----------------------------
The table draws **no rows at all** on 68.2% of germline and 72.1% of somatic variants — which is
the other half of why it is collapsed: on most variants it would be a heading over an apology.
The emptiness is **variant-level, not a missing file**: 144,172 of 144,253 empty germline rows
(99.9%) and 78,577 of 83,939 somatic ones (93.6%) are files that *carry* the columns and hold no
value on that row, and 82–85% of those are introns and silent changes. Where the columns are
there, the table scores 93.8% of germline and 85.9% of somatic missense rows.

So the empty state names which of the two it is, and never says *"no predictions available"* flat:
that reads as a gap in knowledge, when the usual truth is that the variant is outside what these
predictors score. It stops at what was measured — the columns the file carries and the
classification the row records — and leaves the inference to the reader, because CADD is scored
genome-wide while REVEL and VEST4 are missense-only, so a single sentence about *why* the cells
are blank would be false for one of them.

Contamination
-------------
A column ANNOVAR did not supply holds the **chromosome** (issue #194). ``CADD_phred`` is
contaminated on 38 real files. Such a row is drawn with its value and a warning, and takes no
place in the tally's numerator or denominator — a chromosome must not vote on pathogenicity. The
check is per column per file and arrives via ``untrustworthy``, measured at load time.
"""

from typing import Optional

import pandas as pd
import streamlit as st

from config.contaminated_columns import held_value
from config.missing_values import says_nothing

# ACMG PP3/BP4 in silico predictor thresholds — ClinGen SVI calibration.
#
# Source: Pejaver et al. 2022 (Am J Hum Genet 109:2163-2177), transcribed from
# the machine-readable reproduction in Nelson et al. 2025 (medRxiv
# 2025.07.29.25331316, Table 1). ClinGen recommends designating ONE predictor
# in advance (one that reaches Strong for pathogenicity / Moderate for
# benignity): BayesDel, MutPred2, REVEL or VEST4; CADD is included for context.
# Weaker Pejaver tools (GERP++, phyloP, FATHMM, SIFT, PolyPhen-2, PrimateAI...)
# are intentionally omitted as inappropriate single designated tools.
#
# That omission list is worth reading beside the label: four of the tools ClinGen judged too
# weak to designate — GERP++, FATHMM, SIFT, PolyPhen-2 — are exactly the ones CancerVar's CBP10
# counts, and GERP++ is one of InterVar's three. The two scales disagree about which predictors
# are fit for the question, which is why this table's answer is not a second opinion on the
# classifier's. It is not grounds for calling either wrong, and the label does not.

#: ClinGen SVI (Pejaver 2022) PP3/BP4 thresholds.
#:
#: Each entry: ``(display_name, dbnsfp_column, direction, pp3, bp4)``.
#: ``direction "higher"`` -> larger score = more pathogenic: PP3 met when score >= cutoff, BP4
#: met when score <= cutoff; ``"lower"`` reverses both. ``pp3``/``bp4`` map strength -> cutoff;
#: a strength absent from the map was not reached in the calibration.
#:
#: Named for the calibration rather than for the panel (it was ``_PREDICTOR_TABLE_CONFIG`` in
#: ``variant_detail``): since issue #190 the panel draws **three** predictor tables, and only
#: this one is ClinGen's.
CLINGEN_SVI_TABLE = [
    ("REVEL", "REVEL_score", "higher",
     {"strong": 0.932, "moderate": 0.773, "supporting": 0.644},
     {"very_strong": 0.003, "strong": 0.016, "moderate": 0.183, "supporting": 0.290}),
    ("BayesDel (no AF)", "BayesDel_noAF_score", "higher",
     {"strong": 0.50, "moderate": 0.270, "supporting": 0.130},
     {"moderate": -0.360, "supporting": -0.180}),
    ("VEST4", "VEST4_score", "higher",
     {"strong": 0.965, "moderate": 0.861, "supporting": 0.764},
     {"moderate": 0.302, "supporting": 0.449}),
    # dbNSFP ships MutPred v1 as 'MutPred_score'; MutPred2 (the calibrated tool)
    # requires separate scoring, so no dbNSFP column is wired by default.
    ("MutPred2", None, "higher",
     {"strong": 0.932, "moderate": 0.829, "supporting": 0.737},
     {"strong": 0.010, "moderate": 0.197, "supporting": 0.391}),
    ("CADD", "CADD_phred", "higher",
     {"moderate": 28.1, "supporting": 25.3},
     {"strong": 0.15, "moderate": 17.3, "supporting": 22.7}),
]

# Not calibrated here — populate from Stenton/Nadeau et al. 2024 (Genet Med
# 26:101213) before use. All three reach PP3 Strong / >=BP4 Moderate under
# calibration, but the AlphaMissense 0.564 and ESM1b -7.5 defaults do NOT even
# reach Supporting, so they are intentionally left out until verified:
#   AlphaMissense (AlphaMissense_score), ESM1b (ESM1b_score), VARITY_R (VARITY_R_score)

_STRENGTH_DISPLAY = {
    "supporting": "Supporting",
    "moderate": "Moderate",
    "strong": "Strong",
    "very_strong": "Very Strong",
}
_PP3_ORDER = ["very_strong", "strong", "moderate", "supporting"]  # strongest first
_BP4_ORDER = ["very_strong", "strong", "moderate", "supporting"]

#: What the two strength cells of a warned row hold. An em dash rather than an empty cell,
#: which in this table means *this predictor reached no threshold* — a claim about the score,
#: and one this row is in no position to make.
_UNSCORED = "—"

#: The inputs each classifier's own criterion reads — ``(display name, MAF column)``, spelled and
#: ordered exactly as :mod:`components.predictor_context` draws them.
#:
#: Held here only so the label can name them. ``tests/test_reference_scales.py`` asserts this is
#: **exactly** what that module reads, both ways: a predictor added there cannot leave this
#: sentence quietly naming four inputs out of five, and — the claim that matters — it cannot
#: silently overlap :data:`CLINGEN_SVI_TABLE` and make *"read none of these predictors"* false.
#: Columns are carried beside the names for that second check, which is about columns.
CLASSIFIER_INPUTS = {
    "InterVar": (
        ("MetaSVM", "MetaSVM_score"),
        ("GERP++ RS", "GERP++_RS"),
        ("dbscSNV ADA", "dbscSNV_ADA_SCORE"),
        ("dbscSNV RF", "dbscSNV_RF_SCORE"),
    ),
    "CancerVar": (
        ("SIFT", "SIFT_score"),
        ("PolyPhen-2 HDIV", "Polyphen2_HDIV_pred"),
        ("FATHMM", "FATHMM_pred"),
        ("MutationAssessor", "MutationAssessor_pred"),
        ("GERP++ RS", "GERP++_RS"),
    ),
}

#: Which criterion of each tool the disjointness claim is about, for the label's last sentence.
_CRITERION_OF = {"InterVar": "PP3 and BP4", "CancerVar": "CBP10"}


def classify_score(score: float, direction: str, pp3: dict, bp4: dict) -> tuple:
    """Classify a predictor score into PP3/BP4 evidence strength.

    Args:
        score: the predictor's value for this variant.
        direction: ``"higher"`` when a larger score is more pathogenic, else ``"lower"``.
        pp3: ``{strength: cutoff}`` for the pathogenic end.
        bp4: ``{strength: cutoff}`` for the benign end.

    Returns:
        tuple: ``(pp3_strength, bp4_strength)`` — at most one is non-None. Strengths are the
        lowercase keys (``'supporting'``..``'very_strong'``); see :data:`_STRENGTH_DISPLAY`.
        PP3 (pathogenic) takes precedence over BP4 (benign).
    """
    higher = direction == "higher"

    # PP3 (pathogenic end), strongest first
    for strength in _PP3_ORDER:
        cutoff = pp3.get(strength)
        if cutoff is None:
            continue
        if (higher and score >= cutoff) or (not higher and score <= cutoff):
            return strength, None

    # BP4 (benign end) only if no PP3 match
    for strength in _BP4_ORDER:
        cutoff = bp4.get(strength)
        if cutoff is None:
            continue
        if (higher and score <= cutoff) or (not higher and score >= cutoff):
            return None, strength

    return None, None


def parse_score(row: pd.Series, col: str) -> Optional[float]:
    """Parse a numeric score from a MAF row, handling ';'-separated and missing values.

    Splits on ``;`` for the dbNSFP multi-transcript convention. Issue #186 read 210,792 rows
    cell-by-cell and found **no** multi-valued cell on disk, so the split is defensive and has
    never fired; it is kept because the convention is real, not because the corpus needs it.
    """
    if col not in row.index:
        return None
    val = row[col]
    if says_nothing(val):
        return None
    val_str = str(val).strip()
    val_str = val_str.split(";")[0].strip()
    try:
        return float(val_str)
    except (ValueError, TypeError):
        return None


def _untrustworthy_note(name: str, column: str) -> str:
    """The sentence beneath the table explaining why one row did not score.

    Args:
        name: the predictor's display name, as the table's first column spells it.
        column: the MAF column it was read from.

    Returns:
        str: one caption, naming both — the name so the reader can find the row, and the column
        so they can check the claim against their own file.

    Says what was measured and no more. It does not say *why* the column holds a chromosome:
    that explanation is the annotator resolving the name positionally, which is true of the
    corpus issue #186 measured but is a fact about the annotator's release rather than about
    this file, and the check that fired here verified none of it.
    """
    return (
        f"⚠ {name} is excluded: this file's {column} column holds the chromosome, "
        f"not a {name} score."
    )


def _provenance_captions(tool: Optional[str]) -> list:
    """The label: why these scores are here and what they are not.

    Args:
        tool: ``"InterVar"``, ``"CancerVar"``, or None when the row carries neither classifier.

    Returns:
        list: the caption lines, in order, leading with **provenance** — that the classifier read
        none of these predictors and none of these scores fed its verdict.

    Leading with provenance is issue #191's decision, and it is the only one of the three true
    framings that takes no side. *"A different calibration"* invites the reader to pick the
    better one; *"the correct thresholds"* is a claim this repo has no authority to make, and
    map #184 rules disagreeing with a classifier's thresholds out of scope as upstream work.
    That the four columns fed nothing is neither — it is a fact about wiring, and it is the fact
    that blocks the error the placement was changed to prevent.
    """
    scale = (
        "A separate published calibration: ClinGen SVI (Pejaver et al. 2022) thresholds "
        "for ACMG/AMP PP3 and BP4."
    )

    if tool is None:
        # No classifier drew above, so there is no verdict to disclaim against and naming one
        # would be false. The scale still has to introduce itself.
        return [scale, "No guideline classifier scored this variant, so nothing here is compared."]

    names = [name for name, _ in CLASSIFIER_INPUTS[tool]]
    read = f"{', '.join(names[:-1])} and {names[-1]}"
    criterion = _CRITERION_OF[tool]
    verdict = "tier" if tool == "CancerVar" else "classification"

    captions = [
        f"{tool} read none of these predictors. Nothing here contributed to the {verdict} above.",
        scale,
    ]
    if tool == "CancerVar":
        # The wrong-arm fact, said outright. Second, not first: issue #191 chose provenance-led
        # copy, and on a somatic file the guideline mismatch is the *reason* the calibration is
        # foreign, not a separate complaint.
        captions.append(
            "Those thresholds are calibrated for **germline** ACMG/AMP classification. This "
            "variant was classified under AMP/ASCO/CAP, a different guideline."
        )
    captions.append(
        f"{tool}'s own {criterion} read {read} — shown under the evidence above."
    )
    return captions


def _empty_caption(row: pd.Series, carried: list) -> str:
    """What the section says when no row drew.

    Args:
        row: the variant.
        carried: the wired columns this file actually has.

    Returns:
        str: one sentence, naming which of the two absences this is.

    Never *"no in silico predictions available"*, which is what this said before issue #191 and
    reads as a gap in knowledge. On 99.9% of empty germline rows and 93.6% of empty somatic ones
    the file carries the columns and simply holds no value here — most often an intron or a
    silent change. It stops at the measurement rather than explaining it: CADD is scored
    genome-wide while REVEL and VEST4 are missense-only, so one sentence about *why* the cells
    are blank would be false for one of them.
    """
    if not carried:
        return (
            "This file carries none of the ClinGen SVI predictors "
            "(REVEL, BayesDel (no AF), VEST4, CADD)."
        )

    names = ", ".join(carried)
    classification = row.get("Variant_Classification")
    if says_nothing(classification):
        return f"This file carries {names}, but records no value for this variant."
    return (
        f"This file carries {names}, but records no value on this "
        f"{str(classification).strip()} variant."
    )


def _table_html(table_rows: list) -> str:
    """``Source | Score | ACMG/AMP Benign (BP4) | ACMG/AMP Pathogenic (PP3)``.

    Explicit dark cell text for the same reason the other two predictor tables carry it: the
    light row backgrounds would render white-on-white in Streamlit's dark theme.
    """
    th = (
        "background-color:#4472C4;color:white;padding:8px 12px;"
        "text-align:center;font-weight:600;"
    )
    html_parts = [
        '<table style="width:100%;border-collapse:collapse;font-size:0.9em;'
        'table-layout:fixed;">',
        "<colgroup>",
        '<col style="width:20%;">',
        '<col style="width:15%;">',
        '<col style="width:32.5%;">',
        '<col style="width:32.5%;">',
        "</colgroup>",
        "<thead><tr>",
        f'<th style="{th}">Source</th>',
        f'<th style="{th}">Score</th>',
        f'<th style="{th}">ClinGen Benign (BP4)</th>',
        f'<th style="{th}">ClinGen Pathogenic (PP3)</th>',
        "</tr></thead><tbody>",
    ]

    for i, (src, score_str, bp4_str, pp3_str) in enumerate(table_rows):
        bg = "#f0f4fa" if i % 2 == 0 else "#ffffff"
        cell = "padding:6px 12px;text-align:center;color:#1a1a1a;"
        bp4_color = "color:#4472C4;font-weight:600;" if bp4_str else ""
        pp3_color = "color:#C00000;font-weight:600;" if pp3_str else ""
        html_parts.append(
            f'<tr style="background-color:{bg};">'
            f'<td style="{cell}color:#4472C4;font-weight:500;">{src}</td>'
            f'<td style="{cell}">{score_str}</td>'
            f'<td style="{cell}{bp4_color}">{bp4_str}</td>'
            f'<td style="{cell}{pp3_color}">{pp3_str}</td>'
            f"</tr>"
        )

    html_parts.append("</tbody></table>")
    return "\n".join(html_parts)


def _rows_for(row: pd.Series, untrustworthy: frozenset) -> tuple:
    """Score the wired predictors for one variant.

    Returns:
        tuple: ``(table_rows, n_pp3, n_bp4, n_scored, warnings, carried)``.
    """
    table_rows = []
    warnings = []
    carried = []
    n_pp3 = n_bp4 = n_scored = 0

    for name, score_col, direction, pp3_thresh, bp4_thresh in CLINGEN_SVI_TABLE:
        # Tools without a wired dbNSFP column (e.g. MutPred2) cannot be sourced.
        if score_col is None:
            continue
        if score_col in row.index:
            carried.append(name)

        # A column holding the chromosome instead of a score (issue #194). The value is still
        # shown — the reader is entitled to see what their file holds, and a row that silently
        # vanished would look like a variant no predictor scored. What it does *not* get is a
        # place on the ClinGen scale: a chromosome must not cast a vote on pathogenicity, so
        # both strength cells stay empty and the row is left out of the tally's numerator and
        # denominator alike. `says_nothing` cannot reach this case, because the cell is not
        # blank — 22 is a plausible CADD phred, which is exactly the danger.
        #
        # Asked *before* the value is parsed, and that order is the whole of one bug. On a
        # `chrX` variant the contaminated cell holds the letter `X`, which `float()` refuses —
        # so parsing first drops the row through the `score is None` branch below and the
        # warning is never reached. The column is contaminated on every row of the file, and a
        # panel that warned on the autosomes and went quiet on X and Y would be at its least
        # trustworthy exactly where the reader had least reason to suspect it.
        if score_col in untrustworthy:
            held = held_value(row, score_col)
            if held is None:
                continue
            table_rows.append((name, f"⚠ {held}", _UNSCORED, _UNSCORED))
            warnings.append(_untrustworthy_note(name, score_col))
            continue

        score = parse_score(row, score_col)
        if score is None:
            continue

        score_fmt = f"{score:.3g}" if abs(score) < 100 else f"{score:.1f}"
        pp3, bp4 = classify_score(score, direction, pp3_thresh, bp4_thresh)

        table_rows.append(
            (
                name,
                score_fmt,
                _STRENGTH_DISPLAY.get(bp4, "") if bp4 else "",
                _STRENGTH_DISPLAY.get(pp3, "") if pp3 else "",
            )
        )
        n_scored += 1
        if pp3:
            n_pp3 += 1
        if bp4:
            n_bp4 += 1

    return table_rows, n_pp3, n_bp4, n_scored, warnings, carried


def render_reference_scale(
    row: pd.Series, tool: Optional[str], untrustworthy: frozenset
) -> None:
    """Draw the ClinGen SVI PP3/BP4 scores for one variant, in their own collapsed section.

    Args:
        row: the variant.
        tool: the classifier whose verdict the panel drew above — ``"InterVar"``,
            ``"CancerVar"``, or None when the row carries neither. Names the tool in the label;
            it does **not** gate the table, which is neither arm's scale and draws on both.
        untrustworthy: the columns this file holds the chromosome in (issue #194).

    Collapsed by default, and below both evidence sections rather than beneath the verdict.
    Issue #191 settled both: adjacency was doing the arguing, and no label under the verdict can
    outrun the reading that the table directly beneath a classification is its evidence. The
    count rides in the expander label so a collapsed section still says whether there is
    anything inside — necessary when it is empty on ~70% of variants.
    """
    table_rows, n_pp3, n_bp4, n_scored, warnings, carried = _rows_for(row, untrustworthy)

    count = len(table_rows)
    label = (
        f"Other scales — ClinGen SVI PP3/BP4 reference scores ({count})"
        if count
        else "Other scales — ClinGen SVI PP3/BP4 reference scores (none for this variant)"
    )

    with st.expander(label, expanded=False):
        for caption in _provenance_captions(tool):
            st.caption(caption)

        if not table_rows:
            st.caption(_empty_caption(row, carried))
            return

        st.markdown(_table_html(table_rows), unsafe_allow_html=True)

        # `n_scored`, not `len(table_rows)`: a warned row is on screen but did not score, and a
        # denominator that counted it would report a larger evidence base than actually spoke.
        # Skipped altogether when nothing scored, rather than printing `0/0 predictors` — the
        # warnings below say what happened, and a tally over an empty base says nothing.
        if n_scored:
            st.caption(
                f"PP3 (pathogenic): {n_pp3}/{n_scored} predictors "
                f"| BP4 (benign): {n_bp4}/{n_scored} predictors"
            )

        for warning in warnings:
            st.caption(warning)
