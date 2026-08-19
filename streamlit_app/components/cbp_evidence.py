"""AMP/ASCO/CAP evidence, as CancerVar reports it.

The somatic counterpart of :mod:`components.acmg_evidence`, built to the same seam:
:func:`parse_cancervar` takes CancerVar's ``CancerVar and Evidence`` string and returns the
twelve CBP scores, and is pure. Only :func:`render_cbp_evidence` touches Streamlit.

The input, on every somatic file measured::

     CancerVar: 10#Tier_II_potential EVS=[2, 2, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2]

Three things this module says that the germline one does not have to, all of them from issue #185's
read of ``resources/CancerVar/CancerVar.py`` and issue #188's measurements over 109,416 real rows:

1. **A criterion's strength is its value, not a fixed tier.** ``_ACMG_CRITERIA`` carries one
   tier per code because an ACMG criterion means one thing when it fires. A CBP criterion
   carries its strength in the number CancerVar assigned it, and the value sets differ per
   criterion — only CBP1/2/3/9/12 reach 2, only CBP8 and CBP10 go negative, CBP4/7/11 cap at 1.
2. **"Met / not met" is the wrong axis.** 59.0% of real rows carry at least one *negative*
   criterion, and 9.9% have no positive one at all, so the table shows what **fired** — meaning
   non-zero, in either direction — and not what was met.
3. **The description knows the value.** ``CBP8 = 1`` and ``CBP8 = -1`` are opposite claims about
   the same evidence source, so one sentence per code cannot serve both.

Who is responsible for what
---------------------------
The twelve **categories** and the four **tier names** are AMP/ASCO/CAP's. Every number — the
per-criterion grading and the sum thresholds — is CancerVar's own operationalisation; the guideline
assigns no weights at all (#185 §6). So nothing here may be captioned "AMP/ASCO/CAP score", and the
tier is presented as CancerVar's answer on CancerVar's scale. The dev confirmed the institute holds
no term table for these tiers, unlike ClinVar's (issue #114).

The tier lives in two places, and they agree
--------------------------------------------
The ``CancerVar`` column holds the tier, and the tier is *also* printed inside the evidence
string. They have different provenance — the column is what
``bin/add_guidelines_escat_to_funcotator.py:194`` extracted with the regex
``[\\d]#(Tier_[\\w\\W]+) E`` — so issue #189 measured them against each other rather than
assuming: **110,064 comparable rows across 108 real somatic files, 0 disagreements**, and no tier
name outside CancerVar's four. That is what licenses this section to badge the tier from the
string on the **16 of the 124** files that carry the evidence vector with no ``CancerVar`` column
at all. (121 of those 124 spell the column ``' CancerVar: CancerVar and Evidence '`` and 3 spell
it ``'CancerVar..CancerVar.and.Evidence'``; both reach here, and both are fixtures on disk.)
Reading a printed value is not re-deriving one, so map #184's *read the tool's answer* stands
either way.
"""

import html
import re
from typing import Optional

import streamlit as st

from config.missing_values import says_nothing

# The markers behind CBP1/CBP2/CBP3, resolved from CancerVar's own table — issue #198.
# The resolution lives in its own module because it reads a file, which nothing else in
# `components/` does; this module keeps its single Streamlit seam.
from .cancervar_markers import MarkerSet

# The germline palette, read rather than re-typed. #188 settled that a CBP score cell reuses
# the ACMG tier-cell colours — warm for clinical significance, cool for the benign side — so a
# clinician reading a germline file and a somatic file reads one colour language. Deriving the
# scores from that map instead of copying its hex values is what stops the two arms drifting
# apart: there is no second definition to update.
from .acmg_evidence import _TIER_CELL_COLORS as _GERMLINE_CELL_COLORS


# ============================================================================
# CancerVar AMP/ASCO/CAP (CBP) evidence
# ============================================================================

#: The twelve CBP criteria: code, AMP/ASCO/CAP category name, the values CancerVar can assign,
#: and what each of those values *says*.
#:
#: ``values`` is per-criterion because the sets genuinely differ, and a renderer that assumed a
#: uniform -1..2 range would offer cells the tool cannot produce. ``says`` is keyed by value —
#: see the module docstring. ``generic`` is the value-independent phrasing, kept as the fallback
#: for a value outside ``values``: CancerVar has never emitted one in 119,194 measured vectors,
#: and an unexpected number must still draw a row rather than an empty cell.
_CBP_CRITERIA = (
    {
        "code": "CBP1",
        "name": "Therapeutic",
        "values": (0, 1, 2),
        "generic": "therapy marker for this variant",
        "says": {
            2: "therapy marker at evidence level A or B, in this sample's tumour type",
            1: "therapy marker at level C or D, or level A/B in a different tumour type",
            0: "no therapeutic marker matched",
        },
    },
    {
        "code": "CBP2",
        "name": "Diagnostic",
        "values": (0, 1, 2),
        "generic": "diagnostic marker for this variant",
        "says": {
            2: "diagnostic marker in this sample's tumour type",
            1: "diagnostic marker in a different tumour type",
            0: "no diagnostic marker matched",
        },
    },
    {
        "code": "CBP3",
        "name": "Prognostic",
        "values": (0, 1, 2),
        "generic": "prognostic marker for this variant",
        "says": {
            2: "prognostic marker in this sample's tumour type",
            1: "prognostic marker in a different tumour type",
            0: "no prognostic marker matched",
        },
    },
    {
        "code": "CBP4",
        "name": "Mutation type",
        "values": (0, 1),
        "generic": "consequence type, in a gene where loss of function is known",
        "says": {
            1: (
                "a damaging consequence in a listed loss-of-function cancer gene, or predicted "
                "to alter splicing"
            ),
            0: "not a damaging consequence in a listed loss-of-function gene",
        },
    },
    {
        "code": "CBP5",
        "name": "Variant frequency",
        "values": (0,),
        "generic": "mosaic variant allele fraction",
        "says": {0: "CancerVar does not evaluate this — it needs manual input"},
        "never_evaluated": True,
    },
    {
        "code": "CBP6",
        "name": "Potential germline",
        "values": (0,),
        "generic": "non-mosaic variant allele fraction",
        "says": {0: "CancerVar does not evaluate this — it needs manual input"},
        "never_evaluated": True,
    },
    {
        "code": "CBP7",
        "name": "Population database",
        "values": (0, 1),
        "generic": "frequency in the population databases",
        "says": {
            1: "absent, or rarer than 0.01%, in the population frequency databases",
            0: "present in the population databases above CancerVar's rarity cutoff",
        },
    },
    {
        "code": "CBP8",
        "name": "Germline database",
        "values": (-1, 0, 1),
        "generic": "ClinVar's germline significance",
        "says": {
            1: "ClinVar calls this pathogenic",
            0: "ClinVar is silent, uncertain, or conflicting",
            -1: "ClinVar calls this benign or likely benign",
        },
    },
    {
        "code": "CBP9",
        "name": "Somatic database",
        "values": (0, 1, 2),
        "generic": "presence in the somatic databases",
        "says": {
            2: "reported in both COSMIC and ICGC",
            1: "reported in COSMIC or ICGC",
            0: "reported in neither COSMIC nor ICGC",
        },
    },
    {
        "code": "CBP10",
        "name": "Predictive software",
        "values": (-1, 0, 1, 2),
        "generic": "CancerVar's predictor consensus",
        # Not "predictors unavailable" for 0. #185 found an upstream typo (``Prep`` for ``PreP``)
        # that leaves CancerVar's "mostly unknown" branch assigning to a dead local, so a 0
        # licenses no claim at all about how many predictors resolved.
        "says": {
            2: "all five of CancerVar's predictors call this damaging",
            1: "most of CancerVar's predictors call this damaging",
            0: "no predictor consensus either way",
            -1: "most of CancerVar's predictors call this benign",
        },
    },
    {
        "code": "CBP11",
        "name": "Pathway involvement",
        "values": (0, 1),
        "generic": "the gene's standing in the cancer gene lists",
        "says": {
            1: "the gene is in a KEGG cancer pathway or the Cancer Gene Census",
            0: "the gene is in neither cancer gene list",
        },
    },
    {
        "code": "CBP12",
        "name": "Publications",
        "values": (0, 1, 2),
        "generic": "published evidence for this variant",
        "says": {
            2: "published evidence in this sample's tumour type",
            1: "published evidence in a different tumour type",
            0: "no publication matched",
        },
    },
)

#: How many scores a CancerVar evidence vector has. Not a magic number in three places: a
#: vector of any other length is a truncated line, and :func:`parse_cancervar` refuses it.
CBP_COUNT = 12

#: CancerVar's own word for each score. Not a bar or a meter: the value sets are four discrete
#: points and they are not uniform across criteria, so a continuous axis would imply cells the
#: tool cannot produce.
_VALUE_WORD = {2: "Strong", 1: "Supporting", 0: "No support", -1: "Benign / neutral"}

#: Score-cell backgrounds, derived from the germline tier-cell palette — see the import above.
_SCORE_COLORS = {
    2: _GERMLINE_CELL_COLORS[("P", "Strong")],
    1: _GERMLINE_CELL_COLORS[("P", "Moderate")],
    -1: _GERMLINE_CELL_COLORS[("B", "Supporting")],
}

#: CancerVar's four tiers, **most specific first**, with the badge colour and the readable name.
#:
#: The order is load-bearing and this is the table issue #187 fixed: ``"Tier_I"`` is a substring
#: of all four names, so a chain that tested it first matched every tier and painted 99.9% of
#: somatic variants the red reserved for maximal clinical significance. Matched by substring
#: rather than equality because real cells carry the whole evidence string where the tier should
#: be, and the tier inside one of those is still the tool's answer.
#:
#: It lives here rather than in :mod:`components.variant_detail` because these are CancerVar's
#: vocabulary and its thresholds — the same rule that put the ACMG criteria in
#: :mod:`components.acmg_evidence`. One definition also means the tier badge this section draws
#: and the guideline-row badge two rows above it cannot draw one tier in two colours.
_TIER_COLORS = (
    ("Tier_IV", "#2ca02c", "Tier IV — benign"),
    ("Tier_III", "#f0c420", "Tier III — uncertain"),
    ("Tier_II", "#ff7f0e", "Tier II — potential"),
    ("Tier_I", "#d62728", "Tier I — strong"),
)

#: What a tier name outside :data:`_TIER_COLORS` is drawn as. Not green: green is this table's
#: *benign* colour, and a value CancerVar's vocabulary does not contain is unknown rather than
#: reassuring.
_UNKNOWN_TIER_COLOR = "#7f7f7f"

#: The sum thresholds, from ``CancerVar.py:958``. Stated on screen so the arithmetic is legible,
#: and attributed to CancerVar because the guideline assigns no weights.
_TIER_THRESHOLDS = "≤2 benign · 3–7 uncertain · 8–10 potential · ≥11 strong"

_EVS_RE = re.compile(r"EVS\s*=\s*\[([^\]]*)\]")
_SUM_TIER_RE = re.compile(r"(-?\d+)\s*#\s*(Tier_\w+)")
_LABEL_RE = re.compile(r"^\s*CancerVar\s*:\s*", flags=re.IGNORECASE)


def tier_color(value: str) -> str:
    """The badge colour for a CancerVar tier, matched most-specific-name-first."""
    for name, color, _label in _TIER_COLORS:
        if name in str(value):
            return color
    return _UNKNOWN_TIER_COLOR


def tier_label(value: str) -> str:
    """A CancerVar tier as a readable phrase, or the raw value if it is not one of the four.

    An unrecognised name is echoed rather than smoothed into one of CancerVar's own: the file
    said something this module does not know, and inventing a tier for it would be the failure
    :data:`_UNKNOWN_TIER_COLOR` exists to colour.
    """
    text = str(value).strip()
    for name, _color, label in _TIER_COLORS:
        if name in text:
            return label
    return text


def parse_cancervar(result) -> Optional[dict]:
    """Parse a CancerVar 'CancerVar and Evidence' value into a structured dict.

    Returns ``{"scores": [12 ints], "sum": int, "printed_sum": int | None, "tier": str}``, or
    ``None`` when there is no parseable content. An already-parsed dict of the same shape is
    returned unchanged, mirroring :func:`components.acmg_evidence.parse_intervar`.

    Tolerates whitespace, an optional leading ``CancerVar:`` label, and the leading space real
    files carry in both the column name and the value.

    **It returns ``None`` rather than a half-answer**, and this is stricter than
    ``parse_intervar`` on purpose. ``parse_intervar`` zeroes an unreadable list entry because an
    ACMG criterion it cannot read is simply one it cannot claim was met. Here every score feeds
    a **sum**, and the sum is the tier — so zeroing one unreadable entry would quietly move the
    variant down CancerVar's scale and badge a tier the tool never assigned. Refused, therefore:

    * no ``EVS=[...]`` section — **73 real cells**, all in one file, holding a tier and a sum
      with no vector behind them (``'1#Tier_IV_benign'``). What is missing is the reasoning this
      section exists to show, so this section draws nothing rather than a tier with no evidence.

      This refusal was justified here on the grounds that *"the guideline row still badges it
      from the ``CancerVar`` column"*, and issue #210 measured that premise **false**: the one
      file holding these 73 cells has no ``CancerVar`` column at all, so the tier CancerVar
      printed is on screen nowhere. The dev re-read the choice with that in hand and kept the
      refusal — a tier with no evidence behind it is still not what this section is for — so what
      changed is the caption above, which had been promising the tier this function declines to
      draw. See :func:`components.variant_detail._verdict_from_evidence`.
    * a non-integer entry, or any vector that is not exactly :data:`CBP_COUNT` long.

    The tier and the printed sum are read from *inside the string* and are **not** required: a
    vector with no ``<sum>#<Tier_>`` prefix still parses, and :func:`render_cbp_evidence` then
    says the string carries no tier rather than banding the sum into one CancerVar never wrote.
    """
    # Already parsed -> don't re-parse.
    if isinstance(result, dict):
        if "scores" in result and "sum" in result:
            return result
        return None

    if says_nothing(result):
        return None
    text = _LABEL_RE.sub("", str(result).strip())

    m = _EVS_RE.search(text)
    if not m:
        return None
    entries = [entry.strip() for entry in m.group(1).split(",") if entry.strip() != ""]
    try:
        scores = [int(entry) for entry in entries]
    except ValueError:
        return None
    if len(scores) != CBP_COUNT:
        return None

    printed_sum = None
    tier = ""
    sm = _SUM_TIER_RE.search(text)
    if sm:
        printed_sum = int(sm.group(1))
        tier = sm.group(2)

    return {
        "scores": scores,
        "sum": sum(scores),
        "printed_sum": printed_sum,
        "tier": tier,
    }


def _criteria_with_scores(parsed: dict) -> list:
    """Pair each criterion with its score, in CancerVar's own order."""
    return list(zip(_CBP_CRITERIA, parsed["scores"]))


def _describe(criterion: dict, score: int) -> str:
    """What this criterion says *at this score*."""
    return criterion["says"].get(score, criterion["generic"])


def _score_cell(score: int, muted: bool) -> str:
    """A score cell: the number and CancerVar's own word for it."""
    text = f"{score:+d}" if score else "0"
    word = _VALUE_WORD.get(score, "")
    cell = "padding:5px 10px;vertical-align:top;text-align:center;"
    if muted or score == 0:
        return f'<td style="{cell}color:#999999;">{text} · {word}</td>'
    color = _SCORE_COLORS.get(score, _UNKNOWN_TIER_COLOR)
    return (
        f'<td style="{cell}background:{color};color:#fff;font-weight:600;">'
        f"{text} · {word}</td>"
    )


# ============================================================================
# The markers behind a criterion (issue #198)
# ============================================================================
#
# CBP1/CBP2/CBP3 are the three criteria CancerVar backs with a list of database markers,
# and the disclosure below is where those markers appear. Four measurements from #198
# shape it, and none of them is cosmetic:
#
# 1. **A marker list is not this patient's therapy list.** `CancerVar.py:1092` appends
#    every hit *before* the tumour-type test at `:1104`, which only demotes the score.
#    Over the 94 files carrying `tumor_tissue`: 10.6% of therapy markers match the
#    sample's own tissue, and 80.1% of populated cells have none that do. Level A runs
#    781 in-tissue against 6,542 off. So the groups are named and the sample's tissue is
#    stated — a flat list would read as "the therapies for this variant" and be wrong on
#    four rows in five.
# 2. **The response is load-bearing.** ~11k of ~101k resolved rows say `Resistance`,
#    `Resistant`, `no benefit` or `Poor Outcome`. A drug named without its response
#    inverts the claim, so every line carries one and the adverse ones are marked.
# 3. **Uncapped, on purpose.** Groups reach 72 in-tissue and 143 off-tissue. The summary
#    states the count *before* the reader opens it, so opening is informed — and a
#    truncated list inside static HTML is a dead end with no "show all" to offer.
# 4. **`<details>`, not `st.expander`.** The table is one `st.markdown` string, so a
#    Streamlit expander cannot sit inside a cell without breaking the row alignment #188
#    chose. Streamlit sanitises with DOMPurify 3.2.6 under `USE_PROFILES: {html: true}`,
#    whose tag set includes `details` and `summary`.

#: PubMed links per marker before the rest become a count. Three keeps a line readable;
#: the tail is long (p90 23 PMIDs, max 129) and is not what the reader came for.
_PMID_LINKS = 3

#: What CancerVar's evidence levels are worth to the score. `CancerVar.py:1095-1102`
#: grades A and B alike and C and D alike, so these are two strengths spelled four ways —
#: said once in the disclosure footer rather than implied per line.
_LEVEL_NOTE = (
    "Evidence levels are CancerVar's: A and B both score 2, C and D both score 1."
)


def _marker_line(marker) -> str:
    """One marker as a list item: level, drug, tumour type, response, PMIDs.

    Everything here is text out of a vendored data file, so it is escaped. The level is
    first because it is what the score was graded on, and the response is never dropped —
    see note 2 above.
    """
    level = html.escape(marker.level) or "—"
    name = html.escape(marker.drug) if marker.drug else "<em>no drug named</em>"
    where = html.escape(marker.cancer) if marker.cancer else "tumour type not stated"

    response = ""
    if marker.response:
        text = html.escape(marker.response)
        if marker.is_adverse:
            # Not a warning icon: this is the record's own finding, not a fault. Coloured
            # because a reader skimming drug names must not read a resistance marker as
            # an option.
            response = f' · <span style="color:#b3261e;font-weight:600;">{text}</span>'
        else:
            response = f" · {text}"

    pmids = ""
    if marker.pmids:
        shown = [
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/">{p}</a>'
            for p in marker.pmids[:_PMID_LINKS]
        ]
        rest = len(marker.pmids) - len(shown)
        more = f" +{rest}" if rest > 0 else ""
        pmids = (
            f' <span style="font-size:0.9em;color:#666;">'
            f'PMID {", ".join(shown)}{more}</span>'
        )

    return (
        f'<li style="margin:2px 0;">'
        f'<span style="display:inline-block;min-width:1.4em;font-weight:700;">{level}</span>'
        f"{name} — {where}{response}{pmids}</li>"
    )


def _marker_group(heading: str, markers: tuple) -> str:
    """A named group of markers, or the heading alone where the group is empty.

    An empty in-tissue group is drawn rather than skipped: "none in this tumour type" is
    the answer on 80.1% of populated cells, and #187 settled that an absence which is a
    finding gets named instead of leaving a gap the reader fills in wrongly.
    """
    label = (
        f'<div style="margin:6px 0 2px 0;font-weight:600;font-size:0.95em;">'
        f"{heading}</div>"
    )
    if not markers:
        return label + '<div style="color:#888;margin-left:2px;">none</div>'
    items = "".join(_marker_line(m) for m in markers)
    return label + f'<ul style="margin:0 0 0 18px;padding:0;">{items}</ul>'


def _marker_summary(mset: MarkerSet, tissue: str) -> str:
    """The one line the disclosure shows while closed.

    Carries the count, the strongest evidence level, and how many of them are in this
    sample's tumour type — which is the fact that decides whether the rest is worth
    opening. The count comes from the MAF cell, so it is stated even when the table
    cannot be read.
    """
    parts = [f"{mset.cited} {mset.noun}{'s' if mset.cited != 1 else ''}"]

    if mset.table_missing:
        parts.append("cannot be named")
    else:
        if mset.best_level:
            parts.append(f"best level {html.escape(mset.best_level)}")
        if mset.tissue_known and tissue:
            count = len(mset.in_tissue)
            where = html.escape(tissue)
            parts.append(f"{count} in {where}" if count else f"none in {where}")
        if mset.unresolved:
            parts.append(f"{mset.unresolved} unresolved")

    return " · ".join(parts)


def _marker_disclosure(mset: MarkerSet, tissue: str, score: int) -> str:
    """A ``<details>`` block for one criterion's markers, closed by default.

    Closed is what keeps the section the height #188 measured — median 2 rows — while
    still admitting the tail exists.
    """
    body = []

    if mset.table_missing:
        body.append(
            f'<div style="color:#666;">CancerVar recorded {mset.cited} '
            f"{html.escape(mset.noun)}{'s' if mset.cited != 1 else ''} for this variant. "
            f"MAFigate's copy of CancerVar's marker table is unavailable, so it cannot "
            f"name them.</div>"
        )
    elif not mset.markers:
        # Every cited index failed to resolve. Distinct from the table being absent: the
        # table is there and disagrees with the indices, which is what a mismatched
        # vintage looks like.
        body.append(
            f'<div style="color:#666;">None of the {mset.cited} indices CancerVar '
            f"recorded resolved to a {html.escape(mset.noun)} in MAFigate's copy of its "
            f"marker table.</div>"
        )
    else:
        if mset.tissue_known and tissue:
            body.append(
                _marker_group(
                    f"In this sample's tumour type ({html.escape(tissue)})", mset.in_tissue
                )
            )
            body.append(_marker_group("Other tumour types", mset.other_tissue))
        else:
            # No `tumor_tissue` column — 15 of 109 CancerVar-bearing files. Grouping
            # would assert a comparison the app cannot make, so the list is flat and the
            # footer says why rather than implying every marker is off-tissue.
            body.append(_marker_group("Markers recorded", mset.markers))

        if mset.unresolved:
            body.append(
                f'<div style="color:#666;margin-top:4px;">{mset.unresolved} of '
                f"{mset.cited} recorded indices did not resolve against MAFigate's copy "
                f"of the marker table.</div>"
            )

        footer = [
            "CancerVar records a marker whatever the tumour type; only the score is "
            "adjusted for this sample's."
        ]
        if not (mset.tissue_known and tissue):
            footer.append(
                "This file carries no tumour type, so MAFigate cannot say which of these "
                "match the sample's."
            )
        footer.append(_LEVEL_NOTE)
        body.append(
            f'<div style="color:#666;font-size:0.9em;margin-top:6px;">'
            f"{' '.join(footer)}</div>"
        )

    if score == 0:
        # #185 §4: `level` is reset inside the per-transcript loop while the list
        # accumulates across all of them, so the last transcript parsed decides the
        # score. 45 of 9,758 Therap_list cells and 103 of 3,905 Prog_list cells look like
        # this. Named where it happens, because a populated list under a zero score reads
        # as a contradiction otherwise.
        body.append(
            '<div style="color:#666;font-size:0.9em;margin-top:6px;">CancerVar scored '
            "this criterion 0 even so: it grades the last transcript it parsed, while "
            "the marker list accumulates over all of them.</div>"
        )

    return (
        f'<details style="margin-top:6px;">'
        f'<summary style="cursor:pointer;color:#2f5597;">'
        f"{_marker_summary(mset, tissue)}</summary>"
        f'<div style="margin:4px 0 2px 4px;font-weight:400;">{"".join(body)}</div>'
        f"</details>"
    )


def _cbp_table_html(
    pairs: list, muted: bool = False, markers: Optional[dict] = None, tissue: str = ""
) -> str:
    """``Criterion | CancerVar score | What it says`` over (criterion, score) pairs.

    Mirrors :func:`components.acmg_evidence._acmg_table_html`'s three-column layout and its
    explicit dark text — the light row backgrounds would otherwise render white-on-white in
    Streamlit's dark theme.
    """
    th = (
        "background-color:#4472C4;color:#fff;padding:6px 10px;"
        "text-align:left;font-weight:600;"
    )
    parts = [
        '<table style="width:100%;border-collapse:collapse;font-size:0.88em;'
        'table-layout:fixed;">',
        '<colgroup><col style="width:14%;"><col style="width:22%;">'
        '<col style="width:64%;"></colgroup>',
        f'<thead><tr><th style="{th}">Criterion</th>'
        f'<th style="{th}">CancerVar score</th>'
        f'<th style="{th}">What it says</th></tr></thead><tbody>',
    ]
    for i, (criterion, score) in enumerate(pairs):
        bg = "#f7f9fc" if i % 2 == 0 else "#ffffff"
        text_color = "#999999" if (muted or score == 0) else "#1a1a1a"
        cell = f"padding:5px 10px;vertical-align:top;color:{text_color};"
        marker = ""
        if criterion.get("never_evaluated"):
            marker = (
                ' <span title="CancerVar never evaluates this criterion" '
                'style="color:#b8860b;">&#9888;</span>'
            )
        # The markers behind this criterion, in the same cell as the sentence they
        # qualify (#198). Third column rather than a fourth: `table-layout:fixed` splits
        # 14/22/64, and a fourth column would take width from the sentences #188 chose.
        disclosure = ""
        if markers:
            mset = markers.get(criterion["code"])
            if mset is not None and mset.cited:
                disclosure = _marker_disclosure(mset, tissue, score)

        parts.append(
            f'<tr style="background-color:{bg};">'
            f'<td style="{cell}font-weight:600;">{criterion["code"]}{marker}'
            f'<br><span style="font-weight:400;font-size:0.85em;">'
            f'{criterion["name"]}</span></td>'
            f"{_score_cell(score, muted)}"
            f'<td style="{cell}">{_describe(criterion, score)}{disclosure}</td></tr>'
        )
    parts.append("</tbody></table>")
    return "\n".join(parts)


def _render_tier_badge(parsed: dict) -> None:
    """The tier, read from the evidence string rather than banded from the sum.

    On the **15** real files that carry the vector with no ``CancerVar`` column, this is the only
    place the tier appears at all — and it is the tool's own printed answer, matched against the
    column on 110,064 rows with no disagreement. (16 when issue #189 counted; issue #210 re-read
    the corpus with the repo's own fixtures excluded from it.)

    With one exception, which is the whole of issue #210's defect: on **73 rows of one such
    file** the cell holds a tier and a sum with no vector, :func:`parse_cancervar` refuses it,
    and this function is never reached — so there the tier appears in *no* place at all.
    """
    tier = parsed.get("tier") or ""
    caption = (
        '<span style="font-size:0.78em;color:#666;text-transform:uppercase;'
        'letter-spacing:0.5px;">CancerVar AMP/ASCO/CAP tier</span>'
    )
    if not tier:
        st.markdown(f'<div style="margin:2px 0 6px 0;">{caption}</div>', unsafe_allow_html=True)
        st.caption(
            "CancerVar recorded no tier in this variant's evidence string. MAFigate does not "
            "add the twelve scores up: the tier is a sum CancerVar itself computes, on "
            "thresholds of its own."
        )
        return
    st.markdown(
        f'<div style="margin:2px 0 12px 0;">{caption}<br>'
        f'<span style="display:inline-block;margin-top:5px;padding:8px 22px;'
        f'border-radius:8px;background:{tier_color(tier)};color:#ffffff;'
        f'font-size:1.2em;font-weight:700;">{tier_label(tier)}</span></div>',
        unsafe_allow_html=True,
    )


def _render_arithmetic(parsed: dict) -> None:
    """One caption: what CancerVar summed, to what, and whose thresholds those are.

    No distance-to-next-tier (#188). "2 points from Tier II" implies the app could move the
    variant, and it invites reading CancerVar's thresholds as the guideline's.
    """
    total = parsed["sum"]
    tier = parsed.get("tier") or ""
    where = f", which it recorded as **{tier_label(tier)}**" if tier else ""
    st.caption(
        f"CancerVar scored each of the twelve AMP/ASCO/CAP evidence categories and summed them "
        f"to **{total:+d}**{where} on its own scale ({_TIER_THRESHOLDS}). The categories and "
        f"the four tier names are the guideline's; the scores and the thresholds are "
        f"CancerVar's."
    )

    # Never seen: the printed sum matched the vector's own sum on all 119,194 vectors measured,
    # which is structural — CancerVar's ``classfy()`` sums the same list ``assign()`` printed.
    # Named rather than hidden if it ever happens, because the badge above shows what CancerVar
    # printed and this caption shows what the vector adds to.
    printed = parsed.get("printed_sum")
    if printed is not None and printed != total:
        st.caption(
            f"⚠️ CancerVar printed a sum of {printed:+d} for a vector of twelve scores that "
            f"adds to {total:+d}. The tier above is the one CancerVar printed."
        )


def render_cbp_evidence(
    cancervar_result, markers: Optional[dict] = None, tissue: str = ""
) -> None:
    """Render CancerVar's AMP/ASCO/CAP evidence: the tier, then the criteria that fired.

    Shape B of issue #188, which the dev drove and picked: **one signed table of the criteria that
    fired**, strongest claim first, with the remainder behind two expanders — *scored zero* and
    *never evaluated* — that are different facts and so are not one list. Median 2 rows over 109,416
    real rows; no row in the corpus has more than 10 of the twelve non-zero.

    Args:
        cancervar_result: the row's ``CancerVar and Evidence`` value.
        markers: ``{criterion_code: MarkerSet}`` from
            :func:`components.cancervar_markers.resolve_markers`, or ``None`` to draw the
            table with no marker disclosures. Optional so that the section still renders
            from the evidence string alone — the markers are a detail of a criterion, never
            a precondition for showing it.
        tissue: the sample's own tumour type, used to group the markers and named on
            screen. Empty where the file carries no ``tumor_tissue`` column, which makes
            the disclosure say so rather than group on a comparison it cannot make.

    The disclosures are drawn in the *scored zero* expander too (#198), because 45
    ``Therap_list`` and 103 ``Prog_list`` cells in the corpus are populated behind a
    criterion scored 0 — the rows where the score understates what CancerVar found, and so
    exactly the rows where hiding the markers would hide the discrepancy.

    Renders directly into the current Streamlit container.
    """
    parsed = parse_cancervar(cancervar_result)
    if parsed is None:
        st.caption("No CancerVar AMP/ASCO/CAP evidence available for this variant.")
        return

    _render_tier_badge(parsed)

    pairs = _criteria_with_scores(parsed)
    fired = [(c, s) for c, s in pairs if s != 0]
    # Strongest claim first in either direction, and a negative after a positive of equal
    # magnitude — so a dissent reads as a qualification of what precedes it rather than as the
    # headline. Ties broken by code so the order is stable.
    fired.sort(key=lambda pair: (-abs(pair[1]), pair[1] < 0, pair[0]["code"]))

    if fired:
        st.markdown(
            _cbp_table_html(fired, markers=markers, tissue=tissue), unsafe_allow_html=True
        )
    else:
        # 7.4% of real rows score zero on all twelve. Silence reads as a missing annotation —
        # the exact misreading issue #187 fixed for the badge row — so it is named.
        st.info(
            "CancerVar scored all twelve evidence categories at zero for this variant — it "
            "found nothing supporting and nothing opposing. That is a result, not a missing "
            "annotation."
        )

    scored_zero = [
        (c, s) for c, s in pairs if s == 0 and not c.get("never_evaluated")
    ]
    never_evaluated = [(c, s) for c, s in pairs if c.get("never_evaluated")]

    if scored_zero:
        with st.expander(f"Scored zero — no support either way ({len(scored_zero)})"):
            st.markdown(
                _cbp_table_html(scored_zero, muted=True, markers=markers, tissue=tissue),
                unsafe_allow_html=True,
            )
    if never_evaluated:
        with st.expander(f"Never evaluated by CancerVar ({len(never_evaluated)})"):
            st.caption(
                "CancerVar automates ten of the twelve categories. These two need manual input "
                "and are always zero, on every variant in every file — a zero here is not a "
                "finding about this variant."
            )
            st.markdown(_cbp_table_html(never_evaluated, muted=True), unsafe_allow_html=True)

    _render_arithmetic(parsed)
