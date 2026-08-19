"""ACMG/AMP evidence, as InterVar reports it.

The 28 criteria of the ACMG/AMP variant-classification guideline — their codes, tiers
and descriptions, including the two ClinGen SVI has retired — plus the parser for
InterVar's own encoding of which were met, and the tables that render them.

A published standard rather than a screen: :func:`parse_intervar` takes InterVar's
``InterVar and Evidence`` string and returns all 28 codes as 0/1, and is pure. Only
:func:`render_acmg_evidence` touches Streamlit.
"""

import re
from typing import Optional

import streamlit as st

from config.missing_values import says_nothing


# ============================================================================
# InterVar ACMG/AMP evidence
# ============================================================================

# Number of *real* criteria in each InterVar list. InterVar's lists are
# 0-indexed (criterion N is at index N-1) and carry one unused trailing slot.
_INTERVAR_LIST_SPEC = {"PS": 4, "PM": 6, "PP": 5, "BS": 4, "BP": 7}

_INTERVAR_CLASSES = [
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
]

# (code, tier, side, description). side: 'P' pathogenic / 'B' benign.
_ACMG_CRITERIA = [
    ("PVS1", "VeryStrong", "P", "null variant in a gene where LoF is a known disease mechanism"),
    ("PS1", "Strong", "P", "same amino-acid change as an established pathogenic variant"),
    ("PS2", "Strong", "P", "de novo (maternity & paternity confirmed)"),
    ("PS3", "Strong", "P", "well-established functional studies show a damaging effect"),
    ("PS4", "Strong", "P", "prevalence in affecteds significantly increased vs controls"),
    ("PM1", "Moderate", "P", "in a mutational hot spot / critical functional domain"),
    ("PM2", "Moderate", "P", "absent / ultra-rare in population databases"),
    ("PM3", "Moderate", "P", "detected in trans with a pathogenic variant (recessive)"),
    ("PM4", "Moderate", "P", "protein-length change (in-frame indel / stop-loss), non-repeat"),
    ("PM5", "Moderate", "P", "novel missense at a residue with a different known pathogenic missense"),
    ("PM6", "Moderate", "P", "assumed de novo (maternity/paternity not confirmed)"),
    ("PP1", "Supporting", "P", "cosegregation with disease in multiple affected family members"),
    ("PP2", "Supporting", "P", "missense in a gene with low benign-missense rate, missense is a common mechanism"),
    ("PP3", "Supporting", "P", "multiple computational lines support a deleterious effect"),
    ("PP4", "Supporting", "P", "patient phenotype / family history highly specific for the gene"),
    ("PP5", "Supporting", "P", "reputable source reports pathogenic [RETIRED by ClinGen SVI]"),
    ("BA1", "StandAlone", "B", "allele frequency >5% in a population database"),
    ("BS1", "Strong", "B", "allele frequency greater than expected for the disorder"),
    ("BS2", "Strong", "B", "observed in a healthy adult with full penetrance expected early"),
    ("BS3", "Strong", "B", "well-established functional studies show no damaging effect"),
    ("BS4", "Strong", "B", "lack of segregation in affected family members"),
    ("BP1", "Supporting", "B", "missense in a gene where primarily truncating variants cause disease"),
    ("BP2", "Supporting", "B", "in trans with a pathogenic variant (dominant) or in cis with any pathogenic variant"),
    ("BP3", "Supporting", "B", "in-frame indel in a repetitive region of unknown function"),
    ("BP4", "Supporting", "B", "multiple computational lines suggest no impact"),
    ("BP5", "Supporting", "B", "found in a case with an alternate molecular basis for disease"),
    ("BP6", "Supporting", "B", "reputable source reports benign [RETIRED by ClinGen SVI]"),
    ("BP7", "Supporting", "B", "synonymous, no predicted splice impact, nucleotide not highly conserved"),
]
_RETIRED_CRITERIA = {"PP5", "BP6"}

# Strongest-first ordering within a side (StandAlone only on benign side,
# VeryStrong only on pathogenic side, so a single rank works for both).
_TIER_RANK = {"StandAlone": 0, "VeryStrong": 1, "Strong": 2, "Moderate": 3, "Supporting": 4}

# Tier-cell backgrounds, graded by strength. All dark enough for white text
# (WCAG AA): pathogenic warm/red gradient, benign cool/blue-green gradient.
_TIER_CELL_COLORS = {
    ("P", "VeryStrong"): "#7f1d1d",
    ("P", "Strong"): "#a01b2e",
    ("P", "Moderate"): "#b5451b",
    ("P", "Supporting"): "#8a6d0b",
    ("B", "StandAlone"): "#0b3d66",
    ("B", "Strong"): "#1a5a8a",
    ("B", "Supporting"): "#0f6e60",
}

# Classification badge backgrounds (white text → WCAG AA).
_INTERVAR_BADGE_COLORS = {
    "Pathogenic": "#7f1d1d",
    "Likely pathogenic": "#c0392b",
    "Uncertain significance": "#595959",
    "Likely benign": "#117a65",
    "Benign": "#1f4e79",
}


def _all_acmg_codes():
    """All 28 criterion codes in canonical order."""
    return [c[0] for c in _ACMG_CRITERIA]


def _canonical_intervar_class(raw: str) -> str:
    """Normalize a classification string to InterVar's canonical phrasing."""
    if not raw:
        return ""
    low = raw.strip().lower()
    for canon in _INTERVAR_CLASSES:
        if low == canon.lower():
            return canon
    if low in ("vus", "uncertain", "uncertain_significance"):
        return "Uncertain significance"
    return raw.strip()


def parse_intervar(result) -> Optional[dict]:
    """Parse an InterVar 'InterVar and Evidence' value into a structured dict.

    Returns ``{"classification": str, "criteria": {code: 0/1, ...}}`` with all
    28 criterion codes present, or ``None`` when there is no parseable content.
    An already-parsed dict of the same shape is returned unchanged.

    Tolerates whitespace variation, an optional leading ``InterVar:`` label, and
    the bracketed 0/1 lists (PS/PM/PP/BS/BP), which are 0-indexed with one unused
    trailing slot (criterion N at index N-1).
    """
    # Already parsed -> don't re-parse.
    if isinstance(result, dict):
        if "classification" in result and "criteria" in result:
            return result
        return None

    if says_nothing(result):
        return None
    text = str(result).strip()

    # Drop an optional leading "InterVar:" label, then take the classification
    # as everything before the first criterion token (PVS1).
    body = re.sub(r"^\s*InterVar\s*:\s*", "", text, flags=re.IGNORECASE)
    m = re.search(r"\bPVS1\b", body)
    raw_class = (body[: m.start()] if m else body).strip()
    classification = _canonical_intervar_class(raw_class)

    criteria = {code: 0 for code in _all_acmg_codes()}

    # Scalar flags.
    for scalar in ("PVS1", "BA1"):
        sm = re.search(rf"\b{scalar}\s*=\s*(-?\d+)", text)
        if sm:
            criteria[scalar] = 1 if int(sm.group(1)) != 0 else 0

    # 0/1 lists: criterion N at 0-based index N-1; ignore the trailing slot.
    for prefix, count in _INTERVAR_LIST_SPEC.items():
        lm = re.search(rf"\b{prefix}\s*=\s*\[([^\]]*)\]", text)
        if not lm:
            continue
        values = [v.strip() for v in lm.group(1).split(",") if v.strip() != ""]
        for i in range(count):
            if i < len(values):
                try:
                    criteria[f"{prefix}{i + 1}"] = 1 if int(values[i]) != 0 else 0
                except ValueError:
                    criteria[f"{prefix}{i + 1}"] = 0

    if not classification and not any(criteria.values()):
        return None

    return {"classification": classification, "criteria": criteria}


def _acmg_table_html(rows: list, muted: bool = False) -> str:
    """Build an HTML table of ACMG criteria rows: Code | Tier | Description."""
    th = (
        "background-color:#4472C4;color:#fff;padding:6px 10px;"
        "text-align:left;font-weight:600;"
    )
    parts = [
        '<table style="width:100%;border-collapse:collapse;font-size:0.88em;'
        'table-layout:fixed;">',
        '<colgroup><col style="width:16%;"><col style="width:24%;">'
        '<col style="width:60%;"></colgroup>',
        f'<thead><tr><th style="{th}">Code</th><th style="{th}">Tier</th>'
        f'<th style="{th}">Description</th></tr></thead><tbody>',
    ]
    for i, (code, tier, side, desc) in enumerate(rows):
        bg = "#f7f9fc" if i % 2 == 0 else "#ffffff"
        # Explicit dark text so cells stay readable regardless of the Streamlit
        # theme (light backgrounds would otherwise show white text in dark mode).
        text_color = "#999999" if muted else "#1a1a1a"
        cell = (
            "padding:5px 10px;text-align:left;vertical-align:top;"
            f"color:{text_color};"
        )
        marker = ""
        if code in _RETIRED_CRITERIA:
            marker = (
                ' <span title="Retired by ClinGen SVI" '
                'style="color:#b8860b;">&#9888;</span>'
            )
        if muted:
            tier_cell = f'<td style="{cell}">{tier}</td>'
        else:
            color = _TIER_CELL_COLORS.get((side, tier), "#666")
            tier_cell = (
                f'<td style="{cell}text-align:center;background:{color};'
                f'color:#fff;font-weight:600;">{tier}</td>'
            )
        parts.append(
            f'<tr style="background-color:{bg};">'
            f'<td style="{cell}font-weight:600;">{code}{marker}</td>'
            f"{tier_cell}"
            f'<td style="{cell}">{desc}</td></tr>'
        )
    parts.append("</tbody></table>")
    return "\n".join(parts)


def render_acmg_evidence(intervar_result, *, show_unmet: bool = False) -> None:
    """Render InterVar ACMG/AMP evidence: a classification badge + colored tables.

    Renders directly into the current Streamlit container.
    """
    parsed = parse_intervar(intervar_result)
    if parsed is None:
        st.caption("No InterVar ACMG/AMP evidence available for this variant.")
        return

    classification = parsed.get("classification") or ""
    criteria = parsed.get("criteria") or {}

    # 1) Classification badge.
    badge_bg = _INTERVAR_BADGE_COLORS.get(classification, "#595959")
    label = classification or "Not classified"
    st.markdown(
        '<div style="margin:2px 0 12px 0;">'
        '<span style="font-size:0.78em;color:#666;text-transform:uppercase;'
        'letter-spacing:0.5px;">InterVar ACMG/AMP classification</span><br>'
        f'<span style="display:inline-block;margin-top:5px;padding:8px 22px;'
        f'border-radius:8px;background:{badge_bg};color:#ffffff;'
        f'font-size:1.2em;font-weight:700;">{label}</span></div>',
        unsafe_allow_html=True,
    )

    # 2) Met criteria, split by side and sorted strongest-first.
    def _met(side):
        rows = [c for c in _ACMG_CRITERIA if c[2] == side and criteria.get(c[0], 0) == 1]
        rows.sort(key=lambda c: _TIER_RANK.get(c[1], 99))
        return rows

    path_met = _met("P")
    benign_met = _met("B")

    col_p, col_b = st.columns(2)
    with col_p:
        st.markdown("**🔴 Pathogenic evidence**")
        if path_met:
            st.markdown(_acmg_table_html(path_met), unsafe_allow_html=True)
        else:
            st.caption("No pathogenic criteria met.")
    with col_b:
        st.markdown("**🔵 Benign evidence**")
        if benign_met:
            st.markdown(_acmg_table_html(benign_met), unsafe_allow_html=True)
        else:
            st.caption("No benign criteria met.")

    # 3) Optionally show the unmet criteria, muted and collapsed.
    if show_unmet:
        unmet = [c for c in _ACMG_CRITERIA if criteria.get(c[0], 0) == 0]
        with st.expander(f"Unmet criteria ({len(unmet)})"):
            st.markdown(_acmg_table_html(unmet, muted=True), unsafe_allow_html=True)
