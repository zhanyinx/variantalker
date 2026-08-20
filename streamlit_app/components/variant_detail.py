"""One variant, read out in full: the panel behind a row of the grid.

What a clinician sees after clicking a variant — its identity and read counts, the
annotations other parties and this pipeline's own predictors recorded about it, the
guideline verdict its own annotator reached and the evidence behind it, its population
frequencies, links out to UCSC, gnomAD, ClinVar and dbSNP, and — last, collapsed, and
explicitly not part of the verdict — AlphaMissense's own band and the ClinGen SVI PP3/BP4
reference scores.

Two rows, two different claims (issues #187 and #200)
-----------------------------------------------------
:func:`_render_guideline_classifications` states **the guideline verdict** — InterVar's
ACMG/AMP classification on a germline file, CancerVar's AMP/ASCO/CAP tier on a somatic
one. :func:`_render_clinical_badges` draws what is *not* that.

Which row a badge belongs in is decided by a **three-leg test** (issue #200). It draws in
the guideline row only if all three hold: it is a verdict on a *published clinical
guideline's* scale; *this pipeline's own annotator computed it* over this file's
annotations; and the reasoning behind it is *in the MAF*, so the criteria can sit beneath
it. Fail any leg and it is somebody else's claim, and it draws in the clinical row. Two
corollaries: a column that only *qualifies* another column's claim travels with that claim
(``ClinVar_VCF_CLNREVSTAT``, ``RENOVO_pls``), and a URL that *navigates* asserts nothing
about the variant, so the external links are never duplicates of a badge.

``CancerVar`` and ``InterVar`` used to be drawn in *both*, so a germline file showed
``InterVar: Pathogenic`` twice, in two colours, under two headings. The guideline row
owns them, because they pass all three legs and it is the row the evidence attaches to: a
badge in the clinical row is a string with no route to the criteria behind it.

``ClinVar_VCF_CLNSIG`` and ``am_class`` were drawn twice for longer, deliberately, because
the obvious fix emptied the clinical row on somatic data. Map #199 settled both, and the
de-duplication ran **the opposite way** to the one that objection was about. ClinVar fails
all three legs — another laboratory's aggregate, not computed here, no evidence vector
anywhere in the MAF — so it moved **down** into the clinical row, which is the move that
*keeps* that row populated (0.5% → 54.6% of somatic variants). AlphaMissense fails all
three too and is neither a database opinion nor a guideline verdict, so it left both rows
for a section of its own: :mod:`components.alphamissense`.

What the clinical row holds and what it says when it holds nothing is
:mod:`components.clinical_badges`. Its three members — ClinVar, RENOVO and ESCAT — are not
all "database annotations": RENOVO is an in-silico predictor *this pipeline ran*, which is
exactly the leg a third-party database fails. The row has no heading, so nothing above it
makes that claim; each badge's own label says what its claim is (issue #217).

Neither guideline's criteria table is here: each is a published clinical standard with its
own vocabulary and its own parser. ACMG/AMP as InterVar reports it lives in
:mod:`components.acmg_evidence`; AMP/ASCO/CAP as CancerVar reports it lives in
:mod:`components.cbp_evidence`, which also owns CancerVar's four tier names and their
badge colours — so the tier this panel badges in the guideline row and the tier the
evidence section badges one row below it cannot come apart (issue #189).

The scale that is nobody's verdict (issue #191)
-----------------------------------------------
By the same rule, the ClinGen SVI PP3/BP4 thresholds now live in
:mod:`components.reference_scales` — a published clinical calibration with its own
vocabulary. It is the one predictor surface on this page that **no classifier read**: its
four columns are disjoint from the nine InterVar's PP3/BP4 and CancerVar's CBP10 consult,
so nothing it shows contributed to the verdict above it.

It used to draw directly beneath the guideline badges, inside the same function, which made
it look like that verdict's evidence. It is now the panel's last classification row and
collapsed by default, for three measured reasons: it draws nothing at all on ~70% of
variants, it reaches strengths (Moderate, Strong) that InterVar's PP3 cannot express, and
it points the opposite way from the classifier on 810 germline and 2,632 somatic rows.
"""

import html
import streamlit as st
import pandas as pd
import numpy as np
from typing import NamedTuple, Optional

try:
    import plotly.express as px

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

from config.columns import spelled_in
from config.contaminated_columns import NOTHING_CONTAMINATED
from config.contaminated_columns import SESSION_KEY as _CONTAMINATED_COLUMNS
from config.missing_values import says_nothing, shown
from config.param_labels import label_of
from filters.arm_detection import read_arm_evidence
from filters.attribution import IN_REPORT_BECAUSE, explain_variant
from filters.numeric_columns import UnreadableNumericColumns

from .acmg_evidence import parse_intervar, render_acmg_evidence
from .alphamissense import render_alphamissense
from .cancervar_markers import resolve_markers
from .cbp_evidence import parse_cancervar, render_cbp_evidence, tier_color
from .clinical_badges import clinical_absences, clinical_badges
from .predictor_context import render_cbp10_inputs, render_pp3_bp4_inputs
from .reference_scales import render_reference_scale


# ---------------------------------------------------------------------------
# Variant Detail Panel
# ---------------------------------------------------------------------------

# The clinical row's membership, its vocabularies and its empty state live in
# :mod:`components.clinical_badges` — see :data:`~components.clinical_badges.CLINICAL_ROW`. This
# module renders what that one computes, which is what lets a guard drive the vocabularies
# directly instead of through a Streamlit call.
#
# Four things that used to be in the list there are not, and each absence is a decision:
#
# * ``CancerVar`` and ``InterVar`` (issue #187) — guideline verdicts, drawn three rows down where
#   their evidence can sit beneath them.
# * ``CIViC_Evidence_Level`` (issue #202) — an evidence level attaches to an evidence *item*, not
#   to a variant, and 133 of the 140 cells that would draw one hold a Python list repr, the worst
#   375 characters long. ``CIViC_Variant_URL`` already links out of the same panel, which is the
#   honest surface for a per-evidence-item scale. The column is still read by the filters.
# * ``am_class`` (issue #203) — neither a database opinion nor a guideline verdict, so it left
#   both badge rows for :mod:`components.alphamissense`, which draws its score's own band.
#
# ``RENOVO``'s column is ``RENOVO_Class``. It read ``filter_renovo`` until issue #108 — the name
# of a filter *parameter*, and a column no MAF carries (0 of 297 byte-distinct MAFs on this
# machine) — so the row walked it against ``row.index``, never matched, and **the RENOVO badge
# had never appeared in this panel**. It is the third surface to hold that one mistake: #95 found
# it in the overview circles, #108 in ``Clinical_Summary``, and this one was unknown to both,
# found by review while resolving #108. What made it survive on all three is that the name is
# plausible and nothing raises — a badge that is never drawn looks exactly like a variant with no
# RENOVO call. ``tests/test_clinical_summary.py`` now checks every per-row render table for the
# same category error.


class _Classifier(NamedTuple):
    """One guideline classifier, and the names its absence has to be spoken about by.

    ``evidence_column`` is written out as a literal because a sentence about a column the
    file does **not** have cannot read that column's spelling off the header. The real
    spellings carry leading and trailing spaces — ``' CancerVar: CancerVar and Evidence '``
    on every somatic file measured — and three files spell it dot-mangled as well, so a
    column that *is* present is looked up from this canonical name through the app's own
    resolver; see :func:`_evidence_column`.

    ``parsed_key`` is the key each tool's parser puts its verdict under, so
    :func:`_verdict_from_evidence` can ask *"will the section below actually show one?"* without
    a conditional on the tool's name. Issue #210 added it: the sentence about a file with no
    verdict column promised a verdict the section below had already refused to draw, on 73 real
    rows, and nothing in the code connected the promise to the refusal.
    """

    arm: str
    tool: str
    column: str
    evidence_column: str
    scale: str
    verdict: str
    parsed_key: str


#: The two classifiers this panel speaks for, and what strict arm-gating gates on (issue #187).
#:
#: **The gate is the row's own columns, not the arm the user selected.** Measured over 177
#: byte-distinct real MAFs on this machine: ``InterVar`` appears on 0 of 95 somatic files and
#: ``CancerVar`` on 0 of 57 germline ones, so column presence *is* strict arm-gating — with no
#: exception in the corpus and no new machinery. The two alternatives both cost something real:
#:
#: * ``filter_params["sample_type"]`` would hide a genuinely present InterVar verdict whenever
#:   the user has the somatic arm selected — the exact wrong-arm session issue #135 exists for,
#:   whose settled rule is *detect and warn, never override*. Hiding the file's own verdict is a
#:   heavier override than the notice that session already gets.
#: * :func:`~filters.arm_detection.read_arm_evidence` alone would be worse still: it cannot place
#:   24 of the 177, and 13 of those carry a CancerVar evidence vector this panel would then
#:   refuse to draw. It cannot be widened to read the evidence column either — that column is
#:   not in ``REQUIRED_INPUTS``, and ``tests/test_arm_detection.py`` requires the markers stay a
#:   subset of it.
#:
#: So ``read_arm_evidence`` is asked for one thing only: the *arm word*, used to say which
#: classification is missing on a file that carries neither column — see
#: :func:`_classification_absences`. And the mismatch between file and setting is not restated
#: here; :mod:`components.arm_notice` reports it once, on the data page.
_CLASSIFIERS = (
    _Classifier(
        arm="germline",
        tool="InterVar",
        column="InterVar",
        evidence_column="InterVar: InterVar and Evidence",
        scale="ACMG/AMP",
        verdict="classification",
        parsed_key="classification",
    ),
    _Classifier(
        arm="somatic",
        tool="CancerVar",
        column="CancerVar",
        evidence_column="CancerVar: CancerVar and Evidence",
        scale="AMP/ASCO/CAP",
        verdict="tier",
        parsed_key="tier",
    ),
)

# The ClinGen SVI PP3/BP4 thresholds and the table that draws them moved to
# :mod:`components.reference_scales` in issue #191, with the section itself. They are a
# published clinical calibration with their own vocabulary — the same reason the two guideline
# criteria tables are not here either — and, unlike either of those, they are *not* this page's
# answer: nothing the classifier above computed read any of their four columns.
#
# `_STRENGTH_POINTS` went with them and then went entirely: issue #194 found it defined and used
# nowhere, and moving a dead constant into a new module would have given it a second home
# instead of an obituary.

# Population frequency columns with display labels
_POP_FREQ_COLUMNS = [
    ("Freq_ExAC_ALL", "ExAC"),
    ("Freq_esp6500siv2_all", "ESP6500"),
    ("Freq_1000g2015aug_all", "1000G"),
    ("gnomAD_exome_AF", "gnomAD Exome"),
    ("gnomAD_genome_AF", "gnomAD Genome"),
]


def _badge_html(label: str, value: str, color: str, tooltip: Optional[str] = None) -> str:
    """One badge span, with everything read from the row escaped.

    Args:
        label: what the badge is — a tool or scale name, written here rather than read from data.
        value: what it says for this variant. **Comes from the MAF**, so it is escaped.
        color: the badge's background.
        tooltip: the hover text, or ``None`` for a badge with no hover. Also from the MAF, and
            escaped with ``quote=True`` because it lands inside an HTML *attribute*.

    Returns:
        str: the span, for ``st.markdown(..., unsafe_allow_html=True)``.

    **The escaping is here rather than at the call sites, and that is the point.** This
    interpolates into raw HTML, and until issue #200 it escaped nothing at all: the ESCAT list
    repr reached the page verbatim that way. A hover is worse than a cosmetic glitch — an
    unescaped quote inside ``title="…"`` breaks out of the attribute — so no future caller has to
    remember, and ``tests/test_clinical_badges.py`` asserts a caller cannot get an unescaped
    value onto the page by handing one in.

    ``label`` is escaped too, on the same principle: it is written by this app today, and a
    function that escapes two of its three text arguments is a trap for whoever adds the fourth.
    """
    attributes = f'style="background-color:{color};color:white;' \
                 f'padding:3px 10px;border-radius:4px;margin:2px 4px;' \
                 f'display:inline-block;font-size:0.85em;"'
    if tooltip is not None:
        attributes = f'title="{html.escape(tooltip, quote=True)}" {attributes}'
    return (
        f"<span {attributes}>"
        f"<b>{html.escape(label)}</b>: {html.escape(value)}</span>"
    )


def _untrustworthy_columns() -> frozenset:
    """The predictor columns this file holds the chromosome in, from the load-time verdict.

    Returns:
        frozenset: the column names ``page_modules.data_loading`` measured and stashed, or an
        empty set when there is no verdict to read.

    Empty covers two different situations on purpose. A file with no contaminated column and a
    render reached without one of the app's load paths — a test driving the panel directly, or a
    session predating the key — take the same branch, because the panel's job when it does not
    know is to draw what the file says, exactly as it did before issue #194. The measurement is
    a whole-column question and the panel holds one row, so it cannot answer it itself.
    """
    return st.session_state.get(_CONTAMINATED_COLUMNS, NOTHING_CONTAMINATED)


def _render_clinical_badges(row: pd.Series) -> None:
    """ClinVar, RENOVO and ESCAT for this variant — or what each of them is silent about.

    No heading, by decision (issue #217): the words *"clinical annotations"* used to reach a user
    from exactly one rendered string — the caption below — so on 100% of germline rows and 54.6%
    of somatic ones this row drew with no name at all, and one heading true of another
    laboratory's conclusion, a model this pipeline ran and an ESMO actionability tier at once was
    live and was not taken.

    The old *"No clinical annotations available for this variant."* is **deleted**, not reworded.
    It was false on all 47,095 somatic rows where this row is empty — 45.4% of the arm, 0.0% of
    germline — because on **100%** of them the guideline row below drew a badge. A variant reading
    ``AMP/ASCO/CAP (CancerVar): Tier III`` was being told it had no clinical annotation. What
    replaces it makes no global claim at all: one sentence per member state, computed off the row
    by :func:`~components.clinical_badges.clinical_absences` rather than hand-written, because a
    hand-written list goes stale (issue #203 found exactly that one section down).
    """
    badges = clinical_badges(row)
    if badges:
        st.markdown(
            " ".join(
                _badge_html(badge.label, badge.value, badge.color, badge.tooltip)
                for badge in badges
            ),
            unsafe_allow_html=True,
        )
        return

    for note in clinical_absences(row):
        st.caption(note)


def _render_population_frequencies(row: pd.Series) -> None:
    """Horizontal bar chart comparing population frequencies."""
    freq_data = []
    for col, label in _POP_FREQ_COLUMNS:
        if col in row.index:
            val = row[col]
            try:
                freq = float(val)
                if not np.isnan(freq):
                    freq_data.append({"Database": label, "Frequency": freq})
            except (ValueError, TypeError):
                continue

    if not freq_data:
        st.caption("No population frequency data available.")
        return

    freq_df = pd.DataFrame(freq_data)

    if not PLOTLY_AVAILABLE:
        st.dataframe(freq_df)
        return

    fig = px.bar(
        freq_df,
        x="Frequency",
        y="Database",
        orientation="h",
        title="Population Frequencies",
        color="Frequency",
        color_continuous_scale=["#2ca02c", "#f0c420", "#d62728"],
    )
    # Add 1% threshold line
    fig.add_vline(
        x=0.01, line_dash="dash", line_color="red",
        annotation_text="1%", annotation_position="top right",
    )
    fig.update_layout(
        margin=dict(l=10, r=10, t=35, b=10),
        height=250,
        coloraxis_showscale=False,
        xaxis=dict(range=[0, max(freq_df["Frequency"].max() * 1.2, 0.05)]),
    )
    st.plotly_chart(fig, use_container_width=True)


def _get_genome_build(row: pd.Series) -> dict:
    """
    Normalise NCBI_Build from the MAF row into link-ready identifiers.

    Returns a dict with:
        ucsc_db     : 'hg38' or 'hg19'          (for UCSC ?db=)
        gnomad_ds   : 'gnomad_r4' or 'gnomad_r2_1'  (for gnomAD ?dataset=)
    """
    raw = str(row.get("NCBI_Build", "")).strip().lower()
    if raw in ("grch38", "hg38", "38"):
        return {"ucsc_db": "hg38", "gnomad_ds": "gnomad_r4"}
    if raw in ("grch37", "hg19", "37", "b37"):
        return {"ucsc_db": "hg19", "gnomad_ds": "gnomad_r2_1"}
    # Default to hg38 when unknown
    return {"ucsc_db": "hg38", "gnomad_ds": "gnomad_r4"}


def _render_external_links(row: pd.Series) -> None:
    """Generate clickable external database links."""
    links = []
    build = _get_genome_build(row)

    # CIViC links from data
    for col, label in [
        ("CIViC_Variant_URL", "CIViC Variant"),
        ("CIViC_Entity_URL", "CIViC Entity"),
    ]:
        if col in row.index:
            url = str(row[col]).strip()
            if not says_nothing(url):
                links.append(f"- [{label}]({url})")

    # UCSC Genome Browser. UCSC wants `chr17`, and the value already is one: the prefix is
    # settled once when the file is opened (`config/chromosome_spelling.py`), so this used
    # to carry the app's only copy of that rule and no longer needs one. Nor does it want a
    # second: a guard here would be a claim that a row can reach this function unsettled,
    # which is the kind of claim issue #149 was opened to stop the app making.
    chrom = row.get("Chromosome")
    start = row.get("Start_Position")
    end = row.get("End_Position")
    if not says_nothing(chrom) and not says_nothing(start):
        end_str = start if says_nothing(end) else end
        links.append(
            f"- [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks"
            f"?db={build['ucsc_db']}&position={chrom}:{start}-{end_str})"
        )

    # gnomAD — prefer ID-based lookup, fall back to coordinate-based
    _gnomad_linked = False
    for id_col in ("gnomAD_exome_ID", "gnomAD_genome_ID"):
        gnomad_id = row.get(id_col)
        if not says_nothing(gnomad_id):
            links.append(
                f"- [gnomAD](https://gnomad.broadinstitute.org/variant/"
                f"{str(gnomad_id).strip()}?dataset={build['gnomad_ds']})"
            )
            _gnomad_linked = True
            break
    if not _gnomad_linked:
        ref = row.get("Reference_Allele")
        alt = row.get("Tumor_Seq_Allele2")
        if not any(says_nothing(v) for v in (chrom, start, ref, alt)):
            chrom_clean = str(chrom).replace("chr", "")
            links.append(
                f"- [gnomAD](https://gnomad.broadinstitute.org/variant/"
                f"{chrom_clean}-{start}-{ref}-{alt}?dataset={build['gnomad_ds']})"
            )

    # ClinVar — prefer direct ID link, fall back to search
    _clinvar_linked = False
    clinvar_id = row.get("ClinVar_VCF_ID")
    if not says_nothing(clinvar_id):
        try:
            cid = int(float(clinvar_id))
            links.append(
                f"- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/variation/"
                f"{cid}/)"
            )
            _clinvar_linked = True
        except (ValueError, TypeError):
            pass
    if not _clinvar_linked:
        allele_id = row.get("ClinVar_VCF_ALLELEID")
        if not says_nothing(allele_id):
            try:
                aid = int(float(allele_id))
                links.append(
                    f"- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/"
                    f"?term={aid}[alleleid])"
                )
                _clinvar_linked = True
            except (ValueError, TypeError):
                pass
    if not _clinvar_linked:
        gene = row.get("Hugo_Symbol")
        # The same fallback as the detail panel's Protein Change field, and it was dead in the
        # same way: ``a or b`` never reached ``b``, because the NaN pandas leaves in an empty
        # cell is truthy. Written out rather than shortened, so the two sites cannot drift.
        protein = row.get("Protein_Change")
        if says_nothing(protein):
            protein = row.get("AAChange.refGene")
        if not says_nothing(gene):
            search_term = f"{gene}"
            if not says_nothing(protein):
                search_term += f"+{protein}"
            links.append(
                f"- [ClinVar Search](https://www.ncbi.nlm.nih.gov/clinvar/"
                f"?term={search_term})"
            )

    # dbSNP
    for id_col in ("dbSNP_ID", "avsnp151"):
        rsid = row.get(id_col)
        if not says_nothing(rsid):
            rs = str(rsid).strip()
            links.append(
                f"- [dbSNP ({rs})](https://www.ncbi.nlm.nih.gov/snp/{rs})"
            )
            break

    if links:
        st.markdown("**External Links:**")
        st.markdown("\n".join(links))
    else:
        st.caption("No external links available.")


#: The classifiers by tool name, so a read can get at the canonical spelling of the column it
#: wants without a second copy of it. :data:`_CLASSIFIERS` is the one place either name is written.
_CLASSIFIER_BY_TOOL = {classifier.tool: classifier for classifier in _CLASSIFIERS}


def _evidence_column(row: pd.Series, tool: str) -> Optional[str]:
    """The ``'<tool> and Evidence'`` column's name as *this file* spells it, or ``None``.

    Resolved through :func:`~config.columns.spelled_in`, the app's one column-name resolver, from
    the canonical spelling :data:`_CLASSIFIERS` already records.

    **This used to be a substring match** — *any* column with the tool's name and the word
    ``Evidence`` in it — which issue #189 chose because the real spelling is space-padded
    (``' CancerVar: CancerVar and Evidence '``) and #189's own note called it survival "by luck
    rather than by design". Issue #212 measured both halves of that.

    The luck holds: across 167 real MAFs the substring matches **more than one** column on **0** of
    the 96 files carrying CancerVar's evidence column and **0** of the 56 carrying InterVar's, so
    it never had a tie to resolve arbitrarily. And the resolver now agrees with it exactly — same
    column on all 152 files, 0 disagreements, neither finding one the other misses.

    So this is not a bug fix; it is the panel keeping **one** rule instead of two. The substring
    match could not have served ``GERP++_RS`` (``GERP`` alone also matches ``GERP++_NR``, a
    different score) and the resolver could not have served this column until #212 taught it to
    strip padding. Two mechanisms, each blind where the other sees, was the thing worth removing:
    ``tests/test_column_spelling.py`` can hold the panel to a rule about *every* exposed column
    name only if there is one way for such a name to reach a row.

    The resolver is also the stricter of the two. It refuses an ambiguous match rather than taking
    the first, which on a number a clinician reads is the answer this map keeps choosing.
    """
    return spelled_in(row.index, _CLASSIFIER_BY_TOOL[tool].evidence_column)


#: Each tool's own parser, keyed by tool name — the same two functions the evidence sections
#: below are given, so :func:`_verdict_from_evidence` cannot come apart from what they draw.
_EVIDENCE_PARSERS = {"CancerVar": parse_cancervar, "InterVar": parse_intervar}


def _verdict_from_evidence(row: pd.Series, classifier: _Classifier) -> str:
    """The verdict the evidence section below will show for ``classifier``, or ``""`` for none.

    Asked of the same cell, through the same gate and the same parser that section is given —
    :func:`_evidence_value` and :data:`_EVIDENCE_PARSERS` — because the alternative is what issue
    #210 measured: a caption reading *"the AMP/ASCO/CAP tier below is the one CancerVar printed
    inside its own evidence string"* above a section that had already declined to draw one, on 73
    real rows of a file with no ``CancerVar`` column to fall back on. The tier was in the file
    and on screen nowhere.

    Three states return ``""``, and the corpus reaches exactly one of them:

    * the cell says nothing, so ``_evidence_value`` returns ``None`` and the section is not
      drawn at all — **0 real rows** on a file that also lacks the verdict column;
    * the parser refuses the cell — **73 real rows in 1 file**, every one holding a tier and a
      sum with no ``EVS=[...]`` vector behind it (``'1#Tier_IV_benign'``);
    * the vector parses but printed no verdict — **0 real rows**, and there the section says so
      in its own words.

    Not folded into one sentence per state: the caption's job is to say *no verdict is shown*,
    which is true of all three, and guessing *why* is how the sentence it replaces went wrong.
    """
    value = _evidence_value(row, classifier.tool)
    if value is None:
        return ""
    parsed = _EVIDENCE_PARSERS[classifier.tool](value)
    if parsed is None:
        return ""
    return str(parsed.get(classifier.parsed_key) or "").strip()


def _classifier_shown(row: pd.Series) -> Optional[str]:
    """Which classifier this panel spoke for on this row, for the reference scale's label.

    Args:
        row: the variant.

    Returns:
        str: ``"InterVar"`` or ``"CancerVar"``, or None when the row carries neither.

    Reads **the row's own columns**, in ``_CLASSIFIERS`` order, exactly as issue #187 settled
    arm-gating: either the verdict column or the evidence column counts as present, because **15
    of the 109** files carrying a CancerVar evidence vector have no ``CancerVar`` column and 1 of
    the 59 carrying ``InterVar`` has the reverse (16 of 124 when #189 counted; issue #210 re-read
    the corpus with the repo's own fixtures excluded from it). Never ``filter_params``, which
    would name the wrong tool on
    a wrong-arm session — and here that would produce a sentence claiming a tool that drew
    nothing on this page did not read these predictors, which is true but unreadable.

    Returns the first match rather than complaining about two. Issue #187 measured **zero** real
    MAFs carrying both — the three that appeared to were fixture and scratch copies — so a row
    with both is not a case the corpus has, and the label needs one name, not an argument. Issue
    #210 re-walked 184 files and found exactly one carrying both, which is **this repository's own
    ``tests/fixtures/test_sample.maf``**, reached through a copy of the tree dropped in
    ``~/Downloads``. So #187's finding stands, and now names the shape of its own exception.
    """
    for classifier in _CLASSIFIERS:
        if not says_nothing(row.get(classifier.column)):
            return classifier.tool
        if _evidence_column(row, classifier.tool) is not None:
            return classifier.tool
    return None


def _classification_absences(row: pd.Series) -> list:
    """What this panel has to *say* rather than draw, per classifier (issues #187 and #210).

    A silent absence reads exactly like "this variant has no classification", and with
    ClinVar's badge sitting beside it that misreading is the likely one. Each state below is a
    different fact about the file, so each gets its own sentence.

    Re-counted for issue #210 over **184 byte-distinct real MAFs and 336,170 rows** — every row
    of every file, with *this function* called on each rather than its branches re-implemented,
    so the branch order is inside the measurement instead of beside it. Two exclusions are
    load-bearing: a VS Code editor-history snapshot and the repo's own
    ``tests/fixtures/test_sample.maf`` share a first-2 MB hash, and that fixture is the only file
    on this machine carrying **both** verdict columns — counting it would have overturned #187's
    *"zero real MAFs carry both arms' markers"* on the strength of a test file, and it supplies
    two of the three empty verdict cells as well.

    * **a germline marker, no verdict column** — ``RENOVO_Class`` and no ``InterVar`` at all:
      **6,597 rows in 2 files**. The sentence names **the marker, not the arm**, and issue #210
      changed it for a measured reason: on **6,524 of those 6,597 rows** a CancerVar tier badge
      draws one row below, so *"this file was annotated for germline analysis"* was contradicted
      by the page itself. :func:`~filters.arm_detection.read_arm_evidence` calls those files
      germline because they carry ``RENOVO_Class`` and no ``CancerVar`` *column* — they do carry
      CancerVar's evidence vector — so the arm was a conclusion where the marker is a fact, and
      ``ArmEvidence.present`` exists to be quoted back.

      **The somatic form of this sentence cannot fire at all.**
      ``ANNOTATOR_MARKERS["somatic"]`` is ``("CancerVar",)`` — the very column this branch has
      just found absent — so ``read_arm_evidence`` can never answer ``"somatic"`` here. Quoting
      the markers actually found is what makes that safe rather than latent: on the somatic arm
      there are none to quote, and the empty tuple is the whole of the gate.
    * **verdict, no evidence column** — **4,125 rows in 5 files**: 4,044 somatic rows in 4 of the
      98 files carrying ``CancerVar``, and 81 germline rows in 1 of the 59 carrying ``InterVar``.
      The tier or class is real; what is missing is the reasoning, which is why the evidence
      section below draws nothing. Checked rather than assumed — on those five files the only
      other columns with *Evidence* in the name are CIViC's, which is another tool's opinion and
      not this one's reasoning.
    * **the column, an empty cell** — **1 row of the 107,296 carrying ``CancerVar``, in 1 of 98
      files, and 0 of the 211,634 carrying ``InterVar``**. The count here previously read *53 of
      24,330 somatic and 55 of 25,957 germline rows, across 53 of 96 and 55 of 56 files*, which
      issue #210 could not reproduce at any scale: the state is orders of magnitude rarer, and
      "present in nearly every file" was the reverse of true. Kept regardless — it is a
      *variant*-level absence rather than a file-level one (the tool ran and reached no verdict
      on this row), and a state with one row is still a state the sentence must be right about.
    * **evidence, no verdict column** — **9,121 rows in 15 of the 109 files carrying a CancerVar
      evidence vector**, and 0 germline rows. **Two sentences**, because one was false on 73 of
      them. Where the section below does badge a verdict (**9,048 rows**) the caption points at
      it, which is #189's wording after #188 found its predecessor's *"holds no AMP/ASCO/CAP
      tier"* false. Where it does not (**73 rows**, all in one file, cells holding a tier and a
      sum with no vector) the caption used to promise a tier that appeared nowhere on the page —
      the same defect one iteration later — and now says none is shown. Which of the two fires
      is decided by :func:`_verdict_from_evidence`, i.e. by asking the section rather than by
      predicting it.

      MAFigate still does not *add the vector up*: the tier is a sum CancerVar itself computes,
      and map #184's standing decision is to read the tool's answer rather than re-derive it —
      reading a printed value is not re-deriving one. Issue #189 measured the two against each
      other before trusting the string: **110,064 comparable rows across 108 real somatic files,
      0 disagreements**, and no tier name outside CancerVar's four.

    A file the header cannot place *and* that carries neither column says nothing here: there
    is no classifier to name, and the escalated fill warning is already on that screen.
    """
    evidence = read_arm_evidence(row.index)
    notes = []

    for classifier in _CLASSIFIERS:
        column_present = classifier.column in row.index
        has_verdict = column_present and not says_nothing(row.get(classifier.column))
        evidence_column = _evidence_column(row, classifier.tool)

        if not column_present and evidence_column is None:
            # The markers that place the file on this arm, less the one just found absent — so
            # the sentence quotes a column the reader can look for. Empty on the somatic arm by
            # construction; see the docstring.
            markers = tuple(
                marker
                for marker in evidence.evidence_for(classifier.arm)
                if marker != classifier.column
            )
            if evidence.detected == classifier.arm and markers:
                held = ", ".join(f"`{marker}`" for marker in markers)
                notes.append(
                    f"This file carries {held} but no `{classifier.column}` column, so "
                    f"MAFigate has no {classifier.scale} {classifier.verdict} to show for "
                    "this variant."
                )
            continue

        if has_verdict:
            if evidence_column is None:
                notes.append(
                    f"This file carries no `{classifier.evidence_column}` column, so the "
                    f"evidence {classifier.tool} recorded behind this {classifier.verdict} "
                    "is not in it."
                )
            continue

        if column_present:
            notes.append(
                f"{classifier.tool} recorded no {classifier.verdict} for this variant."
            )
        elif _verdict_from_evidence(row, classifier):
            notes.append(
                f"This file carries no `{classifier.column}` column, so the {classifier.scale} "
                f"{classifier.verdict} below is the one {classifier.tool} printed inside its "
                "own evidence string."
            )
        else:
            notes.append(
                f"This file carries no `{classifier.column}` column, so the only "
                f"{classifier.scale} {classifier.verdict} it could show is the one "
                f"{classifier.tool} printed inside this variant's evidence string — and "
                "MAFigate could not read one there, so it shows none."
            )

    return notes


def _render_guideline_classifications(row: pd.Series) -> None:
    """Render the guideline verdict for this variant, and name any that is absent.

    Renamed in issue #191, which took the ClinGen PP3/BP4 table out of this function. It used to
    be ``_render_variant_effect_predictions`` and drew two claims under one name: the guideline
    verdict, and a table of predictors that fed no guideline at all. The table now has its own
    section and its own module (:mod:`components.reference_scales`); what is left is the verdict,
    so the name says that.

    **The heading needed no rename — the removals were the fix** (issue #217). It drew a ClinVar
    aggregate and an AlphaMissense band beside the two verdicts, which is what made *"Guideline
    Classifications"* false of its own membership. Both left under the three-leg test, so the two
    that remain are exactly what the heading claims. On the 24 files the header cannot place, the
    row drew a badge on 1,390 rows and **every one of them was ClinVar** — so the heading's entire
    content there was a third-party database aggregate, on files where this pipeline reached no
    guideline verdict at all.
    """

    # --- Guideline classifications ---
    st.markdown("**Guideline Classifications:**")

    badges = []

    # ACMG (germline) — from InterVar. Arm-gated by the column's own presence, which is what
    # ``_CLASSIFIERS`` records the measurement for: no real MAF carries the other arm's
    # classification column, so the off-arm classifier is not drawn without a rule that says so.
    intervar = row.get("InterVar")
    if not says_nothing(intervar):
        val = str(intervar).strip()
        color = "#d62728" if "athogenic" in val else "#2ca02c" if "enign" in val else "#f0c420"
        badges.append(_badge_html("ACMG (InterVar)", val, color))

    # AMP/ASCO/CAP (somatic) — from CancerVar
    cancervar = row.get("CancerVar")
    if not says_nothing(cancervar):
        val = str(cancervar).strip()
        color = tier_color(val)
        badges.append(_badge_html("AMP/ASCO/CAP (CancerVar)", val, color))

    # ClinVar and AlphaMissense used to draw here as well, and both moved out under the three-leg
    # test rather than being deleted — ClinVar down to :func:`_render_clinical_badges` with its
    # review status, AlphaMissense to :mod:`components.alphamissense`. Neither is drawn twice now,
    # and neither is drawn in two colours: this row colours by tier, the clinical row by source.

    if badges:
        st.markdown(" ".join(badges), unsafe_allow_html=True)

    # Beneath the badges rather than instead of them, and issue #204 sharpened what the common
    # case is. It used to be a file carrying ClinVar and AlphaMissense but no tier for *this*
    # variant; both have since left this row, so what remains is the file that carries a verdict
    # column with no evidence column behind it, or an evidence vector with no verdict column —
    # 4,125 and 9,121 real rows. Either way there is a badge to draw and an absence to name at the
    # same time, which the old caption could not do because it fired only when nothing drew at all.
    notes = _classification_absences(row)
    for note in notes:
        st.caption(note)

    if not badges and not notes:
        st.caption("No guideline classifications available.")

def _evidence_value(row: pd.Series, tool: str):
    """Return the ``'<tool> and Evidence'`` column value for a row, or None.

    One function for both classifiers, because the two arms' evidence sections gate on exactly
    the same thing: the column being present and the cell saying something. The column is
    located by :func:`_evidence_column`, so this function and the sentence
    :func:`_classification_absences` writes about a missing evidence column cannot come apart on
    how the column is recognised.

    This *is* the arm gate (issue #187). ``InterVar`` appears on 0 of 95 real somatic files and
    ``CancerVar`` on 0 of 57 germline ones, so column presence gates strictly with no new
    machinery — and no real MAF carries both arms' evidence, so nothing needs to enforce mutual
    exclusion between the two sections below.
    """
    column = _evidence_column(row, tool)
    if column is None:
        return None
    val = row.get(column)
    return None if says_nothing(val) else val


#: The columns that say *which variant this is*. Used to confirm that the row the filter is
#: about to be asked about is the row on screen — see :func:`_row_as_the_filter_saw_it`.
_IDENTITY_COLUMNS = (
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
)


def _row_as_the_filter_saw_it(row: pd.Series) -> "pd.DataFrame | None":
    """The one-row frame this variant came from, or ``None`` if it cannot be identified.

    Two reasons this is not simply ``row.to_frame().T``.

    **Dtypes.** ``create_data_table`` builds the frame the dialog is handed by
    ``data.fillna("")`` and then ``astype(str)`` on every object column, and the vendored
    comparisons raise on object dtypes — so a frame rebuilt from this ``Series`` would either
    raise or answer about stringified data. Slicing ``maf_data`` by index preserves both the
    dtypes and the values the filter actually read.

    **Identity.** ``row.name`` is trusted here as a label of ``maf_data``, and whether it is
    one depends on how the caller recovered the row. This paragraph used to describe the AgGrid
    path recovering it from a browser-built frame whose index is a **position** rather than the
    MAF's own, so that a position which happens to exist in the report's index looks up **the
    wrong row**. That account is stale as of issue #310: ``variant_table._full_row`` no longer
    reads the component's identifier as a label at all, and both of its resolving branches now
    return a row of the report itself — ``full_data.iloc[position]`` or the sole column match —
    whose ``.name`` *is* the real label. The unresolved branch returns the browser's partial
    row, and on the ``🔍 View details`` route that row is still named by a position, as a
    string; against a MAF's int64 labels it is simply absent, so this function returns ``None``
    and draws nothing, which is the documented degradation.

    The cross-check below stays regardless, because the reason for it never depended on that
    account being current: the row arrives from a round trip we do not own, and a label is a
    claim about identity that this surface is not entitled to take on trust. It earned its keep
    — while the dialog was opening the wrong variant, this section drew nothing rather than
    explaining a variant nobody picked. A mismatch draws nothing: an explanation of a variant
    the user did not select is worse than no explanation, and this is the one surface where
    being confidently wrong costs the most.

    Both blanks count as agreement, because the two frames spell a missing value differently
    — ``fillna("")`` on one side and ``NaN`` on the other — and that is a rendering
    difference rather than a different variant.
    """
    maf = st.session_state.get("maf_data")
    if maf is None or row.name is None:
        return None
    try:
        if row.name not in maf.index:
            return None
        candidate = maf.loc[[row.name]]
    except (KeyError, TypeError):
        return None
    if len(candidate) != 1:
        return None

    seen = candidate.iloc[0]
    for column in _IDENTITY_COLUMNS:
        if column not in candidate.columns or column not in row.index:
            continue
        mine, theirs = seen[column], row[column]
        if says_nothing(mine) and says_nothing(theirs):
            continue
        if str(mine).strip() != str(theirs).strip():
            return None
    return candidate


def _render_report_standing(row: pd.Series) -> None:
    """Why this variant is, or is not, in the report (issue #147).

    Drawn from :func:`filters.attribution.explain_variant`, which re-runs the shipped filter
    over this one row with each setting isolated — so every claim here is the filter's own
    verdict rather than a comparison written in this module. Seven single-row filter runs,
    measured at a median of 9.6ms and a maximum of 16.0ms and flat in file size, which is why
    it is affordable on the open of a dialog.

    Silent in three states, each for its own reason:

    * **no run to describe** — ``filter_run`` is absent, so there are no parameters that
      produced a report and nothing to be outside of;
    * **the row cannot be identified** — see :func:`_row_as_the_filter_saw_it`;
    * **the filter refuses this row** — a MAF whose depth or VAF columns hold unreadable
      values has no verdict at all, and the refusal banner has already said so next door.

    The failing settings are listed **all** of them, never just one. That is the ticket's
    measured finding rather than a stylistic choice: 70–95% of excluded variants fall outside
    more than one setting, so naming a single one would be true of that setting and
    misleading about the variant.
    """
    run = st.session_state.get("filter_run")
    if run is None or not getattr(run, "params", None):
        return
    one_row = _row_as_the_filter_saw_it(row)
    if one_row is None:
        return

    try:
        explanation = explain_variant(one_row, run.params)
    except UnreadableNumericColumns:
        return

    st.markdown("---")

    if explanation.in_report:
        because = IN_REPORT_BECAUSE.get(explanation.reason)
        if because is None:
            # A reason column carrying a cell this table does not know. Report the standing
            # without inventing an account of it.
            st.markdown("**✅ In this report.**")
            return
        st.markdown(f"**✅ In this report** — {because}.")
        return

    if explanation.unreachable:
        st.markdown(
            "**❌ Not in this report.** MAFigate cannot attribute this to any of your "
            "settings: this variant carries no usable read-depth or VAF value, so it is "
            "left out whatever those are set to."
        )
        return

    failing = explanation.failing
    if not failing:
        # Unreachable: a rejected row fails at least one clause, or the criteria would have
        # admitted it — asserted over 322,238 real rows, where the count of rejected rows
        # failing nothing was 0 on both arms at both presets. Handled anyway, because the
        # alternative is a dialog reading "It falls outside 0 of your settings".
        st.markdown(
            "**❌ Not in this report.** MAFigate cannot attribute this to any one of your "
            "settings."
        )
        return

    st.markdown(
        f"**❌ Not in this report.** It falls outside "
        f"{'one of your settings' if len(failing) == 1 else f'{len(failing)} of your settings'}:"
    )

    # One markdown block, so the sources render as a real nested list rather than as
    # non-breaking spaces behind `unsafe_allow_html`. Streamlit needs the whole list in a
    # single call for the nesting to survive — a `st.markdown` per line renders each as its
    # own top-level list.
    lines = []
    for outcome in failing:
        lines.append(f"- **{outcome.label}**")
        lines.extend(f"    - {name} — {value}" for name, value in outcome.values)
    st.markdown("\n".join(lines))

    if explanation.retention_declined:
        st.caption(
            f"{label_of('skip_pathogenic')} is on and did not apply: nothing in your file "
            "calls this variant pathogenic."
        )

    recoverable = explanation.recoverable
    if len(recoverable) == 1:
        st.caption(f"Changing **{recoverable[0].label}** alone would bring it back.")
    elif recoverable:
        st.caption(
            "Changing any one of "
            + ", ".join(f"**{outcome.label}**" for outcome in recoverable)
            + " would bring it back."
        )
    else:
        st.caption(
            "No single change brings it back — every setting above would have to change."
        )

    if explanation.filled_inputs:
        st.warning(
            "Your file carries no "
            + ", ".join(f"`{column}`" for column in explanation.filled_inputs)
            + ", so MAFigate filtered on a stand-in value. Part of this answer is about a "
            "column your file does not have rather than a setting you chose."
        )


def render_variant_detail_panel(row: pd.Series) -> None:
    """Render a detail panel for one selected variant."""

    # Row 1: Basic variant info + sequencing metrics
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown(f"**Gene:** `{shown(row.get('Hugo_Symbol'))}`")
        st.markdown(f"**Classification:** {shown(row.get('Variant_Classification'))}")
        st.markdown(f"**Type:** {shown(row.get('Variant_Type'))}")

        protein = row.get("Protein_Change")
        if says_nothing(protein):
            protein = row.get("AAChange.refGene")
        if not says_nothing(protein):
            st.markdown(f"**Protein Change:** `{str(protein).strip()}`")

    with col2:
        chrom = shown(row.get("Chromosome"))
        start = shown(row.get("Start_Position"))
        end = shown(row.get("End_Position"))
        # No `chr` literal: the column carries the prefix by the time any of this draws
        # (`config/chromosome_spelling.py`). Written into the string here, it printed
        # `chrchr17` on 174 of the 187 real MAFs measured for issue #149 — on the one line
        # a clinician reads to identify the variant in front of them.
        st.markdown(f"**Position:** {chrom}:{start}-{end}")

        ref = shown(row.get("Reference_Allele"))
        alt = shown(row.get("Tumor_Seq_Allele2"))
        st.markdown(f"**Alleles:** {ref} > {alt}")

    with col3:
        vaf = row.get("tumor_f")
        if not says_nothing(vaf):
            try:
                vaf_val = float(vaf)
                st.metric("VAF", f"{vaf_val:.3f}")
            except (ValueError, TypeError):
                st.metric("VAF", str(vaf))

        # Named for the column it comes from, because the number directly beneath it is the
        # *other* depth — the tumour's alt and ref reads, which is what the minimum-read-depth
        # setting gates on. Labelling this one plainly "Depth" beside that pair put the two
        # numbers on one screen under one word, and they differ on most rows (issue #127).
        dp = row.get("DP")
        if not says_nothing(dp):
            st.metric("Depth (DP)", str(dp))

        # Asked of the values rather than of their truthiness, because a read count of zero is
        # a measurement and ``if t_alt and t_ref`` discarded it: ``t_ref_count`` is 0 on 67,797
        # rows across 157 of the 303 MAFs measured for issue #131 — every variant called at 100%
        # VAF — and the caption was simply absent on all of them.
        t_alt = row.get("t_alt_count")
        t_ref = row.get("t_ref_count")
        if not says_nothing(t_alt) and not says_nothing(t_ref):
            st.caption(f"Alt/Ref reads: {t_alt}/{t_ref}")

    # Row 1b: why this variant is, or is not, in the report (issue #147). Above the
    # annotations because for a variant opened from the Failed Filters tab it *is* the
    # question, and below the identity block because an answer wants its subject named first.
    _render_report_standing(row)

    # Row 2: Clinical annotations
    st.markdown("---")
    _render_clinical_badges(row)

    # Row 3: the guideline verdict (ACMG / AMP-ASCO-CAP)
    st.markdown("---")
    _render_guideline_classifications(row)

    # Row 3b: InterVar ACMG/AMP evidence (germline variants only), then the scores behind the two
    # criteria that have any — issue #190. Attached to the evidence section rather than given a
    # rule of its own: it is the reasoning behind two of the rows above it, not a fifth claim.
    intervar_value = _evidence_value(row, "InterVar")
    if intervar_value is not None:
        st.markdown("---")
        render_acmg_evidence(intervar_value, show_unmet=False)
        render_pp3_bp4_inputs(row, intervar_value, _untrustworthy_columns())

    # Row 3c: CancerVar AMP/ASCO/CAP (CBP) evidence (somatic variants only) — issue #189, then
    # the scores behind CBP10 (issue #190).
    #
    # Beneath the guideline row for the same reason row 3b is: issue #187 settled that the
    # guideline row owns the classification, so the evidence behind one attaches to the row that
    # states it. Gated on the file's own evidence column, which is what makes this strictly
    # arm-gated — see :func:`_evidence_value`.
    #
    # Three depths, outermost first: the twelve criteria, then the *markers* behind three of
    # them (#198), then the *scores* behind CBP10 (#190). The markers ride inside the criteria
    # table because a marker belongs to one criterion; CBP10's scores follow the table because
    # they are one criterion's inputs and get their own section.
    #
    # The markers are resolved here rather than inside the section, because the section takes one
    # value while this needs four more columns off the row — and because `resolve_markers` reads a
    # file, which keeps that read at the row boundary rather than inside a renderer.
    #
    # `tumor_tissue` is passed, not looked up in `filter_params`: issue #187 settled that
    # arm-gating and row facts come from the row's own columns. 15 of 109 CancerVar-bearing
    # files carry no such column, and there the disclosure says so instead of grouping on
    # a tissue it does not have.
    cancervar_value = _evidence_value(row, "CancerVar")
    if cancervar_value is not None:
        st.markdown("---")
        tissue = row.get("tumor_tissue")
        tissue = "" if says_nothing(tissue) else str(tissue).strip()
        render_cbp_evidence(
            cancervar_value,
            markers=resolve_markers(row.get, tissue),
            tissue=tissue,
        )
        render_cbp10_inputs(row, cancervar_value, _untrustworthy_columns())

    # Row 3d: AlphaMissense's own band — issue #203.
    #
    # In the *other scales* region rather than beside a verdict, and that position is the whole of
    # the decision: mounted under the guideline badges the prototype interrupted the
    # verdict→evidence flow and read as part of the classification. The shape survived, the
    # position did not — the second time on this panel that placement rather than copy decided a
    # question (issue #191 was the first).
    #
    # A sibling of the ClinGen table, never a row inside it: that table is calibrated strictly to
    # ClinGen SVI, and an uncalibrated fifth row made its own tally report AlphaMissense as having
    # failed a threshold that does not exist for it. Drawn *above* it so the ClinGen table stays
    # the panel's last classification row.
    #
    # Given the same `_classifier_shown` answer as the table below, for the same reason: the
    # provenance caption names whose verdict this did not contribute to, and it is read from the
    # row's own columns per issue #187 — never from `filter_params`.
    render_alphamissense(row, _classifier_shown(row), _untrustworthy_columns())

    # Row 3e: the ClinGen SVI PP3/BP4 reference scale — issue #191.
    #
    # Last of the classification rows, collapsed, and deliberately *after* the two evidence
    # sections rather than under the verdict where it used to draw. Everything above this line
    # is the classifier's own answer and the classifier's own reasoning; this is a different
    # scale, computed by the app, that fed none of it. Position carries that as much as the
    # label does — a table directly beneath a verdict reads as its evidence whatever the
    # heading says, and the two point opposite ways on 13.0% of the somatic rows where both
    # take a side.
    #
    # Not arm-gated: it is neither arm's scale, so it draws on both. The tool is passed only so
    # the label can name whose verdict it did not contribute to, and is read from the row's own
    # columns per issue #187 — never from `filter_params`.
    st.markdown("---")
    render_reference_scale(row, _classifier_shown(row), _untrustworthy_columns())

    # Row 4: Population frequencies + External links
    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        _render_population_frequencies(row)
    with col2:
        _render_external_links(row)
