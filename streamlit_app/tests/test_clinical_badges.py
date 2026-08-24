"""The clinical row, after map #199 re-homed every annotation the panel draws (issue #204).

Five decisions land in one place, so five guards do. Each was **made to fail before being
trusted**, per this repo's standing rule; issue #204's notes record what each mutation was and
which assertion caught it.

1. **Drawn exactly once.** Read off a recording ``st`` planted in *every* module the panel draws
   through, so the claim is about what a clinician sees rather than about one function's source.
   A text grep over ``variant_detail.py`` is satisfied by that module's own comments — it explains
   at length which columns are *no longer* drawn where — and this file is the reason the ticket
   asked for something else.
2. **ClinVar's star vocabulary is total, and its two edges do not render alike.** The mapping has
   to reach a value no MAF on this machine carries (``practice_guideline``), an unmapped value must
   fall back to the full string and **never** to zero stars, and an earned zero must be
   distinguishable from *no review status at all*.
3. **A qualifier that is absent must degrade, not render empty.** ``RENOVO_pls`` is missing on
   6,597 real rows across 2 files, and those are the only germline rows where the guideline row
   draws no verdict — so RENOVO's badge is the whole of the pathogenicity claim on screen there and
   must not read ``(nan)``, ``()`` or a bare ``(``.
4. **ESCAT's three tiering-table vintages, and a cell holding several tiers.** The tumour type
   comes from whichever column the file's vintage filled, no tumour type may be invented for the
   2 files carrying neither, a tier the scale does not define must not borrow a gloss, and
   grouping must keep every tumour type rather than de-duplicate tiers.
5. **The empty state fires exactly two captions.** One ClinVar state of three and one ESCAT state
   of two, on every empty row — never one, never three — and RENOVO is never named.

The two rarest states are reachable from the corpus but thin (the 16-row ``CLNREVSTAT`` orphan;
ClinVar-absent, which co-occurs *only* with ESCAT-absent), so they are named here rather than
sampled. ``practice_guideline`` is not reachable at all, which is exactly why it is written down.
"""

from __future__ import annotations

import ast
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from components.clinical_badges import (  # noqa: E402
    ABSENCE_NEVER_NAMES,
    ABSENCE_SPEAKS_FOR,
    CLINICAL_ROW,
    CLINVAR_REVIEW_STARS,
    RENOVO_MENTION,
    clinical_absences,
    clinical_badges,
    clinvar_badge,
    clinvar_stars,
    escat_badges,
    escat_facts,
    renovo_badge,
)

# ---------------------------------------------------------------------------
# A panel that records what it drew, and which module drew it
# ---------------------------------------------------------------------------

#: Every module the detail panel draws through. A section that draws through a module missing from
#: this tuple is invisible to :func:`panel_records`, so *"drawn exactly once"* would be measured
#: over half the page — the failure shape this repo keeps meeting. The list is asserted against the
#: panel's own imports by :func:`test_every_module_the_panel_draws_through_is_recorded`.
PANEL_DRAWING_MODULES = (
    "components.variant_detail",
    "components.acmg_evidence",
    "components.cbp_evidence",
    "components.cancervar_markers",
    "components.predictor_context",
    "components.reference_scales",
    "components.alphamissense",
)


class _FakeSessionState(dict):
    """Enough of ``st.session_state`` for the panel: ``.get`` and attribute reads."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:  # pragma: no cover - a panel read of a key we did not seed
            raise AttributeError(name) from exc


class _Recorder:
    """A stand-in ``st`` that appends every drawn string, tagged with the module that drew it."""

    def __init__(self, module: str, sink: list):
        self._module = module
        self._sink = sink
        self.session_state = _FakeSessionState()

    # --- the calls that put text on screen ---------------------------------
    def markdown(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def caption(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def write(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def warning(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def info(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def error(self, body="", *args, **kwargs):
        self._sink.append((self._module, str(body)))

    def metric(self, label="", value="", *args, **kwargs):
        self._sink.append((self._module, f"{label}: {value}"))

    # --- containers ---------------------------------------------------------
    def columns(self, spec, *args, **kwargs):
        count = spec if isinstance(spec, int) else len(spec)
        return [self for _ in range(count)]

    def expander(self, label="", *args, **kwargs):
        # The label is on screen even while the expander is collapsed, which is the whole reason
        # `reference_scales` puts its count there — so it is recorded like any other text.
        self._sink.append((self._module, str(label)))
        return self

    def container(self, *args, **kwargs):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    # --- drawn, but not as text ---------------------------------------------
    def plotly_chart(self, *args, **kwargs):
        return None

    def dataframe(self, *args, **kwargs):
        return None

    def __getattr__(self, name):  # pragma: no cover - anything else the panel reaches for
        def _noop(*args, **kwargs):
            return None

        return _noop


def panel_records(monkeypatch, row: pd.Series) -> list:
    """Every ``(module, text)`` the whole detail panel draws for one variant, in document order.

    Document order, and every module — which is what separates this from
    ``tests/test_components.py``'s ``_panel_text``, whose fake ``st`` is planted in
    ``variant_detail`` alone. AlphaMissense's section and the ClinGen table draw through their own
    modules, so a duplication between the badge row and either of them would be invisible to a
    single-module recorder: it would see one of the two draws and report the annotation drawn once.
    """
    from components import variant_detail

    sink: list = []
    for module_name in PANEL_DRAWING_MODULES:
        module = __import__(module_name, fromlist=["st"])
        if hasattr(module, "st"):
            monkeypatch.setattr(module, "st", _Recorder(module_name, sink))

    monkeypatch.setattr(variant_detail.st, "session_state", _FakeSessionState())
    variant_detail.render_variant_detail_panel(row)
    return sink


#: Values chosen so that counting them counts renders rather than coincidences — the convention
#: ``tests/test_components.py`` set for the same question one map earlier.
MARKERS = {
    "ClinVar_VCF_CLNSIG": "Pathogenic_ZZCLINVAR",
    "RENOVO_Class": "HP Benign_ZZRENOVO",
    "ESCAT_CANCER": "ZZTUMOUR",
    "am_pathogenicity": "0.777",
    "am_class": "likely_pathogenic_ZZALPHA",
    "CIViC_Evidence_Level": "B_ZZCIVIC",
}


def a_variant(**overrides) -> pd.Series:
    """A row carrying every annotation map #199 re-homed, plus enough identity to render."""
    row = {
        "Hugo_Symbol": "TP53",
        "Variant_Classification": "Missense_Mutation",
        "Variant_Type": "SNP",
        "Chromosome": "chr17",
        "Start_Position": "7577120",
        "End_Position": "7577120",
        "Reference_Allele": "C",
        "Tumor_Seq_Allele2": "T",
        "NCBI_Build": "GRCh38",
        "ClinVar_VCF_CLNSIG": MARKERS["ClinVar_VCF_CLNSIG"],
        "ClinVar_VCF_CLNREVSTAT": "criteria_provided,_single_submitter",
        "RENOVO_Class": MARKERS["RENOVO_Class"],
        "RENOVO_pls": "0.00432",
        "ESCAT": "IA",
        "ESCAT_TISSUE": "BONE_MARROW",
        "ESCAT_CANCER": MARKERS["ESCAT_CANCER"],
        "CIViC_Evidence_Level": MARKERS["CIViC_Evidence_Level"],
        "am_class": MARKERS["am_class"],
        "am_pathogenicity": MARKERS["am_pathogenicity"],
    }
    row.update(overrides)
    return pd.Series(row)


# ---------------------------------------------------------------------------
# 1. Drawn exactly once, in the home its ticket named
# ---------------------------------------------------------------------------

#: ``column -> (times the value may appear on the page, the modules allowed to draw it)``.
#:
#: Written as an expectation per annotation rather than as one blanket rule, because the five
#: answers are genuinely different: three moved home, one left the badge rows for a section of its
#: own, and one left the panel's *claims* altogether while keeping its link. A single "appears
#: once" rule would have to make an exception for CIViC on the day it was written.
DRAWN_ONCE = {
    "ClinVar_VCF_CLNSIG": (1, {"components.variant_detail"}),
    "RENOVO_Class": (1, {"components.variant_detail"}),
    "ESCAT_CANCER": (1, {"components.variant_detail"}),
    "am_pathogenicity": (1, {"components.alphamissense"}),
    # Issue #203: the class label is not drawn at all. The band is named from the score instead,
    # because 100 of 139 real files spell `am_class` in a vocabulary AlphaMissense does not
    # publish, so echoing it would reproduce an overstatement.
    "am_class": (0, set()),
    # Issue #202: an evidence level attaches to an evidence item, not to a variant. The panel
    # still links `CIViC_Variant_URL`, which is a route rather than a claim.
    "CIViC_Evidence_Level": (0, set()),
}


@pytest.mark.parametrize("column", sorted(DRAWN_ONCE))
def test_each_re_homed_annotation_is_drawn_exactly_once(monkeypatch, column):
    """The ticket's central rule, measured on the page rather than read off the source.

    Every one of these was drawn twice, in two colours, at some point in this map's history:
    ``ClinVar_VCF_CLNSIG`` and ``am_class`` in both badge rows for as long as the panel has had
    two, and ``RENOVO_Class`` never at all — it was walked against a filter *parameter* name until
    issue #108, so a badge that never appeared looked exactly like a variant with no RENOVO call.

    The assertion is over the *rendered text of the whole panel*, so a duplicate cannot hide in a
    section drawn by another module.
    """
    records = panel_records(monkeypatch, a_variant())
    marker = MARKERS[column]
    expected_count, expected_modules = DRAWN_ONCE[column]

    drawn = [(module, text) for module, text in records if marker in text]
    assert len(drawn) == expected_count, (
        f"{column} is drawn {len(drawn)} times, expected {expected_count}: "
        f"{[module for module, _text in drawn]}"
    )
    assert {module for module, _text in drawn} == expected_modules


def test_the_recorder_can_see_a_duplicate_when_there_is_one(monkeypatch):
    """The guard above is worth its green only if it can go red.

    A recorder planted in the wrong modules, or a marker no section ever echoes, would pass
    :func:`test_each_re_homed_annotation_is_drawn_exactly_once` forever. So the duplicate is
    reintroduced deliberately here — a second badge drawn from the guideline row, exactly the shape
    the ticket removed — and the count has to notice.
    """
    from components import variant_detail

    real = variant_detail._render_guideline_classifications

    def also_draws_clinvar(row):
        real(row)
        variant_detail.st.markdown(
            variant_detail._badge_html("ClinVar", str(row.get("ClinVar_VCF_CLNSIG")), "#006400"),
            unsafe_allow_html=True,
        )

    monkeypatch.setattr(
        variant_detail, "_render_guideline_classifications", also_draws_clinvar
    )
    records = panel_records(monkeypatch, a_variant())
    drawn = [text for _module, text in records if MARKERS["ClinVar_VCF_CLNSIG"] in text]
    assert len(drawn) == 2, "the recorder cannot see a second draw, so its green means nothing"


def test_every_module_the_panel_draws_through_is_recorded():
    """:data:`PANEL_DRAWING_MODULES` against the panel's own imports, so it cannot fall behind.

    A section moved into a new module — which is exactly what issue #204 did twice — would
    otherwise draw outside the recorder's view, and every count above would keep passing while
    measuring less of the page each time.
    """
    source = (STREAMLIT_APP / "components" / "variant_detail.py").read_text()
    imported = {
        f"components.{node.module.lstrip('.')}"
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module
    }

    # A module with no `st` of its own cannot draw, so it needs no recorder. That is not a
    # convenience: `components.clinical_badges` is deliberately Streamlit-free — it turns a row
    # into values and `variant_detail` renders them — and demanding a recorder for it would be
    # demanding it grow one.
    draws = set()
    for name in imported:
        module = __import__(name, fromlist=["st"])
        if hasattr(module, "st"):
            draws.add(name)

    missing = sorted(draws - set(PANEL_DRAWING_MODULES))
    assert not missing, (
        f"the panel draws through {missing}, which the recorder does not patch — every "
        "drawn-exactly-once count above is measured over less than the whole page"
    )
    assert "components.alphamissense" in draws, (
        "the check has stopped seeing the section issue #204 added, so it would not see the next"
    )


def test_the_guideline_row_no_longer_reads_either_re_homed_column():
    """The source half of the same claim, on the one function that used to draw both.

    Asked of the parsed function rather than of the module's text, because
    ``variant_detail.py``'s own comments name both columns at length in order to say where they
    went — a text guard over this well-commented file is satisfied by the sentence explaining the
    rule it is meant to enforce.
    """
    from components import variant_detail

    tree = ast.parse((STREAMLIT_APP / "components" / "variant_detail.py").read_text())
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef)
        and node.name == "_render_guideline_classifications"
    )
    literals = {
        node.value
        for node in ast.walk(function)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    for gone in ("ClinVar_VCF_CLNSIG", "ClinVar_VCF_CLNREVSTAT", "am_class", "am_pathogenicity"):
        assert gone not in literals, f"{gone} is still read in the guideline row"

    assert variant_detail._render_guideline_classifications is not None


# ---------------------------------------------------------------------------
# 2. ClinVar's stars: total, and with two edges that must not look alike
# ---------------------------------------------------------------------------

#: Every review status the corpus holds, with its star level and the rows carrying it, re-walked
#: for issue #204. **Ten** values, not the nine issue #200 recorded:
#: ``no_interpretation_for_the_single_variant`` is on 50 real rows and was
#: falling through to the unrecognised branch, printing its raw string where the stars belong.
CORPUS_REVIEW_STATUSES = {
    "reviewed_by_expert_panel": (3, 4831),
    "criteria_provided,_multiple_submitters,_no_conflicts": (2, 101612),
    "criteria_provided,_single_submitter": (1, 13298),
    "criteria_provided,_conflicting_classifications": (1, 1301),
    "criteria_provided,_conflicting_interpretations": (1, 243),
    "no_assertion_criteria_provided": (0, 3901),
    "no_classification_for_the_single_variant": (0, 804),
    "no_interpretation_for_the_single_variant": (0, 50),
    "no_classification_provided": (0, 62),
    "no_assertion_provided": (0, 21),
}


@pytest.mark.parametrize("status", sorted(CORPUS_REVIEW_STATUSES))
def test_every_review_status_the_corpus_holds_has_a_star_level(status):
    """The mapping is total over what is on disk, and at the level that was measured.

    Named per value rather than asserted as a count, because a count would have been satisfied by
    the nine issue #200 found — and the tenth, on 50 real rows, was rendering its raw string on the
    badge where two stars' worth of scale belongs.
    """
    stars, _rows = CORPUS_REVIEW_STATUSES[status]
    assert CLINVAR_REVIEW_STARS[status] == stars
    assert clinvar_stars(status).count("★") == stars


def test_the_older_spelling_of_each_level_folds_onto_the_newer_one():
    """Two pairs, one rule: ClinVar renamed a level and the corpus holds both spellings.

    Issue #200 measured the *conflicting* pair. The *single variant* pair was found re-deriving
    the vocabulary for this ticket, and it is the same shape one level down — so it is folded by
    the same rule rather than treated as a new decision.
    """
    pairs = [
        (
            "criteria_provided,_conflicting_classifications",
            "criteria_provided,_conflicting_interpretations",
        ),
        (
            "no_classification_for_the_single_variant",
            "no_interpretation_for_the_single_variant",
        ),
    ]
    for newer, older in pairs:
        assert clinvar_stars(newer) == clinvar_stars(older), (newer, older)
        assert clinvar_stars(older) is not None, f"{older} still falls through to the raw string"


def test_the_unobserved_four_star_level_is_still_handled():
    """``practice_guideline`` is on 0 of 226,020 rows, so a corpus-driven guard cannot reach it.

    These MAFs were annotated 2022–2024. Unobserved is not impossible, and a mapping that only
    covers what one corpus happened to contain is a mapping that falls back to *unrecognised* on
    ClinVar's strongest level — which would print the raw string where four stars belong.
    """
    assert clinvar_stars("practice_guideline") == "★★★★"


def test_the_two_spellings_of_conflicting_are_one_level():
    """Annotation-vintage drift, not two levels.

    ``conflicting_classifications`` (1,301 rows) and ``conflicting_interpretations`` (243) are
    ClinVar's new and old names for the same one-star level. Anything ordering or colouring by
    these strings has to fold them together, or the same level sorts two ways.
    """
    new = clinvar_stars("criteria_provided,_conflicting_classifications")
    old = clinvar_stars("criteria_provided,_conflicting_interpretations")
    assert new == old == "★☆☆☆"


def test_an_unrecognised_status_renders_in_full_and_never_as_zero_stars():
    """Zero stars is a claim; not recognising a value is not that claim.

    Zero stars asserts *ClinVar assessed this with no assertion criteria*. An unmapped value means
    MAFigate did not recognise the level, and guessing the lower of the two promises a reading the
    page cannot support — the category error issue #210 fixed one row up.
    """
    assert clinvar_stars("some_status_ClinVar_added_later") is None

    badge = clinvar_badge(
        pd.Series(
            {
                "ClinVar_VCF_CLNSIG": "Pathogenic",
                "ClinVar_VCF_CLNREVSTAT": "some_status_ClinVar_added_later",
            }
        )
    )
    assert "some_status_ClinVar_added_later" in badge.value
    assert "★" not in badge.value and "☆" not in badge.value


def test_an_earned_zero_and_no_review_status_do_not_render_alike():
    """4,788 rows hold a zero-star status; 1,637 carry a significance with no status at all.

    One says *ClinVar assessed this with no assertion criteria*, the other says *MAFigate does not
    know how this was reviewed*. A badge showing nothing where the stars go says both at once, so
    the earned zero draws the empty scale and the absence draws no scale.
    """
    earned = clinvar_badge(
        pd.Series(
            {
                "ClinVar_VCF_CLNSIG": "Pathogenic",
                "ClinVar_VCF_CLNREVSTAT": "no_assertion_criteria_provided",
            }
        )
    )
    absent = clinvar_badge(pd.Series({"ClinVar_VCF_CLNSIG": "Pathogenic"}))
    missing_column = clinvar_badge(pd.Series({"ClinVar_VCF_CLNSIG": "Pathogenic"}))

    assert earned.value == "Pathogenic ☆☆☆☆"
    assert absent.value == "Pathogenic"
    assert earned.value != absent.value
    assert absent.tooltip == missing_column.tooltip
    assert "does not know" not in (earned.tooltip or "")


def test_the_full_review_status_travels_on_hover():
    """Abbreviated on the badge, in full on hover — issue #200's decision, both halves.

    It matters most on the 887 rows whose status says *no classification was provided* while the
    badge beside it shows a significance anyway: consistent for a zero-star ClinVar aggregate, and
    the pairing most likely to read as a bug.
    """
    badge = clinvar_badge(
        pd.Series(
            {
                "ClinVar_VCF_CLNSIG": "Pathogenic",
                "ClinVar_VCF_CLNREVSTAT": "no_classification_provided",
            }
        )
    )
    assert badge.value == "Pathogenic ☆☆☆☆"
    assert "no classification provided" in badge.tooltip


def test_a_hover_cannot_break_out_of_its_own_attribute(monkeypatch):
    """The hover puts MAF data inside an HTML attribute, which ``_badge_html`` was not safe for.

    It interpolated ``label`` and ``value`` into ``st.markdown(unsafe_allow_html=True)`` with no
    escaping at all — which is how issue #202's ESCAT list repr reached the page verbatim. Inside
    ``title="…"`` an unescaped quote is an attribute break rather than a cosmetic glitch.

    Asserted as *"the string the caller passed does not appear verbatim"* rather than as
    *"something is escaped somewhere"*. The weaker form was written first and mutation testing
    caught it: escaping only the tooltip still leaves ``&quot;`` in the output, so a check for the
    entity passed while both of the other two arguments went through raw.
    """
    from components.variant_detail import _badge_html

    hostile_value = 'Pathogenic" onmouseover="alert(1)'
    hostile_tooltip = 'reviewed" onfocus="alert(2)'
    hostile_label = "<i>ClinVar</i>"

    rendered = _badge_html(hostile_label, hostile_value, "#006400", hostile_tooltip)

    for raw in (hostile_value, hostile_tooltip, hostile_label):
        assert raw not in rendered, f"{raw!r} reached the page unescaped"
    # The words survive as inert text, which is right — what must not survive is the quote that
    # would end an attribute and start a new one.
    assert 'onmouseover="' not in rendered
    assert 'onfocus="' not in rendered
    assert "<i>" not in rendered

    # And with no tooltip, so the value's own escaping is not answered for by the hover's.
    assert hostile_value not in _badge_html("ClinVar", hostile_value, "#006400")


# ---------------------------------------------------------------------------
# 3. RENOVO: a qualifier that is absent must degrade, not render empty
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "score",
    [None, "", ".", "nan", "NaN", "None", float("nan")],
    ids=["absent", "empty", "annotator-dot", "nan-text", "NaN-text", "None-text", "nan-float"],
)
def test_renovo_degrades_to_class_only_when_the_score_says_nothing(score):
    """6,604 real rows across 3 files take this branch, and they are the sharpest rows on the arm.

    The two large files are a 2022 vintage from before the pipeline carried ``PL_score``. They are
    also the **only** germline rows in the corpus with no ``InterVar`` column at all — so on 100%
    of the germline rows where the guideline row can draw no verdict, RENOVO's badge is the only
    pathogenicity claim in the panel, and on exactly those rows its score is missing too. (Issue
    #201 counted 6,597 in 2 files; a third has arrived on disk since, and issue #204
    re-derives the figure through this function.)
    """
    row = {"RENOVO_Class": "HP Benign"}
    if score is not None:
        row["RENOVO_pls"] = score

    badge = renovo_badge(pd.Series(row))
    assert badge.value == "HP Benign"
    for empty in ("(", ")", "nan", "None"):
        assert empty not in badge.value


def test_renovo_draws_its_score_inline_and_rounded():
    """Issue #201 put the score on the badge; issue #214 found it arrives at full float precision.

    A real row rendered as ``RENOVO: IP Benign (0.0246031746031746)``. Three significant figures is
    the panel's own format, and it keeps every class boundary issue #201 measured distinct — the
    rounding never puts a score on the wrong side of the band its own class names.
    """
    badge = renovo_badge(
        pd.Series({"RENOVO_Class": "IP Benign", "RENOVO_pls": "0.0246031746031746"})
    )
    assert badge.value == "IP Benign (0.0246)"


@pytest.mark.parametrize(
    "low,high",
    [("0.00918", "0.00952"), ("0.49830", "0.50254"), ("0.78412", "0.79413")],
    ids=["HP/IP Benign", "LP Benign/LP Pathogenic", "LP/IP Pathogenic"],
)
def test_the_rounding_keeps_the_class_boundaries_apart(low, high):
    """The two sides of a real boundary must not round to the same printed number.

    Otherwise a score and the class beside it would appear to disagree, on the badge that carries
    both.
    """
    below = renovo_badge(pd.Series({"RENOVO_Class": "X", "RENOVO_pls": low})).value
    above = renovo_badge(pd.Series({"RENOVO_Class": "X", "RENOVO_pls": high})).value
    assert below != above


def test_renovo_never_draws_a_score_without_its_class():
    """A qualifier is not homed on its own merits (issue #200, Corollary A).

    The corpus already guarantees it — 0 rows carry a score and no class — so this asserts the
    code cannot invent the case the data does not have.
    """
    assert renovo_badge(pd.Series({"RENOVO_pls": "0.5"})) is None
    assert renovo_badge(pd.Series({"RENOVO_Class": ".", "RENOVO_pls": "0.5"})) is None


# ---------------------------------------------------------------------------
# 4. ESCAT: three tiering-table vintages, and a cell holding several tiers
# ---------------------------------------------------------------------------


def test_the_current_vintage_names_the_acronym_not_the_organ():
    """``ESCAT_TISSUE`` holds the organ and ``ESCAT_CANCER`` the acronym on 2,216 rows.

    *Most specific present* draws the acronym, and the organ still reaches the reader on hover
    because it is a different word from the one drawn.
    """
    badges = escat_badges(
        pd.Series({"ESCAT": "IA", "ESCAT_TISSUE": "BONE_MARROW", "ESCAT_CANCER": "AML"})
    )
    assert [badge.value for badge in badges] == ["IA (AML)"]
    assert "Tissue: Bone marrow." in badges[0].tooltip


def test_the_older_vintage_puts_the_acronym_in_the_tissue_column():
    """329 rows wrote the **acronym** in ``ESCAT_TISSUE`` and left ``ESCAT_CANCER`` at ``.``.

    The rule draws the same label — ``AML`` — under both live vintages, which is what makes *most
    specific present* a rule rather than a guess about which column to trust.
    """
    badges = escat_badges(
        pd.Series({"ESCAT": "IA", "ESCAT_TISSUE": "AML", "ESCAT_CANCER": "."})
    )
    assert [badge.value for badge in badges] == ["IA (AML)"]
    assert "Tissue:" not in badges[0].tooltip


def test_the_oldest_vintage_invents_no_tumour_type_and_borrows_no_gloss():
    """7 somatic rows in 2 files carry an ``ESCAT`` column and **neither** qualifier column.

    They are exactly issue #202's bare ``I``/``II``/``III`` rows, so the oldest vintage has neither
    a graded tier nor a tumour type — and ``ESCAT_DEFINITIONS`` has no entry for those tiers. A
    borrowed gloss is the mistake issue #89 corrected once already, and an invented tumour type
    would be a claim about a disease this file never named.
    """
    badges = escat_badges(pd.Series({"ESCAT": "III"}))
    assert [badge.value for badge in badges] == ["III"]
    assert "predates the graded scale" not in badges[0].tooltip
    assert "predates the ESCAT_TISSUE and ESCAT_CANCER columns" in badges[0].tooltip

    # And where the tier is undefined but a tumour type does exist, the gloss says so rather than
    # reaching for a neighbouring level's words.
    graded = escat_badges(pd.Series({"ESCAT": "III", "ESCAT_CANCER": "AML"}))
    assert "predates the graded scale" in graded[0].tooltip
    for borrowed in ("Hypothetical target", "Ready for routine use", "Investigational"):
        assert borrowed not in graded[0].tooltip


#: The widest real cell in the corpus, from one 2022 germline sample (two files, filtered and
#: unfiltered). Today's ``assign_escat`` cannot write a list at all — it picks one best-scoring
#: match — so this is a legacy shape to survive rather than one the pipeline still produces.
WIDEST_CELL = {
    "ESCAT": "['IIIA', 'IA', 'IIIA', 'IIIA', 'IIIA', 'IIIA', 'IIIA', 'IIIA', 'IIIB']",
    "ESCAT_TISSUE": (
        "['NSCLC', 'MBC', 'MCRC', 'PC', 'MGC', 'PDAC', 'HCC', 'CC', 'HNSCC']"
    ),
    "ESCAT_CANCER": ".",
}


def test_a_multi_valued_cell_groups_by_tier_and_keeps_every_tumour_type():
    """The repeats are **not** duplicates — each ``IIIA`` names a different tumour type.

    Which is what ``IIIA`` means in ``ESCAT_DEFINITIONS``: *benefit is established for this
    alteration, as tier I or II, but in a different tumour type*. De-duplicating the tiers alone
    would silently drop six facts. Grouping costs 3 badges here against 9 for one-per-pair, with
    every fact kept.
    """
    badges = escat_badges(pd.Series(WIDEST_CELL))

    assert len(badges) == 3
    drawn = " ".join(badge.value for badge in badges)
    for tumour in ("NSCLC", "MBC", "MCRC", "PC", "MGC", "PDAC", "HCC", "CC", "HNSCC"):
        assert tumour in drawn, f"{tumour} was dropped by grouping"
    assert drawn.count("IIIA") == 1, "the tiers were not grouped"
    assert "[" not in drawn, "the list repr reached the badge"


def test_the_tiers_are_ordered_by_escat_s_own_ladder():
    """``ESCAT_OPTIONS`` needs no new ordering and does not disagree with the annotator.

    It picks the same winner as ``assign_escat``'s own ``min()`` on 48 of 48 multi-valued cells.
    """
    badges = escat_badges(pd.Series(WIDEST_CELL))
    assert [badge.value.split(" ")[0] for badge in badges] == ["IA", "IIIA", "IIIB"]


def test_a_tier_pairs_with_the_tissue_positionally_and_never_with_the_cancer():
    """Same-length parallel lists on all 48 list cells; ``ESCAT_CANCER`` is a scalar on all 48."""
    facts = escat_facts(pd.Series(WIDEST_CELL))
    assert ("IA", "MBC") in facts
    assert ("IIIB", "HNSCC") in facts
    assert all(tumour != "." for _tier, tumour in facts)


def test_a_single_item_list_is_not_a_special_case():
    """22 of the 48 list cells hold one item, so once normalised they render like a plain cell."""
    listed = escat_badges(pd.Series({"ESCAT": "['IA']", "ESCAT_TISSUE": "['MBC']"}))
    plain = escat_badges(pd.Series({"ESCAT": "IA", "ESCAT_TISSUE": "MBC"}))
    assert [badge.value for badge in listed] == [badge.value for badge in plain] == ["IA (MBC)"]


def test_a_blank_escat_draws_nothing_and_is_never_called_no_actionable_target():
    """``.`` is the annotator's own *no-match* marker, and ESCAT's tier X is a claim never made.

    ``assign_escat`` writes ``.`` on every row before matching and there is no skip flag, so a
    blank cell means *not in ESMO's table* — a much weaker claim than *no evidence that the
    alteration is therapeutically actionable*, which is what tier X asserts (issue #202).
    """
    assert escat_badges(pd.Series({"ESCAT": "."})) == []
    assert escat_badges(pd.Series({"ESCAT": ""})) == []

    sentences = " ".join(clinical_absences(pd.Series({"ESCAT": "."})))
    assert "actionab" not in sentences.lower()
    assert " X" not in sentences


# ---------------------------------------------------------------------------
# 5. The empty state: exactly two captions, and never RENOVO
# ---------------------------------------------------------------------------

#: The five member states issue #217 measured, with the row counts over 180 files / 238,266 rows.
#: The two rarest are named rather than sampled: state 3 is 16 rows, and state 1 co-occurs *only*
#: with state 4.
EMPTY_STATES = {
    "clinvar-column-absent": ({"ESCAT": "."}, 18633),
    "clinvar-blank-no-status": ({"ClinVar_VCF_CLNSIG": "", "ESCAT": "."}, 50061),
    "clinvar-blank-with-status": (
        {
            "ClinVar_VCF_CLNSIG": "",
            "ClinVar_VCF_CLNREVSTAT": "criteria_provided,_single_submitter",
            "ESCAT": ".",
        },
        16,
    ),
    "escat-column-absent": ({"ClinVar_VCF_CLNSIG": ""}, 20990),
    "the-somatic-shape": (
        {"ClinVar_VCF_CLNSIG": "", "ESCAT": ".", "CIViC_Evidence_Level": "."},
        47095,
    ),
    "an-unplaceable-file": ({"Hugo_Symbol": "TP53"}, 21615),
}


@pytest.mark.parametrize("state", sorted(EMPTY_STATES))
def test_exactly_two_captions_fire_on_every_empty_row(state):
    """The invariant issue #217 handed this ticket: never one, never three.

    ClinVar is in one of three states and ESCAT in one of two on all 68,710 empty rows, and
    ClinVar-absent co-occurs *only* with ESCAT-absent. Two sentences is what makes the row never
    silent and never repetitive — and it is what replaces the single global claim that was false
    on 47,095 somatic variants.
    """
    cells, _rows = EMPTY_STATES[state]
    row = pd.Series(cells)

    assert clinical_badges(row) == [], "this row is supposed to be an empty one"
    notes = clinical_absences(row)
    assert len(notes) == 2, f"{state} fired {len(notes)} captions: {notes}"


@pytest.mark.parametrize("state", sorted(EMPTY_STATES))
def test_the_empty_row_never_names_renovo(state):
    """On somatic RENOVO is inapplicable rather than missing; on an unplaceable file, unknowable.

    ``run_germline_renovo`` only ever runs on germline, so listing it among what is absent would
    report an inapplicable annotation as a failed one (issue #201). On a file the header cannot
    place the arm is unknown, so neither *missing* nor *germline-only* is honest. On germline the
    row never empties, so its case is unreachable — which is why silence is right in all three.
    """
    cells, _rows = EMPTY_STATES[state]
    sentences = " ".join(clinical_absences(pd.Series(cells)))
    assert not RENOVO_MENTION.search(sentences), sentences
    assert ABSENCE_NEVER_NAMES[0] == "RENOVO_Class"


@pytest.mark.parametrize("state", sorted(EMPTY_STATES))
def test_each_empty_row_says_one_thing_about_each_member_it_speaks_for(state):
    """One ClinVar sentence and one ESCAT sentence, driven off the same list they are written from.

    Counting to two would pass on two ClinVar sentences and no ESCAT one, which is the failure the
    invariant exists to exclude.
    """
    cells, _rows = EMPTY_STATES[state]
    notes = clinical_absences(pd.Series(cells))

    assert set(ABSENCE_SPEAKS_FOR) == {"ClinVar_VCF_CLNSIG", "ESCAT"}
    assert sum(1 for note in notes if "ClinVar" in note) == 1, notes
    assert sum(1 for note in notes if "ESCAT" in note) == 1, notes


def test_the_two_blanks_are_said_as_different_kinds_of_nothing():
    """ClinVar wrote nothing at all; ESCAT's annotator ran and recorded a no-match.

    Saying both in the shipped *"{tool} recorded no {verdict} for this variant"* shape would read
    identically and flatten exactly the distinction these sentences exist to draw.
    """
    clinvar, escat = clinical_absences(pd.Series({"ClinVar_VCF_CLNSIG": "", "ESCAT": "."}))
    assert clinvar == "ClinVar holds no classification for this variant."
    assert escat == "ESMO's ESCAT table has no entry for this variant."
    assert clinvar != escat


def test_a_review_status_can_outlive_its_significance_and_is_quoted_verbatim():
    """16 rows carry a drawable ``CLNREVSTAT`` with a blank ``CLNSIG``, all inside the empty case.

    So *"this variant has no ClinVar record"* would be false there. And the status is quoted in
    full rather than as stars: the star hierarchy abbreviates a status that *qualifies a
    significance*, and here there is no significance for it to qualify.
    """
    notes = clinical_absences(
        pd.Series(
            {
                "ClinVar_VCF_CLNSIG": "",
                "ClinVar_VCF_CLNREVSTAT": "criteria_provided,_single_submitter",
                "ESCAT": ".",
            }
        )
    )
    assert "criteria provided, single submitter" in notes[0]
    assert "★" not in notes[0] and "☆" not in notes[0]
    assert "no ClinVar record" not in " ".join(notes)


def test_the_no_match_sentence_cannot_name_a_tumour_type():
    """On every one of the 47,095 empty somatic rows the two qualifier columns are ``.`` too.

    2 distinct pairs corpus-wide. There is no tumour type recorded to name, which is what keeps
    this sentence clear of the tumour type issue #214 draws where a tier *does*.
    """
    notes = clinical_absences(
        pd.Series(
            {
                "ClinVar_VCF_CLNSIG": "",
                "ESCAT": ".",
                "ESCAT_TISSUE": ".",
                "ESCAT_CANCER": ".",
            }
        )
    )
    assert notes[1] == "ESMO's ESCAT table has no entry for this variant."


def test_the_deleted_caption_is_gone_from_the_panel(monkeypatch):
    """*"No clinical annotations available for this variant."* is deleted, not reworded.

    It was false on all 47,095 somatic empty rows, where the guideline row drew a badge on 100% of
    them — a variant reading ``AMP/ASCO/CAP (CancerVar): Tier III`` was being told it had no
    clinical annotation. The fix is not a better global sentence; it is making no global claim.

    Asked of the arguments to the calls that put text on screen, **not** of the modules' text.
    Every module here explains at length what the sentence was and why it went, so a text guard
    would be failed by the record of its own fix — the mirror of the failure mode where a text
    guard is *satisfied* by a file's own header comment.
    """
    for path in sorted((STREAMLIT_APP / "components").glob("*.py")):
        drawn_literals = {
            argument.value
            for node in ast.walk(ast.parse(path.read_text()))
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in ("caption", "markdown", "write", "info", "warning")
            for argument in node.args
            if isinstance(argument, ast.Constant) and isinstance(argument.value, str)
        }
        assert not any(
            "No clinical annotations" in literal for literal in drawn_literals
        ), path

    cells, _rows = EMPTY_STATES["the-somatic-shape"]
    row = a_variant(**{key: value for key, value in cells.items()})
    row = row.drop(labels=["RENOVO_Class", "RENOVO_pls"])
    records = panel_records(monkeypatch, row)
    assert not any("No clinical annotations" in text for _module, text in records)


def test_the_row_draws_badges_or_captions_and_never_both(monkeypatch):
    """The captions are the *empty* state, not a running commentary beside the badges.

    Issue #187's shape — name absence per state, beneath the badges rather than instead of them —
    is the guideline row's, where the common case is badges and an absence at the same time. Here
    it is not: ClinVar is blank on 37.4% of germline rows where RENOVO always draws, so a per-member
    absence beside every badge would fire on a third of the arm with nothing missing.
    """
    populated = a_variant()
    assert clinical_badges(populated), "the fixture is supposed to draw badges"

    from components import variant_detail

    drawn: list = []
    monkeypatch.setattr(variant_detail, "st", _Recorder("panel", drawn))
    variant_detail._render_clinical_badges(populated)

    assert len(drawn) == 1, drawn
    assert "No " not in drawn[0][1]


# ---------------------------------------------------------------------------
# The membership itself
# ---------------------------------------------------------------------------


def test_the_clinical_row_holds_exactly_the_three_map_199_left_it():
    """The row's membership, asserted so a fourth annotation cannot arrive unremarked.

    Two left it during this map and each for its own reason: ``CIViC_Evidence_Level`` because a
    level attaches to an evidence item rather than a variant, and ``am_class`` because it is
    neither a database opinion nor a guideline verdict and got a section of its own.
    """
    assert [column for column, _label, _color in CLINICAL_ROW] == [
        "ClinVar_VCF_CLNSIG",
        "RENOVO_Class",
        "ESCAT",
    ]
    assert [label for _column, label, _color in CLINICAL_ROW] == ["ClinVar", "RENOVO", "ESCAT"]


def test_the_row_has_no_heading():
    """Issue #217: it stays unnamed, and that was a decision rather than an omission.

    A heading here would have to be true of another laboratory's conclusion, a model this pipeline
    ran and an ESMO actionability tier at once. If a future ticket wants one, it is a new decision.
    """
    from components import variant_detail

    tree = ast.parse((STREAMLIT_APP / "components" / "variant_detail.py").read_text())
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_render_clinical_badges"
    )
    headings = [
        node.args[0].value
        for node in ast.walk(function)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "markdown"
        and node.args
        and isinstance(node.args[0], ast.Constant)
    ]
    assert headings == [], f"the clinical row grew a heading: {headings}"
    assert variant_detail._render_clinical_badges is not None


def test_the_module_docstring_no_longer_calls_the_row_a_list_of_databases():
    """Issue #201: *"what third-party databases say"* is provably false of RENOVO.

    RENOVO passes the three-leg test's *computed here* leg — ``run_germline_renovo`` runs in this
    pipeline — which is exactly what a third-party database fails. There is no heading to carry the
    correction, so it rides in the module docstring, and the correction must not repeat the
    falsehood while denying it.
    """
    from components import variant_detail

    docstring = variant_detail.__doc__
    assert "what third-party databases say" not in docstring
    assert "moving ClinVar empties the database row" not in docstring
    assert "still drawn twice" not in docstring
    # The correction is made positively rather than as a denial of the old sentence.
    assert "three-leg test" in docstring
