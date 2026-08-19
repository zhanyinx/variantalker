"""The clinical row: three annotations, each drawn once, and what it says when none draws.

Map #199 sorted every annotation the variant panel draws into one home apiece. This module holds the
three that landed in :func:`components.variant_detail._render_clinical_badges` —
``ClinVar_VCF_CLNSIG``, ``RENOVO_Class`` and ``ESCAT`` — as *values*: what each badge says, what its
hover carries, and what the row says on a variant where none of them speaks. The rendering stays in
``variant_detail``; the vocabularies are here, so a guard can drive them directly rather than
through a Streamlit call.

The row has no heading, and that is a decision
----------------------------------------------
#217 measured that the words *"clinical annotations"* reached a user from exactly **one** rendered
string — the empty-state caption — so on 100% of germline rows and 54.6% of somatic ones the row
drew with no name at all. Naming it would have meant one heading true of another laboratory's
conclusion, a model this pipeline ran, and an ESMO actionability tier at once. It stays unnamed, and
each badge's own label says what its claim is.

What is *not* a database opinion, and why the three sit together anyway
-----------------------------------------------------------------------
The row used to be described as *"what third-party databases say"*. That is false of two of its
three members, and the correction is #201's:

* **ClinVar** is a third-party aggregate of other laboratories' submissions — no guideline scale,
  not computed here, no evidence vector in the MAF.
* **RENOVO** is an in-silico predictor *this pipeline ran* (``workflows/variantalker.nf`` invokes
  ``run_germline_renovo``), on its own six-value scale, with one number behind it. Passing
  *computed here* is exactly what a third-party database fails, so it is a third category rather
  than either side of that binary.
* **ESCAT** is ESMO's published actionability scale, assigned by a lookup table this pipeline
  carries — a published scale, but not a verdict this pipeline's classifier reached.

What they share is the negative that decides the home: #200's three-leg test — a guideline scale,
computed here, with the evidence in the MAF — is failed by all three, so none of them belongs under
*"Guideline Classifications"*. Colour follows from the same fact: in this row a badge's colour
identifies its **source**, so each of the three is flat, while in the guideline row colour carries
the tier.

``CIViC_Evidence_Level`` and ``am_class`` are deliberately absent. #202 dropped CIViC's badge — an
evidence level attaches to an evidence *item*, not to a variant, and 133 of the 140 firing cells
held a Python list drawn verbatim — leaving ``CIViC_Variant_URL`` in external links as the honest
surface. #203 moved AlphaMissense to :mod:`components.alphamissense`, an in-silico predictor's own
section rather than a badge. Both columns are still read by the filters, which this map does not
touch.
"""

import ast
import re
from typing import NamedTuple, Optional

import pandas as pd

from config.missing_values import says_nothing
from config.vocabularies import ESCAT_DEFINITIONS, ESCAT_OPTIONS


class Badge(NamedTuple):
    """One badge of the clinical row.

    ``tooltip`` is the ``title=`` attribute's text, or ``None`` for a badge that carries no
    hover. It is escaped where the badge is built into HTML — see
    :func:`components.variant_detail._badge_html` — never here, so a caller cannot receive
    half-escaped text and escape it twice.
    """

    label: str
    value: str
    color: str
    tooltip: Optional[str] = None


#: ``(column, label, colour)`` for the three members, in the order the row draws them.
#:
#: One flat colour each, because in this row colour says *which source this is* (issue #201).
#: The table is the row's membership: ``tests/test_clinical_badges.py`` drives the
#: drawn-exactly-once guard off it, so a fourth annotation cannot be added to the row without
#: that guard being told which one it is.
CLINICAL_ROW = (
    ("ClinVar_VCF_CLNSIG", "ClinVar", "#006400"),
    ("RENOVO_Class", "RENOVO", "#2F4F4F"),
    ("ESCAT", "ESCAT", "#1E90FF"),
)

#: The three colours, unpacked from :data:`CLINICAL_ROW` rather than looked up by name at each
#: use. A dict subscripted with a column-name literal is indistinguishable from a *row* lookup to
#: ``tests/test_column_spelling.py``'s AST walk, and a guard about how a column reaches a row
#: should not have to tell the two apart.
#: Read by position rather than unpacked, so a fourth member added to the row is caught by
#: ``tests/test_clinical_badges.py``'s membership assertion rather than by an unpacking error at
#: import — a module that refuses to load says nothing about what the row should hold.
_CLINVAR_COLOR = CLINICAL_ROW[0][2]
_RENOVO_COLOR = CLINICAL_ROW[1][2]
_ESCAT_COLOR = CLINICAL_ROW[2][2]


# ---------------------------------------------------------------------------
# ClinVar — the significance, and the review status that qualifies it
# ---------------------------------------------------------------------------

#: ClinVar's own star hierarchy, as row counts over the real corpus found it (issue #200).
#:
#: The dev confirmed the levels are **ordered**, on ClinVar's published star scale, and that that
#: hierarchy is this panel's authority. Four facts about this table are measured rather than
#: guessed, and each is a way a shorter one would have been wrong:
#:
#: * **Vintage drift spells one level two ways, and it happens twice.**
#:   ``conflicting_classifications`` (1,301 rows) and ``conflicting_interpretations`` (243) are
#:   ClinVar's new and old names for the same one-star level, which issue #200 measured; and
#:   ``no_classification_for_the_single_variant`` (804) and
#:   ``no_interpretation_for_the_single_variant`` (**50**) are the same pair one level down, which
#:   it did not — that tenth value was found re-deriving the vocabulary for issue #204
#:   (``docs/wayfinder/issue-204/``) and was falling through to :data:`_UNRECOGNISED_STATUS`,
#:   printing its raw string where the stars belong. Both pairs are folded here, because anything
#:   that orders or colours by these strings would otherwise sort one level two ways.
#: * **``practice_guideline`` is unobserved and still handled.** It appears on **0** rows; these
#:   MAFs were annotated 2022–2024, and unobserved is not impossible. A corpus-driven guard cannot
#:   reach it, which is why ``tests/test_clinical_badges.py`` names it.
#: * **The mapping is total over what the corpus holds** — asserted rather than believed, since
#:   believing it is what missed the tenth value: ``docs/wayfinder/issue-204/measure204.py``
#:   re-walks every ``.maf`` on disk and reports anything unmapped.
#: * **An unmapped value is still reachable and must stay handled.** Two of these ten were
#:   themselves unmapped a moment ago, which is the argument for :data:`_UNRECOGNISED_STATUS`
#:   rather than against it.
CLINVAR_REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_for_the_single_variant": 0,
    "no_interpretation_for_the_single_variant": 0,
    "no_classification_provided": 0,
    "no_assertion_provided": 0,
}

#: How many stars the scale has, so a zero-star badge can draw the scale it scored nothing on.
_STARS_ON_THE_SCALE = 4

#: What an unmapped review status is said to be. **Never zero stars** (issue #200): zero stars is
#: ClinVar asserting *no assertion criteria were provided*, and an unrecognised value is MAFigate
#: not knowing which level this is. Guessing the lower of the two promises a reading this panel
#: cannot support — the same category error issue #210 fixed one row up.
_UNRECOGNISED_STATUS = (
    "MAFigate does not recognise this ClinVar review status, so it is shown in full "
    "rather than as stars."
)


def _text(row: pd.Series, column: str) -> str:
    """The cell as a stripped string, or ``""`` when the column or the value is not there."""
    if column not in row.index:
        return ""
    value = row[column]
    return "" if value is None else str(value).strip()


def clinvar_stars(review_status: str) -> Optional[str]:
    """The star glyphs for a review status, or ``None`` when the scale does not know it.

    Args:
        review_status: the ``ClinVar_VCF_CLNREVSTAT`` cell, stripped.

    Returns:
        str: four glyphs — filled stars earned, hollow stars not — or ``None`` for a value
        outside :data:`CLINVAR_REVIEW_STARS`.

    **Zero stars draws ``☆☆☆☆``, and that is the point of drawing the whole scale.** Issue #200
    measured 4,788 rows holding a status that earns no stars against 1,637 rows carrying a
    significance with *no review status at all*, and required that the two not render alike: one
    says *ClinVar assessed this with no assertion criteria*, the other says *MAFigate does not
    know how this was reviewed*. A badge showing nothing where the stars go says both at once, so
    an earned zero shows the empty scale and an absent status shows no scale.
    """
    stars = CLINVAR_REVIEW_STARS.get(review_status)
    if stars is None:
        return None
    return "★" * stars + "☆" * (_STARS_ON_THE_SCALE - stars)


def _review_status_words(review_status: str) -> str:
    """The status as a sentence reads it — ClinVar's underscores are word joins, not spelling."""
    return review_status.replace("_", " ").strip()


def clinvar_badge(row: pd.Series) -> Optional[Badge]:
    """ClinVar's significance for this variant, with its review status, or ``None``.

    The review status **travels with the significance** rather than staying in the guideline row
    (issue #200, Corollary A): it qualifies ClinVar's claim and makes none of its own, so stranding
    it beside a significance that had moved would leave it qualifying nothing.

    It is abbreviated to stars on the badge and carried in full on hover. Three of the nine
    statuses the corpus holds (887 rows) say *no classification was provided* while the badge
    beside them shows a significance anyway — which is exactly what a zero-star ClinVar aggregate
    is, and is the pairing most likely to read as a bug, so the full sentence on hover matters most
    there.
    """
    significance = _text(row, "ClinVar_VCF_CLNSIG")
    if says_nothing(significance):
        return None

    status = _text(row, "ClinVar_VCF_CLNREVSTAT")
    if says_nothing(status):
        return Badge(
            "ClinVar",
            significance,
            _CLINVAR_COLOR,
            "This file records no ClinVar review status for this variant, so MAFigate cannot "
            "say how this classification was reviewed.",
        )

    stars = clinvar_stars(status)
    if stars is None:
        return Badge(
            "ClinVar",
            f"{significance} ({status})",
            _CLINVAR_COLOR,
            f"ClinVar review status: {_review_status_words(status)}. {_UNRECOGNISED_STATUS}",
        )
    return Badge(
        "ClinVar",
        f"{significance} {stars}",
        _CLINVAR_COLOR,
        f"ClinVar review status: {_review_status_words(status)} — "
        f"{stars.count('★')} of {_STARS_ON_THE_SCALE} stars on ClinVar's own scale.",
    )


# ---------------------------------------------------------------------------
# RENOVO — the class, and the score it is a band of
# ---------------------------------------------------------------------------

#: How a score is rounded onto the badge. Three significant figures, the same as
#: :func:`components.predictor_context._format_score` and the ClinGen table, so the panel reads as
#: one panel.
#:
#: Issue #214 surfaced the need for a decision: ``RENOVO_pls`` carries full float precision in the
#: MAF and a real row rendered as ``RENOVO: IP Benign (0.0246031746031746)``. Three significant
#: figures is enough to keep every class boundary issue #201 measured distinct — ``HP Benign``
#: tops out at 0.00918 where ``IP Benign`` starts at 0.00952, and ``LP Benign`` at 0.49830 where
#: ``LP Pathogenic`` starts at 0.50254 — so the rounding never puts a score on the wrong side of
#: the band its own class names.
def _format_score(number: float) -> str:
    return f"{number:.3g}" if abs(number) < 100 else f"{number:.1f}"


def renovo_badge(row: pd.Series) -> Optional[Badge]:
    """RENOVO's class for this variant, with ``RENOVO_pls`` inline, or ``None``.

    The class **stays in this row** (issue #201): it fails the guideline row's scale and evidence
    legs, and a predictor table was rejected because all three the panel draws are PP3/BP4
    surfaces and RENOVO has no such calibration.

    **The score degrades rather than renders empty.** ``RENOVO_pls`` is absent on **6,604 real
    rows across 3 files** — issue #201 measured 6,597 in 2, and a third file has arrived on disk
    since; re-derived through this function in ``docs/wayfinder/issue-204/measure204.py``, which
    found **0** badges rendering an empty score across the whole corpus. The two large files are a
    2022 annotation vintage from before the pipeline carried ``PL_score``, and they are the *only*
    germline rows where the guideline row can draw no verdict at all — so RENOVO's badge is the
    whole of the pathogenicity claim on screen there. It draws ``RENOVO: HP Benign`` on them, never
    ``(nan)``, ``()`` or a bare ``(``. The corpus holds no row with a score and no class, so the
    qualifier never strands.

    A score the file spells in a way ``float`` refuses is drawn as the file spells it rather than
    dropped: unreadable is a fact about the cell, and no such cell exists in the corpus.
    """
    classification = _text(row, "RENOVO_Class")
    if says_nothing(classification):
        return None

    color = _RENOVO_COLOR
    score = _text(row, "RENOVO_pls")
    if says_nothing(score):
        return Badge("RENOVO", classification, color)
    try:
        return Badge("RENOVO", f"{classification} ({_format_score(float(score))})", color)
    except (TypeError, ValueError):
        return Badge("RENOVO", f"{classification} ({score})", color)


# ---------------------------------------------------------------------------
# ESCAT — one badge per tier, carrying the tumour types it was established in
# ---------------------------------------------------------------------------
#
# Ported from ``docs/wayfinder/issue-214/candidates214.py``, the rendering the dev chose from five
# drawn against eight real MAF rows. Issue #214's spec says port rather than re-derive, and the
# names below are that prototype's.

#: Where a tier sits on ESCAT's own ladder. ``ESCAT_OPTIONS`` needs no new ordering and does not
#: disagree with the annotator: it picks the same winner as ``assign_escat``'s own ``min()`` on
#: **48 of 48** multi-valued cells in the corpus.
_TIER_RANK = {tier: index for index, tier in enumerate(ESCAT_OPTIONS)}


def _as_list(text: str) -> Optional[list]:
    """A Python list repr read back as a list, or ``None`` when the cell is not one.

    45 germline cells and 3 somatic ones hold one — a 2022 sample, two files, filtered and
    unfiltered. Today's ``assign_escat`` (``bin/utils.py``) cannot write a list at all: it picks
    one best-scoring match with ``.values[best_score_idx]``. So this is a legacy shape to survive
    rather than one the pipeline still produces, and it is the shape that reached the page
    verbatim before issue #214 — 375 characters on the worst cell.
    """
    if not text.startswith("["):
        return None
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError, MemoryError, RecursionError):
        return None
    return [str(item).strip() for item in parsed] if isinstance(parsed, list) else None


def _rank(tier: str) -> int:
    """A tier the ladder does not know sorts last, and still draws.

    7 somatic rows in 2 files carry a bare ``I``/``II``/``III`` from an older tiering table, and
    those are exactly the rows with no qualifier columns at all.
    """
    return _TIER_RANK.get(tier, len(_TIER_RANK) + 1)


def _pretty(label: str) -> str:
    """``BONE_MARROW`` reads as *Bone marrow*; an acronym is left exactly as written.

    Only a label carrying an underscore is touched. The acronyms and the organ names are not
    separable by shape — ``LUNG`` and ``GIST`` are the same shape and different kinds of word —
    and mangling ``HNSCC`` into *Hnscc* would be worse than shouting.
    """
    if "_" not in label:
        return label
    return label.replace("_", " ").lower().capitalize()


def escat_facts(row: pd.Series) -> list:
    """Every ``(tier, tumour type)`` this row asserts, strongest tier first.

    The tumour type is **the most specific spelling this row carries** — ``ESCAT_CANCER`` when it
    says something, else ``ESCAT_TISSUE``, else ``None``. That last case is a real state rather
    than a default: 2 files in the corpus carry an ``ESCAT`` column and neither qualifier column,
    and no tumour type may be invented for them.

    The rule exists because the pair arrives in **three tiering-table vintages** (issue #214):
    today ``ESCAT_TISSUE`` holds the organ and ``ESCAT_CANCER`` the acronym (2,216 rows); an older
    table put the **acronym in ``ESCAT_TISSUE``** and left ``ESCAT_CANCER`` at ``.`` (329 rows); an
    older one still wrote neither column (7 rows). *Most specific present* draws the same label —
    ``AML`` — under both live vintages.

    **Pairs de-duplicate as pairs, never as tiers.** The worst cell holds
    ``['IIIA', 'IA', 'IIIA', ...]`` beside nine *different* tissues, so each ``IIIA`` names a
    different tumour type — which is what ``IIIA`` means in :data:`ESCAT_DEFINITIONS`: benefit
    established for this alteration, as tier I or II, but in a different tumour type. Collapsing
    the repeated tiers alone would silently drop six facts. Tiers pair with ``ESCAT_TISSUE``
    positionally — same-length parallel lists on all 48 list cells — and **never** with
    ``ESCAT_CANCER``, which is a scalar ``.`` on all 48.
    """
    raw = _text(row, "ESCAT")
    if says_nothing(raw):
        return []
    tiers = _as_list(raw) or [raw]
    cancer, tissue = _text(row, "ESCAT_CANCER"), _text(row, "ESCAT_TISSUE")
    cancers, tissues = _as_list(cancer), _as_list(tissue)

    def tumour_type(index):
        for parsed, scalar in ((cancers, cancer), (tissues, tissue)):
            if parsed is not None:
                if len(parsed) == len(tiers) and not says_nothing(parsed[index]):
                    return parsed[index]
            elif len(tiers) == 1 and not says_nothing(scalar):
                return scalar
        return None

    pairs = []
    for index, tier in enumerate(tiers):
        pair = (tier, tumour_type(index))
        if pair not in pairs:
            pairs.append(pair)
    pairs.sort(key=lambda pair: (_rank(pair[0]), pair[1] or ""))
    return pairs


def _gloss(tier: str, organ: Optional[str] = None) -> str:
    """The hover text: ESCAT's own words for this tier, and the organ if it adds anything.

    A tier the scale does not define gets a sentence saying so rather than a borrowed gloss — the
    mistake issue #89 corrected once already on the help page. :data:`ESCAT_DEFINITIONS` has no
    entry for ``I``/``II``/``III``, so this branch is reachable from the corpus rather than
    defensive.
    """
    level = ESCAT_DEFINITIONS.get(tier)
    if level is None:
        parts = [
            f"ESCAT {tier}: this file's tiering table predates the graded scale, so the level "
            f"carries no A/B/C grade and no definition here."
        ]
    else:
        parts = [f"ESCAT {tier} — {level.group}. {level.evidence} {level.implication}"]
    if organ:
        parts.append(f"Tissue: {_pretty(organ)}.")
    return " ".join(parts)


def _organ_if_different(row: pd.Series, drawn: str) -> Optional[str]:
    """``ESCAT_TISSUE`` when it is a different word from the label already on the badge."""
    tissue = _text(row, "ESCAT_TISSUE")
    if _as_list(tissue) is not None or says_nothing(tissue) or tissue == drawn:
        return None
    return tissue


def _no_tumour_type_note(tier: str) -> str:
    return (
        f"ESCAT {tier}. This file carries no tumour type for the tier — its tiering table "
        f"predates the ESCAT_TISSUE and ESCAT_CANCER columns."
    )


def escat_badges(row: pd.Series) -> list:
    """One badge per distinct tier, with the tumour types it was established in inside it.

    Inline rather than on hover or on a line beneath: an ESCAT tier is defined *per tumour type*
    and is unreadable without one (issue #202), and this is the habit issue #201 chose for
    RENOVO's score. Grouping is what makes the shape affordable — the widest real cell costs 3
    badges here against 9 for one-per-pair, with every fact kept — and **22 of the 48 list cells
    hold a single item**, so a list stops being a special case at all.
    """
    grouped = {}
    for tier, tumour in escat_facts(row):
        grouped.setdefault(tier, []).append(tumour)

    color = _ESCAT_COLOR
    badges = []
    for tier, tumours in grouped.items():
        named = [_pretty(tumour) for tumour in tumours if tumour]
        if not named:
            badges.append(Badge("ESCAT", tier, color, _no_tumour_type_note(tier)))
            continue
        organ = _organ_if_different(row, tumours[0]) if len(named) == 1 else None
        badges.append(
            Badge("ESCAT", f"{tier} ({', '.join(named)})", color, _gloss(tier, organ))
        )
    return badges


# ---------------------------------------------------------------------------
# The row
# ---------------------------------------------------------------------------


def clinical_badges(row: pd.Series) -> list:
    """Every badge this row draws, in :data:`CLINICAL_ROW` order."""
    badges = []
    clinvar = clinvar_badge(row)
    if clinvar is not None:
        badges.append(clinvar)
    renovo = renovo_badge(row)
    if renovo is not None:
        badges.append(renovo)
    badges.extend(escat_badges(row))
    return badges


def clinical_absences(row: pd.Series) -> list:
    """What the row says on a variant where it draws nothing (issue #217).

    **Exactly two sentences fire on every empty row — never one, never three.** ClinVar is in one
    of three states and ESCAT in one of two, and each pair is a different fact, so the row makes
    no global claim at all. That is the whole of the fix: the sentence this replaces,
    *"No clinical annotations available for this variant."*, was false on all **47,095** somatic
    empty rows (45.4% of the arm; 0.0% of germline), where the guideline row drew a badge on
    **100%** of them — a variant reading ``AMP/ASCO/CAP (CancerVar): Tier III`` was being told it
    had no clinical annotation. A better global sentence was not the fix; making no global claim
    was.

    The states, with the row counts issue #217 measured over 180 files and 238,266 rows:

    ==  ==========================================================  ========
    #   state                                                          rows
    ==  ==========================================================  ========
    1   ``ClinVar_VCF_CLNSIG`` column absent                          18,633
    2   present, blank, no review status                              50,061
    3   present, blank, **review status present**                         16
    4   ``ESCAT`` column absent                                       20,990
    5   present and blank                                             47,720
    ==  ==========================================================  ========

    Four things this is right about that a hand-written sentence was not:

    * **The two blanks are different kinds of nothing.** ClinVar wrote nothing at all; ESCAT's
      annotator ran and recorded a no-match — ``assign_escat`` sets ``ESCAT = "."`` on every row
      before matching and there is no skip flag, so the ``.`` is the annotator's own marker
      (issue #202, confirmed cell by cell). Saying both in the shipped *"{tool} recorded no
      {verdict}"* shape would have flattened exactly the distinction these sentences exist to draw.
    * **A review status can outlive its significance.** 16 rows carry a drawable
      ``ClinVar_VCF_CLNREVSTAT`` with a blank ``CLNSIG``, and on all 16 the row is empty — so
      *"this variant has no ClinVar record"* would be false there. State 3 quotes the status
      **verbatim, not as stars**: the star hierarchy abbreviates a status that qualifies a
      significance, and here there is no significance for it to qualify.
    * **RENOVO is never named.** On somatic it is inapplicable rather than missing —
      ``run_germline_renovo`` only ever runs on germline, so listing it among what is absent would
      report an inapplicable annotation as a failed one (issue #201). On a file the header cannot
      place the arm is unknown, so neither *missing* nor *germline-only* is honest, and it stays
      silent there too. On germline the row never empties, so its case is unreachable.
    * **The sentence cannot name the tumour type.** On every one of the 47,095 empty somatic rows
      ``ESCAT_TISSUE`` and ``ESCAT_CANCER`` are ``.`` too, or their columns are absent — 2 distinct
      pairs corpus-wide. There is no tumour type recorded to name, which is what keeps this clear
      of the tumour type issue #214 draws where a tier *does*.

    On the 28 files the header cannot place, both columns are absent and this says so — states 1
    and 4, no arm needed. The sibling guideline row is silent there because it has no classifier to
    name; these two absences are properties of the file and need no classifier, so the silence that
    row earned is not inherited here.
    """
    notes = []

    if "ClinVar_VCF_CLNSIG" not in row.index:
        notes.append(
            "This file carries no `ClinVar_VCF_CLNSIG` column, so no ClinVar classification "
            "was read for this variant."
        )
    else:
        status = _text(row, "ClinVar_VCF_CLNREVSTAT")
        if says_nothing(status):
            notes.append("ClinVar holds no classification for this variant.")
        else:
            notes.append(
                f"ClinVar records a review status ({_review_status_words(status)}) for this "
                "variant but no classification."
            )

    if "ESCAT" not in row.index:
        notes.append("This file carries no `ESCAT` column.")
    else:
        notes.append("ESMO's ESCAT table has no entry for this variant.")

    return notes


#: The three columns whose *absence* :func:`clinical_absences` may speak about, and the one it may
#: not. Held as data so ``tests/test_clinical_badges.py`` can assert the invariant — exactly one
#: ClinVar sentence and exactly one ESCAT sentence per empty row — off the same list the sentences
#: are written from, rather than by counting a hand-kept number.
ABSENCE_SPEAKS_FOR = ("ClinVar_VCF_CLNSIG", "ESCAT")
ABSENCE_NEVER_NAMES = ("RENOVO_Class", "RENOVO_pls", "RENOVO")

#: Matches the word RENOVO however a sentence might reach for it, for that guard.
RENOVO_MENTION = re.compile(r"renovo", re.IGNORECASE)
